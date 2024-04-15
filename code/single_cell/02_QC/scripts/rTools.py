import sys
sys.path.append("/share/home/yzwl_zhangchao/Project/soybean_sn/02_QC/scripts")
import pandas as pd
import numpy as np
import scipy
from scipy import sparse as sp
import loguru
from loguru import logger
import basic
import functools
import anndata
import rpy2
import rpy2.robjects as ro
import rpy2.robjects.packages
from rpy2 import rinterface_lib
from rpy2.robjects.packages import importr
import inspect
from inspect import signature
from typing import (
    Dict,
    List,
    Optional,
    Union,
    Sequence,
    Literal,
    Any,
    Tuple,
    Iterator,
    Mapping,
    Callable,
)
R = ro.r

def rpy2_check(func):
    """Decorator to check whether rpy2 is installed at runtime"""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            import rpy2
        except ImportError:
            raise ImportError("Please install rpy2 package.")
        return func(*args, **kwargs)

    return wrapper

def anndata2ri_check(func):
    """Decorator to check whether anndata2ri is installed at runtime"""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            import anndata2ri
        except ImportError:
            raise ImportError("Please install anndata2ri package.")
        return func(*args, **kwargs)

    return wrapper

def rcontext(func):
    """
    A decorator to run a function in an R context.
    `rEnv` parameter will be auto updated
    """

    @functools.wraps(func)
    def wrapper(*args, **kargs):
        dt_parsedKargs = signature(func).bind_partial(*args, **kargs).arguments
        if not "rEnv" in dt_parsedKargs:
            kargs["rEnv"] = None

        rEnv = kargs["rEnv"]
        if rEnv is None:
            clearEnv = True
            rEnv = ro.Environment()
        else:
            clearEnv = False
        kargs["rEnv"] = rEnv

        if not "rEnv" in signature(func).parameters:
            kargs.pop("rEnv")
        try:
            with ro.local_context(rEnv) as rlc:
                result = func(*args, **kargs)
        except rinterface_lib.embedded.RRuntimeError as e:
            ro.r.traceback()
            raise e
        if clearEnv:
            rEnv.clear()
        ro.r.gc()
        return result
    
    return wrapper

@rpy2_check
def r_set_seed(seed):
    """Set the seed of R random number generator"""
    from rpy2.robjects import r

    set_seed = r("set.seed")
    set_seed(seed)

@rpy2_check
@anndata2ri_check
def r2py(x, name=None, verbose=0):
    """Convert an rpy2 (R)  object to a Python object"""
    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri, pandas2ri
    from rpy2.robjects.conversion import localconverter
    import anndata2ri
    import time
    from tempfile import NamedTemporaryFile

    def _dataframe(objR):
        tpFile = NamedTemporaryFile(suffix=".feather")
        with ro.local_context() as rlc:
            rlc["objR"] = objR
            rlc["filePath"] = tpFile.name
            R(
                """
            library('arrow')
            objR$`index_r2py` <- rownames(objR)
            rownames(objR) <- NULL
            write_feather(objR, filePath)
            """
            )
        obj = pd.read_feather(tpFile.name)
        obj = obj.set_index("index_r2py").rename_axis(None)
        return obj

    if not name:
        name = ""

    try:
        objType = list(x.rclass)[0]
    except:
        objType = "unknown type"
    if verbose:
        print(f"transfer `{objType}` to python: {name} start", end="")
    timeStart = time.time()
    if ro.r("class")(x)[0] == "data.frame":
        x = _dataframe(x)
    else:
        try:
            with localconverter(
                ro.default_converter
                + numpy2ri.converter
                + pandas2ri.converter
                + anndata2ri.scipy2ri.converter
                + anndata2ri.converter
            ):
                x = ro.conversion.rpy2py(x)

        except TypeError:
            # workaround for: https://github.com/theislab/anndata2ri/issues/47
            x = anndata2ri.scipy2ri.rpy2py(x)
    timeEnd = time.time()
    timePass = timeEnd - timeStart
    if verbose:
        print(
            "\r"
            + f"transfer `{objType}` to python: {name} End. Elapsed time: {timePass:.0f}",
            flush=True,
        )
    return x

def py2r_disk(obj, check=False, *args, **kwargs):
    """Convert a Python object to R on disk"""
    from tempfile import NamedTemporaryFile
    import scanpy as sc

    def _adata(obj, X_layer="X"):
        zellkonverter = importr("zellkonverter")
        sce = importr("SingleCellExperiment")
        tpFile = NamedTemporaryFile(suffix=".h5ad")
        obj.var["temp_featureName"] = obj.var.index
        obj.obs["temp_barcodeName"] = obj.obs.index
        obj.write_h5ad(tpFile.name)
        objR = zellkonverter.readH5AD(tpFile.name, X_layer, reader="R")
        dfR_obs = py2r(obj.obs)
        dfR_var = py2r(obj.var)
        with ro.local_context() as rlc:
            rlc["objR"] = objR
            rlc["dfR_obs"] = dfR_obs
            rlc["dfR_var"] = dfR_var
            R(
                """
            objR@rowRanges@partitioning@NAMES <- rowData(objR)$temp_featureName
            objR@colData@rownames <- colData(objR)$temp_barcodeName
            objR@colData <- dfR_obs %>% DataFrame
            objR@rowRanges@elementMetadata <- dfR_var %>% DataFrame
            """
            )
            objR = R("objR")

        tpFile.close()
        return objR

    def _dataframe(obj):
        arrow = importr("arrow")
        rBase = importr("base")
        tpFile = NamedTemporaryFile(suffix=".feather")
        obj = obj.rename(columns=str)
        #  bypass error: `Object of type bool_ is not JSON serializable`
        for colName, colType in obj.dtypes.items():
            if isinstance(colType, pd.CategoricalDtype):
                ls_category = colType.categories
                obj[colName] = obj[colName].astype('object').astype('category').cat.set_categories(ls_category)
        
        if (obj.index == obj.reset_index().index).all():
            obj.to_feather(tpFile.name)
            needSetIndex = False
        else:
            obj.rename_axis("_index_py2r_").reset_index().to_feather(tpFile.name)
            needSetIndex = True

        dfR = arrow.read_feather(tpFile.name, as_data_frame=True)
        dfR = rBase.as_data_frame(dfR)
        if needSetIndex:
            with ro.local_context() as rlc:
                rlc["dfR"] = dfR
                R(
                    """
                srR_index <- dfR$`_index_py2r_`
                dfR$`_index_py2r_` <- NULL
                rownames(dfR) <- srR_index
                """
                )
                dfR = rlc["dfR"]
        return dfR

    def _array(obj):
        obj = pd.DataFrame(obj)
        obj = obj.rename(columns=str)
        dfR = py2r(obj)
        arR = rBase.as_matrix(dfR)
        return arR

    dt_config = {sc.AnnData: _adata, pd.DataFrame: _dataframe, np.ndarray: _array}
    if check:
        for _class in dt_config.keys():
            if isinstance(obj, _class):
                if _class == np.ndarray:  # _array only worked for 2D arrays
                    if len(obj.shape) == 2:
                        return True
                    else:
                        return False
                else:
                    return True
        else:
            return False
        # if type(obj) in dt_config:
        #     if type(obj) == np.asaarray:
        #         if len(obj.shape) == 2: # _array only worked for 2D arrays
        #             return True
        #         else:
        #             return False
        #     else:
        #         return True
        # else:
        #     return False
    for _class in dt_config.keys():
        if isinstance(obj, _class):
            _type = _class
            break
    func = dt_config[_type]
    objR = func(obj, *args, **kwargs)
    return objR

@rpy2_check
@anndata2ri_check
def py2r(x, name=None, on_disk=None, verbose=0):
    """Convert a Python object to an R object using rpy2"""
    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri, pandas2ri
    from rpy2.robjects.conversion import localconverter
    import anndata2ri
    import time

    if not name:
        name = ""
    objType = type(x)

    if on_disk == None:
        on_disk = True if py2r_disk(x, check=True) else False
    if verbose:
        print(
            f"on disk mode: {on_disk}, transfer `{objType}` to R: {name} start.", end=""
        )
    timeStart = time.time()

    if on_disk:
        x = py2r_disk(x)

    else:
        if sp.issparse(x):
            # workaround for: https://github.com/theislab/anndata2ri/issues/47
            x = anndata2ri.scipy2ri.py2rpy(x)

        with localconverter(
            ro.default_converter
            + numpy2ri.converter
            + pandas2ri.converter
            + anndata2ri.converter
        ):
            x = ro.conversion.py2rpy(x)

    timeEnd = time.time()
    timePass = timeEnd - timeStart
    if verbose:
        print(
            "\r"
            + f"on disk mode: {on_disk}, transfer `{objType}` to R: {name} End. Elapsed time: {timePass:.0f}",
            flush=True,
        )
    return x

@rcontext
def ad2so(
    ad,
    layer="raw",
    dataLayer=None,
    scaleLayer=None,
    scaleLayerInObsm=False,
    assay="RNA",
    rEnv=None,
    verbose=0,
    **kwargs,
):
    '''`ad2so` converts an AnnData object to a SeuratObject object
    
    Parameters
    ----------
    ad
        AnnData object
    layer, optional
        the layer to use for the count slot.
    dataLayer
        the name of the layer to use for the data slot.
    scaleLayer
        the name of the layer to use for scaleData slot.
    scaleLayerInObsm, optional
        if True, then the scaleLayer is assumed to be in the obsm of the adata object.
    assay, optional
        The assay to use.
    rEnv
        R environment to use. If None, then a new one is created.
    verbose, optional
        0, 1, 2, 3, 4
    '''
    import scipy.sparse as ss
    import rpy2
    import rpy2.robjects as ro
    import rpy2.robjects.packages
    from rpy2.robjects.packages import importr
    
    importr("Seurat")
    R = ro.r
    mt_count = ad.layers[layer]
    for x in ad.obs.index:
        assert '_' not in x, f'`_` is not allowed in obs.index. {x}'
    for x in ad.var.index:
        assert '_' not in x, f'`_` is not allowed in var.index. {x}'
    
    if ad.var.empty:
        ad.var["project_ad2so"] = "temp"
    if ad.obs.empty:
        ad.obs["project_ad2so"] = "temp"
    rEnv["mtR_count"] = py2r(mt_count.T)
    rEnv["arR_obsName"] = R.unlist(R.c(ad.obs.index.to_list()))
    rEnv["arR_varName"] = R.unlist(R.c(ad.var.index.to_list()))
    rEnv["assay"] = assay

    R(
        """
    colnames(mtR_count) <- arR_obsName
    rownames(mtR_count) <- arR_varName
    so <- CreateSeuratObject(mtR_count, assay=assay)
    """
    )
    if "highly_variable" in ad.var.columns:
        ls_hvgGene = (
            ad.var["highly_variable"]
            .loc[lambda x: x]
            .index.str.replace("_", "-")
            .to_list()
        )
        rEnv["arR_hvgGene"] = R.unlist(R.c(ls_hvgGene))
        R(
            """
        VariableFeatures(so) <- arR_hvgGene
        """
        )
    rEnv["dfR_obs"] = py2r(ad.obs, verbose=verbose)
    rEnv["dfR_var"] = py2r(ad.var, verbose=verbose)
    R(
        """
    so <- AddMetaData(so, dfR_obs)
    so[['RNA']] <- AddMetaData(so[[assay]], dfR_var)
    """
    )

    if dataLayer is None:
        R(
            "NormalizeData(so, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = F)"
        )
    else:
        mt_data = ad.layers[dataLayer]
        rEnv["mtR_data"] = py2r(mt_data.T, verbose=verbose)
        R(
            """
        colnames(mtR_data) <- arR_obsName
        rownames(mtR_data) <- arR_varName
        so <- SetAssayData(so, slot = "data",mtR_data, assay = assay)
        """
        )

    if scaleLayer is None:
        pass
    else:
        if not scaleLayerInObsm:
            mt_scaleData = ad.layers[scaleLayer]
            rEnv["mtR_scaleData"] = py2r(mt_scaleData.T, verbose=verbose)
            R(
                """
            colnames(mtR_scaleData) <- arR_obsName
            rownames(mtR_scaleData) <- arR_varName
            so <- SetAssayData(so, slot = "scale.data", mtR_scaleData, assay = assay)
            """
            )
        else:
            rEnv["dfR_scaleData"] = py2r(
                ad.obsm[scaleLayer].loc[:, lambda df: df.columns.isin(ad.var.index)].T,
                verbose=verbose,
            )
            R(
                """
            mtR_scaleData <- dfR_scaleData %>% as.matrix
            so <- SetAssayData(so, slot = "scale.data", mtR_scaleData, assay = assay)
            """
            )
    ls_obsm = [x for x in ad.obsm.keys() if x.startswith("X_")]
    for obsm in ls_obsm:
        obsm_ = obsm.split("X_", 1)[1]
        df_obsm = pd.DataFrame(
            ad.obsm[obsm],
            index=ad.obs.index,
            columns=[f"{obsm_}_{x}" for x in range(1, 1 + ad.obsm[obsm].shape[1])],
        )
        rEnv["dfR_obsm"] = py2r(df_obsm, verbose=verbose)
        rEnv["obsm"] = obsm_
        R(
            """
        mtR_obsm <- dfR_obsm %>% as.matrix
        so[[obsm]] <- CreateDimReducObject(mtR_obsm, assay=assay, key=paste0(obsm, '_'))
        """
        )

    for obsp in ad.obsp.keys():
        rEnv["mtR_obsp"] = py2r(ss.csc_matrix(ad.obsp[obsp]), verbose=verbose)
        rEnv["obsp"] = obsp
        R(
            """
        colnames(mtR_obsp) <- arR_obsName
        rownames(mtR_obsp) <- arR_obsName
        so[[obsp]] <- as.Graph(x = mtR_obsp)
        """
        )
    return rEnv["so"]
