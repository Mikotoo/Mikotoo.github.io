import basic
import anndata
from rTools import r_set_seed
from rTools import ad2so
from rTools import r2py
import rpy2
import rpy2.robjects as ro
import rpy2.robjects.packages
from rpy2 import rinterface_lib
from rpy2.robjects.packages import importr
import loguru
from loguru import logger
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

def run_ScDblFinder(
    adata: anndata.AnnData,
    layer: str = "X",
    copy: bool = False,
    batch_key: Optional[str] = None,
    doubletRatio: Optional[float] = None,
    skipCheck: bool = False,
    dropDoublet: bool = True,
    BPPARAM=None
) -> Optional[anndata.AnnData]:
    """
    use ScDblFinder detect doublets.

    Parameters
    ----------
    adata : anndata.AnnData
        anndata
    layer : str, optional
        use this layer. must be raw counts. Defaults to X
    copy : bool, optional
        copy adata or not. Defaults to False.
    doubletRatio : float, optional
        expected doublet ratio. Defaults to 0.1

    Returns
    -------
    Optional[anndata.AnnData]
        anndata if copy
    """

    r_set_seed(39)
    R = ro.r
    Seurat = importr('Seurat')
    rBase = importr("base")
    BiocParallel = importr("BiocParallel")
    scDblFinder = importr("scDblFinder")

    ls_obsInfo = []
    if not batch_key:
        batch_key = R("NULL")
    else:
        ls_obsInfo.append(batch_key)
    if not doubletRatio:
        doubletRatio = R("NULL")
    if not BPPARAM:
        BPPARAM = BiocParallel.SerialParam()
    else:
        BPPARAM = BiocParallel.MulticoreParam(BPPARAM)

    if not skipCheck:
        basic.testAllCountIsInt(adata, layer)

    tempAd = basic.getPartialLayersAdata(adata, layer, obsInfoLs=ls_obsInfo)
    tempAd.layers["counts"] = tempAd.X

    logger.info("start to transfer adata to R")
    so = ad2so(tempAd, layer='counts')
    tempAdr = Seurat.as_SingleCellExperiment(so)
    del tempAd
    del so

    logger.info("start to calculate doublet score")

    tempAdr = scDblFinder.scDblFinder(tempAdr, samples=batch_key, dbr=doubletRatio, BPPARAM=BPPARAM)

    logger.info("start to intergrate result with adata")
    # scDblFinderResultDf = r2py(R.as_data_frame(tempAdr.slots["colData"]))
    scDblFinderResultDfR = rBase.as_data_frame(tempAdr.slots["colData"], row_names = R('row.names')(tempAdr.slots["colData"]), stringsAsFactors=False)
    scDblFinderResultDf = r2py(scDblFinderResultDfR)
    scDblFinderResultDf["scDblFinder.class"] = list(rBase.as_character(scDblFinderResultDfR.rx2['scDblFinder.class']))
    # scDblFinderResultDf = r2py(rBase.as_data_frame(tempAdr.slots["colData"], row_names = R('row.names')(tempAdr.slots["colData"]), stringsAsFactors=False))
    # import pdb; pdb.set_trace()

    adata.obsm["scDblFinder"] = (
        scDblFinderResultDf.reindex(adata.obs.index)
        .filter(regex=r"^scDblFinder[\w\W]*")
        .copy(deep=True)
    )
    adata.obsm["scDblFinder"].columns = adata.obsm["scDblFinder"].columns.astype(
        str
    )

    if dropDoublet:
        logger.info(f"before filter: {len(adata)}")
        adata._inplace_subset_obs(
            adata.obsm["scDblFinder"]["scDblFinder.class"] == "singlet"
        )
        logger.info(f"after filter: {len(adata)}")

    if copy:
        return adata
