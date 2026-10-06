#!/usr/bin/env node
/**
 * 从同级仓库聚合分析步骤文档到 Starlight 内容目录。
 *
 * 数据来源（默认 ../Lotus_genome、../Rscript，可用环境变量覆盖）：
 *   - Lotus_genome/README.md、scripts/**、demo/**、third_party/**、environment/
 *   - Rscript/README.md、Singlecell-RNAseq analysis/Script/README
 *
 * 产物写入 src/content/docs/lotus/ 与 src/content/docs/rscript/，这两处已被 .gitignore 忽略，
 * 属于「每次构建重新生成」的内容，仓库里只保留脚本本身。
 *
 * 约定：
 *   - README 是唯一事实来源，脚本不修改上游仓库。
 *   - 相对链接改写为 GitHub 链接（原文件所在位置），站内链接改写为站内路由。
 *   - 图片一律指向 raw.githubusercontent.com，避免大文件进入构建产物。
 */
import { readFile, writeFile, mkdir, rm, readdir, stat } from 'node:fs/promises';
import { existsSync, statSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const SITE_ROOT = path.resolve(HERE, '..');
const DOCS_ROOT = path.join(SITE_ROOT, 'src', 'content', 'docs');

/**
 * 定位上游仓库目录。本地开发时它们与本站同级；CI 或子模块场景下也可能放在本站内部。
 * 依次尝试候选位置，返回第一个存在的；都不存在则返回首选路径（随后会被跳过并给出提示）。
 */
function resolveRepoDir(envVar, siteRelCandidates, siblingName) {
  const explicit = process.env[envVar];
  if (explicit) return path.resolve(explicit);
  const candidates = [
    ...siteRelCandidates.map((c) => path.resolve(SITE_ROOT, c)),
    path.resolve(SITE_ROOT, '..', siblingName),
  ];
  return candidates.find((c) => existsSync(c)) ?? candidates[candidates.length - 1];
}

const REPOS = {
  lotus: {
    dir: resolveRepoDir('LOTUS_GENOME_DIR', ['Lotus_genome', 'vendor/Lotus_genome'], 'Lotus_genome'),
    owner: 'Mikotoo',
    name: 'Lotus_genome',
    branch: 'main',
  },
  rscript: {
    dir: resolveRepoDir('RSCRIPT_DIR', ['Rscript', 'vendor/Rscript'], 'Rscript'),
    owner: 'Mikotoo',
    name: 'Rscript',
    branch: 'main',
  },
};

const stats = { pages: 0, skipped: [], linked: 0, missingDirs: [] };

/* ------------------------------------------------------------------ utils */

const posix = (p) => p.split(path.sep).join('/');

function slugifySegment(seg) {
  return seg
    .toLowerCase()
    .replace(/[\s_]+/g, '-')
    .replace(/[^\p{L}\p{N}-]+/gu, '-')
    .replace(/^-+|-+$/g, '')
    .replace(/-{2,}/g, '-');
}

/**
 * 导航用的中文标签。上游目录名是英文/拼音缩写，直接展示可读性差，
 * 这里按「仓库相对目录 → 中文标签」做一次映射；未列出的目录自动生成标签。
 * 新增步骤后可以不改这里，脚本会退回 humanize() 的结果。
 */
const LABELS = {
  'Lotus_genome': {
    demo: '可运行示例',
    environment: '软件版本记录',
    scripts: '分析流程总览',
    third_party: '第三方代码',
    'demo/03_telomere': '示例：端粒识别',
    'demo/05_centromere': '示例：着丝粒串联重复',
    'scripts/01_assembly': '01 组装',
    'scripts/02_annotation': '02 注释',
    'scripts/03_telomere': '03 端粒',
    'scripts/04_rDNA': '04 rDNA',
    'scripts/05_centromere': '05 着丝粒',
    'scripts/06_genome_comparison': '06 比较基因组',
    'scripts/07_transcriptome': '07 转录组',
    'scripts/08_single_nucleus': '08 单核转录组',
    'scripts/01_assembly/01_hifiasm': 'hifiasm 组装',
    'scripts/01_assembly/02_purge_haplotigs': '去冗余单倍型',
    'scripts/01_assembly/03_remove_organelle': '去除细胞器序列',
    'scripts/01_assembly/04_telomere_check': '端粒检查',
    'scripts/01_assembly/05_contig_ordering': 'contig 排序',
    'scripts/01_assembly/06_gapclose': '补洞',
    'scripts/01_assembly/07_polish': '打磨校正',
    'scripts/01_assembly/08_assembly_qc': '组装质量评估',
    'scripts/01_assembly/09_hic_scaffold': 'Hi-C 挂载',
    'scripts/01_assembly/09_hic_scaffold/haphic': 'HapHiC 分型挂载',
    'scripts/02_annotation/01_repeat_edta': '重复序列屏蔽（EDTA）',
    'scripts/02_annotation/02_rnaseq_align': 'RNA-seq 比对',
    'scripts/02_annotation/03_isoseq': 'IsoSeq 全长转录本',
    'scripts/02_annotation/04_braker': '基因预测（BRAKER）',
    'scripts/02_annotation/05_geta': '基因模型整合（GETA）',
    'scripts/02_annotation/06_gene_model_filter': '基因模型过滤',
    'scripts/02_annotation/07_gene_classification': '基因分类',
    'scripts/02_annotation/08_functional_annotation': '功能注释',
    'scripts/02_annotation/09_annotation_export': '注释导出',
    'scripts/02_annotation/10_gene_correspondence': '基因对应关系',
    'scripts/04_rDNA/01_locus_and_arrays': 'rDNA 位点与阵列',
    'scripts/04_rDNA/03_copy_number': 'rDNA 拷贝数',
    'scripts/05_centromere/01_supercluster': '超级簇合并',
    'scripts/05_centromere/02_repeat_units': '重复单元',
    'scripts/05_centromere/03_consensus': '一致性序列',
    'scripts/05_centromere/04_stainglass': 'StainedGlass 分析',
    'scripts/06_genome_comparison/01_published_assemblies': '已发表组装',
    'scripts/06_genome_comparison/01_published_assemblies/gap_analysis': 'gap 区域分析',
    'scripts/06_genome_comparison/02_accession_alignment': '品种间比对',
    'scripts/06_genome_comparison/03_translocation_check': '易位核查',
    'scripts/06_genome_comparison/05_t2t_vs_previous': 'T2T 与前版本比较',
    'scripts/06_genome_comparison/06_annotation_comparison': '注释比较',
    'scripts/07_transcriptome/00_prepare': '数据准备',
    'scripts/07_transcriptome/01_tissue_atlas': '组织表达图谱',
    'scripts/07_transcriptome/02_mutant_and_public_data': '突变体与公共数据',
    'scripts/07_transcriptome/05_wgcna': 'WGCNA 共表达模块',
    'scripts/07_transcriptome/06_mfuzz': 'Mfuzz 时序聚类',
    'scripts/07_transcriptome/07_coexpression_network': '共表达网络',
    'scripts/07_transcriptome/08_module_consistency': '模块一致性',
    'scripts/08_single_nucleus/Script': '单核分析脚本',
  },
  Rscript: {
    'Bulk-RNAseq analysis': 'Bulk RNA-seq 分析',
    'Homologous genes': '同源基因',
    'New genes finder': '新基因识别',
    'Singlecell-RNAseq analysis': '单核 RNA-seq 分析',
    'Singlecell-RNAseq analysis/Script': '单核分析脚本明细',
    'Structural variation': '结构变异',
  },
};

/** 目录名 → 展示标题：优先取中文标签映射，否则去掉数字前缀后用目录名。 */
function humanize(name) {
  let s = name.replace(/^\d+[._\-\s]*/, '').trim();
  if (!s) s = name;
  const parts = s.split(' / ');
  let title = parts[0].trim();
  title = title.replace(/_/g, ' ');
  return title.charAt(0).toUpperCase() + title.slice(1);
}

/** 按仓库相对目录取导航标签。 */
function labelFor(repoName, dirRel, fallbackName) {
  const map = LABELS[repoName] ?? {};
  if (map[dirRel]) return map[dirRel];
  return humanize(fallbackName);
}

/** 从目录名的数字前缀取侧边栏排序值；没有前缀则排在后面。 */
function orderOf(name) {
  const m = name.match(/^(\d+)/);
  return m ? Number(m[1]) : 999;
}

function yamlString(s) {
  return `'${String(s).replace(/'/g, "''")}'`;
}

function rawUrl(repo, relPath) {
  return `https://raw.githubusercontent.com/${repo.owner}/${repo.name}/${repo.branch}/${posix(relPath)}`;
}

function blobUrl(repo, relPath, anchor) {
  const base = `https://github.com/${repo.owner}/${repo.name}/blob/${repo.branch}/${posix(relPath)}`;
  return anchor ? `${base}#${anchor}` : base;
}

/** 去掉行内的 Markdown 装饰，得到纯文本。 */
function stripInline(text) {
  return text
    .replace(/`([^`]*)`/g, '$1')
    .replace(/\*\*([^*]*)\*\*/g, '$1')
    .replace(/\*([^*]*)\*/g, '$1')
    .replace(/\[([^\]]*)\]\([^)]*\)/g, '$1')
    .replace(/\s+/g, ' ')
    .trim();
}

/** 收集仓库内所有以 README.md（或 README）结尾的文件，相对仓库根。 */
async function collectReadmes(repoDir) {
  const out = [];
  async function walk(dir, rel) {
    let entries;
    try {
      entries = await readdir(dir, { withFileTypes: true });
    } catch {
      return;
    }
    for (const e of entries) {
      if (e.name === '.git' || e.name === 'node_modules') continue;
      const abs = path.join(dir, e.name);
      const relPath = rel ? `${rel}/${e.name}` : e.name;
      if (e.isDirectory()) {
        await walk(abs, relPath);
      } else if (/^README(\.md)?$/i.test(e.name)) {
        out.push(relPath);
      }
    }
  }
  await walk(repoDir, '');
  return out.sort();
}

/* -------------------------------------------------------------- markdown */

/** 相对路径解析：返回 { kind, path } */
function resolveTarget(currentRelPath, target) {
  const hashIdx = target.indexOf('#');
  const anchor = hashIdx >= 0 ? target.slice(hashIdx + 1) : '';
  let clean = (hashIdx >= 0 ? target.slice(0, hashIdx) : target).trim();
  if (!clean) return { kind: 'anchor', anchor };
  if (/^(https?:|mailto:|#)/i.test(clean)) return { kind: 'absolute', url: target };
  if (clean.startsWith('/')) return { kind: 'absolute', url: target };
  const baseDir = path.posix.dirname(posix(currentRelPath));
  const resolved = path.posix.normalize(path.posix.join(baseDir, clean));
  if (resolved.startsWith('..')) return { kind: 'outside', rel: resolved, anchor };
  const trailingSlash = clean.endsWith('/');
  return { kind: 'relative', rel: resolved, anchor, trailingSlash };
}

function makeLinkRewriter(repo, currentRelPath, siteIndex) {
  return function rewrite(raw, isImage) {
    const r = resolveTarget(currentRelPath, raw);
    if (r.kind === 'absolute') return raw;
    if (r.kind === 'anchor') return raw;
    if (r.kind === 'outside') {
      // 指向上游仓库之外（例如本站），退回 GitHub 链接不可靠，保留原文
      stats.missingDirs.push(`${currentRelPath} -> ${raw} (越出仓库)`);
      return raw;
    }
    const upstreamDir = path.posix.dirname(posix(currentRelPath));
    const abs = path.join(repo.dir, r.rel.split('/').join(path.sep));

    if (isImage) return rawUrl(repo, r.rel);

    // 站内已有对应页面 → 转成站内路由，并尽量保留锚点
    const siteSlug = siteIndex.get(r.rel.replace(/\/$/, ''));
    if (siteSlug && !isImage) {
      const anchor = r.anchor ? `#${slugifySegment(r.anchor)}` : '';
      return `/${siteSlug}/${anchor}`;
    }
    // 目录没有独立页面（例如指向某个步骤目录）→ 若能落到其 README，则用其页面
    const dirSlug = siteIndex.get(`${r.rel.replace(/\/$/, '')}/README.md`);
    if (dirSlug) {
      const anchor = r.anchor ? `#${slugifySegment(r.anchor)}` : '';
      return `/${dirSlug}/${anchor}`;
    }
    if (existsSync(abs)) {
      const st = statSyncSafe(abs);
      if (st && st.isDirectory()) {
        const idx = siteIndex.get(`${r.rel.replace(/\/$/, '')}/README.md`);
        if (idx) return `/${idx}/`;
      }
    }
    return blobUrl(repo, r.rel, r.anchor);
  };
}

function statSyncSafe(p) {
  try {
    return existsSync(p) ? statSync(p) : null;
  } catch {
    return null;
  }
}

/** 把英文段落折叠进 <details>，正文默认展示中文。 */
function collapseEnglish(body) {
  const re = /^(##\s+English\s*)$/m;
  const m = body.match(re);
  if (!m) return body;
  const idx = m.index;
  const head = body.slice(0, idx);
  const tail = body.slice(idx + m[0].length);
  const collapsed = `<details>\n<summary>English</summary>\n${tail.trim()}\n\n</details>\n`;
  return `${head.trimEnd()}\n\n${collapsed}`;
}

/**
 * 转换一个 README：去 h1、改写链接、可选折叠英文、补 frontmatter。
 * sidebarLabel / sidebarOrder 会写进 frontmatter，让自动生成的侧边栏显示中文友好名，
 * 并按目录编号排序（而不是按 slug 字母序）。
 */
function convertReadme({ markdown, repo, relPath, title, siteIndex, collapse, sidebarLabel, sidebarOrder }) {
  let body = markdown.replace(/^\uFEFF/, '');
  let docTitle = title;

  // 取第一个 h1 作为标题并移除
  const h1 = body.match(/^#\s+(.+)$/m);
  if (h1) {
    if (!docTitle) docTitle = stripInline(h1[1]).split(' / ')[0].trim();
    body = body.replace(h1[0], '').trimStart();
  }
  if (collapse) body = collapseEnglish(body);

  const rewrite = makeLinkRewriter(repo, relPath, siteIndex);

  // 图片 + 链接（避免二次处理已生成的结果，用占位符保护）
  body = body.replace(/!\[([^\]]*)\]\(([^)\s]+)(\s+"[^"]*")?\)/g, (full, alt, url, t) => {
    const next = rewrite(url, true);
    if (next !== url) stats.linked++;
    return `![${alt}](${next}${t ?? ''})`;
  });
  body = body.replace(/(?<!!)\[([^\]]*)\]\(([^)\s]+)(\s+"[^"]*")?\)/g, (full, text, url, t) => {
    const next = rewrite(url, false);
    if (next !== url) stats.linked++;
    return `[${text}](${next}${t ?? ''})`;
  });

  // MDX 安全：占位符里的花括号、HTML 注释
  body = body.replace(/<!--[\s\S]*?-->/g, '');
  body = body.replace(/\{(PROJ_[A-Z_]+|CONDA_[A-Z0-9_]+|SOFTWARE_[A-Z_]+|HOME_[A-Z_]+|TOKEN)\}/g, '`{$1}`');

  return `${frontmatter({ title: docTitle, sidebarLabel, sidebarOrder })}\n${body.trim()}\n`;
}

/** 生成页面 frontmatter。 */
function frontmatter({ title, sidebarLabel, sidebarOrder }) {
  const lines = ['---', `title: ${yamlString(title)}`];
  if (sidebarLabel) lines.push(`sidebar:\n  label: ${yamlString(sidebarLabel)}`);
  if (Number.isFinite(sidebarOrder)) {
    if (sidebarLabel) lines[lines.length - 1] += `\n  order: ${sidebarOrder}`;
    else lines.push(`sidebar:\n  order: ${sidebarOrder}`);
  }
  lines.push('---', '');
  return `${lines.join('\n')}\n`;
}

/* ------------------------------------------------------------ page build */

/** 预扫描：决定哪些目录有页面，并建立「仓库相对路径 → 站内 slug」索引。 */
function buildPlan(repo, readmes) {
  const pages = new Map(); // repoRelDir（'' 表示根） -> { slug, readmeRel, synthetic }
  const siteIndex = new Map(); // 仓库相对路径（文件或目录） -> slug

  for (const rel of readmes) {
    const dir = path.posix.dirname(rel) === '.' ? '' : path.posix.dirname(rel);
    if (pages.has(dir)) continue;
    const slug = pageSlugFor(repo, dir);
    pages.set(dir, { slug, readmeRel: rel });
    siteIndex.set(rel, slug);
    siteIndex.set(dir, slug);
  }

  // 上游有些中间目录本身没有 README（例如 scripts/01_assembly/ 的 README 在各步骤里）。
  // 为这些目录补一个索引页，侧边栏才会形成「阶段 → 步骤」的层级，而不是摊平。
  const dirsWithReadme = new Set(pages.keys());
  const allDirs = new Set();
  for (const dir of dirsWithReadme) {
    if (dir === '') continue;
    const parts = dir.split('/');
    for (let i = 1; i <= parts.length; i++) allDirs.add(parts.slice(0, i).join('/'));
  }
  for (const dir of [...allDirs].sort()) {
    if (dirsWithReadme.has(dir)) continue;
    const slug = pageSlugFor(repo, dir);
    pages.set(dir, { slug, readmeRel: null, synthetic: true });
    siteIndex.set(dir, slug);
  }

  return { pages, siteIndex };
}

function pageSlugFor(repo, repoRelDir) {
  if (repo === REPOS.lotus) {
    if (repoRelDir === '') return 'lotus';
    return `lotus/${repoRelDir.split('/').map(slugifySegment).join('/')}`;
  }
  if (repoRelDir === '') return 'rscript';
  return `rscript/${repoRelDir.split('/').map(slugifySegment).join('/')}`;
}

async function writePage(slug, content) {
  const file = path.join(DOCS_ROOT, `${slug}.md`);
  await mkdir(path.dirname(file), { recursive: true });
  await writeFile(file, content, 'utf8');
  stats.pages++;
}

/** 生成阶段页（目录 README + 子步骤索引表）。 */
function buildStageIndex({ repo, dirRel, childEntries }) {
  const rows = childEntries
    .map(({ name, slug, desc }) => {
      const label = labelFor(repo.name, `${dirRel}/${name}`, name);
      const link = slug ? `[${label}](/${slug}/)` : label;
      return `| ${link} | ${desc || '—'} |`;
    })
    .join('\n');
  const table = childEntries.length
    ? `\n### 本阶段步骤\n\n| 步骤 | 内容 |\n| --- | --- |\n${rows}\n`
    : '';
  return table;
}

/** 从 README 的第一个表格里取「步骤 → 内容」摘要，用于索引表。 */
function extractSummaryTable(markdown) {
  const map = new Map();
  const lines = markdown.split('\n');
  for (const line of lines) {
    const m = line.match(/^\|\s*\[?`?([^`\]|]+?)\/`?\]?\(?[^|]*\)?\s*\|\s*([^|]*)\|/);
    if (m) {
      const key = m[1].trim();
      const desc = stripInline(m[2]);
      if (key && desc && !/^-+$/.test(key) && key !== '阶段' && key !== 'Stage') map.set(key, desc);
    }
  }
  return map;
}

/* ------------------------------------------------------------------ main */

async function syncRepo(key, repo) {
  if (!existsSync(repo.dir)) {
    console.warn(`[sync] 跳过 ${key}：目录不存在 ${repo.dir}`);
    return;
  }
  const readmes = await collectReadmes(repo.dir);
  const plan = buildPlan(repo, readmes);
  const outDir = path.join(DOCS_ROOT, key);
  if (existsSync(outDir)) await rm(outDir, { recursive: true, force: true });

  let summaryTables = new Map();

  // 先读根 README 的索引表，用于给阶段/步骤补「内容」摘要
  const rootPage = plan.pages.get('');
  if (rootPage?.readmeRel) {
    const rootMd = await readFile(path.join(repo.dir, rootPage.readmeRel.split('/').join(path.sep)), 'utf8');
    summaryTables = extractSummaryTable(rootMd);
  }

  for (const [dirRel, page] of plan.pages) {
    const { slug, readmeRel, synthetic } = page;

    // 子目录：本目录下直接子目录中有页面的
    const children = [...plan.pages.keys()]
      .filter((d) => d !== dirRel && path.posix.dirname(d) === (dirRel || '.') && d !== '')
      .sort()
      .map((d) => {
        const name = path.posix.basename(d);
        const desc = summaryTables.get(name) ?? summaryTables.get(`${name}/`) ?? '';
        return { name, slug: plan.pages.get(d).slug, desc };
      });

    const isRoot = dirRel === '';
    const baseName = isRoot ? '' : path.posix.basename(dirRel);
    const label = isRoot
      ? key === 'lotus'
        ? 'Lotus T2T 分析流程'
        : '论文分析脚本'
      : labelFor(repo.name, dirRel, baseName);
    // 目录自己的索引页要排在它的子页前面，否则会按 slug 字母序掉到中间
    const order = isRoot ? 0 : children.length ? -1 : orderOf(baseName);
    const title = label;

    let content;
    if (synthetic || !readmeRel) {
      // 上游没有这个目录的 README：生成一个纯索引页
      const rows = children.map(
        ({ name, slug: s, desc }) => `| [${labelFor(repo.name, `${dirRel}/${name}`, name)}](/${s}/) | ${desc || '—'} |`
      );
      const table = rows.length
        ? `| 步骤 | 内容 |\n| --- | --- |\n${rows.join('\n')}`
        : '_该目录下暂无可展示的子页面。_';
      content = `${frontmatter({ title, sidebarLabel: label, sidebarOrder: order })}\n本阶段包含以下步骤，内容取自各步骤目录的 README。\n\n${table}\n`;
    } else {
      const markdown = await readFile(path.join(repo.dir, readmeRel.split('/').join(path.sep)), 'utf8');
      content = convertReadme({
        markdown,
        repo,
        relPath: readmeRel,
        title,
        siteIndex: plan.siteIndex,
        collapse: !isRoot,
        sidebarLabel: label,
        sidebarOrder: order,
      });
      if (children.length) {
        const extra = buildStageIndex({ repo, dirRel, childEntries: children });
        content = `${content.trimEnd()}\n${extra}\n`;
      }
    }
    await writePage(slug, `${content}\n`);
  }

  console.log(
    `[sync] ${key}: ${plan.pages.size} 页（${readmes.length} 个 README，补 ${plan.pages.size - readmes.length} 个目录索引）`
  );
  return plan;
}

async function main() {
  await mkdir(DOCS_ROOT, { recursive: true });
  for (const [key, repo] of Object.entries(REPOS)) {
    await syncRepo(key, repo);
  }
  if (stats.missingDirs.length) {
    console.log(`[sync] 提示：${stats.missingDirs.length} 个链接越出仓库被保留，示例：`);
    stats.missingDirs.slice(0, 5).forEach((s) => console.log(`   - ${s}`));
  }
  console.log(`[sync] 完成：生成 ${stats.pages} 个页面，改写 ${stats.linked} 个链接`);
}

await main();
