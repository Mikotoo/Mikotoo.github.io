#!/usr/bin/env node
/**
 * 项目文档同步：扫描「带 content.config.ts 的仓库」，把它们的 README 聚合成站点文档区。
 *
 * 为什么要这样设计：
 *   本站在持续成长，以后可能新增任意主题的仓库/工具项目，不应该每次都在站点里手写一份介绍。
 *   约定很简单——**只要一个目录里有 content.config.ts，它就是一个可被收录的项目**。
 *   出处分两类：
 *     A. 同级目录（默认）：../Lotus_genome、../Rscript —— 你本地的代码仓库
 *     B. 本站内 vendor/ 下的目录（可选）：适合把外部小工具作为 git submodule 收进来
 *
 * 每个项目在站点里的结构：
 *   /_projects/<key>/                    项目首页（来自仓库根 README）
 *   /_projects/<key>/<子目录>/            子目录 README 自动成页（保留层级）
 * 站点对外入口是 /projects/（见 src/content/docs/projects.md）。
 */
import { readFile, writeFile, mkdir, rm, readdir, stat, lstat, cp } from 'node:fs/promises';
import { existsSync, readFileSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const SITE_ROOT = path.resolve(HERE, '..');
const DOCS_ROOT = path.join(SITE_ROOT, 'src', 'content', 'docs');
/** 项目文档统一放在带下划线前缀的目录下，避免与手写内容混在一起 */
const PROJECTS_DIR = '_projects';

/**
 * 导航用的中文标签（按「项目 key → 仓库内相对目录」映射）。
 * 未列出的目录会自动用目录名生成标签，因此新增项目不改这里也能跑。
 */
const LABELS = {
  lotus: {
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
  rscript: {
    'Bulk-RNAseq analysis': 'Bulk RNA-seq 分析',
    'Homologous genes': '同源基因',
    'New genes finder': '新基因识别',
    'Singlecell-RNAseq analysis': '单核 RNA-seq 分析',
    'Singlecell-RNAseq analysis/Script': '单核分析脚本明细',
    'Structural variation': '结构变异',
  },
};

const stats = { pages: 0, projects: [], skipped: [], linked: 0, notes: [] };

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

function humanize(name) {
  let s = name.replace(/^\d+[._\-\s]*/, '').trim();
  if (!s) s = name;
  let title = s.split(' / ')[0].trim().replace(/_/g, ' ');
  return title.charAt(0).toUpperCase() + title.slice(1);
}

function labelFor(projectKey, dirRel, fallbackName) {
  const map = LABELS[projectKey] ?? {};
  return map[dirRel] ?? humanize(fallbackName);
}

function orderOf(name) {
  const m = name.match(/^(\d+)/);
  return m ? Number(m[1]) : 999;
}

function yamlString(s) {
  return `'${String(s).replace(/'/g, "''")}'`;
}

function frontmatter({ title, sidebarLabel, sidebarOrder, description }) {
  const lines = ['---', `title: ${yamlString(title)}`];
  if (description) lines.push(`description: ${yamlString(description)}`);
  const sb = [];
  if (sidebarLabel) sb.push(`  label: ${yamlString(sidebarLabel)}`);
  if (Number.isFinite(sidebarOrder)) sb.push(`  order: ${sidebarOrder}`);
  if (sb.length) lines.push('sidebar:', ...sb);
  lines.push('---', '');
  return `${lines.join('\n')}\n`;
}

function stripInline(text) {
  return text
    .replace(/`([^`]*)`/g, '$1')
    .replace(/\*\*([^*]*)\*\*/g, '$1')
    .replace(/\*([^*]*)\*/g, '$1')
    .replace(/\[([^\]]*)\]\([^)]*\)/g, '$1')
    .replace(/\s+/g, ' ')
    .trim();
}

function rawUrl(repo, relPath) {
  return `https://raw.githubusercontent.com/${repo.owner}/${repo.name}/${repo.branch}/${posix(relPath)}`;
}

function blobUrl(repo, relPath, anchor) {
  const base = `https://github.com/${repo.owner}/${repo.name}/blob/${repo.branch}/${posix(relPath)}`;
  return anchor ? `${base}#${anchor}` : base;
}

function statSyncSafe(p) {
  try {
    return existsSync(p) ? statSync(p) : null;
  } catch {
    return null;
  }
}

/* ------------------------------------------------------- project discovery */

/**
 * 找一个仓库的 remote 信息，用于把相对链接改写回 GitHub。
 * 读 .git/config，失败则回退到目录名 + 默认 owner。
 */
async function readRepoRemote(dir) {
  const fallback = { owner: 'Mikotoo', name: path.basename(dir), branch: 'main' };
  const cfg = path.join(dir, '.git', 'config');
  if (!existsSync(cfg)) return fallback;
  try {
    const text = await readFile(cfg, 'utf8');
    const m = text.match(/\[remote "origin"\][\s\S]*?url\s*=\s*(.+)/);
    if (!m) return fallback;
    const url = m[1].trim();
    const gh = url.match(/github\.com[/:]([^/]+)\/([^/\s]+?)(?:\.git)?$/);
    if (!gh) return fallback;
    return { owner: gh[1], name: gh[2], branch: fallback.branch };
  } catch {
    return fallback;
  }
}

/** 读仓库 content.config.ts 里声明的项目名（可选）。 */
async function readProjectName(dir) {
  const f = path.join(dir, 'content.config.ts');
  if (!existsSync(f)) return null;
  try {
    const text = await readFile(f, 'utf8');
    const m = text.match(/name:\s*['"]([^'"]+)['"]/);
    return m ? m[1] : null;
  } catch {
    return null;
  }
}

/**
 * 发现可收录的项目，两个来源合并：
 *
 *   1. 显式登记：src/data/projects.json 里的 `key` + `dir`（或同名的同级目录）。
 *      这是正常工作方式——你的仓库不动，只在站点里登记一条。
 *   2. 自动发现：同级目录或 vendor/ 下，凡是含 content.config.ts 的目录。
 *      给「以后新加的仓库、外部小工具」用，不用先登记也能出现在文档区。
 *
 * 两种方式都找不到的仓库会被静默跳过（例如 CI 上没有同级仓库），不影响其余内容构建。
 */
async function discoverProjects() {
  const found = new Map();

  const register = async (dir, keyHint, displayHint) => {
    if (!existsSync(dir) || !existsSync(path.join(dir, 'README.md'))) return null;
    if (path.resolve(dir) === SITE_ROOT) return null;
    const base = path.basename(dir);
    const key = (keyHint ?? base).toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '');
    if (!key || found.has(key)) return null;
    const project = {
      key,
      dir,
      repo: await readRepoRemote(dir),
      displayName: displayHint ?? (await readProjectName(dir)) ?? base,
      via: keyHint ? 'projects.json' : 'content.config.ts',
    };
    found.set(key, project);
    return project;
  };

  // 1) 显式登记
  const manifestPath = path.join(SITE_ROOT, 'src', 'data', 'projects.json');
  const manifest = await readJsonFile(manifestPath, { projects: [] });
  for (const entry of manifest.projects ?? []) {
    if (!entry?.key) continue;
    const names = [entry.dir, entry.repo, entry.key].filter(Boolean);
    const bases = [path.resolve(SITE_ROOT, '..'), path.join(SITE_ROOT, 'vendor')];
    for (const base of bases) {
      for (const name of names) {
        if (await register(path.join(base, name), entry.key, entry.title)) break;
      }
      if (found.has(entry.key)) break;
    }
    if (!found.has(entry.key)) {
      stats.skipped.push(`${entry.key}（未找到同级目录：${names.join(' / ')}）`);
    }
  }

  // 2) 自动发现带标记文件的仓库
  const scan = async (baseDir) => {
    if (!existsSync(baseDir)) return;
    let entries;
    try {
      entries = await readdir(baseDir, { withFileTypes: true });
    } catch {
      return;
    }
    for (const e of entries) {
      if (!e.isDirectory()) continue;
      const dir = path.join(baseDir, e.name);
      const looksLikeProject =
        existsSync(path.join(dir, 'content.config.ts')) &&
        existsSync(path.join(dir, 'README.md'));
      if (looksLikeProject) await register(dir);
    }
  };
  await scan(path.resolve(SITE_ROOT, '..'));
  await scan(path.join(SITE_ROOT, 'vendor'));

  return [...found.values()].sort((a, b) => a.key.localeCompare(b.key));
}

async function readJsonFile(file, fallback) {
  if (!existsSync(file)) return fallback;
  try {
    return JSON.parse(await readFile(file, 'utf8'));
  } catch {
    return fallback;
  }
}

/** 从快照里的项目首页取 sidebar.label，作为沿用快照时的导航标题 */
function projectLabelFromSnapshot(rootPageFile) {
  try {
    if (!existsSync(rootPageFile)) return null;
    const text = statSync(rootPageFile) && readFileSync(rootPageFile, 'utf8');
    return text.match(/^\s*label:\s*'?([^'\n]+)'?\s*$/m)?.[1]?.trim() ?? null;
  } catch {
    return null;
  }
}

/* -------------------------------------------------------------- markdown */

function resolveTarget(currentRelPath, target) {
  const hashIdx = target.indexOf('#');
  const anchor = hashIdx >= 0 ? target.slice(hashIdx + 1) : '';
  const clean = (hashIdx >= 0 ? target.slice(0, hashIdx) : target).trim();
  if (!clean) return { kind: 'anchor', anchor };
  if (/^(https?:|mailto:|#)/i.test(clean) || clean.startsWith('/')) {
    return { kind: 'absolute', url: target };
  }
  const baseDir = path.posix.dirname(posix(currentRelPath));
  const resolved = path.posix.normalize(path.posix.join(baseDir, clean));
  if (resolved.startsWith('..')) return { kind: 'outside', rel: resolved, anchor };
  return { kind: 'relative', rel: resolved, anchor };
}

function makeLinkRewriter(project, currentRelPath, siteIndex) {
  return function rewrite(raw, isImage) {
    const r = resolveTarget(currentRelPath, raw);
    if (r.kind === 'absolute' || r.kind === 'anchor') return raw;
    if (r.kind === 'outside') {
      stats.notes.push(`${currentRelPath} -> ${raw}（越出仓库，保留原文）`);
      return raw;
    }
    if (isImage) return rawUrl(project.repo, r.rel);

    const normalized = r.rel.replace(/\/$/, '');
    const siteSlug = siteIndex.get(normalized) ?? siteIndex.get(`${normalized}/README.md`);
    if (siteSlug) {
      const anchor = r.anchor ? `#${slugifySegment(r.anchor)}` : '';
      return `/${siteSlug}/${anchor}`;
    }
    return blobUrl(project.repo, r.rel, r.anchor);
  };
}

function collapseEnglish(body) {
  const m = body.match(/^(##\s+English\s*)$/m);
  if (!m) return body;
  const head = body.slice(0, m.index);
  const tail = body.slice(m.index + m[0].length);
  return `${head.trimEnd()}\n\n<details>\n<summary>English</summary>\n${tail.trim()}\n\n</details>\n`;
}

function convertReadme({ markdown, project, relPath, title, siteIndex, collapse, sidebarLabel, sidebarOrder }) {
  let body = markdown.replace(/^\uFEFF/, '');
  let docTitle = title;

  const h1 = body.match(/^#\s+(.+)$/m);
  if (h1) {
    if (!docTitle) docTitle = stripInline(h1[1]).split(' / ')[0].trim();
    body = body.replace(h1[0], '').trimStart();
  }
  if (collapse) body = collapseEnglish(body);

  const rewrite = makeLinkRewriter(project, relPath, siteIndex);
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

  body = body.replace(/<!--[\s\S]*?-->/g, '');
  body = body.replace(/\{(PROJ_[A-Z_]+|CONDA_[A-Z0-9_]+|SOFTWARE_[A-Z_]+|HOME_[A-Z_]+|TOKEN)\}/g, '`{$1}`');

  return `${frontmatter({ title: docTitle, sidebarLabel, sidebarOrder })}\n${body.trim()}\n`;
}

/* ------------------------------------------------------------ page build */

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
      if (e.isDirectory()) await walk(abs, relPath);
      else if (/^README(\.md)?$/i.test(e.name)) out.push(relPath);
    }
  }
  await walk(repoDir, '');
  return out.sort();
}

function buildPlan(project, readmes) {
  const pages = new Map();
  const siteIndex = new Map();
  const slugFor = (dirRel) =>
    dirRel === ''
      ? `${PROJECTS_DIR}/${project.key}`
      : `${PROJECTS_DIR}/${project.key}/${dirRel.split('/').map(slugifySegment).join('/')}`;

  for (const rel of readmes) {
    const dir = path.posix.dirname(rel) === '.' ? '' : path.posix.dirname(rel);
    if (pages.has(dir)) continue;
    const slug = slugFor(dir);
    pages.set(dir, { slug, readmeRel: rel });
    siteIndex.set(rel, slug);
    siteIndex.set(dir, slug);
  }

  // 上游有些中间目录没有 README（例如 scripts/01_assembly/ 的说明都在各步骤里），
  // 补一个索引页，侧边栏才能形成「阶段 → 步骤」的层级。
  const withReadme = new Set(pages.keys());
  const allDirs = new Set();
  for (const dir of withReadme) {
    if (dir === '') continue;
    const parts = dir.split('/');
    for (let i = 1; i <= parts.length; i++) allDirs.add(parts.slice(0, i).join('/'));
  }
  for (const dir of [...allDirs].sort()) {
    if (withReadme.has(dir)) continue;
    const slug = slugFor(dir);
    pages.set(dir, { slug, readmeRel: null, synthetic: true });
    siteIndex.set(dir, slug);
  }

  return { pages, siteIndex };
}

async function writePage(slug, content) {
  const file = path.join(DOCS_ROOT, `${slug}.md`);
  await mkdir(path.dirname(file), { recursive: true });
  await writeFile(file, content, 'utf8');
  stats.pages++;
}

function extractSummaryTable(markdown) {
  const map = new Map();
  for (const line of markdown.split('\n')) {
    const m = line.match(/^\|\s*\[?`?([^`\]|]+?)\/`?\]?\(?[^|]*\)?\s*\|\s*([^|]*)\|/);
    if (m) {
      const key = m[1].trim();
      const desc = stripInline(m[2]);
      if (key && desc && !/^-+$/.test(key) && key !== '阶段' && key !== 'Stage') map.set(key, desc);
    }
  }
  return map;
}

function buildStageIndex({ project, dirRel, childEntries }) {
  if (!childEntries.length) return '';
  const rows = childEntries
    .map(({ name, slug, desc }) => {
      const label = labelFor(project.key, `${dirRel}/${name}`, name);
      return `| [${label}](/${slug}/) | ${desc || '—'} |`;
    })
    .join('\n');
  return `\n### 本阶段步骤\n\n| 步骤 | 内容 |\n| --- | --- |\n${rows}\n`;
}

async function syncProject(project) {
  const readmes = await collectReadmes(project.dir);
  if (!readmes.length) {
    console.warn(`[sync] 跳过 ${project.key}：目录里没有 README`);
    return null;
  }
  const plan = buildPlan(project, readmes);
  const outDir = path.join(DOCS_ROOT, PROJECTS_DIR, project.key);
  if (existsSync(outDir)) await rm(outDir, { recursive: true, force: true });

  let summaryTables = new Map();
  const rootPage = plan.pages.get('');
  if (rootPage?.readmeRel) {
    const rootMd = await readFile(path.join(project.dir, rootPage.readmeRel.split('/').join(path.sep)), 'utf8');
    summaryTables = extractSummaryTable(rootMd);
  }

  const sidebarEntry = { label: project.displayName, items: [] };

  for (const dirRel of [...plan.pages.keys()].sort()) {
    const { slug, readmeRel, synthetic } = plan.pages.get(dirRel);
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
    const label = isRoot ? project.displayName : labelFor(project.key, dirRel, baseName);
    // 有子页的目录（含自动补的索引页）：自身索引排到最前（-1）。
    // 没有子页的目录（例如独立的 03_telomere）：沿用目录编号，否则它会掉到列表末尾。
    const order = isRoot ? -1 : children.length ? -1 : orderOf(baseName);
    const title = label;

    let content;
    if (synthetic || !readmeRel) {
      const rows = children.map(
        ({ name, slug: s, desc }) => `| [${labelFor(project.key, `${dirRel}/${name}`, name)}](/${s}/) | ${desc || '—'} |`
      );
      const table = rows.length ? `| 步骤 | 内容 |\n| --- | --- |\n${rows.join('\n')}` : '_该目录下暂无可展示的子页面。_';
      content = `${frontmatter({ title, sidebarLabel: label, sidebarOrder: order })}\n本阶段包含以下步骤，内容取自各步骤目录的 README。\n\n${table}\n`;
    } else {
      const markdown = await readFile(path.join(project.dir, readmeRel.split('/').join(path.sep)), 'utf8');
      content = convertReadme({
        markdown,
        project,
        relPath: readmeRel,
        title,
        siteIndex: plan.siteIndex,
        collapse: !isRoot,
        sidebarLabel: isRoot ? project.displayName : label,
        sidebarOrder: order,
      });
      content = `${content.trimEnd()}\n${buildStageIndex({ project, dirRel, childEntries: children })}\n`;
    }
    await writePage(slug, `${content}\n`);

    sidebarEntry.items.push({ slug, label, order, dirRel, childCount: children.length });
  }

  // 侧边栏要体现层级：先放项目首页，再按目录结构嵌套
  const byDir = new Map(sidebarEntry.items.map((i) => [i.dirRel, i]));
  const toItem = (entry) => ({ label: entry.label, slug: entry.slug });
  const buildTree = (dirRel) => {
    const kids = [...byDir.values()]
      .filter((i) => i.dirRel !== dirRel && path.posix.dirname(i.dirRel || '.') === (dirRel || '.'))
      .sort((a, b) => a.order - b.order || a.dirRel.localeCompare(b.dirRel));
    return kids.map((k) => {
      const sub = buildTree(k.dirRel);
      return sub.length ? { label: k.label, items: sub } : toItem(k);
    });
  };
  const root = byDir.get('');
  const tree = root ? [{ label: root.label, slug: root.slug }, ...buildTree('')] : buildTree('');

  stats.projects.push({
    key: project.key,
    label: project.displayName,
    pages: sidebarEntry.items.length,
    sidebar: { label: project.displayName, items: tree },
    rootSlug: root?.slug ?? null,
  });

  console.log(`[sync] ${project.key}: ${sidebarEntry.items.length} 页（${readmes.length} 个 README）`);
  return stats.projects.at(-1);
}

/* ------------------------------------------------------------------ main */

async function main() {
  await mkdir(DOCS_ROOT, { recursive: true });

  // 上一轮生成的清单：决定「快照兜底」时要还原哪几个项目。
  // 不能只看本轮 discoverProjects 的结果——CI 上没有同级仓库，那个列表是空的。
  const prevManifest = await readJsonFile(path.join(SITE_ROOT, 'src', 'data', 'projects.generated.json'), {
    projects: [],
  });
  const prevKeys = (prevManifest.projects ?? []).map((p) => p.key).filter(Boolean);

  const projects = await discoverProjects();
  if (!projects.length) {
    console.info(
      '[sync] 没有发现项目仓库（CI 上属正常：同级目录不存在）。\n' +
        '        将沿用上一次生成的文档快照，保证线上项目文档仍然可用。'
    );
  }
  for (const s of stats.skipped) console.info(`[sync] 未找到本地仓库，沿用快照：${s}`);

  const outRoot = path.join(DOCS_ROOT, PROJECTS_DIR);
  // 上一次生成的结果先备份：本地拿不到仓库时（例如 CI）用它兜底，
  // 这样线上项目文档不会因为「没有同级仓库」而整体消失。
  const backupRoot = `${outRoot}.bak`;
  if (existsSync(backupRoot)) await rm(backupRoot, { recursive: true, force: true });
  const baselineKeys = existsSync(outRoot)
    ? (await readdir(outRoot, { withFileTypes: true })).filter((e) => e.isDirectory()).map((e) => e.name)
    : [];
  if (existsSync(outRoot)) {
    await cp(outRoot, backupRoot, { recursive: true });
    await rm(outRoot, { recursive: true, force: true });
  }

  const syncedKeys = new Set();
  for (const project of projects) {
    if (await syncProject(project)) syncedKeys.add(project.key);
  }

  // 本轮没能从仓库生成的项目：从快照补齐。
  // 注意每个项目在站点上对应两部分：`_projects/<key>/`（子目录页）与 `_projects/<key>.md`（首页），
  // 两者都要补，且只补缺失的项目——「先整体还原再删」会把本轮刚生成的内容覆盖掉。
  const cachedKeys = new Set();
  const restoreCandidates = [...new Set([...prevKeys, ...baselineKeys])];
  for (const key of restoreCandidates) {
    if (syncedKeys.has(key)) continue;
    const backupDir = path.join(backupRoot, key);
    const backupRootPage = path.join(backupRoot, `${key}.md`);
    if (!existsSync(backupDir) && !existsSync(backupRootPage)) continue;

    if (existsSync(backupDir)) {
      const dest = path.join(outRoot, key);
      if (existsSync(dest)) await rm(dest, { recursive: true, force: true });
      await cp(backupDir, dest, { recursive: true });
    }
    if (existsSync(backupRootPage) && !existsSync(path.join(outRoot, `${key}.md`))) {
      await mkdir(outRoot, { recursive: true });
      await cp(backupRootPage, path.join(outRoot, `${key}.md`));
    }

    cachedKeys.add(key);
    const known = (prevManifest.projects ?? []).find((p) => p.key === key);
    stats.projects.push({
      key,
      label: known?.label ?? projectLabelFromSnapshot(backupRootPage) ?? key,
      pages: known?.pages ?? null,
      rootSlug: `${PROJECTS_DIR}/${key}`,
    });
    console.info(`[sync] ${key}: 沿用上一次生成的文档快照`);
  }

  if (existsSync(backupRoot)) await rm(backupRoot, { recursive: true, force: true });

  const manifest = {
    generatedAt: new Date().toISOString(),
    /** 本轮真正从仓库重新生成的 key */
    synced: [...syncedKeys],
    /** 沿用快照的 key */
    cached: [...cachedKeys],
    /** 站点上实际存在文档目录的 key → 这些链接是有效的 */
    available: [...syncedKeys, ...cachedKeys],
    projects: stats.projects
      .map((p) => ({
        key: p.key,
        label: p.label,
        pages: p.pages,
        rootSlug: p.rootSlug,
        from: syncedKeys.has(p.key) ? 'repo' : 'snapshot',
      }))
      .sort((a, b) => a.key.localeCompare(b.key)),
  };
  await writeFile(path.join(SITE_ROOT, 'src', 'data', 'projects.generated.json'), `${JSON.stringify(manifest, null, 2)}\n`, 'utf8');

  if (stats.notes.length) {
    console.log(`[sync] 提示：${stats.notes.length} 个链接越出仓库、已保留原文，示例：`);
    stats.notes.slice(0, 5).forEach((s) => console.log(`   - ${s}`));
  }
  const parts = [`${syncedKeys.size} 个仓库重新生成`];
  if (cachedKeys.size) parts.push(`${cachedKeys.size} 个沿用快照`);
  console.log(`[sync] 完成：${stats.projects.length} 个项目（${parts.join('、')}）、共 ${stats.pages} 页，改写 ${stats.linked} 个链接`);
}

void lstat;
void stat;

await main();
