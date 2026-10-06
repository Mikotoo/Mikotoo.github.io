#!/usr/bin/env node
/**
 * 生成站点里「按现状自动成形」的部分：
 *   src/content/docs/projects.md        项目文档总览
 *   src/content/docs/index-all.md       全部页面一览
 *   src/content/docs/code.md            代码与资源
 *   src/sidebar.generated.mjs           完整侧边栏（项目文档树 + 笔记）
 *
 * 数据来源：
 *   src/data/projects.json           手工维护的项目清单（展示顺序、首页文案）
 *   src/data/projects.generated.json 上一轮 sync-docs 实际生成的项目文档
 *   _posts/*.md                      笔记
 *
 * 目标：新增一个项目或一篇笔记时，都不需要手改这个脚本。
 */
import { readdir, readFile, writeFile } from 'node:fs/promises';
import { existsSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const SITE_ROOT = path.resolve(HERE, '..');
const DOCS = path.join(SITE_ROOT, 'src', 'content', 'docs');
const PROJECTS_DIR = '_projects';
const REPO_BLOB = 'https://github.com/Mikotoo/Mikotoo.github.io/blob/main/';
const REPO_TREE = 'https://github.com/Mikotoo/Mikotoo.github.io/tree/main/';

const readJson = async (f, fallback) => {
  if (!existsSync(f)) return fallback;
  try {
    return JSON.parse(await readFile(f, 'utf8'));
  } catch {
    return fallback;
  }
};

async function walk(dir, base = '') {
  const out = [];
  if (!existsSync(dir)) return out;
  for (const e of await readdir(dir, { withFileTypes: true })) {
    const rel = base ? `${base}/${e.name}` : e.name;
    if (e.isDirectory()) {
      // 跳过以 _ 开头的目录（项目文档单独处理）
      if (base === '' && e.name.startsWith('_')) continue;
      out.push(...(await walk(path.join(dir, e.name), rel)));
    } else if (e.name.endsWith('.md')) out.push(rel);
  }
  return out;
}

async function readFrontmatter(file) {
  const raw = await readFile(file, 'utf8');
  const fm = raw.match(/^---\r?\n([\s\S]*?)\r?\n---/);
  const get = (k) => fm?.[1].match(new RegExp(`^${k}:\\s*(.+)$`, 'm'))?.[1];
  return {
    title: (get('title') ?? path.basename(file, '.md')).replace(/^['"]|['"]$/g, '').trim(),
    date: get('date')?.trim(),
  };
}

/** 递归列出项目文档页，返回目录相对路径（相对项目根）等元信息，便于拼侧边栏 */
const PROJECT_ROOT = '.'; // 项目首页的哨兵值：避免与同名子目录（如 _projects/lotus/）撞键

async function listProjectPages(key) {
  const root = path.join(DOCS, PROJECTS_DIR, key);
  const pages = [];

  const readMeta = async (file) => {
    const raw = await readFile(file, 'utf8');
    const fm = raw.match(/^---\r?\n([\s\S]*?)\r?\n---/);
    return {
      label: fm?.[1].match(/^\s*label:\s*(.+)$/m)?.[1]?.replace(/^['"]|['"]$/g, '').trim(),
      order: Number(fm?.[1].match(/^\s*order:\s*(-?\d+)/m)?.[1] ?? 999),
      title: fm?.[1].match(/^title:\s*(.+)$/m)?.[1]?.replace(/^['"]|['"]$/g, '').trim(),
    };
  };

  // 项目首页：_projects/<key>.md（仓库根 README）
  const rootFile = path.join(DOCS, `${PROJECTS_DIR}/${key}.md`);
  if (existsSync(rootFile)) {
    const meta = await readMeta(rootFile);
    pages.push({
      dirRel: PROJECT_ROOT,
      slug: `${PROJECTS_DIR}/${key}`,
      label: meta.label ?? meta.title ?? key,
      order: Number.isFinite(meta.order) ? meta.order : 999,
    });
  }

  const walkPages = async (dir, rel) => {
    let entries;
    try {
      entries = await readdir(dir, { withFileTypes: true });
    } catch {
      return;
    }
    for (const e of entries) {
      const r = rel ? `${rel}/${e.name}` : e.name;
      if (e.isDirectory()) {
        await walkPages(path.join(dir, e.name), r);
        continue;
      }
      if (!e.name.endsWith('.md')) continue;

      const meta = await readMeta(path.join(dir, e.name));
      const dirRel = `${rel ? `${rel}/` : ''}${e.name}`.replace(/\/?index\.md$/, '').replace(/\.md$/, '');
      const slug = `${PROJECTS_DIR}/${key}/${r.replace(/\.md$/, '')}`.replace(/\/index$/, '');
      pages.push({
        dirRel,
        slug,
        label: meta.label ?? meta.title ?? r,
        order: Number.isFinite(meta.order) ? meta.order : 999,
      });
    }
  };
  if (existsSync(root)) await walkPages(root, '');
  // 去重（极端情况下同一 slug 出现两次）
  const seen = new Set();
  return pages.filter((p) => (seen.has(p.slug) ? false : seen.add(p.slug)));
}

function buildSidebarTree(pages) {
  const byDir = new Map(pages.map((p) => [p.dirRel, p]));
  const parentOf = (dirRel) => path.posix.dirname(dirRel);

  const childrenOf = (dirRel) =>
    [...byDir.values()].filter((p) => p.dirRel !== dirRel && parentOf(p.dirRel) === dirRel);

  const build = (dirRel) => {
    const kids = childrenOf(dirRel);
    /**
     * 项目文档的侧边栏条目分两类：
     *   0 = 目录级条目（分组标题 order:-1，或独占一个阶段的目录如 03_telomere）→ 按阶段编号排
     *   1 = 具体步骤页 → 按自身编号排
     * 同属 0 时，分组标题（order:-1）排在同编号的独立阶段之前，
     * 这样「Hi-C 挂载」这类分组会待在自己阶段内部的正确位置。
     */
    const kind = (p) => {
      const base = path.posix.basename(p.dirRel);
      if (p.order === -1) return 0;
      const n = Number(base.match(/^(\d+)/)?.[1]);
      return Number.isFinite(n) && n >= 1 && n <= 99 ? 0 : 1;
    };
    const stageNo = (p) => {
      const n = Number(path.posix.basename(p.dirRel).match(/^(\d+)/)?.[1]);
      return Number.isFinite(n) ? n : 999;
    };
    const isGroupHeader = (p) => p.order === -1;
    kids.sort((a, b) => {
      const ka = kind(a);
      const kb = kind(b);
      if (ka !== kb) return ka - kb;
      if (ka === 0) {
        const d = stageNo(a) - stageNo(b);
        if (d !== 0) return d;
        if (isGroupHeader(a) !== isGroupHeader(b)) return isGroupHeader(a) ? -1 : 1;
        return a.dirRel.localeCompare(b.dirRel);
      }
      const d = a.order - b.order;
      if (d !== 0) return d;
      const depth = (p) => p.dirRel.split('/').length;
      if (depth(a) !== depth(b)) return depth(a) - depth(b);
      return a.dirRel.localeCompare(b.dirRel);
    });
    return kids.map((k) => {
      const sub = build(k.dirRel);
      return sub.length ? { label: k.label, items: [toLink(k), ...sub] } : toLink(k);
    });
  };
  return build(PROJECT_ROOT);
}

/** 侧边栏叶子项 */
function toLink(page) {
  return { label: page.label, slug: page.slug };
}

async function listScripts(repoDir, exts, limit = 400) {
  if (!existsSync(repoDir)) return [];
  const out = [];
  const walkAll = async (dir) => {
    if (out.length >= limit) return;
    let entries;
    try {
      entries = await readdir(dir, { withFileTypes: true });
    } catch {
      return;
    }
    for (const e of entries) {
      if (out.length >= limit) return;
      if (e.name === '.git' || e.name === 'node_modules') continue;
      const abs = path.join(dir, e.name);
      if (e.isDirectory()) await walkAll(abs);
      else if (exts.some((x) => e.name.endsWith(x))) {
        out.push(path.relative(repoDir, abs).split(path.sep).join('/'));
      }
    }
  };
  await walkAll(repoDir);
  return out.sort();
}

async function main() {
  const manifest = await readJson(path.join(SITE_ROOT, 'src', 'data', 'projects.json'), { projects: [] });
  const generated = await readJson(path.join(SITE_ROOT, 'src', 'data', 'projects.generated.json'), { projects: [] });
  const generatedByKey = new Map((generated.projects ?? []).map((p) => [p.key, p]));

  // 合并：以手工清单的顺序/文案为准，页数等运行信息取自生成结果
  const projects = (manifest.projects ?? [])
    .map((p) => ({ ...p, ...(generatedByKey.get(p.key) ?? {}), summary: p.summary, href: p.href, repo: p.repo }))
    .filter((p) => p.key)
    .sort((a, b) => (a.order ?? 999) - (b.order ?? 999));
  const orphan = (generated.projects ?? []).filter((g) => !projects.some((p) => p.key === g.key));

  /* ---------------- 项目文档总览页 ---------------- */
  const card = (p) =>
    `  <a href="/${PROJECTS_DIR}/${p.key}/">\n` +    `    <strong>${p.label ?? p.key}</strong>\n` +
    `    <span>${p.summary ?? ''}</span>\n` +
    `  </a>`;
  const cards = projects.length ? projects.map(card).join('\n') : '  <p>暂无项目文档。</p>';
  const chips = projects
    .filter((p) => p.tags?.length)
    .map((p) => `| ${p.label ?? p.key} | ${p.tags.join(' · ')} | [\`${p.repo ?? p.key}\`](${p.href ?? '#'}) |`)
    .join('\n');

  const projectsMd = `---
title: 项目文档
description: 已收录的项目仓库与其分析流程文档
---

<div class="hero-cards">
${cards}
</div>

## 收录约定

本站的项目文档**不在站点里手写**：构建时自动抓取项目仓库里的所有 \`README.md\`，按目录结构生成带导航和检索的文档区。收录一个项目有两种方式：

| 方式 | 做法 | 适用场景 |
| --- | --- | --- |
| 登记（推荐） | 在 \`src/data/projects.json\` 里加一条，用 \`dir\`/\`repo\` 指向本地仓库 | 自己的项目，需要定制标题、排序、关键词 |
| 自动发现 | 在仓库根目录放一个 \`content.config.ts\`（内容随意，只是标记） | 外部工具/临时克隆的仓库，不想先登记 |

无论哪种方式，都要求：

1. 仓库根目录有 \`README.md\`（每个步骤目录也建议一份，讲清做法、输入、参数）；
2. 仓库放在本站**同级目录**（本地开发）或 \`vendor/\` 下（想作为 submodule 收进站点仓库）。

找不到本地仓库时会被跳过并给出提示，不影响其它内容构建。

${chips ? `## 一览\n\n| 项目 | 关键词 | 代码仓库 |\n| --- | --- | --- |\n${chips}\n` : ''}
${orphan.length ? `\n## 已同步但未登记\n\n以下项目已生成文档，但还没写进 \`src/data/projects.json\`（因此不会出现在首页卡片里）：\n\n${orphan.map((o) => `- \`${o.key}\`（${o.pages} 页）`).join('\n')}\n` : ''}
`;

  await writeFile(path.join(DOCS, 'projects.md'), projectsMd, 'utf8');

  /* ---------------- 全部页面 ---------------- */
  const files = await walk(DOCS);
  const groups = { notes: [], other: [] };
  for (const rel of files) {
    const slug = rel.replace(/\.md$/, '');
    if (['index-all', 'code', 'projects'].includes(slug)) continue;
    const { title, date } = await readFrontmatter(path.join(DOCS, rel));
    const item = { slug, title, date };
    if (slug.startsWith('notes/')) groups.notes.push(item);
    else groups.other.push(item);
  }
  const notesSorted = groups.notes.sort((a, b) => String(b.date).localeCompare(String(a.date)));

  const fmt = (item) => `| [${item.title.replace(/\|/g, '\\|')}](/${item.slug}/) | \`${item.slug}\` |${item.date ? ` ${item.date} |` : ' — |'}`;
  const table = (items) => [`| 页面 | 路径 | 日期 |`, '| --- | --- | --- |', ...items.map(fmt)].join('\n');

  const projectSections = [];
  for (const p of projects) {
    const pages = await listProjectPages(p.key);
    if (!pages.length) continue;
    projectSections.push(
      `## 项目文档 · ${p.label ?? p.key}（${pages.length} 页）\n\n${table(
        pages
          .sort((a, b) => a.slug.localeCompare(b.slug))
          .map((x) => ({ slug: x.slug, title: x.label, date: null }))
      )}`
    );
  }

  const indexAll = `---
title: 内容索引
description: 本站全部页面的清单
---

这是自动生成的页面清单，每次构建时按实际内容重新扫描，不需要手动维护。

${projectSections.join('\n\n')}

## 笔记与复现（${notesSorted.length} 篇）

${table(notesSorted)}

## 其他页面（${groups.other.length} 页）

${table(groups.other)}
`;

  await writeFile(path.join(DOCS, 'index-all.md'), indexAll, 'utf8');

  /* ---------------- 侧边栏 ---------------- */
  const projectItems = [];
  for (const p of projects) {
    const pages = await listProjectPages(p.key);
    if (!pages.length) continue;
    const tree = buildSidebarTree(pages);
    // 项目首页（仓库根 README）放在该项目的第一个条目，方便回到总览
    const root = pages.find((x) => x.dirRel === PROJECT_ROOT);
    projectItems.push({
      label: p.label ?? p.key,
      items: root ? [{ label: root.label, slug: root.slug }, ...tree] : tree,
    });
  }
  const noteItems = notesSorted.map((n) => ({ label: n.title, slug: n.slug }));

  const sidebarFile = `/**
 * 由 scripts/gen-index.mjs 自动生成，请勿手工编辑。
 *
 * 项目文档：${projectItems.map((p) => p.label).join('、') || '（无）'}
 * 笔记：${noteItems.length} 篇，按日期倒序
 *
 * 新增项目请改 src/data/projects.json；新增笔记直接放 _posts/。
 */
export const projectsSidebar = ${JSON.stringify(projectItems, null, 2)};

export const notesSidebar = ${JSON.stringify(noteItems, null, 2)};
`;
  await writeFile(path.join(SITE_ROOT, 'src', 'sidebar.generated.mjs'), sidebarFile, 'utf8');

  /* ---------------- 代码与资源 ---------------- */
  const repoScripts = await listScripts(path.join(SITE_ROOT, 'code'), ['.py', '.sh', '.R', '.yaml', '.yml', '.csv', '.html', '.md']);
  const localList = repoScripts.length
    ? repoScripts
        .map((rel) => `| \`code/${rel}\` | [查看](${REPO_BLOB}code/${rel.split('/').map(encodeURIComponent).join('/')}) |`)
        .join('\n')
    : '| — | — |';

  const repoBlocks = [];
  for (const p of projects) {
    const dir = existsSync(path.resolve(SITE_ROOT, '..', p.repo ?? ''))
      ? path.resolve(SITE_ROOT, '..', p.repo)
      : path.resolve(SITE_ROOT, 'vendor', p.repo ?? '');
    const scripts = await listScripts(dir, ['.sh', '.py', '.R', '.smk'], 60);
    if (!scripts.length) continue;
    const rows = scripts
      .map(
        (rel) =>
          `| \`${rel}\` | [GitHub](https://github.com/Mikotoo/${p.repo}/blob/main/${rel.split('/').map(encodeURIComponent).join('/')}) |`
      )
      .join('\n');
    repoBlocks.push(`## ${p.label ?? p.key} 脚本（列出前 ${scripts.length} 个）\n\n| 脚本 | 链接 |\n| --- | --- |\n${rows}`);
  }

  const codeMd = `---
title: 代码与资源
description: 随文脚本、分析图件，以及可供下载的资料
---

<div class="hero-cards">
  <a href="/resources/">
    <strong>下载资源</strong>
    <span>书籍 PDF 等可下载文件</span>
  </a>
  <a href="/index-all/">
    <strong>内容索引</strong>
    <span>本站全部页面清单</span>
  </a>
  <a href="/projects/">
    <strong>项目文档</strong>
    <span>已收录的代码仓库与分析流程</span>
  </a>
</div>

## 本仓库中的随文代码（${repoScripts.length} 个文件）

| 文件 | 链接 |
| --- | --- |
${localList}

${repoBlocks.join('\n\n')}

## 关于大文件

\`code/single_cell/\` 目录下曾有分析缓存（\`.h5ad\`、稀疏矩阵导出等，单个文件可达数百 MB）。这些文件不会进入站点构建产物，需要的话直接从 GitHub 仓库获取。

如果以后要缩小仓库体积，可把这类数据移到外部存储（Zenodo / 网盘）并在文中给链接，再用 \`git filter-repo\` 重写历史——这会改写提交历史，需单独规划。
`;

  await writeFile(path.join(DOCS, 'code.md'), codeMd, 'utf8');

  console.log(
    `[index] 项目文档 ${projectItems.length} 个（${projectItems.map((p) => p.label).join('、') || '无'}）；` +
      `笔记 ${noteItems.length} 篇；其他页面 ${groups.other.length}；code 页 ${repoScripts.length} 个本地文件`
  );
  void REPO_TREE;
}

await main();
