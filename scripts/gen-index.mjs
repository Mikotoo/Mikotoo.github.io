#!/usr/bin/env node
/**
 * 生成两个索引页：
 *   1. src/content/docs/index-all.md   全部页面一览（按来源分组）
 *   2. src/content/docs/code.md        代码与资源（自动列出仓库里可公开的脚本/图件）
 *
 * 目的是「新增笔记不用改配置」：这两个页面每次构建时按实际内容重新生成，
 * 因此不会出现链接失效或漏列。
 */
import { readdir, readFile, writeFile, stat } from 'node:fs/promises';
import { existsSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const SITE_ROOT = path.resolve(HERE, '..');
const DOCS = path.join(SITE_ROOT, 'src', 'content', 'docs');
const REPO_BLOB = 'https://github.com/Mikotoo/Mikotoo.github.io/blob/main/';
const REPO_TREE = 'https://github.com/Mikotoo/Mikotoo.github.io/tree/main/';

function yamlString(s) {
  return `'${String(s).replace(/'/g, "''")}'`;
}

async function walk(dir, base = '') {
  const out = [];
  for (const e of await readdir(dir, { withFileTypes: true })) {
    const rel = base ? `${base}/${e.name}` : e.name;
    if (e.isDirectory()) out.push(...(await walk(path.join(dir, e.name), rel)));
    else if (e.name.endsWith('.md')) out.push(rel);
  }
  return out;
}

async function readTitle(file) {
  const raw = await readFile(file, 'utf8');
  const fm = raw.match(/^---\r?\n([\s\S]*?)\r?\n---/);
  const t = fm?.[1].match(/^title:\s*(.+)$/m);
  const d = fm?.[1].match(/^date:\s*(.+)$/m);
  const title = (t?.[1] ?? path.basename(file, '.md')).replace(/^['"]|['"]$/g, '').trim();
  return { title, date: d?.[1]?.trim() };
}

/** 在 Lotus_genome / Rscript 中列出可展示的脚本文件 */
async function listScripts(repoDir, exts, limit = 400) {
  if (!existsSync(repoDir)) return [];
  const out = [];
  async function walk(dir) {
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
      if (e.isDirectory()) await walk(abs);
      else if (exts.some((x) => e.name.endsWith(x))) {
        out.push(path.relative(repoDir, abs).split(path.sep).join('/'));
      }
    }
  }
  await walk(repoDir);
  return out.sort();
}

async function main() {
  const files = await walk(DOCS);
  const groups = {
    lotus: [],
    rscript: [],
    notes: [],
    other: [],
  };
  for (const rel of files) {
    const slug = rel.replace(/\.md$/, '');
    if (slug === 'index-all' || slug === 'code') continue;
    const { title, date } = await readTitle(path.join(DOCS, rel));
    const item = { slug, title, date, rel };
    if (slug.startsWith('lotus/')) groups.lotus.push(item);
    else if (slug.startsWith('rscript/')) groups.rscript.push(item);
    else if (slug.startsWith('notes/')) groups.notes.push(item);
    else groups.other.push(item);
  }

  const fmt = (item) => {
    const href = `/${item.slug}/`;
    const label = item.title.replace(/\|/g, '\\|');
    return `| [${label}](${href}) | \`${item.slug}\` |${item.date ? ` ${item.date} |` : ' — |'}`;
  };
  const table = (items, dateHeader = '日期') =>
    [`| 页面 | 路径 | ${dateHeader} |`, '| --- | --- | --- |', ...items.map(fmt)].join('\n');

  const notesSorted = groups.notes.sort((a, b) => String(b.date).localeCompare(String(a.date)));

  const indexAll = `---
title: 内容索引
description: 本站全部页面的清单
---

这是自动生成的页面清单。每次构建时按实际内容重新扫描，不需要手动维护。

## 分析流程 · Lotus T2T（${groups.lotus.length} 页）

${table(groups.lotus)}

## 分析脚本 · Rscript（${groups.rscript.length} 页）

${table(groups.rscript)}

## 笔记与复现（${notesSorted.length} 篇）

${table(notesSorted)}

## 其他页面（${groups.other.length} 页）

${table(groups.other)}
`;

  await writeFile(path.join(DOCS, 'index-all.md'), indexAll, 'utf8');

  // ---- 侧边栏「笔记与复现」：按 _posts 实际内容生成，新增文章无需改配置 ----
  const sidebarPath = path.join(SITE_ROOT, 'src', 'sidebar.generated.mjs');
  const sidebarItems = notesSorted
    .map((n) => `  { label: ${JSON.stringify(n.title)}, slug: ${JSON.stringify(n.slug)} },`)
    .join('\n');
  const sidebarFile = `/**
 * 由 scripts/gen-index.mjs 自动生成，请勿手工编辑。
 * 来源：_posts/*.md（共 ${notesSorted.length} 篇），按日期倒序。
 * 新增文章后重新执行 pnpm gen / pnpm build 即可。
 */
export const notesSidebar = [
${sidebarItems}
];
`;
  await writeFile(sidebarPath, sidebarFile, 'utf8');

  // ---- code.md ----
  const lotusDir = path.resolve(SITE_ROOT, '..', 'Lotus_genome');
  const rscriptDir = path.resolve(SITE_ROOT, '..', 'Rscript');
  const repoScripts = await listScripts(path.join(SITE_ROOT, 'code'), ['.py', '.sh', '.R', '.yaml', '.yml', '.csv', '.html', '.md']);
  const lotusScripts = await listScripts(path.join(lotusDir, 'scripts'), ['.sh', '.py', '.R', '.smk'], 200);
  const rscriptScripts = await listScripts(rscriptDir, ['.sh', '.py', '.R'], 200);

  const localList = repoScripts.length
    ? repoScripts
        .map((rel) => `| \`code/${rel}\` | [查看](${REPO_BLOB}code/${rel.split('/').map(encodeURIComponent).join('/')}) |`)
        .join('\n')
    : '| — | — |';

  const lotusList = lotusScripts.slice(0, 60)
    .map((rel) => `| \`${rel}\` | [GitHub](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/${rel.split('/').map(encodeURIComponent).join('/')}) |`)
    .join('\n');

  const rscriptList = rscriptScripts.slice(0, 60)
    .map((rel) => `| \`${rel}\` | [GitHub](https://github.com/Mikotoo/Rscript/blob/main/${rel.split('/').map(encodeURIComponent).join('/')}) |`)
    .join('\n');

  const codeMd = `---
title: 代码与资源
description: 随文脚本、分析图件，以及可供下载的资料
---

import { CardGrid, LinkCard } from '@astrojs/starlight/components';

<CardGrid>
  <LinkCard title="下载资源" href="/resources/" description="书籍 PDF 等可下载文件" />
  <LinkCard title="内容索引" href="/index-all/" description="本站全部页面清单" />
  <LinkCard title="Lotus_genome 仓库" href="https://github.com/Mikotoo/Lotus_genome" description="T2T 组装、注释与下游分析全部脚本" />
  <LinkCard title="Rscript 仓库" href="https://github.com/Mikotoo/Rscript" description="论文 Fig.2–4 的分析脚本" />
</CardGrid>

## 本仓库中的随文代码（${repoScripts.length} 个文件）

| 文件 | 链接 |
| --- | --- |
${localList}

## Lotus_genome 脚本（共 ${lotusScripts.length} 个，列出前 60）

| 脚本 | 链接 |
| --- | --- |
${lotusList}

## Rscript 脚本（共 ${rscriptScripts.length} 个，列出前 60）

| 脚本 | 链接 |
| --- | --- |
${rscriptList}

## 关于大文件

\`code/single_cell/02_QC/\` 与 \`code/single_cell/02_QC/cache/\` 下曾提交过分析缓存（\`.h5ad\` 约 220 MB、稀疏矩阵导出等）。这些文件**不会**进入本站构建产物（构建后由 \`scripts/prune-dist.mjs\` 从 \`dist/\` 中剔除），但仍保存在仓库里，可直接从 GitHub 访问。

如果以后要缩小仓库体积，可把这类数据移到外部存储（Zenodo / 网盘）并在文中给链接，再用 \`git filter-repo\` 重写历史——注意这会改写提交历史，需单独规划。
`;

  await writeFile(path.join(DOCS, 'code.md'), codeMd, 'utf8');
  console.log(
    `[index] 内容索引：lotus ${groups.lotus.length} / rscript ${groups.rscript.length} / notes ${notesSorted.length} / other ${groups.other.length}；侧边栏写入 ${notesSorted.length} 条；code 页列出 ${repoScripts.length} 个本地文件、${lotusScripts.length} 个 Lotus 脚本、${rscriptScripts.length} 个 Rscript 脚本`
  );

  void REPO_TREE;
  void stat;
  void yamlString;
}

await main();
