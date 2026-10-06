#!/usr/bin/env node
/**
 * Explicit snapshot update (pnpm sync); ordinary builds consume committed snapshots.
 * Only config/projects.json registers projects. Sources are read-only.
 * Files: snapshots/projects/<key>.md and <key>/...; URLs: /_projects/<key>/...
 * Preflight reads every README before staging; publish only after all projects succeed.
 */
import { readFile, writeFile, mkdir, rm, readdir, rename } from 'node:fs/promises';
import { existsSync, statSync } from 'node:fs';
import { createHash } from 'node:crypto';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const SITE_ROOT = path.resolve(HERE, '..');
const SNAPSHOT_ROOT = path.join(SITE_ROOT, 'snapshots');
const STAGING_ROOT = path.join(SITE_ROOT, '.generated', 'sync-staging');
const PROJECTS_URL_PREFIX = '_projects';
const LABELS = JSON.parse(await readFile(path.join(SITE_ROOT, 'config', 'project-labels.json'), 'utf8'));
const stats = { pages: 0, projects: [], linked: 0, notes: [] };

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

/* --------------------------------------------------------- source preflight */

async function readRepoRemote(dir, entry) {
  const configured = entry.href?.match(/^https:\/\/github\.com\/([^/]+)\/([^/]+?)\/?$/);
  if (configured) return { owner: configured[1], name: configured[2].replace(/\.git$/, ''), branch: entry.branch ?? 'main' };
  const fallback = { owner: entry.owner ?? '', name: entry.repo ?? path.basename(dir), branch: entry.branch ?? 'main' };
  const cfg = path.join(dir, '.git', 'config');
  if (!existsSync(cfg)) {
    if (fallback.owner) return fallback;
    throw new Error(`Project ${entry.key} needs a GitHub href in config/projects.json`);
  }
  const text = await readFile(cfg, 'utf8');
  const m = text.match(/\[remote "origin"\][\s\S]*?url\s*=\s*(.+)/);
  const gh = m?.[1].trim().match(/github\.com[/:]([^/]+)\/([^/\s]+?)(?:\.git)?$/);
  return gh ? { owner: gh[1], name: gh[2], branch: fallback.branch } : fallback;
}

async function collectReadmes(repoDir) {
  const out = [];
  async function walk(dir, rel) {
    const entries = await readdir(dir, { withFileTypes: true });
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

async function loadProjects() {
  const config = JSON.parse(await readFile(path.join(SITE_ROOT, 'config', 'projects.json'), 'utf8'));
  if (!Array.isArray(config.projects)) throw new Error('config/projects.json must contain a projects array');
  const keys = new Set();
  const projects = [];
  for (const entry of config.projects) {
    const key = entry?.key;
    if (typeof key !== 'string' || !/^[a-z0-9]+(?:-[a-z0-9]+)*$/.test(key) || keys.has(key)) {
      throw new Error(`Invalid or duplicate registered project key: ${key}`);
    }
    keys.add(key);
    // dir names are relative to the site's parent, as in the existing registry.
    // An explicit absolute dir can register a source located elsewhere.
    const source = entry.dir ?? entry.repo ?? key;
    if (typeof source !== 'string' || !source.trim()) throw new Error(`Invalid source dir for ${key}`);
    const dir = path.resolve(SITE_ROOT, '..', source);
    if (dir === SITE_ROOT || !statSync(dir).isDirectory()) throw new Error(`Invalid source directory for ${key}: ${dir}`);
    const readmes = await collectReadmes(dir);
    if (!readmes.some((rel) => !rel.includes('/'))) throw new Error(`Registered source ${key} has no root README: ${dir}`);
    const contents = new Map();
    const hash = createHash('sha256');
    for (const rel of readmes) {
      const bytes = await readFile(path.join(dir, ...rel.split('/')));
      contents.set(rel, bytes);
      // Length-prefixed UTF-8 relative names and raw bytes avoid ambiguous boundaries.
      hash.update(`${Buffer.byteLength(rel, 'utf8')}:`).update(rel, 'utf8');
      hash.update(`${bytes.length}:`).update(bytes);
    }
    const project = {
      key, dir, readmes, contents,
      displayName: entry.title ?? path.basename(dir),
      repo: await readRepoRemote(dir, entry),
      sourceSha256: hash.digest('hex'),
    };
    project.plan = buildPlan(project, readmes);
    projects.push(project);
  }
  return projects.sort((a, b) => a.key < b.key ? -1 : a.key > b.key ? 1 : 0);
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

function buildPlan(project, readmes) {
  const pages = new Map();
  const siteIndex = new Map();
  const slugFor = (dirRel) => dirRel === ''
    ? `${PROJECTS_URL_PREFIX}/${project.key}`
    : `${PROJECTS_URL_PREFIX}/${project.key}/${dirRel.split('/').map(slugifySegment).join('/')}`;
  for (const rel of readmes) {
    const dir = path.posix.dirname(rel) === '.' ? '' : path.posix.dirname(rel);
    if (pages.has(dir)) continue;
    const slug = slugFor(dir);
    pages.set(dir, { slug, readmeRel: rel });
    siteIndex.set(rel, slug);
    siteIndex.set(dir, slug);
  }
  // Missing intermediate READMEs get index pages to preserve the sidebar hierarchy.
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
  const slugs = new Set();
  for (const [dir, { slug }] of pages) {
    if (slug.split('/').some((segment) => !segment) || slugs.has(slug)) {
      throw new Error(`Empty or duplicate page slug in ${project.key}: ${dir} -> ${slug}`);
    }
    slugs.add(slug);
  }
  return { pages, siteIndex };
}

async function writePage(slug, content) {
  const relative = slug.slice(`${PROJECTS_URL_PREFIX}/`.length);
  const file = path.join(STAGING_ROOT, 'projects', `${relative}.md`);
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
  const rows = childEntries.map(({ name, slug, desc }) => {
    const label = labelFor(project.key, `${dirRel}/${name}`, name);
    return `| [${label}](/${slug}/) | ${desc || '—'} |`;
  }).join('\n');
  return `\n### 本阶段步骤\n\n| 步骤 | 内容 |\n| --- | --- |\n${rows}\n`;
}

async function syncProject(project) {
  const plan = project.plan;
  const rootPage = plan.pages.get('');
  const summaryTables = extractSummaryTable(project.contents.get(rootPage.readmeRel).toString('utf8'));
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
      const markdown = project.contents.get(readmeRel).toString('utf8');
      content = convertReadme({
        markdown, project, relPath: readmeRel, title, siteIndex: plan.siteIndex,
        collapse: !isRoot, sidebarLabel: isRoot ? project.displayName : label, sidebarOrder: order,
      });
      content = `${content.trimEnd()}\n${buildStageIndex({ project, dirRel, childEntries: children })}\n`;
    }
    await writePage(slug, `${content}\n`);
  }
  stats.projects.push({
    key: project.key,
    label: project.displayName,
    pages: plan.pages.size,
    rootSlug: rootPage.slug,
    sourceSha256: project.sourceSha256,
  });
  console.log(`[sync] ${project.key}: ${plan.pages.size} 页（${project.readmes.length} 个 README）`);
}

/* ------------------------------------------------------------------ main */

async function main() {
  // Missing/unreadable registered sources fail here, before any filesystem mutation.
  const projects = await loadProjects();
  await rm(STAGING_ROOT, { recursive: true, force: true });
  await mkdir(path.join(STAGING_ROOT, 'projects'), { recursive: true });
  for (const project of projects) await syncProject(project);
  await writeFile(path.join(STAGING_ROOT, 'manifest.json'), `${JSON.stringify({ projects: stats.projects }, null, 2)}\n`, 'utf8');

  // All input reads and rendering have succeeded. Replace the complete generated tree
  // so removed projects/pages disappear too. Publication is not a filesystem transaction.
  await mkdir(SNAPSHOT_ROOT, { recursive: true });
  await rm(path.join(SNAPSHOT_ROOT, 'projects'), { recursive: true, force: true });
  await rename(path.join(STAGING_ROOT, 'projects'), path.join(SNAPSHOT_ROOT, 'projects'));
  await rename(path.join(STAGING_ROOT, 'manifest.json'), path.join(SNAPSHOT_ROOT, 'manifest.json'));
  await rm(STAGING_ROOT, { recursive: true, force: true });
  if (stats.notes.length) {
    console.log(`[sync] 提示：${stats.notes.length} 个链接越出仓库、已保留原文，示例：`);
    stats.notes.slice(0, 5).forEach((s) => console.log(`   - ${s}`));
  }
  console.log(`[sync] 完成：${stats.projects.length} 个项目、共 ${stats.pages} 页，改写 ${stats.linked} 个链接`);
}

await main();
