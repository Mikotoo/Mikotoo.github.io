#!/usr/bin/env node
/** Assemble Starlight input from editable content and committed snapshots.
 * Does not read sibling repositories or update any source file.
 */
import { readdir, readFile, writeFile, mkdir, rm, cp } from 'node:fs/promises';
import { existsSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { site } from '../config/site.mjs';

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const DOCS = path.join(ROOT, 'src/content/docs');
const read = (rel) => readFile(path.join(ROOT, rel), 'utf8');
const manifest = JSON.parse(await read('config/projects.json'));
const keys = new Set();
for (const p of manifest.projects) {
  if (!/^[a-z0-9]+(?:-[a-z0-9]+)*$/.test(p.key) || keys.has(p.key)) throw new Error(`Invalid or duplicate project key: ${p.key}`);
  keys.add(p.key);
  if (p.landing && (!/^projects\/[a-z0-9-]+$/.test(p.landing) || !existsSync(path.join(ROOT, `content/pages/${p.landing}.md`)))) throw new Error(`Invalid project landing page: ${p.landing}`);
}
const projects = [...manifest.projects].sort((a,b) => (a.order ?? 999) - (b.order ?? 999));
const course = JSON.parse(await read('config/tutorials.json'));
const lessons = [...course.lessons].sort((a,b) => a.order - b.order);
const lessonSlugs = new Set();
for (const lesson of lessons) {
  if (!/^[a-z0-9]+(?:-[a-z0-9]+)*$/.test(lesson.slug) || lessonSlugs.has(lesson.slug)) throw new Error(`Invalid or duplicate lesson slug: ${lesson.slug}`);
  lessonSlugs.add(lesson.slug);
  if (!existsSync(path.join(ROOT, `snapshots/tutorials/${lesson.slug}.md`))) throw new Error(`Missing tutorial snapshot: ${lesson.slug}. Run pnpm import:tutorials and commit the outputs.`);
}

async function walk(dir, prefix = '') {
  const result = [];
  for (const e of (await readdir(dir, { withFileTypes: true })).sort((a,b) => a.name.localeCompare(b.name))) {
    const rel = prefix + e.name;
    if (e.isDirectory()) result.push(...await walk(path.join(dir,e.name), `${rel}/`));
    else result.push(rel);
  }
  return result;
}
function metadata(raw) {
  const fm = raw.match(/^---\r?\n([\s\S]*?)\r?\n---/)?.[1] ?? '';
  const value = (name) => fm.match(new RegExp(`^${name}:\\s*(.+)$`, 'm'))?.[1]?.trim().replace(/^['"]|['"]$/g, '').replace(/''/g, "'");
  return { title: value('title'), date: value('date') ?? '', label: value('\\s+label'), order: Number(value('\\s+order') ?? 999) };
}
const esc = (s) => String(s).replaceAll('&','&amp;').replaceAll('<','&lt;').replaceAll('>','&gt;').replaceAll('"','&quot;');
const cell = (s) => String(s).replaceAll('|','\\|');
const archiveList = (items) => `<ul class="archive-list">${items.map(p => `<li><a href="/${p.slug}/">${esc(p.title)}</a>${p.date ? `<span>${esc(p.date)}</span>` : ''}</li>`).join('')}</ul>`;

// This entire directory is disposable; authoring files live exclusively in content/.
await rm(DOCS, { recursive: true, force: true });
await mkdir(DOCS, { recursive: true });
await cp(path.join(ROOT,'content/pages'), DOCS, { recursive: true });
await cp(path.join(ROOT,'content/notes'), path.join(DOCS,'notes'), { recursive: true });
await mkdir(path.join(DOCS, 'tutorials'), { recursive: true });
for (const lesson of lessons) await cp(path.join(ROOT, `snapshots/tutorials/${lesson.slug}.md`), path.join(DOCS, `tutorials/${lesson.slug}.md`));
await mkdir(path.join(DOCS,'_projects'), { recursive: true });
for (const p of projects) {
  const root = path.join(ROOT,`snapshots/projects/${p.key}.md`);
  if (!existsSync(root)) throw new Error(`Missing snapshot for ${p.key}. Run pnpm sync and commit snapshots/.`);
  await cp(root, path.join(DOCS,`_projects/${p.key}.md`));
  const dir = path.join(ROOT,`snapshots/projects/${p.key}`);
  if (existsSync(dir)) await cp(dir, path.join(DOCS,`_projects/${p.key}`), { recursive: true });
}
const pages = [];
for (const rel of await walk(DOCS)) {
  if (!rel.endsWith('.md')) continue;
  const meta = metadata(await readFile(path.join(DOCS,rel),'utf8'));
  pages.push({ ...meta, title: meta.title ?? rel, slug: rel.slice(0,-3) });
}
const notes = pages.filter(p => p.slug.startsWith('notes/')).sort((a,b) => b.date.localeCompare(a.date) || a.slug.localeCompare(b.slug));
const others = pages.filter(p => !p.slug.includes('/') && !['index','index-all','intro'].includes(p.slug));
const readingConfig = JSON.parse(await read('config/readings.json')).readings;
const readingByKey = new Map(readingConfig.map(r => [r.key, r]));
if (readingByKey.size !== readingConfig.length) throw new Error('Duplicate reading keys');
for (const r of readingConfig) {
  if (!notes.some(n => n.slug === `notes/${r.key}`)) throw new Error(`Missing note for reading card: ${r.key}`);
}
const readings = notes.map(n => ({ ...n, ...(readingByKey.get(n.slug.slice(6)) ?? {}) }));
const noteCard = n => `<a class="reading-card" href="/${n.slug}/"><span class="note-topic">${esc(n.topic ?? '技术笔记')}</span><h3>${esc(n.title)}</h3>${n.summary ? `<p>${esc(n.summary)}</p>` : ''}<span class="reading-meta"><time datetime="${esc(n.date)}">${esc(n.date)}</time><span aria-hidden="true">↗</span></span></a>`;
const readingGrid = items => `<div class="reading-grid">${items.map(noteCard).join('\n')}</div>`;
const projectCard = p => `<article class="project-card"><p class="eyebrow">${esc(p.eyebrow ?? 'RESEARCH PROJECT')}</p><h3><a href="/${esc(p.landing ?? `_projects/${p.key}`)}/">${esc(p.title)}</a></h3><p>${esc(p.summary ?? '')}</p><div class="tags">${(p.tags ?? []).map(t => `<span>${esc(t)}</span>`).join('')}</div><div class="card-actions"><a href="/${esc(p.landing ?? `_projects/${p.key}`)}/">探索项目 →</a>${p.href && p.showRepository !== false ? `<a href="${esc(p.href)}">GitHub ↗</a>` : ''}</div></article>`;
const projectGrid = items => `<div class="project-grid">${items.map(projectCard).join('\n')}</div>`;
const link = (p) => ({ label: p.label ?? p.title, slug: p.slug });
function projectTree(p) {
  const rootSlug = `_projects/${p.key}`;
  const docs = pages.filter(d => d.slug === rootSlug || d.slug.startsWith(rootSlug+'/'));
  const children = (parent) => docs.filter(d => path.posix.dirname(d.slug) === parent).sort((a,b) => {
    const number = d => Number(path.posix.basename(d.slug).match(/^\d+/)?.[0] ?? 999);
    return number(a)-number(b) || a.order-b.order || a.slug.localeCompare(b.slug);
  }).map(d => {
    const sub = children(d.slug);
    return sub.length ? { label: d.label ?? d.title, items: [link(d), ...sub] } : link(d);
  });
  return { label: p.title, items: [...(p.landing ? [{ label: '项目介绍', slug: p.landing }] : []), link(docs.find(d => d.slug === rootSlug)), ...children(rootSlug)] };
}
async function fill(file, replacements) {
  let text = await readFile(path.join(DOCS,file),'utf8');
  for (const [name,value] of Object.entries(replacements)) {
    const marker = `<!-- AUTO:${name} -->`;
    if (text.split(marker).length !== 2) throw new Error(`${file}: expected exactly one ${marker}`);
    text = text.replace(marker, () => value);
  }
  await writeFile(path.join(DOCS,file), text);
}
await fill('index.md', {
  FEATURED_PROJECTS: projectGrid(projects.filter(p => p.featured)),
  FEATURED_NOTES: readingGrid(readingConfig.filter(r => r.featured).map(r => readings.find(n => n.key === r.key))),
});
await fill('projects.md', { PROJECT_CARDS: projectGrid(projects) });
await fill('notes.md', { NOTE_COLLECTION: readingGrid(readings) });
await fill('tutorials.md', {
  LESSON_CARDS: `<div class="project-grid">${lessons.map(l => `<article class="project-card"><p class="eyebrow">LESSON ${String(l.order).padStart(2, '0')}</p><h3><a href="/tutorials/${l.slug}/">${esc(l.title)}</a></h3><p>${esc(l.description)}</p><div class="tags">${(l.topics ?? []).map(t => `<span>${esc(t)}</span>`).join('')}</div><div class="card-actions"><a href="/tutorials/${l.slug}/">进入本节 →</a><a href="/tutorials/bioinformatics/${l.slug}.md" download>下载讲义 ↓</a></div></article>`).join('\n')}</div>`,
});
await fill('index-all.md', {
  CONTENT_INDEX: `<h2>生物信息学基础教程</h2>\n${archiveList(lessons.map(l => ({ title: l.title, slug: `tutorials/${l.slug}` })))}\n<h2>笔记与复现</h2>\n${archiveList(readings)}\n<h2>项目分析步骤</h2>\n` + projects.map(p => `<details><summary>${esc(p.title)}</summary>${archiveList(pages.filter(d => d.slug === `_projects/${p.key}` || d.slug.startsWith(`_projects/${p.key}/`)))}</details>`).join('\n') + `\n<h2>继续探索</h2>\n${archiveList(others)}`,
});
// Explicit inventory prevents a local untracked file from changing published content.
const inventory = JSON.parse(await read('config/code-files.json'));
for (const rel of inventory.files) {
  if (!rel.startsWith('code/') || rel.includes('..') || !existsSync(path.join(ROOT,rel))) throw new Error(`Invalid code inventory entry: ${rel}`);
}
const blob = `${site.repository}/blob/${site.branch}/`;
await fill('code.md', {
  CODE_FILES: ['| 文件 | 链接 |','| --- | --- |',...inventory.files.map(rel => `| \`${cell(rel)}\` | [查看](${blob}${rel.split('/').map(encodeURIComponent).join('/')}) |`)].join('\n'),
  PROJECT_REPOS: projectGrid(projects),
});
await mkdir(path.join(ROOT,'.generated'), { recursive: true });
await writeFile(path.join(ROOT,'.generated/sidebar.mjs'), `// Generated by pnpm gen. Do not edit.\nexport const projectsSidebar = ${JSON.stringify(projects.map(projectTree),null,2)};\nexport const notesSidebar = ${JSON.stringify(notes.map(link),null,2)};\nexport const tutorialsSidebar = ${JSON.stringify(lessons.map(l => ({ label: l.title, slug: `tutorials/${l.slug}` })),null,2)};\n`);
console.log(`[content] ${pages.length} pages assembled from content/ and snapshots/; ${notes.length} notes, ${projects.length} projects.`);
