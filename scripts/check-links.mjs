#!/usr/bin/env node
/** Check HTML href/src targets, including relative URLs. Anchors and external
 * websites are outside this check. Strict by default; --warn-only is optional.
 * --dir <path> supports isolated fixtures and alternate build directories.
 */
import { readdir, readFile, stat } from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { site } from '../config/site.mjs';
const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const args = process.argv.slice(2);
const dirIndex = args.indexOf('--dir');
const DIST = path.resolve(ROOT, dirIndex < 0 ? 'dist' : args[dirIndex+1]);
const ORIGIN = new URL(site.url).origin;
async function walk(dir) {
  const files = [];
  for (const e of await readdir(dir, { withFileTypes: true })) {
    const full = path.join(dir,e.name);
    if (e.isDirectory()) files.push(...await walk(full));
    else if (e.name.endsWith('.html')) files.push(full);
  }
  return files;
}
const cache = new Map();
async function resolves(urlPath) {
  if (cache.has(urlPath)) return cache.get(urlPath);
  let clean;
  try { clean = decodeURIComponent(urlPath).replace(/^\/+/, ''); } catch { return false; }
  const base = path.resolve(DIST, clean);
  if (path.relative(DIST,base).startsWith('..')) return false;
  const candidates = [base, path.join(base,'index.html')];
  // Extensionless pages may be emitted as foo.html; trailing / still needs an index.
  if (!urlPath.endsWith('/')) candidates.push(base+'.html');
  let found = false;
  for (const file of candidates) {
    try { if ((await stat(file)).isFile()) { found = true; break; } } catch {}
  }
  cache.set(urlPath,found);
  return found;
}
const files = await walk(DIST);
if (!files.length) throw new Error('[links] No HTML files to check.');
const broken = new Map();
let checked = 0;
for (const file of files) {
  const rel = path.relative(DIST,file).split(path.sep).join('/');
  const pagePath = '/'+rel.replace(/index\.html$/, '');
  const html = (await readFile(file,'utf8')).replace(/<script\b[^>]*>[\s\S]*?<\/script>/gi, '').replace(/<!--[\s\S]*?-->/g, '');
  for (const match of html.matchAll(/<[a-z][^>]*>/gi)) {
    const tag = match[0];
    for (const attr of tag.matchAll(/\b(href|src)\s*=\s*(["'])(.*?)\2/gi)) {
      const raw = attr[3].replaceAll('&amp;','&');
      if (!raw || raw.startsWith('#')) continue;
      let url;
      try { url = new URL(raw, ORIGIN+pagePath); } catch { continue; }
      if (!['http:','https:'].includes(url.protocol) || url.origin !== ORIGIN) continue;
      // Starlight emits /404/ as canonical metadata for 404.html. This is not a navigation link.
      if (rel === '404.html' && /^<link\b/i.test(tag) && /\brel=["']canonical["']/i.test(tag) && url.pathname === '/404/') continue;
      checked++;
      if (!await resolves(url.pathname)) {
        if (!broken.has(url.pathname)) broken.set(url.pathname,new Set());
        broken.get(url.pathname).add('/'+rel);
      }
    }
  }
}
console.log(`[links] ${files.length} HTML pages, ${checked} internal href/src targets checked.`);
for (const [target,from] of broken) console.error(`[links] Missing ${target} ← ${[...from].slice(0,3).join(', ')}`);
if (broken.size) {
  console.error(`[links] ${broken.size} missing targets.`);
  process.exitCode = args.includes('--warn-only') ? 0 : 1;
} else console.log('[links] No missing internal files.');
