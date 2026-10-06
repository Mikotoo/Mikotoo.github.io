#!/usr/bin/env node
/**
 * 构建后结构检查：扫 dist 里所有 HTML 的站内链接与图片，确认目标文件真实存在。
 *
 * 为什么需要它：站点内容一半是「构建时生成」的（项目文档、笔记、索引页），
 * 一旦某个生成环节在 CI 上退化（例如没有同级仓库），就会出现指向 404 的链接。
 * 这类问题在本地预览时看不出来，所以放在构建末尾自动检查。
 *
 * 用法：node scripts/check-links.mjs        （由 pnpm build 的 postbuild 调用）
 *     严格模式：CHECK_LINKS_STRICT=1     发现死链直接失败退出
 */
import { readdir, readFile, stat } from 'node:fs/promises';
import { existsSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const SITE_ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const DIST = path.join(SITE_ROOT, 'dist');
/** 站点自身的域名：只有这些链接算站内链接 */
const SITE_HOST = 'mikotoo.github.io';

/**
 * 已知且可接受的目标（不视为死链）：
 *   /404/  —— Starlight 默认 404 页面里的「返回首页」旧链接，属于主题自带页面
 */
const IGNORED = new Set(['/404/']);

async function walk(dir, acc = []) {
  for (const e of await readdir(dir, { withFileTypes: true })) {
    const abs = path.join(dir, e.name);
    if (e.isDirectory()) await walk(abs, acc);
    else if (e.name.endsWith('.html')) acc.push(abs);
  }
  return acc;
}

/** 判断站内路径是否对应 dist 里的真实文件 */
async function resolves(distRel) {
  let clean = decodeURIComponent(distRel.split('#')[0].split('?')[0]);
  if (!clean.startsWith('/')) return true; // 相对链接暂不处理
  if (IGNORED.has(clean)) return true;
  clean = clean.replace(/^\//, '');
  if (clean === '') return existsSync(path.join(DIST, 'index.html'));
  const candidates = [
    path.join(DIST, clean),
    path.join(DIST, `${clean}.html`),
    path.join(DIST, clean, 'index.html'),
  ];
  for (const c of candidates) {
    try {
      const s = await stat(c);
      if (s.isFile()) return true;
    } catch {
      /* 继续试下一个 */
    }
  }
  return false;
}

const files = existsSync(DIST) ? await walk(DIST) : [];
if (!files.length) {
  console.error('[links] dist 为空，跳过检查');
  process.exit(0);
}

const broken = new Map(); // 目标 -> Set(来源页)
let checked = 0;
let media = 0;

for (const file of files) {
  const html = await readFile(file, 'utf8');
  const from = `/${path.relative(DIST, file).split(path.sep).join('/')}`;

  const hrefs = [...html.matchAll(/(?:href|src)="([^"]+)"/g)].map((m) => m[1]);
  for (const raw of hrefs) {
    if (!raw || raw.startsWith('#') || raw.startsWith('mailto:')) continue;
    if (/^(https?:)?\/\//.test(raw)) {
      // 绝对 URL 里的本站链接也要检查（例如 index.md 里写死的 /notes/xxx/）
      try {
        const u = new URL(raw.startsWith('//') ? `https:${raw}` : raw);
        if (u.host !== SITE_HOST) continue;
        checked++;
        if (!(await resolves(u.pathname))) {
          if (!broken.has(u.pathname)) broken.set(u.pathname, new Set());
          broken.get(u.pathname).add(from);
        }
      } catch {
        /* 非法 URL 忽略 */
      }
      continue;
    }
    if (!raw.startsWith('/') || raw.startsWith('//')) continue;
    if (/\.(png|jpe?g|gif|svg|webp|ico|woff2?|ttf|otf|eot|pdf|txt|xml|json|js|css)$/i.test(raw)) media++;
    checked++;
    if (!(await resolves(raw))) {
      if (!broken.has(raw)) broken.set(raw, new Set());
      broken.get(raw).add(from);
    }
  }
}

const strict = process.env.CHECK_LINKS_STRICT === '1';
console.log(`[links] 检查 ${files.length} 个页面、${checked} 个站内链接（其中静态资源 ${media} 个）`);
if (!broken.size) {
  console.log('[links] 未发现死链');
  process.exit(0);
}

console.error(`[links] 发现 ${broken.size} 个死链目标：`);
for (const [target, froms] of [...broken].slice(0, 25)) {
  console.error(`   ${target}   ←  ${[...froms].slice(0, 3).join(', ')}${froms.size > 3 ? ` 等 ${froms.size} 处` : ''}`);
}
if (broken.size > 25) console.error(`   … 以及另外 ${broken.size - 25} 个`);

if (strict) {
  console.error('[links] CHECK_LINKS_STRICT=1，构建失败');
  process.exit(1);
}
console.error('[links] 提示：设 CHECK_LINKS_STRICT=1 可让构建在此失败');
