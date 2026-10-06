#!/usr/bin/env node
/**
 * 把 Jekyll 时代的 _posts/*.md 迁移为 Starlight 笔记页。
 *
 * 设计：**源文件继续留在 _posts/ 下**（保持原有的 git 写作习惯），每次构建前重新生成
 * src/content/docs/notes/*.md。这样你仍然用「新建 _posts/日期-标题.md」的方式写文章，
 * 不需要学新目录结构。
 *
 * 处理内容：
 *   - Jekyll frontmatter（title/date/category/layout）→ Starlight frontmatter
 *   - 正文里指向本仓库的图片地址（raw.githubusercontent.com / github.com/.../raw/）
 *     → 复制到 public/media/<slug>/ 并用站点绝对路径 /media/<slug>/... 引用
 *   - 其余外链原样保留
 */
import { readFile, writeFile, mkdir, rm, readdir, copyFile, stat } from 'node:fs/promises';
import { existsSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const SITE_ROOT = path.resolve(HERE, '..');
const POSTS_DIR = path.join(SITE_ROOT, '_posts');
const NOTES_OUT = path.join(SITE_ROOT, 'src', 'content', 'docs', 'notes');
/** public/ 下的文件原样进入站点产物，用绝对路径引用最稳妥 */
const MEDIA_ROOT = path.join(SITE_ROOT, 'public', 'media');

const RAW_PREFIXES = [
  'https://github.com/Mikotoo/Mikotoo.github.io/raw/main/',
  'https://github.com/Mikotoo/Mikotoo.github.io/raw/master/',
  'https://raw.githubusercontent.com/Mikotoo/Mikotoo.github.io/main/',
  'https://raw.githubusercontent.com/Mikotoo/Mikotoo.github.io/master/',
];

/** 个别文章用稳定的语义化 slug，便于固定链接 */
const SLUG_OVERRIDES = {
  '2024-4-7-steroSeq': 'soybean-snrna',
};

function slugFor(base) {
  if (SLUG_OVERRIDES[base]) return SLUG_OVERRIDES[base];
  return base
    .replace(/^\d{4}-\d{1,2}-\d{1,2}-/, '')
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, '-')
    .replace(/^-+|-+$/g, '');
}

function parseFrontmatter(raw) {
  const m = raw.match(/^---\r?\n([\s\S]*?)\r?\n---\r?\n?/);
  if (!m) return { data: {}, body: raw };
  const data = {};
  for (const line of m[1].split(/\r?\n/)) {
    const kv = line.match(/^([A-Za-z_][\w-]*):\s*(.*)$/);
    if (kv) data[kv[1]] = kv[2].replace(/^["']|["']$/g, '').trim();
  }
  return { data, body: raw.slice(m[0].length) };
}

function yamlString(s) {
  return `'${String(s).replace(/'/g, "''")}'`;
}

/** 把上游图片地址换成本站静态资源，返回改写后的正文。 */
async function localizeImages(body, slug, copied) {
  let out = body;
  for (const prefix of RAW_PREFIXES) {
    const re = new RegExp(prefix.replace(/[.*+?^${}()|[\]\\]/g, '\\$&') + '([^)\\s"]+)', 'g');
    for (const m of out.matchAll(re)) {
      const repoRel = decodeURIComponent(m[1]);
      const src = path.join(SITE_ROOT, repoRel.split('/').join(path.sep));
      if (!existsSync(src)) {
        console.warn(`  [warn] 资源不存在，保留外链：${repoRel}`);
        continue;
      }
      // downloads/image/blog2_sirna/srna_1.png → <slug>/blog2_sirna/srna_1.png
      const relParts = repoRel.split('/');
      const tail = relParts[0] === 'downloads' ? relParts.slice(2) : relParts;
      const destRel = path.join(slug, ...tail);
      const dest = path.join(MEDIA_ROOT, destRel);
      if (!copied.has(destRel)) {
        await mkdir(path.dirname(dest), { recursive: true });
        await copyFile(src, dest);
        copied.add(destRel);
      }
      const url = `/media/${destRel.split(path.sep).map(encodeURIComponent).join('/')}`;
      out = out.split(m[0]).join(url);
    }
  }
  return out;
}

async function main() {
  if (!existsSync(POSTS_DIR)) {
    console.warn('[notes] 没有 _posts 目录，跳过');
    return;
  }

  // 清掉上一次的生成结果（含早期版本遗留的 src/assets/notes）
  const legacyAssets = path.join(SITE_ROOT, 'src', 'assets', 'notes');
  if (existsSync(legacyAssets)) await rm(legacyAssets, { recursive: true, force: true });
  if (existsSync(NOTES_OUT)) await rm(NOTES_OUT, { recursive: true, force: true });
  await mkdir(NOTES_OUT, { recursive: true });
  if (existsSync(MEDIA_ROOT)) await rm(MEDIA_ROOT, { recursive: true, force: true });
  await mkdir(MEDIA_ROOT, { recursive: true });
  // 媒体文件是构建期从 downloads/ 复制出来的派生物，不入库；
  // 这里放一个 .gitkeep，保证 public/media/ 目录本身存在于仓库中。
  await writeFile(path.join(MEDIA_ROOT, '.gitkeep'), '', 'utf8');

  const files = (await readdir(POSTS_DIR)).filter((f) => /\.(md|markdown)$/i.test(f)).sort();
  const copied = new Set();
  const written = [];
  let bytes = 0;

  for (const file of files) {
    const base = path.basename(file, path.extname(file));
    const slug = slugFor(base);
    const raw = await readFile(path.join(POSTS_DIR, file), 'utf8');
    const { data, body } = parseFrontmatter(raw);

    const date = (data.date || base.slice(0, 10)).trim();
    const title = (data.title || base).trim();
    let content = await localizeImages(body.trim(), slug, copied);

    // MDX 安全：HTML 注释在 MDX 里属于 JSX 语法，先移除避免意外
    content = content.replace(/<!--[\s\S]*?-->/g, '');

    // 侧边栏用「日期倒序」：order 越小越靠前，因此取 9999 - YYYYMMDD
    const compact = date.replace(/[^0-9]/g, '').slice(0, 8);
    const order = /^\d{8}$/.test(compact) ? 99999999 - Number(compact) : 99999999;

    const frontmatter = `${[
      '---',
      `title: ${yamlString(title)}`,
      `date: ${date}`,
      data.category ? `category: ${yamlString(data.category)}` : null,
      'sidebar:',
      `  order: ${order}`,
      '---',
      '',
    ]
      .filter((l) => l !== null)
      .join('\n')}\n`;

    await writeFile(path.join(NOTES_OUT, `${slug}.md`), `${frontmatter}\n${content}\n`, 'utf8');
    written.push({ slug, file, title, date });
  }

  for (const rel of copied) {
    try {
      bytes += (await stat(path.join(MEDIA_ROOT, rel))).size;
    } catch {
      /* ignore */
    }
  }

  console.log(
    `[notes] 生成 ${written.length} 篇笔记，本地化 ${copied.size} 个资源（${(bytes / 1024 / 1024).toFixed(1)} MB → public/media/）`
  );
  for (const w of written) console.log(`   ${w.date}  ${w.slug}  ←  ${w.file}`);
}

await main();
