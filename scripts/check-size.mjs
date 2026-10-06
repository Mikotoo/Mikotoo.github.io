#!/usr/bin/env node
/**
 * 构建后体积检查：GitHub Pages 的站点体积软上限为 1 GB，
 * 这里给出明确读数，超过阈值（默认 400 MB）直接失败，避免默默上线一个巨站。
 */
import { readdir, stat } from 'node:fs/promises';
import { existsSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const SITE_ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const DIST = path.join(SITE_ROOT, 'dist');
const THRESHOLD_MB = Number(process.env.SITE_SIZE_LIMIT_MB ?? 400);

async function walk(dir, acc) {
  for (const e of await readdir(dir, { withFileTypes: true })) {
    const abs = path.join(dir, e.name);
    if (e.isDirectory()) await walk(abs, acc);
    else {
      const s = await stat(abs);
      acc.total += s.size;
      acc.files.push({ rel: path.relative(DIST, abs).split(path.sep).join('/'), size: s.size });
    }
  }
}

if (!existsSync(DIST)) {
  console.error('[size] dist 不存在');
  process.exit(1);
}

const acc = { total: 0, files: [] };
await walk(DIST, acc);
const mb = acc.total / 1024 / 1024;
acc.files.sort((a, b) => b.size - a.size);

console.log(`[size] 站点产物：${acc.files.length} 个文件，合计 ${mb.toFixed(1)} MB`);
for (const f of acc.files.slice(0, 5)) {
  console.log(`   ${(f.size / 1024 / 1024).toFixed(2)} MB  ${f.rel}`);
}

if (mb > THRESHOLD_MB) {
  console.error(`[size] 超出阈值 ${THRESHOLD_MB} MB，请检查 public/ 下是否有不该入站的大文件`);
  process.exit(1);
}
