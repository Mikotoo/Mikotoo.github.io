#!/usr/bin/env node
/**
 * 构建后清理 dist：把不适合放进站点的体积型文件剔除。
 *
 * 原因：本仓库历史上把分析缓存（h5ad）和大图件一起提交进了 git，总计约 220 MB+。
 * Astro 会把 public/ 原样复制到 dist/，若不清理，GitHub Pages 产物会接近上限。
 * 这些文件仍可通过 GitHub 原始链接访问，只是不进站点产物。
 */
import { readdir, rm, stat } from 'node:fs/promises';
import { existsSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const SITE_ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const DIST = path.join(SITE_ROOT, 'dist');

/** 从站点产物中排除的规则 */
const RULES = [
  { test: (rel) => rel.endsWith('.h5ad'), why: '单细胞分析缓存' },
  { test: (rel) => rel.endsWith('sparse_matrix.csv'), why: '稀疏矩阵导出' },
  { test: (rel) => rel.endsWith('.pyc'), why: 'Python 字节码' },
  { test: (rel) => /__pycache__|\.pytest_cache/.test(rel), why: 'Python 缓存目录' },
];

const removed = [];
let freed = 0;

async function walk(dir) {
  const entries = await readdir(dir, { withFileTypes: true });
  for (const e of entries) {
    const abs = path.join(dir, e.name);
    const rel = path.relative(DIST, abs).split(path.sep).join('/');
    if (e.isDirectory()) {
      if (RULES.some((r) => r.test(rel) && rel.split('/').length === new URL('file:///').pathname.length)) {
        // 目录级命中在下面统一处理
      }
      await walk(abs);
      // 清理空目录
      try {
        const rest = await readdir(abs);
        if (rest.length === 0) await rm(abs, { recursive: true, force: true });
      } catch {
        /* ignore */
      }
      continue;
    }
    const rule = RULES.find((r) => r.test(rel));
    if (!rule) continue;
    try {
      const s = await stat(abs);
      freed += s.size;
      await rm(abs, { force: true });
      removed.push({ rel, size: s.size, why: rule.why });
    } catch {
      /* ignore */
    }
  }
}

if (!existsSync(DIST)) {
  console.warn('[prune] dist 不存在，跳过');
  process.exit(0);
}

await walk(DIST);

if (removed.length) {
  console.log(`[prune] 从站点产物移除 ${removed.length} 个文件，释放 ${(freed / 1024 / 1024).toFixed(1)} MB`);
  for (const r of removed.slice(0, 10)) {
    console.log(`   - ${(r.size / 1024 / 1024).toFixed(1)} MB  ${r.rel}  (${r.why})`);
  }
  if (removed.length > 10) console.log(`   … 以及另外 ${removed.length - 10} 个`);
} else {
  console.log('[prune] 站点产物中没有需要剔除的大文件');
}
