#!/usr/bin/env node
/** Watch authoring sources; Astro watches the assembled content as usual. */
import { watch } from 'node:fs';
import { spawn, spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
const root = fileURLToPath(new URL('../', import.meta.url));
let timer;
const watchers = ['content', 'config', 'snapshots'].map(dir => watch(new URL(`../${dir}/`, import.meta.url), { recursive: true }, () => {
  clearTimeout(timer);
  timer = setTimeout(() => {
    const result = spawnSync(process.execPath, ['scripts/gen-index.mjs'], { cwd: root, stdio: 'inherit' });
    if (result.status !== 0) console.error('[dev] 内容装配失败，请修正上述问题后保存文件。');
  }, 250);
}));
const astro = spawn(process.execPath, ['node_modules/astro/bin/astro.mjs', 'dev', ...process.argv.slice(2)], { cwd: root, stdio: 'inherit' });
function cleanup() { clearTimeout(timer); watchers.forEach(w => w.close()); }
astro.on('exit', code => { cleanup(); process.exitCode = code ?? 1; });
astro.on('error', error => { cleanup(); console.error(error); process.exitCode = 1; });
for (const signal of ['SIGINT','SIGTERM']) process.on(signal, () => { cleanup(); astro.kill(signal); });
