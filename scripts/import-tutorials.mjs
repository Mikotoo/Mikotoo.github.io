#!/usr/bin/env node
/** Explicitly import teaching originals. Builds use committed snapshots only. */
import { readFile, readdir, mkdir, writeFile, rm, cp } from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { createHash } from 'node:crypto';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const config = JSON.parse(await readFile(path.join(root, 'config/tutorials.json'), 'utf8'));
const sourceRoot = path.resolve(root, config.sourceDir);
const lessons = [...config.lessons].sort((a, b) => a.order - b.order);
const hash = bytes => createHash('sha256').update(bytes).digest('hex');
const yaml = value => JSON.stringify(value);
const assetBase = '/tutorials/bioinformatics';
const slugs = new Set();
for (const lesson of lessons) {
  if (!/^[a-z0-9]+(?:-[a-z0-9]+)*$/.test(lesson.slug) || slugs.has(lesson.slug)) throw new Error(`Invalid or duplicate lesson slug: ${lesson.slug}`);
  slugs.add(lesson.slug);
}

function formatBody(text) {
  const lines = text.replace(/\r\n/g, '\n').split('\n');
  if (!/^#\s+/.test(lines[0])) throw new Error('Expected the original lesson to begin with a title');
  lines.shift();
  let fence = null;
  const result = lines.map(line => {
    const marker = line.match(/^ {0,3}(`{3,}|~{3,})(.*)$/);
    if (marker) {
      if (!fence) fence = marker[1];
      else if (marker[1][0] === fence[0] && marker[1].length >= fence.length && !marker[2].trim()) fence = null;
      return line;
    }
    if (fence) return line;
    const heading = line.match(/^(#{1,6})\s+(.+)$/);
    if (!heading) return line;
    const level = heading[1].length === 1 || /^[一二三四五六七八九十百]+、/.test(heading[2]) ? 2 : Math.min(6, heading[1].length + 1);
    return `${'#'.repeat(level)} ${heading[2]}`;
  });
  if (fence) throw new Error('Unclosed code fence in teaching original');
  return result.join('\n').trim() + '\n';
}

// ZIP STORE is sufficient for these small text tables; fixed dates keep bytes deterministic.
function zipStore(entries) {
  const locals = [], central = [];
  let offset = 0;
  for (const { name, bytes } of entries) {
    const filename = Buffer.from(name);
    let crc = 0xffffffff;
    for (const byte of bytes) {
      crc ^= byte;
      for (let bit = 0; bit < 8; bit++) crc = (crc >>> 1) ^ (0xedb88320 & -(crc & 1));
    }
    crc = (crc ^ 0xffffffff) >>> 0;
    const local = Buffer.alloc(30);
    local.writeUInt32LE(0x04034b50, 0);
    local.writeUInt16LE(20, 4);
    local.writeUInt16LE(0x800, 6);
    local.writeUInt16LE(0x21, 12);
    local.writeUInt32LE(crc, 14);
    local.writeUInt32LE(bytes.length, 18);
    local.writeUInt32LE(bytes.length, 22);
    local.writeUInt16LE(filename.length, 26);
    locals.push(local, filename, bytes);
    const record = Buffer.alloc(46);
    record.writeUInt32LE(0x02014b50, 0);
    record.writeUInt16LE(20, 4);
    record.writeUInt16LE(20, 6);
    record.writeUInt16LE(0x800, 8);
    record.writeUInt16LE(0x21, 14);
    record.writeUInt32LE(crc, 16);
    record.writeUInt32LE(bytes.length, 20);
    record.writeUInt32LE(bytes.length, 24);
    record.writeUInt16LE(filename.length, 28);
    record.writeUInt32LE(offset, 42);
    central.push(record, filename);
    offset += local.length + filename.length + bytes.length;
  }
  const directory = Buffer.concat(central);
  const end = Buffer.alloc(22);
  end.writeUInt32LE(0x06054b50, 0);
  end.writeUInt16LE(entries.length, 8);
  end.writeUInt16LE(entries.length, 10);
  end.writeUInt32LE(directory.length, 12);
  end.writeUInt32LE(offset, 16);
  return Buffer.concat([...locals, directory, end]);
}

// Read and format all inputs before touching any committed output.
const dataSource = path.join(sourceRoot, config.dataDir);
const data = [];
for (const entry of (await readdir(dataSource, { withFileTypes: true })).sort((a,b) => a.name.localeCompare(b.name, 'en'))) {
  if (!entry.isFile() || !/\.(txt|bed)$/.test(entry.name)) throw new Error(`Unexpected practice-data entry: ${entry.name}`);
  data.push({ name: entry.name, bytes: await readFile(path.join(dataSource, entry.name)) });
}
if (!data.length) throw new Error('No practice data found');
const imported = [];
for (const [index, lesson] of lessons.entries()) {
  const raw = await readFile(path.join(sourceRoot, lesson.source));
  const prev = index ? { label: lessons[index - 1].title, link: `/tutorials/${lessons[index - 1].slug}/` } : { label: '课程目录', link: '/tutorials/' };
  const next = index < lessons.length - 1 ? { label: lessons[index + 1].title, link: `/tutorials/${lessons[index + 1].slug}/` } : false;
  const frontmatter = `---\ntitle: ${yaml(lesson.title)}\ndescription: ${yaml(lesson.description)}\nsidebar:\n  label: ${yaml(lesson.title)}\n  order: ${lesson.order}\nprev: ${yaml(prev)}\nnext: ${yaml(next)}\ntableOfContents:\n  minHeadingLevel: 2\n  maxHeadingLevel: 2\n---\n\n`;
  let intro = `[← 课程目录](/tutorials/) · 第 ${index + 1} 节 / 共 ${lessons.length} 节 · [下载本节原始讲义（Markdown）](${assetBase}/${lesson.slug}.md)\n\n`;
  if (lesson.hasPracticeData) {
    intro += `> **练习数据**：[下载 Lotus 圈图数据包（ZIP）](${assetBase}/lotus-circos-data.zip)。解压后得到 \`lotus_circos_data/\` 文件夹，放在本节工作目录中即可。基础绘图仍使用第一节大豆数据生成的统计表；这个数据包用于后面的 \`circlize\` 进阶练习。\n\n`;
    intro += '<details>\n<summary>单独查看或下载圈图数据文件</summary>\n\n';
    intro += data.map(d => `- [${d.name}](${assetBase}/lotus_circos_data/${encodeURIComponent(d.name)})`).join('\n');
    intro += '\n\n</details>\n\n';
  }
  imported.push({ ...lesson, raw, page: frontmatter + intro + formatBody(raw.toString('utf8')) });
}
const zip = zipStore(data.map(d => ({ name: `lotus_circos_data/${d.name}`, bytes: d.bytes })));
const manifest = {
  lessons: imported.map(l => ({ slug: l.slug, source: l.source, sourceSha256: hash(l.raw), snapshotSha256: hash(l.page) })),
  data: data.map(d => ({ name: d.name, bytes: d.bytes.length, sha256: hash(d.bytes) })),
  archive: { name: 'lotus-circos-data.zip', bytes: zip.length, sha256: hash(zip) },
};
const stage = path.join(root, '.generated/tutorial-import');
await rm(stage, { recursive: true, force: true, maxRetries: 3 });
await mkdir(path.join(stage, 'snapshots'), { recursive: true });
await mkdir(path.join(stage, 'assets/lotus_circos_data'), { recursive: true });
for (const lesson of imported) {
  await writeFile(path.join(stage, 'snapshots', `${lesson.slug}.md`), lesson.page);
  await writeFile(path.join(stage, 'assets', `${lesson.slug}.md`), lesson.raw);
}
for (const item of data) await writeFile(path.join(stage, 'assets/lotus_circos_data', item.name), item.bytes);
await writeFile(path.join(stage, 'assets/lotus-circos-data.zip'), zip);
for (const [from, to] of [['snapshots', 'snapshots/tutorials'], ['assets', 'public/tutorials/bioinformatics']]) {
  const destination = path.join(root, to);
  await rm(destination, { recursive: true, force: true, maxRetries: 3 });
  await mkdir(path.dirname(destination), { recursive: true });
  await cp(path.join(stage, from), destination, { recursive: true });
}
await writeFile(path.join(root, 'snapshots/tutorials.manifest.json'), JSON.stringify(manifest, null, 2) + '\n');
console.log(`[tutorials] Imported ${imported.length} complete lessons, ${data.length} practice-data files, and a ${zip.length}-byte ZIP. Original sources unchanged.`);
