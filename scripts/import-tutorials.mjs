#!/usr/bin/env node
/** Explicitly import teaching originals. Builds use committed snapshots and assets only. */
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
const MAX_ASSET_BYTES = 5 * 1024 * 1024;
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

// ZIP STORE is sufficient for these small teaching files; fixed date fields keep bytes deterministic.
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

/* ------------------------------------------------------------------ */
/* Teaching demo reads are generated from the authored FASTA + GTF, so */
/* the practice dataset is reproducible and carries no real data.      */
/* ------------------------------------------------------------------ */
function mulberry32(seed) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
const REV = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N' };
const revcomp = s => s.split('').reverse().map(c => REV[c] ?? 'N').join('');

function parseFasta(text) {
  const seqs = new Map();
  let name = null;
  for (const line of text.replace(/\r\n/g, '\n').split('\n')) {
    if (line.startsWith('>')) { name = line.slice(1).split(/\s+/)[0]; seqs.set(name, ''); }
    else if (name && line.trim()) seqs.set(name, seqs.get(name) + line.trim().toUpperCase());
  }
  if (!seqs.size) throw new Error('Practice FASTA contains no sequences');
  return seqs;
}
function parseGtf(text) {
  const transcripts = new Map();
  for (const line of text.replace(/\r\n/g, '\n').split('\n')) {
    if (!line || line.startsWith('#')) continue;
    const f = line.split('\t');
    if (f.length < 9 || f[2] !== 'exon') continue;
    const gene = f[8].match(/gene_id\s+"([^"]+)"/)?.[1];
    const tx = f[8].match(/transcript_id\s+"([^"]+)"/)?.[1] ?? gene;
    const key = `${f[0]}|${tx}`;
    if (!transcripts.has(key)) transcripts.set(key, { chrom: f[0], gene, tx, strand: f[6], exons: [] });
    transcripts.get(key).exons.push([Number(f[3]), Number(f[4])]);
  }
  if (!transcripts.size) throw new Error('Practice GTF contains no exon records');
  return [...transcripts.values()].map(t => {
    t.exons.sort((a, b) => a[0] - b[0]);
    const ordered = t.strand === '-' ? [...t.exons].reverse() : t.exons;
    t.orderedExons = ordered;
    const lengths = ordered.map(([s, e]) => e - s + 1);
    let acc = 0;
    t.junctions = [];
    for (const length of lengths) { acc += length; t.junctions.push(acc); }
    t.junctions.pop();                                    // the last boundary is the transcript end
    return t;
  });
}
const transcribe = (seqs, tx) => {
  const joined = tx.orderedExons.map(([s, e]) => seqs.get(tx.chrom).slice(s - 1, e)).join('');
  return tx.strand === '-' ? revcomp(joined) : joined;
};

function generateDemoReads(practice, fastaText, gtfText) {
  const r = practice.reads;
  const seqs = parseFasta(fastaText);
  const transcripts = parseGtf(gtfText);
  const random = mulberry32(20240407);
  const int = (lo, hi) => lo + Math.floor(random() * (hi - lo + 1));
  const fragments = [];
  const push = (template, start, fragLen) => fragments.push({
    r1: template.slice(start, start + r.readLength),
    r2: revcomp(template.slice(start + fragLen - r.readLength, start + fragLen)),
  });
  for (const tx of transcripts) {
    const seq = transcribe(seqs, tx);
    const from = fragments.length;
    // Guarantee reads that straddle every exon-exon junction, so CIGAR shows N.
    for (const junction of tx.junctions) {
      for (let i = 0; i < r.junctionPairsPerSite; i++) {
        const start = Math.max(0, junction - Math.floor(r.readLength * 0.4));
        const maxFrag = Math.min(r.maxFragment, seq.length - start);
        if (maxFrag < r.readLength * 2) throw new Error(`${tx.tx}: transcript too short for a junction-spanning pair`);
        push(seq, start, maxFrag >= r.minFragment ? int(r.minFragment, maxFrag) : maxFrag);
      }
    }
    while (fragments.length - from < r.pairsPerGene) {
      const maxFrag = Math.min(r.maxFragment, seq.length);
      if (maxFrag < r.readLength * 2) throw new Error(`${tx.tx}: transcript too short for a fragment`);
      const fragLen = maxFrag >= r.minFragment ? int(r.minFragment, maxFrag) : maxFrag;
      push(seq, int(0, seq.length - fragLen), fragLen);
    }
  }
  // Reads from a block duplicated on two chromosomes: each has two equally good loci.
  const block = seqs.get(r.repeatRegion.chrom).slice(r.repeatRegion.start - 1, r.repeatRegion.end);
  if (block.length !== r.repeatRegion.end - r.repeatRegion.start + 1) throw new Error('Repeat region falls outside the practice genome');
  for (let i = 0; i < r.repeatPairs; i++) {
    const fragLen = int(Math.max(r.minFragment, r.readLength * 2), Math.min(r.maxFragment, block.length));
    push(block, int(0, block.length - fragLen), fragLen);
  }
  // Reads that do not exist in the genome at all, so the BAM also holds unmapped records.
  const noise = mulberry32(987654321);
  for (let i = 0; i < r.noisePairs; i++) {
    const seq = Array.from({ length: r.readLength }, () => 'ACGT'[Math.floor(noise() * 4)]).join('');
    fragments.push({ r1: seq, r2: revcomp(seq) });
  }
  const expected = transcripts.length * r.pairsPerGene + r.repeatPairs + r.noisePairs;
  if (fragments.length !== expected) throw new Error(`Generated ${fragments.length} pairs, expected ${expected}`);
  // Real FASTQ is not grouped by locus; shuffle reproducibly.
  for (let i = fragments.length - 1; i > 0; i--) {
    const j = Math.floor(random() * (i + 1));
    [fragments[i], fragments[j]] = [fragments[j], fragments[i]];
  }
  const quality = () => Array.from({ length: r.readLength }, () => String.fromCharCode(int(r.qualityMin, r.qualityMax) + 33)).join('');
  const r1 = [], r2 = [];
  fragments.forEach((f, i) => {
    const name = `demo_${String(i + 1).padStart(5, '0')}`;
    r1.push(`@${name}`, f.r1, '+', quality());
    r2.push(`@${name}`, f.r2, '+', quality());
  });
  return {
    pairs: fragments.length,
    // Plain FASTQ: some static servers answer `.gz` with Content-Encoding: gzip, which makes a
    // browser save decompressed bytes under a `.gz` name. The lecture explains the `.gz` convention.
    files: [
      { name: 'demo_reads_R1.fastq', bytes: Buffer.from(`${r1.join('\n')}\n`) },
      { name: 'demo_reads_R2.fastq', bytes: Buffer.from(`${r2.join('\n')}\n`) },
    ],
  };
}

/* ---------------------- read every input before writing ---------------------- */
const imported = [];
const practices = new Map();
for (const [index, lesson] of lessons.entries()) {
  const raw = await readFile(path.join(sourceRoot, lesson.source));
  const practice = lesson.practice ?? null;
  let files = null;
  if (practice) {
    if (!/^[a-z0-9_-]+$/.test(practice.assetName)) throw new Error(`Invalid assetName: ${practice.assetName}`);
    const dir = path.join(sourceRoot, practice.dir);
    const collected = [];
    for (const entry of (await readdir(dir, { withFileTypes: true })).sort((a, b) => a.name.localeCompare(b.name, 'en'))) {
      if (!entry.isFile()) throw new Error(`Practice asset directories must be flat: ${practice.dir}/${entry.name}`);
      const bytes = await readFile(path.join(dir, entry.name));
      if (bytes.length > MAX_ASSET_BYTES) throw new Error(`Practice asset too large: ${entry.name}`);
      collected.push({ name: entry.name, bytes, generated: false });
    }
    if (!collected.length) throw new Error(`No practice files found in ${practice.dir}`);
    if (practice.generateReads) {
      const fasta = collected.find(f => /\.(fa|fasta)$/.test(f.name));
      const gtf = collected.find(f => f.name.endsWith('.gtf'));
      if (!fasta || !gtf) throw new Error(`${practice.dir}: generateReads needs one .fa/.fasta and one .gtf`);
      const produced = generateDemoReads(practice, fasta.bytes.toString('utf8'), gtf.bytes.toString('utf8'));
      collected.push(...produced.files.map(f => ({ ...f, generated: true })));
      console.log(`[tutorials] ${lesson.slug}: generated ${produced.pairs} teaching read pairs from the practice reference`);
    }
    files = collected.sort((a, b) => a.name.localeCompare(b.name, 'en'));
    const archive = zipStore(files.map(f => ({ name: `${practice.archiveRoot}/${f.name}`, bytes: f.bytes })));
    practices.set(lesson.slug, { practice, files, archive });
  }
  const prev = index ? { label: lessons[index - 1].title, link: `/tutorials/${lessons[index - 1].slug}/` } : { label: '课程目录', link: '/tutorials/' };
  const next = index < lessons.length - 1 ? { label: lessons[index + 1].title, link: `/tutorials/${lessons[index + 1].slug}/` } : false;
  const frontmatter = `---\ntitle: ${yaml(lesson.title)}\ndescription: ${yaml(lesson.description)}\nsidebar:\n  label: ${yaml(lesson.title)}\n  order: ${lesson.order}\nprev: ${yaml(prev)}\nnext: ${yaml(next)}\ntableOfContents:\n  minHeadingLevel: 2\n  maxHeadingLevel: 2\n---\n\n`;
  let intro = `[← 课程目录](/tutorials/) · 第 ${index + 1} 节 / 共 ${lessons.length} 节 · [下载本节原始讲义（Markdown）](${assetBase}/${lesson.slug}.md)\n\n`;
  if (practice) {
    intro += `> **练习数据**：[下载${practice.label}](${assetBase}/${practice.archive})。${practice.summary}\n\n`;
    intro += '<details>\n<summary>单独查看或下载数据文件</summary>\n\n';
    intro += files.map(f => `- [${f.name}](${assetBase}/${practice.assetName}/${encodeURIComponent(f.name)})`).join('\n');
    intro += '\n\n</details>\n\n';
  }
  imported.push({ ...lesson, raw, page: frontmatter + intro + formatBody(raw.toString('utf8')) });
}

/* --------------------------------- publish --------------------------------- */
const stage = path.join(root, '.generated/tutorial-import');
await rm(stage, { recursive: true, force: true, maxRetries: 3 });
await mkdir(path.join(stage, 'snapshots'), { recursive: true });
await mkdir(path.join(stage, 'assets'), { recursive: true });
for (const lesson of imported) {
  await writeFile(path.join(stage, 'snapshots', `${lesson.slug}.md`), lesson.page);
  await writeFile(path.join(stage, 'assets', `${lesson.slug}.md`), lesson.raw);
}
for (const { practice, files, archive } of practices.values()) {
  await mkdir(path.join(stage, 'assets', practice.assetName), { recursive: true });
  for (const file of files) await writeFile(path.join(stage, 'assets', practice.assetName, file.name), file.bytes);
  await writeFile(path.join(stage, 'assets', practice.archive), archive);
}
for (const [from, to] of [['snapshots', 'snapshots/tutorials'], ['assets', 'public/tutorials/bioinformatics']]) {
  const destination = path.join(root, to);
  await rm(destination, { recursive: true, force: true, maxRetries: 3 });
  await mkdir(path.dirname(destination), { recursive: true });
  await cp(path.join(stage, from), destination, { recursive: true });
}
const manifest = {
  lessons: imported.map(l => {
    const published = practices.get(l.slug);
    return {
      slug: l.slug,
      source: l.source,
      sourceSha256: hash(l.raw),
      snapshotSha256: hash(l.page),
      practice: published ? {
        assetName: published.practice.assetName,
        files: published.files.map(f => ({ name: f.name, bytes: f.bytes.length, sha256: hash(f.bytes), generated: f.generated })),
        archive: { name: published.practice.archive, bytes: published.archive.length, sha256: hash(published.archive) },
      } : null,
    };
  }),
};
await writeFile(path.join(root, 'snapshots/tutorials.manifest.json'), JSON.stringify(manifest, null, 2) + '\n');
const summary = manifest.lessons.map(l => l.practice ? `${l.slug} (${l.practice.files.length} files)` : l.slug);
console.log(`[tutorials] Imported ${imported.length} lessons: ${summary.join(', ')}. Original sources unchanged.`);
