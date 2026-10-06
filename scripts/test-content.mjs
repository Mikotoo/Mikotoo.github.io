#!/usr/bin/env node
/** Isolated regression checks; never runs sync against this checkout's snapshots. */
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { cp, mkdir, mkdtemp, readFile, readdir, rm, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import path from 'node:path';
import { spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { site } from '../config/site.mjs';

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const temporary = await mkdtemp(path.join(tmpdir(), 'site-content-test-'));

async function put(file, content) {
  await mkdir(path.dirname(file), { recursive: true });
  await writeFile(file, content, 'utf8');
}

function run(script, args, cwd, expectedStatus = 0) {
  const result = spawnSync(process.execPath, [script, ...args], { cwd, stdio: 'inherit' });
  assert.ifError(result.error);
  assert.equal(result.signal, null, `${path.basename(script)} was interrupted`);
  assert.equal(result.status, expectedStatus, `${path.basename(script)} exit status`);
}

async function snapshotTree(dir) {
  const files = {};
  async function walk(current, prefix) {
    const entries = await readdir(current, { withFileTypes: true });
    entries.sort((a, b) => a.name < b.name ? -1 : a.name > b.name ? 1 : 0);
    for (const entry of entries) {
      const relative = prefix ? `${prefix}/${entry.name}` : entry.name;
      const absolute = path.join(current, entry.name);
      if (entry.isDirectory()) await walk(absolute, relative);
      else files[relative] = (await readFile(absolute)).toString('base64');
    }
  }
  await walk(dir, '');
  return files;
}

function sourceDigest(readmes) {
  const hash = createHash('sha256');
  for (const name of Object.keys(readmes).sort()) {
    const bytes = Buffer.from(readmes[name], 'utf8');
    hash.update(`${Buffer.byteLength(name, 'utf8')}:`);
    hash.update(name);
    hash.update(`${bytes.length}:`);
    hash.update(bytes);
  }
  return hash.digest('hex');
}

async function testSync() {
  console.log('[test] Isolated snapshot sync');
  const fixture = path.join(temporary, 'site');
  const source = path.join(temporary, 'source');
  const script = path.join(fixture, 'scripts', 'sync-docs.mjs');
  const registry = path.join(fixture, 'config', 'projects.json');
  const snapshots = path.join(fixture, 'snapshots');
  const manifestPath = path.join(snapshots, 'manifest.json');
  await mkdir(path.dirname(script), { recursive: true });
  await mkdir(path.dirname(registry), { recursive: true });
  await cp(path.join(ROOT, 'scripts', 'sync-docs.mjs'), script);
  await cp(path.join(ROOT, 'config', 'project-labels.json'), path.join(fixture, 'config', 'project-labels.json'));
  const project = { key: 'fixture', title: 'Fixture 项目', dir: source, href: 'https://github.com/fixture/source' };
  await put(registry, JSON.stringify({ projects: [project] }));
  const readmes = {
    'README.md': '# Fixture\n\n[Guide](guide/README.md)\n',
    'guide/README.md': '# Guide\n\n[Home](../README.md)\n\n指南。\n',
  };
  for (const [name, content] of Object.entries(readmes)) await put(path.join(source, name), content);
  const sourceBefore = await snapshotTree(source);

  run(script, [], fixture);
  const rootPage = await readFile(path.join(snapshots, 'projects', 'fixture.md'), 'utf8');
  const subPage = await readFile(path.join(snapshots, 'projects', 'fixture', 'guide.md'), 'utf8');
  assert.ok(rootPage.includes('[Guide](/_projects/fixture/guide/)'), 'root page retains the public URL prefix');
  assert.ok(subPage.includes('[Home](/_projects/fixture/)'), 'nested README links to the project root');
  const manifest = JSON.parse(await readFile(manifestPath, 'utf8'));
  assert.deepEqual(manifest, {
    projects: [{ key: 'fixture', label: 'Fixture 项目', pages: 2, rootSlug: '_projects/fixture', sourceSha256: sourceDigest(readmes) }],
  }, 'manifest contains deterministic metadata and the hash of README relative names plus bytes');
  const first = await snapshotTree(snapshots);
  run(script, [], fixture);
  assert.deepEqual(await snapshotTree(snapshots), first, 'identical sync produces byte-identical snapshots and manifest');
  assert.deepEqual(await snapshotTree(source), sourceBefore, 'sync never changes source files');

  readmes['guide/README.md'] += '\nChanged source content.\n';
  await put(path.join(source, 'guide', 'README.md'), readmes['guide/README.md']);
  run(script, [], fixture);
  const changedManifest = JSON.parse(await readFile(manifestPath, 'utf8'));
  assert.notEqual(changedManifest.projects[0].sourceSha256, manifest.projects[0].sourceSha256);
  assert.equal(changedManifest.projects[0].sourceSha256, sourceDigest(readmes));
  assert.ok((await readFile(path.join(snapshots, 'projects', 'fixture', 'guide.md'), 'utf8')).includes('Changed source content.'));

  const beforeFailure = await snapshotTree(snapshots);
  // Change the available project too: a partial run must not publish its new content.
  await put(path.join(source, 'README.md'), `${readmes['README.md']}\nMust not be published.\n`);
  await put(registry, JSON.stringify({ projects: [project, { key: 'missing', dir: path.join(temporary, 'unavailable-source') }] }));
  console.log('[test] Expected failure: unavailable second registered source');
  run(script, [], fixture, 1);
  assert.deepEqual(await snapshotTree(snapshots), beforeFailure, 'missing registered source leaves existing snapshots byte-identical');
}

async function testLinks() {
  console.log('[test] Link checker fixtures');
  const fixture = path.join(temporary, 'html');
  const script = path.join(ROOT, 'scripts', 'check-links.mjs');
  const validHome = '<a href="guide/">Relative guide</a><a href="/guide/">Root guide</a><a href="/standalone">Standalone page</a>';
  await put(path.join(fixture, 'index.html'), validHome);
  await put(path.join(fixture, 'guide', 'index.html'), [
    '<a href="../">Home</a>',
    '<img src="../assets/pixel.svg">',
    '<link rel="stylesheet" href="/assets/site.css">',
    `<a href="${site.url}/guide/?mode=test#section">Same origin</a>`,
  ].join('\n'));
  await put(path.join(fixture, 'standalone.html'), '<a href="/">Home</a>');
  await put(path.join(fixture, 'assets', 'pixel.svg'), '<svg xmlns="http://www.w3.org/2000/svg"/>');
  await put(path.join(fixture, 'assets', 'site.css'), 'body { color: black; }');
  await put(path.join(fixture, '404.html'), `<link rel="canonical" href="${site.url}/404/"><a href="/">Home</a>`);
  run(script, ['--dir', fixture], ROOT);

  await put(path.join(fixture, 'index.html'), `${validHome}<a href="missing-page/">Missing page</a>`);
  console.log('[test] Expected failure: missing navigation target');
  run(script, ['--dir', fixture], ROOT, 1);

  await put(path.join(fixture, 'index.html'), validHome);
  // Put the anchor in 404.html itself to catch an overbroad canonical exemption.
  await put(path.join(fixture, '404.html'), `<link rel="canonical" href="${site.url}/404/"><a href="/404/">Actual navigation</a>`);
  console.log('[test] Expected failure: /404/ anchor is not canonical metadata');
  run(script, ['--dir', fixture], ROOT, 1);
}

try {
  await testSync();
  await testLinks();
  console.log('[test] All content regression checks passed.');
} finally {
  await rm(temporary, { recursive: true, force: true });
}
