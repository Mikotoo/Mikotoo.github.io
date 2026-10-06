import { defineCollection } from 'astro:content';
import { docsLoader } from '@astrojs/starlight/loaders';
import { docsSchema } from '@astrojs/starlight/schema';

/**
 * Starlight 的 docs 集合。
 * 目录结构：
 *   src/content/docs/            手写页面（index / intro / about / code / resources）
 *   src/content/docs/notes/      由 _posts/*.md 构建前自动生成
 *   src/content/docs/lotus/      由 ../Lotus_genome 的 README 自动生成
 *   src/content/docs/rscript/    由 ../Rscript 的 README 自动生成
 */
export const collections = {
  docs: defineCollection({ loader: docsLoader(), schema: docsSchema() }),
};
