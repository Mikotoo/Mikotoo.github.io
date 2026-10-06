import { defineCollection } from 'astro:content';
import { docsLoader } from '@astrojs/starlight/loaders';
import { docsSchema } from '@astrojs/starlight/schema';

// pnpm gen 将 content/ 与 snapshots/ 装配到 src/content/docs/。
// 此目录仅供 Starlight 加载，不是手写内容源。
export const collections = {
  docs: defineCollection({ loader: docsLoader(), schema: docsSchema() }),
};
