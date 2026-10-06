// @ts-check
import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';
import { site } from './src/site.config.mjs';
import { projectsSidebar, notesSidebar } from './src/sidebar.generated.mjs';

// 用户站点仓库（Mikotoo.github.io）→ 站点根路径为 https://mikotoo.github.io/
export default defineConfig({
  site: 'https://mikotoo.github.io',
  base: '/',
  trailingSlash: 'always',
  integrations: [
    starlight({
      title: site.title,
      description: site.description,
      defaultLocale: 'root',
      locales: {
        root: { label: '简体中文', lang: 'zh-CN' },
      },
      social: [{ icon: 'github', label: 'GitHub', href: site.github }],
      customCss: ['./src/styles/custom.css'],
      credits: false,
      lastUpdated: false,
      pagination: true,
      tableOfContents: { minHeadingLevel: 2, maxHeadingLevel: 3 },
      head: [{ tag: 'meta', attrs: { name: 'author', content: site.author } }],
      sidebar: [
        {
          label: '开始',
          items: [
            { label: '首页', slug: 'index' },
            { label: '站点说明', slug: 'intro' },
            { label: '内容索引', slug: 'index-all' },
          ],
        },
        // 以下两块由 scripts/gen-index.mjs 生成，新增项目/笔记都不需要改本文件
        {
          label: '项目文档',
          items: [{ label: '全部项目', slug: 'projects' }, ...projectsSidebar],
        },
        {
          label: '笔记与复现',
          items: notesSidebar,
        },
        {
          label: '代码与资源',
          items: [
            { label: '代码与资源', slug: 'code' },
            { label: '下载资源', slug: 'resources' },
          ],
        },
        {
          label: '关于',
          items: [{ label: '关于与联系', slug: 'about' }],
        },
      ],
    }),
  ],
});
