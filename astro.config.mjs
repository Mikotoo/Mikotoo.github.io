import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';
import { site } from './config/site.mjs';
import { navigation } from './config/navigation.mjs';
import { projectsSidebar, notesSidebar, tutorialsSidebar } from './.generated/sidebar.mjs';

export default defineConfig({
  site: site.url,
  base: '/',
  trailingSlash: 'always',
  integrations: [starlight({
    title: site.title,
    description: site.description,
    defaultLocale: 'root',
    locales: { root: { label: '简体中文', lang: 'zh-CN' } },
    social: [{ icon: 'github', label: 'GitHub', href: site.github }],
    favicon: '/favicon.ico',
    customCss: ['./src/styles/custom.css'],
    components: { Header: './src/components/Header.astro' },
    credits: false,
    lastUpdated: false,
    pagination: true,
    tableOfContents: { minHeadingLevel: 2, maxHeadingLevel: 3 },
    head: [{ tag: 'meta', attrs: { name: 'author', content: site.author } }],
    sidebar: navigation(projectsSidebar, notesSidebar, tutorialsSidebar),
  })],
});
