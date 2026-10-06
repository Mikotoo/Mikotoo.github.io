/** 文章页的阅读导航；展示页通过顶部导航进入。 */
export function navigation(projects, notes, tutorials) {
  return [
    { label: '探索', items: [
      { label: '首页', slug: 'index' },
      { label: '全部内容', slug: 'index-all' },
    ] },
    { label: '生物信息学基础教程', items: [{ label: '课程目录', slug: 'tutorials' }, ...tutorials] },
    { label: '研究项目', items: [{ label: '项目概览', slug: 'projects' }, ...projects] },
    { label: '笔记与复现', items: [{ label: '全部笔记', slug: 'notes' }, ...notes] },
    { label: '代码与资源', items: [
      { label: '分析代码', slug: 'code' },
      { label: '学习资源', slug: 'resources' },
    ] },
    { label: '关于', items: [{ label: '关于 Mikotoo', slug: 'about' }] },
  ];
}
