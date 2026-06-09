// ============================================
// CP2K Tutorial Website — Application Logic
// ============================================

// ---- Chapter Data ----
const chapters = [
  {
    id: 1, folder: '01-hello-cp2k',
    title: '你的第一个CP2K计算', titleEn: 'First CP2K Calculation',
    diff: 1, level: '入门',
    desc: '了解CP2K的输入文件结构，运行你的第一个单原子能量计算。',
    tags: ['GLOBAL', 'FORCE_EVAL', 'KIND', 'RUN_TYPE'],
    file: 'README.md'
  },
  {
    id: 2, folder: '02-water-static',
    title: '水分子的静态DFT计算', titleEn: 'Static DFT Calculation of Water',
    diff: 1, level: '入门',
    desc: '使用DFT计算水分子的总能量和原子力，理解QUICKSTEP模块。',
    tags: ['DFT', 'SCF', 'XC_FUNCTIONAL', 'BASIS_SET'],
    file: 'README.md'
  },
  {
    id: 3, folder: '03-basis-sets',
    title: '基组与赝势详解', titleEn: 'Basis Sets and Pseudopotentials',
    diff: 2, level: '基础',
    desc: '深入了解高斯基组（SZV/DZVP/TZV2P）和GTH赝势的选择与收敛性。',
    tags: ['GTO', 'GTH', 'SZV', 'DZVP', 'TZV2P', 'Pseudopotential'],
    file: 'README.md'
  },
  {
    id: 4, folder: '04-cutoff-convergence',
    title: '收敛CUTOFF与REL_CUTOFF', titleEn: 'Converging CUTOFF & REL_CUTOFF',
    diff: 2, level: '基础',
    desc: '学习QUICKSTEP多网格系统，系统性地收敛平面波截断能。',
    tags: ['MGRID', 'CUTOFF', 'REL_CUTOFF', 'NGRIDS', 'Multi-grid'],
    file: 'README.md'
  },
  {
    id: 5, folder: '05-periodic-silicon',
    title: '周期性体系：硅晶体', titleEn: 'Periodic Systems: Bulk Silicon',
    diff: 2, level: '基础',
    desc: '学习周期性边界条件、晶胞定义、能带理论基础和展宽技术。',
    tags: ['PBC', 'CELL', 'KPOINTS', 'SMEAR', 'Band Structure'],
    file: 'README.md'
  },
  {
    id: 6, folder: '06-geo-opt',
    title: '几何优化', titleEn: 'Geometry Optimization',
    diff: 2, level: '基础',
    desc: '使用GEO_OPT找到分子的平衡构型，理解收敛判据和优化算法。',
    tags: ['GEO_OPT', 'CG', 'BFGS', 'FORCE', 'FIXED_ATOMS'],
    file: 'README.md'
  },
  {
    id: 7, folder: '07-cell-opt',
    title: '晶胞优化', titleEn: 'Cell Optimization',
    diff: 3, level: '中级',
    desc: '同时优化原子位置和晶胞参数，理解应力张量和压强收敛。',
    tags: ['CELL_OPT', 'Stress', 'Pressure', 'Virial'],
    file: 'README.md'
  },
  {
    id: 8, folder: '08-md-nve',
    title: 'NVE分子动力学', titleEn: 'NVE Molecular Dynamics',
    diff: 3, level: '中级',
    desc: '运行微正则系综分子动力学，检验能量守恒作为模拟质量测试。',
    tags: ['MD', 'NVE', 'TIMESTEP', 'EXTRAPOLATION', 'Trajectory'],
    file: 'README.md'
  },
  {
    id: 9, folder: '09-md-thermostats',
    title: 'NVT/NPT分子动力学与热浴', titleEn: 'NVT/NPT with Thermostats',
    diff: 3, level: '中级',
    desc: '使用Nosé-Hoover和CSVR热浴进行恒温/恒压分子动力学模拟。',
    tags: ['NVT', 'NPT', 'NOSE', 'CSVR', 'BAROSTAT', 'Langevin'],
    file: 'README.md'
  },
  {
    id: 10, folder: '10-metadynamics',
    title: '元动力学增强采样', titleEn: 'Metadynamics Enhanced Sampling',
    diff: 4, level: '高级',
    desc: '使用元动力学加速稀有事件采样，重建自由能面。',
    tags: ['METADYN', 'CV', 'HILLS', 'FES', 'Collective Variable'],
    file: 'README.md'
  },
  {
    id: 11, folder: '11-xtb',
    title: 'GFN1-xTB半经验方法', titleEn: 'GFN1-xTB Semi-empirical Method',
    diff: 3, level: '中级',
    desc: '使用xTB快速进行大体系计算，了解半经验方法的理论基础。',
    tags: ['xTB', 'Semi-empirical', 'GFN1', 'OT', 'D3'],
    file: 'README.md'
  },
  {
    id: 12, folder: '12-qmmm-basics',
    title: 'QM/MM混合方法基础', titleEn: 'QM/MM Hybrid Methods',
    diff: 4, level: '高级',
    desc: '将体系分为QM和MM区域，学习多尺度建模的基本概念。',
    tags: ['QMMM', 'QM_KIND', 'LINK', 'Force Field', 'Embedding'],
    file: 'README.md'
  },
  {
    id: 13, folder: '13-vibrational',
    title: '振动分析与红外光谱', titleEn: 'Vibrational Analysis & IR',
    diff: 3, level: '中级',
    desc: '计算分子振动频率和红外光谱，理解Hessian矩阵和简正模式。',
    tags: ['VIBRATIONAL_ANALYSIS', 'Hessian', 'IR', 'Normal Mode'],
    file: 'README.md'
  },
  {
    id: 14, folder: '14-tddft',
    title: 'TDDFT激发态计算', titleEn: 'TDDFT Excited States',
    diff: 4, level: '高级',
    desc: '使用线性响应TDDFT计算激发能和吸收光谱。',
    tags: ['TDDFT', 'Casida', 'Excitation', 'Oscillator Strength'],
    file: 'README.md'
  },
  {
    id: 15, folder: '15-mp2',
    title: '后Hartree-Fock方法：RI-MP2', titleEn: 'Post-HF: RI-MP2',
    diff: 5, level: '专家',
    desc: '使用RI-MP2计算关联能，了解后Hartree-Fock方法的理论框架。',
    tags: ['MP2', 'RI', 'WF_CORRELATION', 'Post-HF', 'Correlation'],
    file: 'README.md'
  },
  {
    id: 16, folder: '16-gw-bse',
    title: 'GW/BSE能带结构计算', titleEn: 'GW/BSE Band Structure',
    diff: 5, level: '专家',
    desc: '使用GW近似计算准粒子能级和带隙，BSE计算激子效应。',
    tags: ['GW', 'BSE', 'Quasiparticle', 'Band Gap', 'Self-energy'],
    file: 'README.md'
  }
];

// ---- Glossary Data ----
const glossary = [
  { term: 'DFT', def: '密度泛函理论 (Density Functional Theory)。基于电子密度而非波函数来求解多电子体系薛定谔方程的理论框架。' },
  { term: 'GPW', def: '高斯和平面波方法 (Gaussian and Plane Waves)。CP2K QUICKSTEP模块的核心方法，用高斯基组表示轨道，平面波/实空间网格表示密度。' },
  { term: 'GAPW', def: '高斯增广平面波方法 (Gaussian Augmented Plane Waves)。GPW的扩展，支持全电子计算和核心性质。' },
  { term: 'QUICKSTEP', def: 'CP2K的DFT计算模块，实现了GPW和GAPW方法。' },
  { term: 'GTO', def: '高斯型轨道 (Gaussian Type Orbital)。形如 exp(-αr²) 的基函数，用于展开分子轨道。' },
  { term: 'GTH', def: 'Goedecker-Teter-Hutter赝势。CP2K中使用的模守恒赝势族。' },
  { term: 'Basis Set (基组)', def: '用于展开分子轨道的基函数集合。CP2K中常用SZV(单ζ)、DZVP(双ζ极化)、TZV2P(三ζ双极化)等。' },
  { term: 'Pseudopotential (赝势)', def: '替代核心电子和核的等效势，只处理价电子，大幅降低计算成本。' },
  { term: 'SCF', def: '自洽场 (Self-Consistent Field)。迭代求解Kohn-Sham方程直到电子密度自洽的过程。' },
  { term: 'XC Functional', def: '交换关联泛函。DFT中近似电子交换和关联效应的泛函，如LDA(PADE)、GGA(PBE)、杂化泛函(B3LYP)等。' },
  { term: 'CUTOFF', def: '平面波截断能(Ry)。控制实空间网格的精细程度，值越大网格越细、计算越精确但越慢。' },
  { term: 'REL_CUTOFF', def: '相对截断能。控制高斯函数到多网格层级的映射，决定每个高斯函数被分配到哪个网格层级。' },
  { term: 'MGRID', def: '多网格(Multi-grid)设置。QUICKSTEP使用多层级网格系统，宽高斯映射到粗网格，窄高斯映射到细网格。' },
  { term: 'PBC', def: '周期性边界条件 (Periodic Boundary Conditions)。模拟无限周期性体系的边界处理方法。' },
  { term: 'K-points (K点)', def: '布里渊区中的采样点。用于周期性体系的电子结构计算，Γ点仅采样布里渊区中心。' },
  { term: 'Smearing (展宽)', def: '电子占据数的平滑处理方法(Fermi-Dirac等)，用于改善金属体系的SCF收敛。' },
  { term: 'GEO_OPT', def: '几何优化。通过迭代调整原子位置找到势能面上的能量极小点(平衡构型)。' },
  { term: 'CELL_OPT', def: '晶胞优化。同时优化原子位置和晶胞参数(晶格常数)。' },
  { term: 'MD', def: '分子动力学 (Molecular Dynamics)。通过求解牛顿运动方程模拟原子随时间的运动。' },
  { term: 'NVE', def: '微正则系综。粒子数N、体积V、能量E恒定。用于检验能量守恒。' },
  { term: 'NVT', def: '正则系综。粒子数N、体积V、温度T恒定。需要热浴(thermostat)控制温度。' },
  { term: 'NPT', def: '等温等压系综。粒子数N、压强P、温度T恒定。需要热浴和压浴。' },
  { term: 'Nosé-Hoover', def: '一种扩展拉格朗日热浴方法，通过引入额外自由度来控制温度。' },
  { term: 'CSVR', def: 'Canonical Sampling through Velocity Rescaling。一种高效的恒温方法，比Nosé-Hoover更快达到平衡。' },
  { term: 'Metadynamics (元动力学)', def: '一种增强采样方法，通过在集合变量空间沉积高斯偏置势来加速稀有事件。' },
  { term: 'Collective Variable (CV)', def: '集合变量。描述体系集体运动的低维变量，如配位数、键长、二面角等。' },
  { term: 'xTB', def: '扩展紧束缚方法 (eXtended Tight Binding)。GFN1-xTB是一种半经验量子化学方法，速度远快于DFT。' },
  { term: 'QM/MM', def: '量子力学/分子力学混合方法。将体系分为QM区域(精确处理)和MM区域(力场处理)。' },
  { term: 'Hessian', def: 'Hessian矩阵。能量对原子坐标的二阶导数矩阵，用于振动分析和频率计算。' },
  { term: 'TDDFT', def: '含时密度泛函理论 (Time-Dependent DFT)。用于计算激发态性质(激发能、吸收光谱)。' },
  { term: 'MP2', def: '二阶Møller-Plesset微扰理论。最简单的后Hartree-Fock关联方法，O(N⁵)标度。' },
  { term: 'RI', def: 'Resolution of Identity。用辅助基组展开四中心积分的近似方法，大幅降低计算成本。' },
  { term: 'GW', def: '一种准粒子方法，通过自能Σ=iGW修正DFT本征值，给出准确的带隙和能带结构。' },
  { term: 'BSE', def: 'Bethe-Salpeter方程。用于计算激子(电子-空穴对)效应，给出准确的光学吸收谱。' },
  { term: 'OT', def: 'Orbital Transformation。CP2K中的一种SCF求解方法，直接最小化能量而非对角化。' },
  { term: 'Broyden Mixing', def: '一种电荷密度混合方法，用于加速SCF收敛。通过历史信息自适应调整混合参数。' },
  { term: 'Pulay Mixing', def: 'DIIS混合方法。利用多步历史信息进行外推，加速SCF收敛。' },
  { term: 'BSSE', def: '基组重叠误差 (Basis Set Superposition Error)。复合物能量因基组重叠而被人为降低的误差。' },
  { term: 'Band Gap (带隙)', def: '价带顶到导带底的能量差。区分导体、半导体和绝缘体的关键参数。' },
  { term: 'FES', def: '自由能面 (Free Energy Surface)。元动力学的目标输出，描述体系在集合变量空间的自由能分布。' },
];

// ---- State ----
let currentChapter = null;
let completed = JSON.parse(localStorage.getItem('cp2k_completed') || '{}');
let theme = localStorage.getItem('cp2k_theme') || 'light';

// ---- Theme ----
function initTheme() {
  document.documentElement.setAttribute('data-theme', theme);
  document.querySelector('.theme-toggle').textContent = theme === 'dark' ? '☀️' : '🌙';
  document.getElementById('hljs-light').disabled = theme === 'dark';
  document.getElementById('hljs-dark').disabled = theme === 'light';
}

function toggleTheme() {
  theme = theme === 'dark' ? 'light' : 'dark';
  localStorage.setItem('cp2k_theme', theme);
  initTheme();
}

// ---- Navigation ----
function showHome() {
  document.getElementById('home-view').classList.remove('hidden');
  document.getElementById('chapter-view').classList.add('hidden');
  document.getElementById('glossary-view').classList.add('hidden');
  document.getElementById('progress-view').classList.add('hidden');
  currentChapter = null;
  renderChapterCards();
  window.scrollTo(0, 0);
}

function showChapter(id) {
  const ch = chapters.find(c => c.id === id);
  if (!ch) return;
  currentChapter = id;

  document.getElementById('home-view').classList.add('hidden');
  document.getElementById('chapter-view').classList.remove('hidden');
  document.getElementById('glossary-view').classList.add('hidden');
  document.getElementById('progress-view').classList.add('hidden');

  document.getElementById('chapter-indicator').textContent = `第 ${ch.id} 章 / 共 ${chapters.length} 章`;
  document.getElementById('prev-chapter').disabled = id <= 1;
  document.getElementById('next-chapter').disabled = id >= chapters.length;
  document.getElementById('prev-chapter-bottom').disabled = id <= 1;
  document.getElementById('next-chapter-bottom').disabled = id >= chapters.length;

  updateCompleteButton();
  loadChapterContent(ch);
  window.scrollTo(0, 0);
}

function showGlossary() {
  document.getElementById('home-view').classList.add('hidden');
  document.getElementById('chapter-view').classList.add('hidden');
  document.getElementById('glossary-view').classList.remove('hidden');
  document.getElementById('progress-view').classList.add('hidden');
  renderGlossary();
  window.scrollTo(0, 0);
}

function showProgress() {
  document.getElementById('home-view').classList.add('hidden');
  document.getElementById('chapter-view').classList.add('hidden');
  document.getElementById('glossary-view').classList.add('hidden');
  document.getElementById('progress-view').classList.remove('hidden');
  renderProgress();
  window.scrollTo(0, 0);
}

function navigateChapter(delta) {
  if (!currentChapter) return;
  const next = currentChapter + delta;
  if (next >= 1 && next <= chapters.length) {
    showChapter(next);
  }
}

// ---- Chapter Cards ----
function renderChapterCards() {
  const grid = document.getElementById('chapters-grid');
  grid.innerHTML = chapters.map(ch => {
    const stars = '⭐'.repeat(ch.diff);
    const isCompleted = completed[ch.id];
    return `
      <div class="chapter-card ${isCompleted ? 'completed' : ''}" onclick="showChapter(${ch.id})">
        <div class="card-header">
          <div class="card-num">${String(ch.id).padStart(2, '0')}</div>
          <div>
            <div class="card-title">${ch.title}</div>
            <div class="card-title-en">${ch.titleEn}</div>
          </div>
        </div>
        <div class="card-difficulty">
          <span class="card-stars">${stars}</span>
          <span class="card-level-text">${ch.level}</span>
        </div>
        <div class="card-desc">${ch.desc}</div>
        <div class="card-tags">
          ${ch.tags.map(t => `<span class="tag">${t}</span>`).join('')}
        </div>
      </div>
    `;
  }).join('');
}

function filterChapters(level) {
  document.querySelectorAll('.filter-btn').forEach(btn => btn.classList.remove('active'));
  event.target.classList.add('active');

  const cards = document.querySelectorAll('.chapter-card');
  cards.forEach((card, i) => {
    const ch = chapters[i];
    if (level === 'all' || ch.diff === level) {
      card.style.display = '';
    } else {
      card.style.display = 'none';
    }
  });
}

// ---- Chapter Content ----
async function loadChapterContent(ch) {
  const contentEl = document.getElementById('chapter-content');
  contentEl.innerHTML = '<p style="text-align:center; color:var(--text-muted); padding:60px 0;">加载中...</p>';

  try {
    const url = `../${ch.folder}/${ch.file}`;
    const resp = await fetch(url);
    if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
    let md = await resp.text();

    // Configure marked
    marked.setOptions({
      highlight: function(code, lang) {
        if (lang && hljs.getLanguage(lang)) {
          return hljs.highlight(code, { language: lang }).value;
        }
        return hljs.highlightAuto(code).value;
      },
      breaks: false,
      gfm: true
    });

    let html = marked.parse(md);

    // Add copy buttons to code blocks
    html = html.replace(/<pre><code(.*?)>([\s\S]*?)<\/code><\/pre>/g, (match, attrs, code) => {
      return `<div class="code-block-wrapper"><button class="copy-btn" onclick="copyCode(this)">复制</button><pre><code${attrs}>${code}</code></pre></div>`;
    });

    contentEl.innerHTML = html;

    // Add glossary tooltips to key terms
    addGlossaryTooltips(contentEl);

  } catch (err) {
    contentEl.innerHTML = `
      <div style="text-align:center; padding:60px 0; color:var(--text-muted);">
        <p style="font-size:2rem; margin-bottom:16px;">📄</p>
        <p>章节内容加载失败</p>
        <p style="font-size:0.85rem; margin-top:8px;">请确保已生成章节文件: ${ch.folder}/${ch.file}</p>
        <p style="font-size:0.8rem; margin-top:4px; color:var(--danger);">${err.message}</p>
      </div>
    `;
  }
}

// ---- Copy Code ----
function copyCode(btn) {
  const code = btn.parentElement.querySelector('code');
  const text = code.textContent;
  navigator.clipboard.writeText(text).then(() => {
    btn.textContent = '已复制 ✓';
    btn.classList.add('copied');
    setTimeout(() => {
      btn.textContent = '复制';
      btn.classList.remove('copied');
    }, 2000);
  });
}

// ---- Progress ----
function toggleComplete() {
  if (!currentChapter) return;
  if (completed[currentChapter]) {
    delete completed[currentChapter];
  } else {
    completed[currentChapter] = Date.now();
  }
  localStorage.setItem('cp2k_completed', JSON.stringify(completed));
  updateCompleteButton();
  renderChapterCards();
}

function updateCompleteButton() {
  const btn = document.getElementById('complete-btn');
  if (completed[currentChapter]) {
    btn.textContent = '✓ 已完成';
    btn.classList.add('done');
  } else {
    btn.textContent = '✓ 标记为已完成';
    btn.classList.remove('done');
  }
}

function renderProgress() {
  const done = Object.keys(completed).length;
  const total = chapters.length;
  const pct = Math.round((done / total) * 100);

  const statsEl = document.getElementById('progress-stats');
  statsEl.innerHTML = `
    <div class="progress-stat-card">
      <div class="progress-stat-num">${done}/${total}</div>
      <div class="progress-stat-label">已完成章节</div>
    </div>
    <div class="progress-stat-card">
      <div class="progress-stat-num">${pct}%</div>
      <div class="progress-stat-label">完成进度</div>
    </div>
    <div class="progress-stat-card">
      <div class="progress-stat-num">${total - done}</div>
      <div class="progress-stat-label">剩余章节</div>
    </div>
  `;

  // Progress bar
  const barContainer = document.querySelector('.progress-bar-container');
  if (!barContainer) {
    const bar = document.createElement('div');
    bar.className = 'progress-bar-container';
    bar.innerHTML = `<div class="progress-bar-bg"><div class="progress-bar-fill" style="width:${pct}%"></div></div>`;
    statsEl.after(bar);
  } else {
    barContainer.querySelector('.progress-bar-fill').style.width = `${pct}%`;
  }

  const listEl = document.getElementById('progress-list');
  listEl.innerHTML = chapters.map(ch => {
    const isDone = completed[ch.id];
    const stars = '⭐'.repeat(ch.diff);
    return `
      <div class="progress-item">
        <div class="progress-check ${isDone ? 'done' : ''}" onclick="toggleProgressItem(${ch.id})">
          ${isDone ? '✓' : ''}
        </div>
        <div class="progress-item-title ${isDone ? 'done' : ''}" onclick="showChapter(${ch.id})" style="cursor:pointer;">
          ${ch.id}. ${ch.title}
        </div>
        <div class="progress-item-diff">${stars}</div>
      </div>
    `;
  }).join('');
}

function toggleProgressItem(id) {
  if (completed[id]) {
    delete completed[id];
  } else {
    completed[id] = Date.now();
  }
  localStorage.setItem('cp2k_completed', JSON.stringify(completed));
  renderProgress();
  renderChapterCards();
}

// ---- Glossary ----
function renderGlossary(filter = '') {
  const grid = document.getElementById('glossary-grid');
  const filtered = filter
    ? glossary.filter(g => g.term.toLowerCase().includes(filter.toLowerCase()) || g.def.includes(filter))
    : glossary;

  grid.innerHTML = filtered.map(g => `
    <div class="glossary-item">
      <div class="glossary-term">${g.term}</div>
      <div class="glossary-def">${g.def}</div>
    </div>
  `).join('');
}

function searchGlossary(query) {
  renderGlossary(query);
}

function addGlossaryTooltips(container) {
  const tooltip = document.getElementById('glossary-tooltip');

  glossary.forEach(g => {
    const termKey = g.term.split('(')[0].trim().split(' ')[0];
    if (!termKey || termKey.length < 2) return;

    const walker = document.createTreeWalker(container, NodeFilter.SHOW_TEXT, null, false);
    const textNodes = [];
    while (walker.nextNode()) textNodes.push(walker.currentNode);

    textNodes.forEach(node => {
      const parent = node.parentNode;
      if (parent.tagName === 'CODE' || parent.tagName === 'PRE' || parent.tagName === 'A') return;

      const regex = new RegExp(`(?<![\\w-])${termKey.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')}(?![\\w-])`, 'g');
      if (!regex.test(node.textContent)) return;

      const span = document.createElement('span');
      span.innerHTML = node.textContent.replace(regex, match =>
        `<span class="glossary-keyword" data-term="${g.term}" data-def="${g.def}" style="border-bottom:1px dashed var(--accent); cursor:help;">${match}</span>`
      );
      parent.replaceChild(span, node);
    });
  });

  // Tooltip events
  container.querySelectorAll('.glossary-keyword').forEach(el => {
    el.addEventListener('mouseenter', (e) => {
      const tooltip = document.getElementById('glossary-tooltip');
      tooltip.innerHTML = `<div class="tooltip-term">${el.dataset.term}</div><div>${el.dataset.def}</div>`;
      tooltip.classList.add('visible');
      const rect = el.getBoundingClientRect();
      tooltip.style.left = Math.min(rect.left, window.innerWidth - 340) + 'px';
      tooltip.style.top = (rect.bottom + 8) + 'px';
    });
    el.addEventListener('mouseleave', () => {
      document.getElementById('glossary-tooltip').classList.remove('visible');
    });
  });
}

// ---- Init ----
function init() {
  initTheme();
  renderChapterCards();
}

document.addEventListener('DOMContentLoaded', init);
