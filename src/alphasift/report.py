"""
The HTML report: one self-contained page per screen that mirrors the companion
website's protein page, drawn in the browser from the `alphasift/1` document
embedded in it. No network, no build step: open the file.

What it reproduces (website src/components, 2026-10-02): the header with its count
tiles and parameter line; the heatmap (RMSF and recurrence tracks, binder rows in
display order, contacted residues as columns, one interface lane, hover read-out);
the page key (Zappo+H amino-acid colours, motif marks, partner tiers, engine
parameters); one card per interface with the motif engine's bait-side view, motif
lanes, relations line and a panel per motif (partner-sequence logo, carriers by
paralog group with reading direction, elements and bipartite pairs), then the
partner list (H and M shown, L on request, motif badges, motif residues
underlined); and the searchable list of all binders. Selecting an interface,
motif, carrier or partner lights it across the page, as on the website.

What it does not: the 3-D viewer. Its place holds the current selection and the
ChimeraX command that shows it (alphasift.cxc names every interface and motif).
"""

from __future__ import annotations

import json

CSS = r"""
:root{--ground:#fbfbfa;--panel:#fff;--panel-2:#f3f4f6;--plate:#fff;--ink:#191c20;--ink-2:#4a5058;--ink-3:#7c828b;
--line:#e2e4e8;--line-2:#d0d3d9;--sel:#111418;--sel-fill:rgba(17,20,24,.06);--gold:#ffc400;
--mono:ui-monospace,"SF Mono",Menlo,"Cascadia Mono","Roboto Mono",monospace;
--sans:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,"Helvetica Neue",Arial,sans-serif}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--ground:#0e1116;--panel:#161a20;--panel-2:#1d222a;
--ink:#e7eaee;--ink-2:#a8afb9;--ink-3:#767d87;--line:#272d36;--line-2:#333a45;--sel:#e7eaee;--sel-fill:rgba(231,234,238,.09)}}
:root[data-theme="dark"]{--ground:#0e1116;--panel:#161a20;--panel-2:#1d222a;--ink:#e7eaee;--ink-2:#a8afb9;--ink-3:#767d87;
--line:#272d36;--line-2:#333a45;--sel:#e7eaee;--sel-fill:rgba(231,234,238,.09)}
*{box-sizing:border-box}html,body{margin:0;padding:0}
body{background:var(--ground);color:var(--ink);font-family:var(--sans);font-size:15px;line-height:1.5;-webkit-font-smoothing:antialiased}
.mono{font-family:var(--mono);font-variant-ligatures:none}.num{font-variant-numeric:tabular-nums}
.eyebrow{font-size:11px;font-weight:600;letter-spacing:.09em;text-transform:uppercase;color:var(--ink-3)}
button{font:inherit;color:inherit}
.wrap{max-width:1280px;margin:0 auto;padding:26px 22px 80px}
.phead{display:flex;align-items:flex-start;gap:26px;flex-wrap:wrap;padding-bottom:18px;border-bottom:1px solid var(--line-2)}
.gene{font-size:34px;font-weight:680;margin:0}.idline{display:flex;gap:10px;color:var(--ink-2);font-size:13.5px;margin-top:6px;flex-wrap:wrap}
.params{margin-top:12px;font-size:12px;color:var(--ink-3);max-width:70ch}.params .mono{color:var(--ink-2)}
.tiles{display:flex;gap:10px;margin-left:auto;flex-wrap:wrap}
.tile{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:11px 15px;min-width:96px}
.tileV{font-size:23px;font-weight:680;letter-spacing:-.02em}.tileK{font-size:11px;color:var(--ink-3);text-transform:uppercase;letter-spacing:.05em}
.grid{display:grid;grid-template-columns:300px minmax(0,1fr);gap:22px;margin-top:22px;align-items:start}
.side{position:sticky;top:16px}@media(max-width:1040px){.grid{grid-template-columns:1fr}.side{position:static}}
.card{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:16px 16px 14px;margin-bottom:18px}
.ctitle{font-size:15px;font-weight:640}.csub{font-size:12px;color:var(--ink-3);margin:2px 0 10px;max-width:68ch;line-height:1.45}
.plate{background:var(--plate);color:#191c20;border-radius:8px;padding:10px;overflow-x:auto}
.selbox{font-size:13px;color:var(--ink-2)}.selbox b{color:var(--ink)}.selbox code{font-family:var(--mono);font-size:12px;background:var(--panel-2);padding:2px 5px;border-radius:4px;color:var(--ink)}
.clear{border:1px solid var(--line-2);background:var(--panel);border-radius:6px;padding:3px 10px;font-size:12px;cursor:pointer}
.clear:disabled{opacity:.4;cursor:default}
.hm{position:relative}.hm canvas{display:block;width:100%}
.lane{position:relative;height:21px;margin-top:2px}.laneLabel{position:absolute;font:10px var(--mono);color:#7c828b;left:0;top:4px}
.band{position:absolute;top:4px;height:12px;cursor:pointer;border-radius:2px;min-width:2px}
.band[data-on="1"]{outline:2px solid #111418;outline-offset:1px}
.tip{position:fixed;z-index:9;pointer-events:none;background:#111418;color:#fff;font-size:12px;padding:5px 8px;border-radius:6px;max-width:360px;white-space:pre-line}
.legends{display:flex;flex-wrap:wrap;gap:16px;font-size:12px;color:var(--ink-2);margin-top:8px;align-items:center}
.ramp{display:inline-block;width:70px;height:9px;border-radius:2px;vertical-align:middle;margin:0 4px}
.sw{display:inline-block;width:14px;height:9px;vertical-align:middle;margin-right:4px;border-radius:2px}
.key{font-size:12px;color:#4a5058;display:flex;flex-wrap:wrap;gap:6px 16px;margin-top:6px}
.aa{font-weight:800;font-family:var(--mono)}
.icard{border:1px solid var(--line);border-radius:10px;margin-bottom:10px;background:var(--panel)}
.ihead{display:flex;flex-wrap:wrap;gap:10px;align-items:center;padding:10px 12px;cursor:pointer}
.ihead[data-on="1"]{background:var(--sel-fill)}
.iid{font-weight:700}.ns{font-size:11px;border:1px solid var(--line-2);border-radius:4px;padding:0 5px;color:var(--ink-3)}
.seg{font-size:12px;color:var(--ink-2)}.imeta{font-size:12px;color:var(--ink-3);margin-left:auto}
.ibody{padding:0 12px 12px}
.mhead{display:flex;flex-wrap:wrap;gap:8px;align-items:baseline;margin:2px 0 6px}
.hbtn{background:none;border:0;padding:0;font-weight:650;font-size:13px;cursor:pointer;color:#191c20}
.hmeta{font-size:11.5px;color:#7c828b}
.bs{display:flex;gap:6px}.gut{position:relative;flex:0 0 112px}.gl{position:absolute;right:4px;font-size:10.5px;color:#7c828b;white-space:nowrap}
.laneBtn{position:absolute;right:4px;font:600 11px var(--mono);background:#eef0f3;border:1px solid #d0d3d9;border-radius:4px;padding:0 5px;cursor:pointer;color:#191c20}
.laneBtn[aria-pressed="true"]{background:#ffc400;border-color:#111418}
.scroll{overflow-x:auto}svg text{font-family:var(--sans)}
.summary{font-size:12.5px;color:#4a5058;margin:6px 0 2px}
details.motif{border-top:1px solid #e2e4e8;margin-top:8px;padding-top:6px}
details.motif summary{cursor:pointer;font-size:12.5px;color:#191c20}
.mmeta{color:#4a5058;margin-left:6px}.badge{font-size:10.5px;border:1px solid #d0d3d9;border-radius:4px;padding:0 5px;margin-left:6px;color:#4a5058}
.mnote{font-size:12px;color:#7c828b;margin:6px 0}
.chips{display:flex;flex-wrap:wrap;gap:6px;margin:6px 0}.pgroup{display:inline-flex;gap:4px;padding:2px;border:1px dashed #b9bfc7;border-radius:6px}
.chip{font-size:12px;border:1px solid #d0d3d9;background:#fff;border-radius:5px;padding:1px 7px;cursor:pointer;color:#191c20}
.chip[aria-pressed="true"]{background:#ffc400;border-color:#111418}.cun{color:#7c828b}
.elToggle{background:none;border:0;padding:0;font-size:12px;color:#4a5058;cursor:pointer;margin-top:4px}
table.el{border-collapse:collapse;font-size:12px;margin:4px 0;color:#191c20}table.el th,table.el td{padding:2px 8px;text-align:right;border-bottom:1px solid #eceef1}
table.el th{color:#7c828b;font-weight:500}table.el .l{text-align:left}
.pWrap{margin-top:10px}.pHead{display:flex;flex-wrap:wrap;gap:12px;align-items:center}
.pListToggle{background:none;border:0;padding:0;font-size:12.5px;font-weight:600;cursor:pointer;color:var(--ink-2)}
.pToggle{font-size:12px;color:var(--ink-3)}
.pRow{border-top:1px solid var(--line);padding:6px 4px;cursor:pointer}.pRow[data-on="1"]{background:var(--sel-fill)}
.pTop{display:flex;flex-wrap:wrap;gap:6px 14px;align-items:baseline}.pName{font-weight:650;font-size:13px}.pSub{font-size:11.5px;color:var(--ink-3)}
.pStats{display:flex;gap:10px;font-size:11.5px;color:var(--ink-3);margin-left:auto}
.mBadge{font:600 10.5px var(--mono);border:1px solid #d0d3d9;background:#eef0f3;border-radius:4px;padding:0 5px;cursor:pointer;color:#191c20}
.mBadge[aria-pressed="true"]{background:#ffc400;border-color:#111418}
.pRanges{display:flex;flex-wrap:wrap;gap:4px;margin-top:4px}
.range,.tier{display:inline-block;border-radius:4px;font-size:11px;line-height:1.5;padding:0 5px;border:1px solid transparent;font-family:var(--mono)}
.tierH{background:#25603c;border-color:#25603c;color:#fff;font-weight:650}.tierM{background:#adc5b4;border-color:#adc5b4;color:#191c20}
.tierL{background:transparent;border:1px dashed var(--ink-3);color:var(--ink-2)}.rPos{opacity:.75}
.mres{text-decoration:underline;text-decoration-thickness:2px;text-underline-offset:2px;font-weight:800}
.allToggle{background:none;border:0;font-size:12px;color:var(--ink-2);cursor:pointer;padding:6px 0}
.controls{display:flex;gap:10px;align-items:center;flex-wrap:wrap;margin-bottom:8px}
.controls input{font:inherit;font-size:13px;padding:4px 8px;border:1px solid var(--line-2);border-radius:6px;background:var(--panel);color:var(--ink)}
.sortBtn{font-size:12px;border:1px solid var(--line-2);background:var(--panel);border-radius:5px;padding:2px 8px;cursor:pointer}
.sortBtn[data-on="1"]{background:var(--sel-fill);font-weight:650}
.blist{max-height:420px;overflow:auto;border-top:1px solid var(--line)}
.brow{display:flex;justify-content:space-between;width:100%;border:0;border-bottom:1px solid var(--line);background:none;padding:5px 4px;cursor:pointer;text-align:left}
.brow[data-on="1"]{background:var(--sel-fill)}.bc{font-size:12px;color:var(--ink-3)}
.foot{font-size:11.5px;color:var(--ink-3);margin-top:24px}
"""

JS = r"""
const D = JSON.parse(document.getElementById('alphasift-data').textContent);
const P = D.params, SEQ = D.bait.sequence || '';
const IFC = ['#2a78d6','#d06f00','#d4609a','#a0662f','#1a9fc0','#d4462c','#7b61d9','#c2566b'];
const AAG = [['aliphatic / hydrophobic','ILVAM','#242c35'],['aromatic','FWY','#a169eb'],['positive','KR','#2d3ccb'],
  ['histidine','H','#0e9f93'],['negative','DE','#ec5560'],['hydrophilic','STNQ','#335309'],
  ['conformationally special','PG','#904680'],['cysteine','C','#847f13']];
const AAC = {}; AAG.forEach(([,l,c]) => [...l].forEach(a => AAC[a] = c));
const aaColor = l => AAC[(l||'').toUpperCase()] || '#8a9099';
const Q_SIG = 0.05;
const el = (tag, attrs = {}, ...kids) => { const e = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) { if (k === 'on') for (const [ev, f] of Object.entries(v)) e.addEventListener(ev, f);
    else if (k === 'style') Object.assign(e.style, v); else if (v !== null && v !== undefined && v !== false) e.setAttribute(k, v === true ? '' : v); }
  for (const k of kids.flat()) if (k !== null && k !== undefined && k !== false) e.append(k.nodeType ? k : String(k)); return e; };
const NS = 'http://www.w3.org/2000/svg';
const sv = (tag, attrs = {}, ...kids) => { const e = document.createElementNS(NS, tag);
  for (const [k, v] of Object.entries(attrs)) if (v !== null && v !== undefined) e.setAttribute(k, v);
  for (const k of kids.flat()) if (k !== null && k !== undefined) e.append(k.nodeType ? k : document.createTextNode(String(k))); return e; };
const fmtQ = q => q < 0.001 ? '< 0.001' : q.toPrecision(2);
const compact = rs => { rs = [...rs].sort((a,b)=>a-b); const o = []; for (let i = 0; i < rs.length;) { let j = i;
  while (j + 1 < rs.length && rs[j+1] === rs[j] + 1) j++; o.push(j > i ? `${rs[i]}–${rs[j]}` : String(rs[i])); i = j + 1; } return o.join(', '); };
const PRED = {}; D.predictions.forEach(p => PRED[p.id] = p);
const idn = id => (PRED[id] && PRED[id].identity) || null;
const nameOf = id => { const i = idn(id); return (i && (i.gene || (i.entry_name || '').split('_')[0])) || id; };
const subOf = id => { const i = idn(id); return i ? [i.uniprot, i.gene ? null : i.description].filter(Boolean).join(' · ') : ''; };
const paeAt = s => { const v = P.d0_signal * Math.sqrt(1 / s - 1); return Math.abs(v - Math.round(v)) < 0.05 ? String(Math.round(v)) : v.toFixed(1); };
const tierLabel = { H: `high: PAE ≤ ${paeAt(P.partner_tau_high)} Å`, M: `medium: PAE ${paeAt(P.partner_tau_high)}–${paeAt(P.partner_tau_in)} Å`, L: `low: PAE > ${paeAt(P.partner_tau_in)} Å` };
const ifaceFill = (i, sig) => { const c = IFC[i % 8]; return sig ? { background: c } :
  { background: `repeating-linear-gradient(135deg, ${c} 0 2px, ${c}38 2px 5px)`, border: `1px solid ${c}` }; };

// ------------------------------------------------------------------- selection state
const S = { iface: null, motif: null, carrier: null, partner: null };
const listeners = [];
const setSel = patch => { Object.assign(S, patch); listeners.forEach(f => f()); };
const onSel = f => { listeners.push(f); };
function litResidues() {
  if (S.iface === null) return new Set();
  const f = D.interfaces[S.iface];
  if (S.motif) { const m = f.motif.motifs.find(x => x.id === S.motif); if (m) return new Set(m.positions); }
  return new Set(f.residues);
}

// ------------------------------------------------------------------- header
function header() {
  const b = D.bait.identity || {};
  const nSig = D.sites.filter(s => s.significant).length;
  const nMot = D.interfaces.reduce((a, i) => a + ((i.motif && i.motif.motifs.length) || 0), 0);
  const binders = D.matrix.row_order.length;
  const tile = (v, k) => el('div', { class: 'tile' }, el('div', { class: 'tileV num' }, v), el('div', { class: 'tileK' }, k));
  return el('header', { class: 'phead' },
    el('div', {}, el('div', { class: 'eyebrow' }, 'Bait protein · AlphaSIFT report'),
      el('h1', { class: 'gene' }, b.gene || D.screen.name),
      el('div', { class: 'idline' }, b.uniprot ? el('span', { class: 'mono' }, b.uniprot) : null,
        el('span', {}, `${D.bait.length} aa`), b.description ? el('span', {}, b.description) : null, b.organism ? el('span', {}, b.organism) : null),
      el('div', { class: 'params' }, 'computed with ',
        el('span', { class: 'mono' }, `alphasift ${D.alphasift_version} · source ${D.screen.source} · d0 ${P.d0_signal} Å (0.5 = PAE ${P.d0_signal} Å) · contact Cβ–Cβ ≤ ${P.d0_contact} Å · pct ${P.pct} · ${P.P} null shifts · strength: ≥${P.min_single} res at signal ≥${P.tau_single} (PAE ≤ ${paeAt(P.tau_single)} Å), gap ≤ ${P.single_gap} · merge ${P.merge_method} / support ${P.merge_support} / d_merge ${P.d_merge} Å · motif engine ${P.motif_engine}`),
        D.screen.single_pair ? el('div', {}, 'Single prediction: recurrence not tested, strength regions only.') : null,
        Object.keys(D.screen.skipped || {}).length ? el('div', {}, `${Object.keys(D.screen.skipped).length} prediction(s) skipped (see summary.json).`) : null)),
    el('div', { class: 'tiles' }, tile(D.screen.n_predictions, 'predictions'), tile(binders, 'binders'),
      tile(nSig, 'sig. sites'), tile(D.interfaces.length, 'interfaces'), tile(nMot, 'motifs')));
}

// ------------------------------------------------------------------- side: selection
function side() {
  const box = el('div', { class: 'selbox' });
  const clear = el('button', { class: 'clear', on: { click: () => setSel({ iface: null, motif: null, carrier: null, partner: null }) } }, 'Clear selection');
  const draw = () => {
    box.replaceChildren();
    clear.disabled = S.iface === null && S.partner === null;
    if (S.iface === null && S.partner === null) {
      box.append(el('p', {}, 'Click an interface band or card, a motif lane, a carrier or a partner: it is lit across the page.'),
        el('p', {}, '3-D: open ', el('code', {}, 'alphasift.cxc'), ' in ChimeraX — every interface and motif is a named selection.'));
      return;
    }
    if (S.iface !== null) {
      const f = D.interfaces[S.iface];
      box.append(el('p', {}, el('b', {}, f.id), ` · bait ${compact(f.residues)} · ${f.significant ? 'significant' : 'not significant'}`));
      if (S.motif) { const m = f.motif.motifs.find(x => x.id === S.motif);
        box.append(el('p', {}, el('b', {}, `${f.id} ${m.id}`), ` · bait ${compact(m.positions)} · ${m.n_groups} groups · ${m.carriers.length} carriers`)); }
      box.append(el('p', {}, 'ChimeraX: ', el('code', {}, `view ${S.motif ? f.id + '_' + S.motif : f.id}`)));
    }
    if (S.partner) {
      box.append(el('p', {}, el('b', {}, nameOf(S.partner)), ' ', el('span', { class: 'mono' }, S.partner)));
      if (S.carrier) { const f = D.interfaces[S.iface]; const al = f.motif.alignments[S.motif];
        const c = al.carriers.find(x => x.id === S.carrier);
        if (c) box.append(el('p', {}, `${S.motif} residues on ${nameOf(c.id)}: ${compact(c.residues)} (${c.letters}), read ${c.dir === 1 ? 'N→C' : 'backwards, C→N'}`)); }
    }
  };
  onSel(draw); draw();
  return el('div', { class: 'card' }, el('div', { class: 'ctitle' }, 'Selection'), el('div', { class: 'csub' }, 'The website shows this in its 3-D viewer.'), box, clear);
}

// ------------------------------------------------------------------- heatmap
function heatmap() {
  const rows = D.matrix.row_order, sig = D.matrix.signal, pae = D.matrix.pae;
  const residues = [];
  D.matrix.residues.forEach((r, c) => { if (rows.some(i => sig[i][c] > 0)) residues.push([r, c]); });
  const nC = residues.length, nB = rows.length;
  const colOf = new Map(residues.map(([r], k) => [r, k]));
  const PAD_L = 52, PAD_R = 12, PAD_T = 6, TRACK_H = 13, GAP = 3, AXIS_H = 18;
  const rowH = Math.max(1, Math.min(3.4, 340 / Math.max(1, nB)));
  const matTop = PAD_T + TRACK_H * 2 + GAP + 8, matH = nB * rowH, H = matTop + matH + AXIS_H;
  const rm = {}; D.tracks.rmsf.residues.forEach((r, k) => rm[r] = D.tracks.rmsf.values[k]);
  const rmax = Math.max(1e-6, ...D.tracks.rmsf.values, 0);
  const ss = D.tracks.summed_signal, smax = Math.max(1e-6, ...residues.map(([r]) => ss[r - 1]));
  const lerp = (lo, hi, t) => { const h = x => [1,3,5].map(i => parseInt(x.slice(i, i+2), 16)); const a = h(lo), b = h(hi);
    t = Math.max(0, Math.min(1, t)); return `rgb(${a.map((v, i) => Math.round(v + (b[i] - v) * t)).join(',')})`; };
  const conf = t => lerp('#f7fcf5', '#00441b', t), rmsf = t => lerp('#efedf5', '#3f007d', t);
  const cv = el('canvas'), wrap = el('div', { class: 'hm' }, cv), lane = el('div', { class: 'lane' });
  const tip = el('div', { class: 'tip', style: { display: 'none' } }); document.body.append(tip);
  let W = 720;
  const draw = () => {
    W = wrap.clientWidth || 720; const dpr = Math.min(2, window.devicePixelRatio || 1);
    cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr); cv.style.height = H + 'px';
    const g = cv.getContext('2d'); g.setTransform(dpr, 0, 0, dpr, 0, 0); g.fillStyle = '#fff'; g.fillRect(0, 0, W, H);
    const cw = (W - PAD_L - PAD_R) / Math.max(1, nC), w = Math.max(1, Math.ceil(cw));
    residues.forEach(([r], k) => { const x = PAD_L + k * cw;
      g.fillStyle = rmsf((rm[r] || 0) / rmax); g.fillRect(x, PAD_T, w, TRACK_H);
      g.fillStyle = conf(ss[r - 1] / smax); g.fillRect(x, PAD_T + TRACK_H + GAP, w, TRACK_H); });
    rows.forEach((i, b) => { const y = matTop + b * rowH;
      residues.forEach(([, c], k) => { const v = sig[i][c]; if (v > 0.02) { g.fillStyle = conf(v);
        g.fillRect(PAD_L + k * cw, y, w, Math.max(1, Math.ceil(rowH))); } }); });
    const lit = litResidues();
    if (lit.size) { g.fillStyle = '#111418'; g.globalAlpha = 0.08;
      lit.forEach(r => { const k = colOf.get(r); if (k !== undefined) g.fillRect(PAD_L + k * cw, PAD_T, Math.max(1, cw), matTop + matH - PAD_T); });
      g.globalAlpha = 1; }
    const focus = S.partner && D.predictions.findIndex(p => p.id === S.partner);
    if (focus !== null && focus >= 0) { const b = rows.indexOf(focus); if (b >= 0) { g.strokeStyle = '#ffc400'; g.lineWidth = 2;
      g.strokeRect(PAD_L - 1, matTop + b * rowH - 1, W - PAD_L - PAD_R + 2, Math.max(2, rowH) + 2); } }
    g.fillStyle = '#7c828b'; g.font = '10px ui-monospace, monospace'; g.textAlign = 'right';
    g.fillText('RMSF', PAD_L - 6, PAD_T + TRACK_H - 2); g.fillText('recur', PAD_L - 6, PAD_T + 2 * TRACK_H + GAP - 2);
    g.save(); g.translate(12, matTop + matH / 2); g.rotate(-Math.PI / 2); g.textAlign = 'center'; g.fillText(`${nB} binders`, 0, 0); g.restore();
    g.textAlign = 'center'; const step = Math.max(1, Math.round(nC / 8));
    for (let k = 0; k < nC; k += step) g.fillText(String(residues[k][0]), PAD_L + (k + 0.5) * cw, matTop + matH + 13);
    // interface lane (HTML, so every band is a button)
    lane.replaceChildren(el('span', { class: 'laneLabel' }, 'interfaces'));
    D.interfaces.forEach((f, i) => f.segments.forEach(([a, b2]) => {
      const ks = []; for (let r = a; r <= b2; r++) if (colOf.has(r)) ks.push(colOf.get(r)); if (!ks.length) return;
      const x0 = PAD_L + Math.min(...ks) * cw, x1 = PAD_L + (Math.max(...ks) + 1) * cw;
      lane.append(el('span', { class: 'band', 'data-on': S.iface === i ? 1 : 0, title: `${f.id}: bait ${compact(f.residues)} · ${f.significant ? 'significant' : 'not significant'}`,
        style: Object.assign({ left: x0 + 'px', width: Math.max(2, x1 - x0) + 'px' }, ifaceFill(i, f.significant)),
        on: { click: () => { setSel({ iface: S.iface === i && !S.motif ? null : i, motif: null, carrier: null }); openCard(i); } } }));
    }));
  };
  cv.addEventListener('mousemove', e => {
    const rect = cv.getBoundingClientRect(), x = e.clientX - rect.left, y = e.clientY - rect.top;
    const cw = (W - PAD_L - PAD_R) / Math.max(1, nC), k = Math.floor((x - PAD_L) / cw), b = Math.floor((y - matTop) / rowH);
    if (k < 0 || k >= nC) { tip.style.display = 'none'; return; }
    const [r, c] = residues[k]; let t = `bait ${SEQ[r - 1] || ''}${r}`;
    if (y < matTop) t += `\nRMSF ${(rm[r] || 0).toFixed(2)} Å · summed signal ${ss[r - 1].toFixed(2)}`;
    else if (b >= 0 && b < nB) { const i = rows[b], p = D.predictions[i];
      t += `\n${nameOf(p.id)} (${p.id})\n` + (sig[i][c] > 0 ? `signal ${sig[i][c].toFixed(2)}` + (pae[i][c] !== null ? ` · PAE ${pae[i][c].toFixed(1)} Å` : '') : 'no contact'); }
    else { tip.style.display = 'none'; return; }
    tip.textContent = t; tip.style.display = 'block'; tip.style.left = (e.clientX + 12) + 'px'; tip.style.top = (e.clientY + 12) + 'px';
  });
  cv.addEventListener('mouseleave', () => tip.style.display = 'none');
  cv.addEventListener('click', e => { const rect = cv.getBoundingClientRect(), b = Math.floor((e.clientY - rect.top - matTop) / rowH);
    if (b >= 0 && b < nB) { const id = D.predictions[rows[b]].id; setSel({ partner: S.partner === id ? null : id, carrier: null }); } });
  onSel(draw); new ResizeObserver(draw).observe(wrap);
  const legends = el('div', { class: 'legends' },
    el('span', {}, 'confidence', el('span', { class: 'ramp', style: { background: 'linear-gradient(90deg,#f7fcf5,#00441b)' } }), `low → high (0.5 = PAE ${P.d0_signal} Å)`),
    el('span', {}, el('span', { class: 'sw', style: { background: IFC[0] } }), 'significant ',
      el('span', { class: 'sw', style: ifaceFill(0, false) }), 'not significant'),
    el('span', {}, 'RMSF', el('span', { class: 'ramp', style: { background: 'linear-gradient(90deg,#efedf5,#3f007d)' } }), 'low → high'));
  return el('section', { class: 'card' }, el('div', { class: 'ctitle' }, 'Confidence matrix'),
    el('div', { class: 'csub' }, 'One row per binder (clustered), one column per bait residue some binder contacts. Hover to read PAE; click a row to select that partner, a band to select an interface.'),
    el('div', { class: 'plate' }, wrap, lane), legends);
}

// ------------------------------------------------------------------- page key
function pageKey() {
  const aa = el('div', { class: 'key' }, ...AAG.map(([n, l, c]) => el('span', {}, el('span', { class: 'aa', style: { color: c } }, l), ' ', n)));
  const marks = el('div', { class: 'key' }, el('span', {}, `✱ q < ${Q_SIG} vs letter shuffles (tint + bar)`),
    el('span', {}, 'hatched = untested: < 3 partner groups, not negative'), el('span', {}, 'column height ∝ √ partner groups (busiest column = full); letters = shares'),
    el('span', {}, 'motif lane — click: light it across the page (gold = selected)'), el('span', {}, '→ carrier read N→C · ← read backwards (antiparallel, flipped)'));
  const tiers = el('div', { class: 'key' }, el('span', {}, 'Partner residues by confidence:'),
    ...['H', 'M', 'L'].map(t => el('span', {}, el('span', { class: `tier tier${t}` }, t), ' ', tierLabel[t], t === 'L' ? ' (hidden unless asked)' : '')),
    el('span', {}, el('span', { class: 'range tierM' }, el('span', { class: 'mres' }, 'LL'), 'ML'), ' underlined = motif residue'));
  const eng = el('div', { class: 'key' }, `Motif engine ${P.motif_engine} · contact Cβ–Cβ ≤ ${P.motif_touch_cb} Å · partners ≥ ${Math.round(P.motif_id_group * 100)}% identical = one group · ${P.motif_n_perm} letter shuffles`);
  return el('div', { class: 'plate', style: { marginBottom: '12px' } }, el('div', { class: 'hbtn' }, 'Key — colours and marks used in every card'), aa, marks, tiers, eng);
}

// ------------------------------------------------------------------- motif engine view
const COL_W = 16, JUMP_W = 5, ROW_GROUPS = 12, ROW_SIG = 12, LOGO_H = 80, ROW_BAIT = 16, ROW_NUM = 26, LANE_H = 13, LANE_GAP = 4;
const ADV = { W: 0.944, M: 0.833 };
function glyph(ch, x, y, w, h, color) {
  if (h < 2.5) return null; const sx = Math.min((w * 0.9) / 7.2, (w * 0.95) / ((ADV[ch] || 0.72) * 10));
  return sv('text', { transform: `translate(${x + w / 2} ${y}) scale(${sx} ${h / 7.2})`, 'font-size': 10, 'font-weight': 700, 'text-anchor': 'middle', fill: color }, ch);
}
function baitSide(fi) {
  const f = D.interfaces[fi], e = f.motif, cols = e.panel1;
  const xs = []; let x = 0; cols.forEach((c, i) => { if (i > 0 && c.pos - cols[i - 1].pos > 1) x += JUMP_W; xs.push(x); x += COL_W; });
  const maxG = Math.max(1, ...cols.map(c => c.n_groups)), width = (xs.length ? xs[xs.length - 1] + COL_W : COL_W) + 2;
  const yLogo = ROW_GROUPS + ROW_SIG, yBait = yLogo + LOGO_H, yNum = yBait + ROW_BAIT, yLanes = yNum + ROW_NUM + 4;
  const height = yLanes + e.motifs.length * (LANE_H + LANE_GAP);
  const gut = el('div', { class: 'gut', style: { height: height + 'px' } });
  [[0, ROW_GROUPS, 'groups'], [ROW_GROUPS, ROW_SIG, `✱ q < ${Q_SIG}`], [yLogo, LOGO_H, 'partner residues'], [yBait, ROW_BAIT, 'bait residue'], [yNum, ROW_NUM, 'position']]
    .forEach(([y, h, t]) => gut.append(el('div', { class: 'gl', style: { top: y + 'px', height: h + 'px', lineHeight: h + 'px' } }, t)));
  const svg = sv('svg', { width, height });
  svg.append(sv('defs', {}, sv('pattern', { id: `hatch${fi}`, width: 5, height: 5, patternUnits: 'userSpaceOnUse', patternTransform: 'rotate(45)' },
    sv('line', { x1: 0, y1: 0, x2: 0, y2: 5, stroke: '#d6d9de', 'stroke-width': 2 }))));
  cols.forEach((c, i) => {
    const X = xs[i], untested = c.q === null, sig = !untested && c.q < Q_SIG, colH = LOGO_H * Math.sqrt(c.n_groups / maxG);
    const bait = SEQ[c.pos - 1] || '?', letters = Object.entries(c.letters).sort((a, b) => b[1] - a[1]);
    const g = sv('g', {}, sv('title', {}, `bait ${bait}${c.pos} · ${c.n_groups} group${c.n_groups === 1 ? '' : 's'} · ` +
      (untested ? 'untested (< 3 partner groups)' : `q ${fmtQ(c.q)}${sig ? ' (significant)' : ''}`) + (c.agree !== null ? ` · agreement ${c.agree.toFixed(2)}` : '') +
      `\npartner residues: ${letters.map(([l, s]) => `${l} ${Math.round(s * 100)}%`).join(', ')}`));
    if (sig) g.append(sv('rect', { x: X, y: yLogo - 1, width: COL_W, height: LOGO_H + 2, fill: '#e6ebf2', rx: 2 }));
    if (untested) g.append(sv('rect', { x: X, y: yLogo, width: COL_W, height: LOGO_H, fill: `url(#hatch${fi})` }));
    if (sig) g.append(sv('rect', { x: X + 1, y: yLogo + LOGO_H, width: COL_W - 2, height: 2, fill: '#191c20' }));
    g.append(sv('text', { x: X + COL_W / 2, y: ROW_GROUPS - 3, 'text-anchor': 'middle', 'font-size': 9, fill: untested ? '#a3a9b1' : '#4a5058' }, c.n_groups));
    if (sig) g.append(sv('text', { x: X + COL_W / 2, y: ROW_GROUPS + ROW_SIG - 2, 'text-anchor': 'middle', 'font-size': 10 }, '✱'));
    let y = yLogo + LOGO_H;
    [...letters].reverse().filter(([, s]) => s >= 0.04).forEach(([l, s]) => { const h = s * colH; y -= h;
      const gl = glyph(l, X, y + h, COL_W, h, untested ? '#b3b8bf' : aaColor(l)); if (gl) g.append(gl); });
    g.append(sv('text', { x: X + COL_W / 2, y: yBait + 12, 'text-anchor': 'middle', 'font-size': 12, 'font-weight': 700, fill: aaColor(bait) }, bait));
    g.append(sv('text', { transform: `translate(${X + COL_W / 2 + 3} ${yNum + 2}) rotate(90)`, 'font-size': 9, fill: '#7c828b' }, c.pos));
    g.append(sv('rect', { x: X, y: 0, width: COL_W, height: yNum + ROW_NUM, fill: 'transparent' }));
    svg.append(g);
  });
  const colOf = new Map(cols.map((c, i) => [c.pos, i]));
  e.motifs.forEach((m, li) => {
    const y = yLanes + li * (LANE_H + LANE_GAP), idx = m.positions.map(p => colOf.get(p)).filter(v => v !== undefined);
    const on = S.iface === fi && S.motif === m.id, toggle = () => setSel(on ? { iface: fi, motif: null, carrier: null } : { iface: fi, motif: m.id, carrier: null });
    gut.append(el('button', { class: 'laneBtn', 'aria-pressed': on ? 'true' : 'false', style: { top: y + 'px', height: LANE_H + 'px', lineHeight: (LANE_H - 2) + 'px' },
      on: { click: toggle } }, m.id));
    if (!idx.length) return;
    const x0 = xs[Math.min(...idx)], x1 = xs[Math.max(...idx)] + COL_W;
    const g = sv('g', { style: 'cursor:pointer' }, sv('title', {}, `${m.id}: bait ${compact(m.positions)}`),
      sv('line', { x1: x0, x2: x1, y1: y + LANE_H / 2, y2: y + LANE_H / 2, stroke: '#9aa1aa', 'stroke-width': 1 }),
      ...idx.map(ci => sv('rect', { x: xs[ci] + 0.5, y, width: COL_W - 1, height: LANE_H, rx: 2, fill: on ? '#ffc400' : '#3b4552', stroke: on ? '#111418' : 'none', 'stroke-width': on ? 1.5 : 0 })));
    g.addEventListener('click', toggle); svg.append(g);
  });
  return el('div', { class: 'bs' }, gut, el('div', { class: 'scroll' }, svg));
}
function partnerLogo(al) {
  const PCOL_W = 18, GAP_W = 30, PLOGO_H = Math.max(48, Math.min(96, al.n * 8));
  const xs = []; let acc = 0; al.columns.forEach(c => { xs.push(acc); acc += c.kind === 'gap' ? GAP_W : PCOL_W; });
  const yLogo = 12, yOff = yLogo + PLOGO_H + 12, height = yOff + 4;
  const gut = el('div', { class: 'gut', style: { height: height + 'px' } },
    el('div', { class: 'gl', style: { top: 0, height: '12px', lineHeight: '12px' } }, `k of ${al.n}`),
    el('div', { class: 'gl', style: { top: yLogo + 'px', height: PLOGO_H + 'px', lineHeight: PLOGO_H + 'px' } }, 'partner residues'),
    el('div', { class: 'gl', style: { top: (yLogo + PLOGO_H) + 'px', height: '14px', lineHeight: '14px' } }, 'offset'));
  const svg = sv('svg', { width: acc + 2, height });
  al.columns.forEach((c, i) => { const x = xs[i];
    if (c.kind === 'gap') { svg.append(sv('g', {}, sv('title', {}, `${c.skipped} partner residues where no carrier has a motif residue (collapsed)`),
      sv('text', { x: x + GAP_W / 2, y: yLogo + PLOGO_H - 4, 'text-anchor': 'middle', fill: '#7c828b' }, '⋯'),
      sv('text', { x: x + GAP_W / 2, y: yOff, 'text-anchor': 'middle', 'font-size': 9, fill: '#7c828b' }, `${c.skipped} aa`))); return; }
    const g = sv('g', {}, sv('title', {}, `offset ${c.offset > 0 ? '+' : ''}${c.offset} · ${c.k}/${al.n} carriers` + (c.counts.length ? `\n${c.counts.map(([l, n]) => `${l} ${n}`).join(', ')}` : '')));
    if (c.offset === 0) g.append(sv('rect', { x, y: yLogo - 1, width: PCOL_W, height: PLOGO_H + 2, fill: '#eef0f3', rx: 2 }));
    g.append(sv('text', { x: x + PCOL_W / 2, y: 9, 'text-anchor': 'middle', 'font-size': 9, fill: c.k ? '#4a5058' : '#c3c7cd' }, c.k || '·'));
    let y = yLogo + PLOGO_H;
    [...c.counts].reverse().forEach(([l, n]) => { const h = (n / Math.max(1, al.n)) * PLOGO_H; y -= h; const gl = glyph(l, x, y + h, PCOL_W, h, aaColor(l)); if (gl) g.append(gl); });
    g.append(sv('text', { x: x + PCOL_W / 2, y: yOff, 'text-anchor': 'middle', 'font-size': 9, 'font-weight': c.offset === 0 ? 700 : 400, fill: c.offset === 0 ? '#191c20' : '#7c828b' }, c.offset > 0 ? `+${c.offset}` : c.offset));
    g.append(sv('rect', { x, y: 0, width: PCOL_W, height, fill: 'transparent' })); svg.append(g); });
  return el('div', { class: 'bs' }, gut, el('div', { class: 'scroll' }, svg));
}
const linkerSummary = ls => { const n = new Map(); ls.forEach(([d, l]) => n.set(`${d} ${l}`, (n.get(`${d} ${l}`) || 0) + 1));
  return [...n].map(([k, c]) => c > 1 ? `${k} ×${c}` : k).join(', ') || '—'; };
function motifPanel(fi, m) {
  const f = D.interfaces[fi], e = f.motif, al = e.alignments[m.id];
  const det = el('details', { class: 'motif' });
  const head = el('summary', {}, el('b', {}, m.id), el('span', { class: 'mmeta' },
    `${m.positions.length} bait res · ${m.n_groups} group${m.n_groups === 1 ? '' : 's'} · ${m.carriers.length} carrier${m.carriers.length === 1 ? '' : 's'} · ${m.elements.length} element${m.elements.length === 1 ? '' : 's'}`),
    m.bipartite ? el('span', { class: 'badge' }, 'bipartite') : null, ...m.linked_by.map(l => el('span', { class: 'badge' }, `linked by ${l}`)));
  const body = el('div');
  det.append(head, body);
  const fill = () => {
    body.replaceChildren();
    const nFl = al.carriers.filter(c => c.dir === -1).length;
    body.append(el('div', { class: 'mnote' }, `Aligned on bait residue ${al.anchor}.`,
      nFl ? ' Carriers marked ← bind antiparallel and are read backwards (C→N, flipped) so every carrier runs along the bait.' : '',
      al.orientation === 'per-carrier' ? ' Direction is set per carrier from its own contacts (the elements disagree or are mixed).' : '',
      al.n < al.carriers.length ? ` ${al.carriers.length - al.n} carrier(s) share no bait position with the anchor holders and are not in the logo.` : '',
      al.missing_letters ? ` ${al.missing_letters} residue(s) are outside the partner's listed ranges and show as ?.` : ''));
    body.append(partnerLogo(al));
    const groups = new Map(); al.carriers.forEach(c => groups.set(c.group, [...(groups.get(c.group) || []), c]));
    const chips = el('div', { class: 'chips' });
    groups.forEach(g => chips.append(el('span', { class: g.length > 1 ? 'pgroup' : '', title: g.length > 1 ? `paralog group: these ${g.length} carriers count as ONE independent group` : null },
      ...g.map(c => { const on = S.iface === fi && S.motif === m.id && S.carrier === c.id;
        return el('button', { class: 'chip', 'aria-pressed': on ? 'true' : 'false', title: `${nameOf(c.id)} · ${c.id} · motif residues ${compact(c.residues)} (${c.letters})`,
          on: { click: () => setSel(on ? { carrier: null, partner: null } : { iface: fi, motif: m.id, carrier: c.id, partner: c.id }) } },
          `${nameOf(c.id)} `, el('b', {}, c.dir === 1 ? '→' : '←'), c.aligned ? null : el('span', { class: 'cun' }, ' (not in logo)')); }))));
    body.append(chips);
    const inM = new Set(m.elements), pairs = e.bipartite.filter(b => inM.has(b.a) && inM.has(b.b));
    const tbl = el('table', { class: 'el' }, el('tr', {}, ...['el.', 'bait residues', 'groups', 'q', 'or.', 'spacing', 'carriers'].map((h, k) => el('th', { class: k === 1 ? 'l' : '' }, h))),
      ...m.elements.map(ei => { const x = e.elements[ei]; return el('tr', {}, el('td', { class: 'l' }, el('b', {}, `E${ei + 1}`)), el('td', { class: 'l mono' }, compact(x.positions)),
        el('td', {}, x.n_groups), el('td', {}, fmtQ(x.q)), el('td', {}, { parallel: '→', antiparallel: '←', mixed: '↔' }[x.orientation] || '?'),
        el('td', {}, x.consistency === null ? '—' : Math.round(x.consistency * 100) + '%'), el('td', {}, x.carriers.length)); }));
    const ptbl = pairs.length ? el('table', { class: 'el' }, el('tr', {}, ...['pair', 'shared groups', 'q', 'linker', 'fixed'].map(h => el('th', {}, h))),
      ...pairs.map(b => el('tr', {}, el('td', { class: 'l mono' }, `E${b.a + 1} + E${b.b + 1}`), el('td', {}, b.k), el('td', {}, fmtQ(b.q)),
        el('td', { class: 'l mono' }, linkerSummary(b.linkers)), el('td', {}, b.fixed_linker ? 'yes' : '')))) : null;
    const box = el('div', { style: { display: 'none' } }, tbl, ptbl);
    const tg = el('button', { class: 'elToggle', on: { click: () => { const o = box.style.display === 'none'; box.style.display = o ? 'block' : 'none';
      tg.firstChild.textContent = (o ? '▾' : '▸') + tg.firstChild.textContent.slice(1); } } },
      `▸ Elements (${m.elements.length})${pairs.length ? ` · ${pairs.length} bipartite pair${pairs.length === 1 ? '' : 's'}` : ''}`);
    body.append(tg, box);
  };
  det.addEventListener('toggle', () => { if (det.open) fill(); });
  onSel(() => { if (det.open) { const st = body.querySelector('div[style*="block"]'); fill(); if (st) body.lastChild.style.display = 'block'; } });
  return det;
}
function motifView(fi) {
  const f = D.interfaces[fi], e = f.motif;
  if (!e) return el('div', { class: 'plate' }, el('div', { class: 'summary' }, 'No motif-engine entry for this interface: none of its partners gave the engine a contact.'));
  const p1 = e.panel1, nT = p1.filter(c => c.q !== null).length, nS = p1.filter(c => c.q !== null && c.q < Q_SIG).length;
  const bsBox = el('div');
  const drawBs = () => bsBox.replaceChildren(p1.length ? baitSide(fi) : el('div', { class: 'summary' }, 'No bait residue of this interface is touched by a partner the engine could use.'));
  drawBs(); onSel(drawBs);
  const rel = e.motifs.length === 0 ? el('div', { class: 'summary' }, 'No motif found in this interface (no element passed the engine\'s tests) — a valid result, not missing data.')
    : e.motifs.length > 1 ? el('div', { class: 'summary' }, e.motif_relations.map(r => `${r.a}–${r.b} ${r.relation === 'same_pocket' ? 'same pocket' : 'separate sites'} · ${r.status === 'different' ? 'different (competing)' : 'unresolved, too few carriers to tell'}`).join(' · ')) : null;
  return el('div', { class: 'plate' }, el('div', { class: 'mhead' }, el('span', { class: 'hbtn' }, 'Motifs — bait side'),
    el('span', { class: 'hmeta' }, `${p1.length} bait residues · ${nT} tested · ${nS} with q < ${Q_SIG} · partners used ${e.partners_used}/${e.partners_total} · shuffle control: ${e.control.false_elements} false element${e.control.false_elements === 1 ? '' : 's'}, columns ${e.control.false_columns.join('/')}`)),
    bsBox, rel, ...e.motifs.map(m => motifPanel(fi, m)));
}

// ------------------------------------------------------------------- partners
function partnerList(fi) {
  const f = D.interfaces[fi], parts = f.partners, wrap = el('div', { class: 'pWrap' });
  if (parts === null) return el('div', { class: 'pWrap summary' }, 'Confidently placed partners n/a: this screen has no pairwise PAE.');
  if (!parts.length) return el('div', { class: 'pWrap summary' }, f.significant ? 'No confidently placed partner (this site is recurrent: many partners converge here with low confidence).' : `No confidently placed partner (no contact at PAE ≤ ${paeAt(P.partner_tau_in)} Å).`);
  const mem = new Map();
  if (f.motif) f.motif.motifs.forEach(mo => f.motif.alignments[mo.id].carriers.forEach(c => mem.set(c.id, [...(mem.get(c.id) || []), { id: mo.id, dir: c.dir, res: c.residues }])));
  const st = { open: false, low: false, all: false };
  const draw = () => {
    wrap.replaceChildren();
    const nIn = parts.filter(p => mem.has(p.id)).length;
    const tg = el('button', { class: 'pListToggle', on: { click: () => { st.open = !st.open; draw(); } } }, `${st.open ? '▾' : '▸'} Confidently placed partners · ${parts.length}${f.motif ? ` · ${nIn} carry a motif` : ''}`);
    if (!st.open) { wrap.append(tg); return; }
    const low = el('input', { type: 'checkbox', on: { change: ev => { st.low = ev.target.checked; draw(); } } }); low.checked = st.low;
    wrap.append(el('div', { class: 'pHead' }, tg, el('label', { class: 'pToggle' }, low, ' show low-confidence residues')));
    (st.all ? parts : parts.slice(0, 10)).forEach(p => {
      const ms = mem.get(p.id) || [], mres = new Set(ms.flatMap(x => x.res));
      const ranges = st.low ? p.ranges : p.ranges.filter(r => r[2] !== 'L');
      wrap.append(el('div', { class: 'pRow', 'data-on': S.partner === p.id ? 1 : 0, on: { click: () => setSel({ partner: S.partner === p.id ? null : p.id, carrier: null }) } },
        el('div', { class: 'pTop' }, el('span', { class: 'pName' }, nameOf(p.id)), el('span', { class: 'pSub mono' }, subOf(p.id)),
          ...ms.map(x => { const on = S.iface === fi && S.motif === x.id && S.carrier === p.id;
            return el('button', { class: 'mBadge', 'aria-pressed': on ? 'true' : 'false', title: `carries ${x.id} (read ${x.dir === 1 ? 'N→C' : 'backwards, C→N'})`,
              on: { click: ev => { ev.stopPropagation(); setSel(on ? { carrier: null } : { iface: fi, motif: x.id, carrier: p.id, partner: p.id }); } } }, `${x.id} ${x.dir === 1 ? '→' : '←'}`); }),
          el('span', { class: 'pStats' }, el('span', { title: tierLabel.H }, 'H ', el('b', {}, p.n_high)), el('span', { title: tierLabel.M }, 'M ', el('b', {}, p.n_medium)),
            st.low ? el('span', { title: tierLabel.L }, 'L ', el('b', {}, p.n_low)) : null,
            el('span', {}, 'best ', el('b', {}, p.best_pae.toFixed(1)), ' Å'), el('span', {}, 'median ', el('b', {}, p.median_pae.toFixed(1)), ' Å'),
            el('span', {}, el('b', {}, p.bait_residues.length), ' bait res'))),
        el('div', { class: 'pRanges' }, ...ranges.map(r => el('span', { class: `range tier${r[2]}`, title: `${nameOf(p.id)} ${r[0] === r[1] ? r[0] : r[0] + '–' + r[1]} ${r[3]} · ${tierLabel[r[2]]}` },
          el('span', { class: 'rPos' }, r[0] === r[1] ? r[0] : `${r[0]}–${r[1]}`), ' ',
          ...[...r[3]].map((ch, k) => mres.has(r[0] + k) ? el('span', { class: 'mres' }, ch) : ch))))));
    });
    if (parts.length > 10) wrap.append(el('button', { class: 'allToggle', on: { click: () => { st.all = !st.all; draw(); } } }, st.all ? '▴ Show first 10' : `▾ Show the other ${parts.length - 10}`));
  };
  draw(); onSel(() => { if (st.open) draw(); });
  return wrap;
}

// ------------------------------------------------------------------- interface cards
const cardOpen = new Map();
function openCard(i) { const c = cardOpen.get(i); if (c) c(true); }
function interfaceCards() {
  const sec = el('section', { class: 'card' }, el('div', { class: 'ctitle' }, `Interfaces (${D.interfaces.length})`),
    el('div', { class: 'csub' }, 'Sites merged in 3-D. Open a card for its motif engine view and its partners; selecting one lights it in the matrix.'), pageKey());
  D.interfaces.forEach((f, i) => {
    const body = el('div', { class: 'ibody', style: { display: 'none' } });
    let built = false;
    const set = open => { if (open && !built) { body.append(motifView(i), partnerList(i)); built = true; } body.style.display = open ? 'block' : 'none'; };
    cardOpen.set(i, set);
    const n = f.partners === null ? 'confidently placed partners n/a' : `${f.partners.length} confidently placed partner${f.partners.length === 1 ? '' : 's'}`;
    const nm = f.motif ? f.motif.motifs.length : 0;
    const head = el('div', { class: 'ihead', on: { click: () => { const o = body.style.display === 'none'; set(o); setSel({ iface: o ? i : (S.iface === i ? null : S.iface), motif: null, carrier: null }); } } },
      el('span', { class: 'sw', style: Object.assign({ width: '14px', height: '14px' }, ifaceFill(i, f.significant)) }), el('span', { class: 'iid' }, f.id),
      f.significant ? null : el('span', { class: 'ns' }, 'not significant'), el('span', { class: 'seg mono' }, f.segments.map(s => `${s[0]}–${s[1]}`).join('  ')),
      el('span', { class: 'imeta' }, `${n} · ${nm} motif${nm === 1 ? '' : 's'}`));
    onSel(() => head.setAttribute('data-on', S.iface === i ? 1 : 0));
    sec.append(el('div', { class: 'icard' }, head, body));
  });
  return sec;
}

// ------------------------------------------------------------------- all binders
function allBinders() {
  const rows = D.matrix.row_order.map(i => { const p = D.predictions[i];
    return { id: p.id, name: nameOf(p.id), sub: subOf(p.id), count: D.matrix.signal[i].filter(v => v > 0).length }; });
  const st = { q: '', sort: 'contacts' }, list = el('div', { class: 'blist' }), count = el('span', { class: 'bc' });
  const draw = () => { const q = st.q.trim().toLowerCase();
    let r = q ? rows.filter(x => (x.name + ' ' + x.sub + ' ' + x.id).toLowerCase().includes(q)) : rows;
    r = [...r].sort((a, b) => st.sort === 'name' ? a.name.localeCompare(b.name) : b.count - a.count);
    count.textContent = `${r.length} shown`;
    list.replaceChildren(...r.map(x => el('button', { class: 'brow', 'data-on': S.partner === x.id ? 1 : 0, on: { click: () => setSel({ partner: S.partner === x.id ? null : x.id, carrier: null }) } },
      el('span', {}, el('b', {}, x.name), ' ', el('span', { class: 'pSub mono' }, x.sub || x.id)), el('span', { class: 'bc mono' }, `${x.count} res`)))); };
  const sb = (k, t) => { const b = el('button', { class: 'sortBtn', on: { click: () => { st.sort = k; draw(); sbs.forEach(([kk, bb]) => bb.setAttribute('data-on', st.sort === kk ? 1 : 0)); } } }, t);
    b.setAttribute('data-on', st.sort === k ? 1 : 0); return [k, b]; };
  const sbs = [sb('contacts', 'contacts'), sb('name', 'name')];
  draw(); onSel(draw);
  return el('section', { class: 'card' }, el('div', { class: 'ctitle' }, `All binders (${rows.length})`),
    el('div', { class: 'csub' }, 'Every candidate partner that produced interface signal. Search by name or UniProt; click one to select it.'),
    el('div', { class: 'controls' }, el('input', { placeholder: 'Search gene / UniProt…', on: { input: ev => { st.q = ev.target.value; draw(); } } }), ...sbs.map(x => x[1]), count), list);
}

document.getElementById('app').append(el('div', { class: 'wrap' }, header(),
  el('div', { class: 'grid' }, el('div', { class: 'side' }, side()), el('div', {}, heatmap(), interfaceCards(), allBinders())),
  el('div', { class: 'foot' }, `AlphaSIFT ${D.alphasift_version} · schema ${D.schema} · ${D.screen.name} · the same data is in alphasift.json`)));
"""


def html(doc):
    """The report for one `alphasift/1` document, as a string."""
    data = json.dumps(doc, separators=(",", ":")).replace("</", "<\\/")
    bait = (doc["bait"].get("identity") or {}).get("gene") or doc["screen"]["name"]
    return (f"<!doctype html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">"
            f"<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">"
            f"<title>{bait} — AlphaSIFT</title><style>{CSS}</style></head>"
            f"<body><div id=\"app\"></div>"
            f"<script type=\"application/json\" id=\"alphasift-data\">{data}</script>"
            f"<script>{JS}</script></body></html>\n")


def write_report(doc, path):
    path.write_text(html(doc), encoding="utf-8")
