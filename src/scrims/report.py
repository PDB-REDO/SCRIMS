"""
The HTML report: one self-contained page per screen that mirrors the companion
website's protein page, drawn in the browser from the `scrims/1` document
embedded in it. No network, no build step: open the file.

What it reproduces (website src/components, 2026-10-02): the header with its count
tiles and parameter line; the heatmap (RMSF and recurrence tracks, binder rows in
display order, contacted residues as columns, one interface lane, hover read-out);
the page key (Jalview Zappo amino-acid colours, motif marks, partner tiers, engine
parameters); one card per interface with the motif engine's bait-side view, motif
lanes, relations line and a panel per motif (partner-sequence logo over a collapsible
per-carrier table in the same columns), then the partner list (H and M shown, L on
request, motif badges, motif residues underlined); and the searchable list of all
binders. Selecting an interface, motif, carrier or partner lights it across the page
(cyan, as on the website).

Motif view (website round 3, 2026-10-07; src/components/MotifEngineView.tsx,
src/lib/logo.ts): both logos have two views -- "Classic logo (bits)", WebLogo
information content with the exact small-sample correction, group-weighted (related
partners count once), y axis scaled to the logo's tallest column; and "Evidence view"
(bait side: height ~ sqrt(groups); partner side: carriers with a residue / all
carriers; letters = shares). Defaults: bait side evidence, partner side classic.
Letters are stacked contiguously from the baseline, largest on top, each fitted to its
slot from its measured ink box; shares too small to draw form one grey band on top.
A motif lane filters the bait-side view to that motif's residues, letters from its
carriers only. The logo maths is LOGO_JS, mirrored in scrims.logo for the figures.

What it does not: the 3-D viewer. Its place holds the current selection and the
ChimeraX command that shows it (scrims.cxc names every interface and motif).
"""

from __future__ import annotations

import json

CSS = r"""
:root{--ground:#fbfbfa;--panel:#fff;--panel-2:#f3f4f6;--plate:#fff;--ink:#191c20;--ink-2:#4a5058;--ink-3:#7c828b;
--line:#e2e4e8;--line-2:#d0d3d9;--sel:#111418;--sel-fill:rgba(17,20,24,.06);--hl:#00ffff;--hl-tint:#d9ffff;
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
.aa{font-weight:800;font-family:var(--mono)}.aak{display:inline-block;min-width:14px;padding:0 3px;text-align:center;font:700 10.5px/15px var(--mono);color:#111418;border-radius:2px;letter-spacing:1px}
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
.laneBtn[aria-pressed="true"]{background:var(--hl);border-color:#111418}
.scroll{overflow-x:auto}svg text{font-family:var(--sans)}
.summary{font-size:12.5px;color:#4a5058;margin:6px 0 2px}
.mplate{--ink:#191c20;--ink-2:#4a5058;--ink-3:#7c828b;--line:#e2e4e8;--line-2:#d0d3d9;--panel-2:#f3f4f6;--sel:#111418;
--sel-fill:rgba(17,20,24,.07);color:var(--ink);border:1px solid var(--line-2);padding:14px;overflow:visible}
.mplate .mhead{align-items:center}
.headBtn{appearance:none;background:none;border:0;padding:0;font:600 11px/1.4 var(--sans);letter-spacing:.09em;text-transform:uppercase;color:var(--ink-3);cursor:pointer}
.headBtn:hover{color:var(--ink)}
.vt{display:inline-flex;margin-left:auto;border:1px solid var(--line-2);border-radius:6px;overflow:hidden}
.vt button{appearance:none;border:0;background:transparent;font:11px/1.6 var(--sans);color:var(--ink-2);padding:1px 9px;cursor:pointer}
.vt button+button{border-left:1px solid var(--line-2)}.vt button[aria-pressed="true"]{background:#111418;color:#fff}
.tick{position:absolute;right:6px;font:8.5px/10px var(--mono);color:var(--ink-3)}
.glAxis{position:absolute;left:0;right:34px;font-size:9.5px;color:var(--ink-3);text-align:right}
svg text.glyph{font-family:Arial,Helvetica,"Liberation Sans",sans-serif}
.filterBanner{display:flex;align-items:center;gap:8px;flex-wrap:wrap;font-size:11.5px;color:var(--ink-2);background:var(--hl-tint);border:1px solid #111418;border-radius:6px;padding:4px 10px;margin-bottom:6px}
.showAll{appearance:none;margin-left:auto;border:1px solid var(--ink-3);background:#fff;border-radius:4px;font:11px/1.5 var(--sans);padding:0 8px;cursor:pointer;color:var(--ink)}
.classicNote{font-size:11.5px;color:#6b4e00;background:#fff8dc;border:1px dashed #d9a400;border-radius:6px;padding:5px 10px;margin:4px 0 6px;line-height:1.45}
.axisNote{font-size:10.5px;color:var(--ink-3);margin:2px 0}
details.motif{border:1px solid var(--line);border-radius:8px;margin-top:8px}
details.motif>summary{display:flex;align-items:center;gap:8px;flex-wrap:wrap;padding:7px 10px;cursor:pointer;font-size:12.5px;list-style:none;color:var(--ink)}
details.motif>summary::-webkit-details-marker{display:none}
details.motif>summary::before{content:'\25B8';color:var(--ink-3);font-size:11px;width:10px}
details.motif[open]>summary::before{content:'\25BE'}
details.motif>summary b{font-family:var(--mono)}
.mbody{padding:4px 10px 10px}
.mmeta{color:var(--ink-2);font-size:12px}.badge{font-size:10.5px;border:1px solid var(--line-2);border-radius:4px;padding:0 6px;line-height:1.6;color:var(--ink-2)}
.mtop{display:flex;align-items:flex-start;gap:10px;flex-wrap:wrap}
.mnote{flex:1;min-width:240px;font-size:11.5px;color:var(--ink-3);line-height:1.45;margin-bottom:6px}.mnote b{color:var(--ink-2)}
.ptab{overflow-x:auto;padding-bottom:4px;margin-top:4px}
.prow{display:flex;align-items:flex-start;min-width:max-content}
.pleft{position:sticky;left:0;z-index:1;width:210px;flex:none;background:var(--plate)}
.pleft .gl{right:6px}
.tableToggle{position:sticky;left:0;display:block;appearance:none;background:transparent;border:1px solid var(--line-2);border-radius:6px;padding:2px 10px;margin:6px 0 4px;font:11.5px/1.5 var(--sans);color:var(--ink-2);cursor:pointer}
.crow{display:flex;align-items:center;min-width:max-content;border-top:1px solid var(--line);cursor:pointer;font-size:11.5px}
.crow:hover,.crow:hover .cleft{background:var(--panel-2)}
.crow[aria-pressed="true"],.crow[aria-pressed="true"] .cleft{background:var(--hl-tint);box-shadow:inset 3px 0 0 var(--hl),inset 4px 0 0 #111418}
.cleft{position:sticky;left:0;z-index:1;flex:none;width:210px;display:flex;align-items:baseline;gap:5px;padding:3px 8px 3px 2px;background:var(--plate);overflow:hidden;white-space:nowrap}
.cname{font-weight:650;overflow:hidden;text-overflow:ellipsis}.cuni{font:10px/1 var(--mono);color:var(--ink-3)}
.cgrp{font:9.5px/1.4 var(--mono);color:var(--ink-2);border:1px dashed var(--ink-3);border-radius:3px;padding:0 3px}
.carrow{font-weight:700;color:var(--ink-2)}.crange{margin-left:auto;font:10px/1 var(--mono);color:var(--ink-2)}
.cun{font-size:10.5px;color:var(--ink-3);font-style:italic}
.sumRow{display:flex;align-items:center;gap:10px;flex-wrap:wrap;margin:4px 0 8px}.consBoxes{display:inline-flex;align-items:center;gap:2px}
.cBox{display:inline-flex;height:20px;border-radius:3px;overflow:hidden;box-shadow:inset 0 0 0 1px rgba(0,0,0,.18)}
.cPart{display:inline-flex;align-items:center;justify-content:center;font:700 11px/20px var(--mono);color:#111418}
.cDot{display:inline-flex;width:12px;justify-content:center;color:var(--ink-3);font:700 13px/20px var(--mono)}
.exChip{display:inline-flex;align-items:center;gap:4px;height:22px;padding:0 9px;border-radius:11px;border:1px solid var(--line-2);font-size:12px}
.howto{margin:4px 0 12px;border:1px solid var(--line-2);border-radius:8px;background:var(--panel)}.howto>summary{cursor:pointer;padding:6px 12px;font-size:12.5px;font-weight:600;color:var(--ink-2)}
.howtoBody{max-height:40vh;overflow-y:auto;padding:2px 14px 10px;font-size:12.5px;line-height:1.5;color:var(--ink-2)}.howtoBody h4{margin:10px 0 2px;font-size:12px;color:var(--ink)}.howtoBody p{margin:0;max-width:80ch}
.sortLab{font-size:11px;color:var(--ink-3);margin-left:6px}
.consRow{display:flex;justify-content:center;align-items:center;gap:10px;flex-wrap:wrap;margin:6px 0 8px}
.cxBtn{appearance:none;height:20px;padding:0 8px;border-radius:10px;border:1px solid var(--line-2);background:transparent;font:500 11px/18px var(--sans);color:var(--ink-2);cursor:pointer;white-space:nowrap}.cxBtn:hover{border-color:var(--ink-3);color:var(--ink)}
.knownRow{display:flex;justify-content:center;align-items:center;gap:8px;flex-wrap:wrap;margin:4px 0 6px}
.kmBadge{display:inline-flex;align-items:center;gap:6px;height:22px;padding:0 9px;border-radius:11px;font-size:12px;font-weight:600}
.kmLikely{background:#111418;color:#fff;border:1px solid #111418}.kmPossible{background:transparent;color:var(--ink);border:1px dashed var(--ink-3)}
.kmTier{font-weight:400;font-size:10.5px;opacity:.75;text-transform:uppercase;letter-spacing:.04em}.kmWarn{color:#b8860b}.kmAlso{font-size:11.5px;color:var(--ink-3)}
.extlink{color:inherit;text-decoration:underline dotted;text-underline-offset:2px}.extlink:hover{text-decoration-style:solid}
a.kmBadge{text-decoration:none}a.kmBadge:hover{filter:brightness(1.15)}.kmExt{font-size:10px;opacity:.7}
.kmRegex{font-size:11px;color:var(--ink-2)}.kmRegex code{font:11px/1.4 var(--mono);color:var(--ink);background:var(--panel-2);border-radius:4px;padding:1px 5px}.kmCount{color:var(--ink-3)}
.segElm{background:rgba(184,134,11,.22);border-radius:2px}.cElm{font:9.5px/1.4 var(--mono);color:#b8860b;border:1px solid #b8860b;border-radius:3px;padding:0 3px;margin-right:6px;letter-spacing:0}
.info{position:relative;display:inline-flex;vertical-align:middle;margin-left:6px;cursor:help}
.info .q{width:16px;height:16px;border-radius:50%;border:1px solid var(--ink-3);color:var(--ink-3);font:700 10.5px/14px var(--sans);display:inline-flex;align-items:center;justify-content:center}
.info:hover .q,.info:focus .q{border-color:var(--ink);color:var(--ink)}
.infoPop{display:none;position:absolute;z-index:50;top:22px;left:-8px;width:max-content;max-width:340px;padding:9px 11px;border:1px solid var(--line-2);border-radius:8px;background:var(--panel);color:var(--ink-2);box-shadow:0 6px 22px rgba(0,0,0,.25);font:400 12.5px/1.5 var(--sans);text-transform:none;letter-spacing:normal;white-space:normal;text-align:left}
.info:hover .infoPop,.info:focus .infoPop,.info.pinned .infoPop{display:block}
.subRow{display:flex;align-items:center;gap:4px;flex-wrap:wrap;margin:2px 0 6px;font-size:11px}.subLabel{color:var(--ink3,#7c828b);margin-right:4px}
.subRow button{appearance:none;border:1px solid #d0d3d9;border-radius:10px;background:transparent;font:11px/1.6 inherit;padding:0 8px;cursor:pointer}
.subRow button[aria-pressed="true"]{background:#111418;border-color:#111418;color:#fff}
.ccells{display:flex;flex:none}.cell{display:inline-flex;justify-content:center;align-items:center;height:18px}
.cdot{color:#c3c7cd;font:11px/1 var(--mono)}
.cgap{display:inline-block;height:18px;background:repeating-linear-gradient(90deg,transparent 0 4px,#eef0f3 4px 5px)}
.aabox{display:inline-block;width:14px;text-align:center;font:700 10.5px/15px var(--mono);color:#111418;border-radius:2px}
.cseg{padding:0 10px;font:11px/1.4 var(--mono);color:var(--ink-3);letter-spacing:.5px;white-space:nowrap}
.segM{color:#111418;font-weight:800;text-decoration:underline;text-underline-offset:2px}.segUnk{color:#c3c7cd}
.kWide{flex-basis:100%;line-height:1.5}.kFlag{color:#8a5a00}
.kOther{display:inline-block;width:12px;height:4px;background:#c3c7cd;border-radius:1px;vertical-align:middle}
.kSig{display:inline-block;width:14px;text-align:center;background:#e6ebf2;border-bottom:2px solid #191c20;font-weight:700;color:#191c20}
.kUnt{display:inline-block;width:14px;height:12px;vertical-align:middle;background:repeating-linear-gradient(45deg,#d6d9de 0 2px,transparent 2px 5px);border:1px solid #d0d3d9}
.kLane{display:inline-block;width:18px;height:9px;border-radius:2px;background:#3b4552;vertical-align:middle}
.kSel{display:inline-block;width:12px;height:9px;border-radius:2px;background:var(--hl);border:1px solid #111418;vertical-align:middle}
details.pkey>summary{cursor:pointer;list-style:none}details.pkey>summary::-webkit-details-marker{display:none}
.pWrap{margin-top:10px}.pHead{display:flex;flex-wrap:wrap;gap:12px;align-items:center}
.pListToggle{background:none;border:0;padding:0;font-size:12.5px;font-weight:600;cursor:pointer;color:var(--ink-2)}
.pToggle{font-size:12px;color:var(--ink-3)}
.pRow{border-top:1px solid var(--line);padding:6px 4px;cursor:pointer}.pRow[data-on="1"]{background:var(--sel-fill)}
.pTop{display:flex;flex-wrap:wrap;gap:6px 14px;align-items:baseline}.pName{font-weight:650;font-size:13px}.pSub{font-size:11.5px;color:var(--ink-3)}
.pStats{display:flex;gap:10px;font-size:11.5px;color:var(--ink-3);margin-left:auto}
.mBadge{font:600 10.5px var(--mono);border:1px solid #d0d3d9;background:#eef0f3;border-radius:4px;padding:0 5px;cursor:pointer;color:#191c20}
.mBadge[aria-pressed="true"]{background:var(--hl);border-color:#111418}
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

LOGO_JS = r"""
// Logo maths -- the website's src/lib/logo.ts (2026-10-07); scrims.logo is the Python twin.
// No DOM at load time (tests run this file under node).
const LOG2_20 = Math.log2(20), MIN_GLYPH_PX = 3, LOGO_FONT = 'Arial, Helvetica, "Liberation Sans", sans-serif';
const _ink = new Map(); let _ctx;
// ink box of one letter at font-size 100, bold; a typical capital when there is no canvas
function inkBox(ch) {
  if (_ink.has(ch)) return _ink.get(ch);
  let ink = { left: -4, right: 68, ascent: 72, descent: 0 };
  try {
    if (_ctx === undefined) _ctx = typeof document !== 'undefined' ? document.createElement('canvas').getContext('2d') : null;
    if (_ctx) { _ctx.font = `700 100px ${LOGO_FONT}`; const m = _ctx.measureText(ch);
      if (m.actualBoundingBoxAscent + m.actualBoundingBoxDescent > 1)
        ink = { left: m.actualBoundingBoxLeft, right: m.actualBoundingBoxRight, ascent: m.actualBoundingBoxAscent, descent: m.actualBoundingBoxDescent }; }
  } catch (e) { /* keep the fallback */ }
  _ink.set(ch, ink); return ink;
}
// transform putting the ink of `ch` (<text font-size=100 x=0 y=0>) exactly into [x, x+w] x [yBottom-h, yBottom]:
// height exact; width 86 % of the column for an H-wide letter, wide letters squeezed to at most 94 %
function glyphTransform(ch, x, yBottom, w, h) {
  const ink = inkBox(ch), ref = inkBox('H'), inkW = Math.max(1, ink.left + ink.right), inkH = Math.max(1, ink.ascent + ink.descent);
  const sx = Math.min((w * 0.86) / Math.max(1, ref.left + ref.right), (w * 0.94) / inkW), sy = h / inkH, x0 = x + (w - inkW * sx) / 2;
  return `translate(${x0.toFixed(2)} ${yBottom.toFixed(2)}) scale(${sx.toFixed(4)} ${sy.toFixed(4)}) translate(${ink.left.toFixed(2)} ${(-ink.descent).toFixed(2)})`;
}
// letters: [letter, height px] summing to the column height -> visible (bottom -> top, smallest first) + one pooled band
function buildStack(letters, minPx = MIN_GLYPH_PX) {
  const vis = letters.filter(([, h]) => h >= minPx).sort((a, b) => a[1] - b[1]), small = letters.filter(([, h]) => h > 0 && h < minPx);
  return { visible: vis.map(([letter, h]) => ({ letter, h })), other: small.reduce((a, [, h]) => a + h, 0), otherLetters: small.map(([l]) => l) };
}
// exact small-sample correction e(n) = log2 20 - E[H_n], E[H_n] = expected entropy of n uniform letters from 20
const _eC = new Map();
function smallSampleCorrection(n) {
  if (!(n >= 1)) return LOG2_20;
  const k = Math.round(n); if (_eC.has(k)) return _eC.get(k);
  const p = 1 / 20; let pmf = Math.pow(1 - p, k), eh = 0;
  for (let c = 1; c <= k; c++) { pmf *= ((k - c + 1) / c) * (p / (1 - p)); const f = c / k; eh -= 20 * pmf * f * Math.log2(f); }
  const e = Math.max(0, LOG2_20 - eh); _eC.set(k, e); return e;
}
// information content (bits) of a frequency vector from n independent sequences, floored at 0
function infoContent(freqs, n) {
  if (!(n > 0)) return 0; const tot = freqs.reduce((a, b) => a + b, 0); if (!(tot > 0)) return 0;
  let H = 0; for (const f0 of freqs) { const f = f0 / tot; if (f > 0) H -= f * Math.log2(f); }
  return Math.max(0, LOG2_20 - (H + smallSampleCorrection(n)));
}
// every group weighs 1 in total, split over its members present; items: [letter, group id]
function groupWeighted(items) {
  const members = new Map(); for (const [, g] of items) members.set(g, (members.get(g) || 0) + 1);
  const n = members.size, w = new Map();
  for (const [l, g] of items) w.set(l, (w.get(l) || 0) + 1 / members.get(g) / n);
  return { freq: [...w.entries()].sort((a, b) => b[1] - a[1]), n };
}
// classic y axis per logo: tallest column rounded UP to the next 0.5 bit (min 0.5), capped at 4.5
function bitsAxisMax(ics) { const m = Math.max(0, ...ics); return Math.min(Math.max(0.5, Math.ceil(m * 2 - 1e-9) / 2), Math.ceil(LOG2_20 * 2) / 2); }
"""

JS = r"""
const D = JSON.parse(document.getElementById('scrims-data').textContent);
const P = D.params, SEQ = D.bait.sequence || '';
const IFC = ['#2a78d6','#d06f00','#d4609a','#a0662f','#1a9fc0','#d4462c','#7b61d9','#c2566b'];
// official Jalview Zappo (jalview.org/help/html/colourSchemes/zappo.html), as the website (lib/colors.ts):
// logo letters are FILLED with these + a thin dark outline; small letters are a coloured box + black letter.
const AAG = [['aliphatic / hydrophobic','ILVAM','#ffafaf'],['aromatic','FWY','#ffc800'],['positive','KRH','#6464ff'],
  ['negative','DE','#ff0000'],['hydrophilic','STNQ','#00ff00'],['conformationally special','PG','#ff00ff'],['cysteine','C','#ffff00']];
const AAC = {}; AAG.forEach(([,l,c]) => [...l].forEach(a => AAC[a] = c));
const aaColor = l => AAC[(l||'').toUpperCase()] || '#c3c7cd';
const LETTER_OUTLINE = '#1f2329';
// the SELECTION highlight: cyan (gold #ffc400 was indistinguishable from Zappo aromatic #ffc800);
// on white it always gets a dark #111418 outline
const HL = '#00ffff';
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
// UniProt accessions and ELM classes link out (new tab); the link never selects the row it sits in (2026-10-12, user)
const extA = (href, text, attrs = {}) => el('a', { href, target: '_blank', rel: 'noopener noreferrer', ...attrs,
  class: ['extlink', attrs.class].filter(Boolean).join(' '), on: { click: ev => ev.stopPropagation(), keydown: ev => ev.stopPropagation() } }, text);
const uniA = (acc, cls) => extA(`https://www.uniprot.org/uniprotkb/${encodeURIComponent(acc)}/entry`, acc, { class: cls, title: `${acc} on UniProt` });
const elmHref = id => `http://elm.eu.org/elms/${encodeURIComponent(id)}`;
// subOf as elements: the accession linked, the description (no gene) as text
const subEl = (id, cls) => { const i = idn(id); if (!i) return el('span', { class: cls }, '');
  return el('span', { class: cls }, i.uniprot ? uniA(i.uniprot) : null, !i.gene && i.description ? `${i.uniprot ? ' · ' : ''}${i.description}` : null); };
const paeAt = s => { const v = P.d0_signal * Math.sqrt(1 / s - 1); return Math.abs(v - Math.round(v)) < 0.05 ? String(Math.round(v)) : v.toFixed(1); };
const tierShort = { H: `≤ ${paeAt(P.partner_tau_high)} Å`, M: `${paeAt(P.partner_tau_high)}–${paeAt(P.partner_tau_in)} Å`, L: `> ${paeAt(P.partner_tau_in)} Å` };
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

// "?" info tip (2026-10-10, user: a ? per section instead of one "How to read" block): hover / focus shows it, click pins it
const info = text => { const w = el('span', { class: 'info', tabindex: 0, role: 'button', 'aria-label': 'What is this?',
  on: { click: ev => { ev.preventDefault(); ev.stopPropagation(); w.classList.toggle('pinned'); } } }, el('span', { class: 'q' }, '?'), el('span', { class: 'infoPop', role: 'tooltip' }, text)); return w; };
// display name of a motif id: "M1" -> "Motif 1" (2026-10-12, user's PI: "M1 is not obvious"); data ids stay M1
const motifLabel = id => { const m = /^M(\d+)$/.exec(id); return m ? `Motif ${m[1]}` : id; };
// the motif's significance: its strongest element's BH q against letter shuffles (2026-10-12, with the bits' rank)
const motifQ = (e, m) => { const qs = m.elements.map(k => e.elements[k] && e.elements[k].q).filter(q => typeof q === 'number' && isFinite(q)); return qs.length ? Math.min(...qs) : null; };
// a button that copies the ChimeraX command showing an interface / motif (2026-10-12, user's PI: scrims.cxc not obvious)
const cxBtn = cmd => { const b = el('button', { type: 'button', class: 'cxBtn', 'data-cx': cmd,
    title: `Copy the ChimeraX command: ${cmd}\n(first run scrims.cxc, in this report's folder, in ChimeraX)` }, 'ChimeraX ⧉');
  b.addEventListener('click', ev => { ev.preventDefault(); ev.stopPropagation();
    const done = () => { b.textContent = 'copied ✓'; setTimeout(() => { b.textContent = 'ChimeraX ⧉'; }, 1400); };
    const fallback = () => { const t = el('textarea', { style: { position: 'fixed', opacity: '0' } }); t.value = cmd; document.body.append(t); t.select();
      try { document.execCommand('copy'); done(); } catch (_) { window.prompt('ChimeraX command:', cmd); } t.remove(); };
    if (navigator.clipboard && window.isSecureContext) navigator.clipboard.writeText(cmd).then(done, fallback); else fallback(); });
  return b; };

// ------------------------------------------------------------------- header
function header() {
  const b = D.bait.identity || {};
  const nMot = D.interfaces.reduce((a, i) => a + ((i.motif && i.motif.motifs.length) || 0), 0);
  const binders = D.matrix.row_order.length;
  const tile = (v, k) => el('div', { class: 'tile' }, el('div', { class: 'tileV num' }, v), el('div', { class: 'tileK' }, k));
  return el('header', { class: 'phead' },
    el('div', {}, el('div', { class: 'eyebrow' }, 'Bait protein · SCRIMS report'),
      el('h1', { class: 'gene' }, b.gene || b.uniprot || D.screen.name),
      el('div', { class: 'idline' }, b.uniprot ? uniA(b.uniprot, 'mono') : null,
        el('span', {}, `${D.bait.length} aa`), b.description ? el('span', {}, b.description) : null, b.organism ? el('span', {}, b.organism) : null),
      el('div', { class: 'params' },
        D.screen.single_pair ? el('div', {}, 'Single prediction: recurrence not tested, strength regions only.') : null,
        Object.keys(D.screen.skipped || {}).length ? el('div', {}, `${Object.keys(D.screen.skipped).length} prediction(s) skipped (see summary.json).`) : null)),
    el('div', { class: 'tiles' }, tile(D.screen.n_predictions, 'predictions'), tile(binders, 'binders'),
      tile(D.interfaces.length, 'interfaces'), tile(nMot, 'motifs')));
}

// ------------------------------------------------------------------- side: selection
function side() {
  const box = el('div', { class: 'selbox' });
  const clear = el('button', { class: 'clear', on: { click: () => setSel({ iface: null, motif: null, carrier: null, partner: null }) } }, 'Clear selection');
  const draw = () => {
    box.replaceChildren();
    clear.disabled = S.iface === null && S.partner === null;
    if (S.iface === null && S.partner === null) {
      box.append(el('p', {}, 'Click anything to select it.'), el('p', {}, '3-D: open ', el('code', {}, 'scrims.cxc'), ' (in this report\'s folder) in ChimeraX, then use the ', el('b', {}, 'ChimeraX ⧉'), ' buttons.'));
      return;
    }
    if (S.iface !== null) {
      const f = D.interfaces[S.iface];
      box.append(el('p', {}, el('b', {}, f.id), ` · bait ${compact(f.residues)} · ${f.significant ? 'significant' : 'not significant'}`));
      if (S.motif) { const m = f.motif.motifs.find(x => x.id === S.motif);
        box.append(el('p', {}, el('b', {}, `${f.id} ${motifLabel(m.id)}`), ` · bait ${compact(m.positions)} · ${m.n_groups} groups · ${m.carriers.length} carriers`)); }
      box.append(el('p', {}, 'ChimeraX: ', el('code', {}, `view ${S.motif ? f.id + '_' + S.motif : f.id}`)));
    }
    if (S.partner) {
      const pu = (idn(S.partner) || {}).uniprot;
      box.append(el('p', {}, el('b', {}, nameOf(S.partner)), ' ', el('span', { class: 'mono' }, S.partner),
        pu && pu !== S.partner ? el('span', {}, ' · ', uniA(pu, 'mono')) : null, pu && pu === S.partner ? el('span', {}, ' ', extA(`https://www.uniprot.org/uniprotkb/${encodeURIComponent(pu)}/entry`, '↗', { title: `${pu} on UniProt` })) : null));
      if (S.carrier) { const f = D.interfaces[S.iface]; const al = f.motif.alignments[S.motif];
        const c = al.carriers.find(x => x.id === S.carrier);
        if (c) box.append(el('p', {}, `${motifLabel(S.motif)} residues on ${nameOf(c.id)}: ${compact(c.residues)} (${c.letters}), read ${c.dir === 1 ? 'N→C' : 'backwards, C→N'}`)); }
    }
  };
  onSel(draw); draw();
  return el('div', { class: 'card' }, el('div', { class: 'ctitle' }, 'Selection', info('Click anything to select it across the page. 3-D: scrims.cxc is a ChimeraX script written next to this report (with a models/ folder). Run “chimerax scrims.cxc”, or “open scrims.cxc” inside ChimeraX: it opens the structures, colours the interfaces and names every interface (I1) and motif (I1_M1) as a selection. Each card’s ChimeraX ⧉ button then copies the command that shows it.')), box, clear);
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
    if (focus !== null && focus >= 0) { const b = rows.indexOf(focus); if (b >= 0) { const r = [PAD_L - 1, matTop + b * rowH - 1, W - PAD_L - PAD_R + 2, Math.max(2, rowH) + 2];
      g.strokeStyle = '#111418'; g.lineWidth = 4; g.strokeRect(...r); g.strokeStyle = HL; g.lineWidth = 2; g.strokeRect(...r); } }
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
  return el('section', { class: 'card' }, el('div', { class: 'ctitle' }, 'Confidence matrix', info('One row per prey (clustered), one column per bait residue that some prey contacts. Colour = confidence that the pair is placed correctly (0.5 ≈ PAE 10 Å). The tracks above show backbone flexibility (RMSF) and recurrence; coloured bands mark the interfaces. Hover for values; click a row to select that prey, a band to select an interface.')),
    el('div', { class: 'plate' }, wrap, lane), legends);
}

// ------------------------------------------------------------------- page key
function pageKey() {
  const aa = el('div', { class: 'key' },
    ...AAG.map(([n, l, c]) => el('span', {}, el('span', { class: 'aak', style: { background: c } }, l), ' ', n)));
  const grouped = D.interfaces.every(f => !f.motif || f.motif.panel1.every(c => c.letters_grouped !== undefined));
  // 2026-10-10: marks only (explanations in "How to read this report"); no engine parameters
  const marks = el('div', { class: 'key' },
    el('span', {}, el('span', { class: 'kSig' }, '✱'), ' significant'),
    el('span', {}, el('span', { class: 'kUnt' }), ' untested'),
    el('span', {}, el('span', { class: 'kOther' }), ' letters too small to draw'),
    el('span', {}, el('span', { class: 'kLane' }), ' motif lane'),
    el('span', {}, el('b', {}, '←'), ' read the other way'),
    grouped ? null : el('span', { class: 'kFlag' }, 'older data: bait-side classic logo uses plain shares'));
  const tiers = el('div', { class: 'key' }, el('span', {}, 'partner PAE'),
    ...['H', 'M', 'L'].map(t => el('span', {}, el('span', { class: `tier tier${t}` }, t), ' ', tierShort[t])),
    el('span', {}, el('span', { class: 'range tierM' }, el('span', { class: 'mres' }, 'LL'), 'ML'), ' underlined = motif'));
  const det = el('details', { class: 'pkey', open: true }, el('summary', { class: 'hbtn' }, 'Key', info('Amino-acid colours: aliphatic ILVAM, aromatic FWY, positive KRH, negative DE, hydrophilic STNQ, special PG, cysteine C. ✱ = significant against letter shuffles (q < 0.05); hatched = untested (fewer than 3 partner groups); a grey band = letters too small to draw; ← = a partner read the other way round.')), aa, marks, tiers);
  return el('div', { class: 'plate', style: { marginBottom: '12px' } }, det);
}

// ------------------------------------------------------------------- motif engine view
const COL_W = 16, JUMP_W = 5, ROW_GROUPS = 12, ROW_SIG = 12, LOGO_H = 80, ROW_BAIT = 16, ROW_NUM = 26, LANE_H = 13, LANE_GAP = 4;
const PCOL_W = 18, GAP_W = 30, LEFT_W = 210, OTHER_FILL = '#c3c7cd', GREY_LETTER = '#d3d7dc';
// the two logo views, one switch per side for the whole page: classic on the partner side, evidence on the bait side
const VIEW = { baitClassic: false, partnerClassic: true };
const setView = patch => { Object.assign(VIEW, patch); listeners.forEach(f => f()); };
// one letter, its measured ink fitted exactly into [x, x+w] x [yBottom-h, yBottom]
const glyph = (ch, x, yBottom, w, h, fill) => sv('text', { class: 'glyph', transform: glyphTransform(ch, x, yBottom, w, h),
  'font-size': 100, 'font-weight': 700, fill, stroke: LETTER_OUTLINE, 'stroke-width': h < 8 ? 0.35 : 0.6, 'vector-effect': 'non-scaling-stroke' }, ch);
// a column: [letter, height px] (heights sum to the column height), baseline at yBase; visible letters
// contiguous from the baseline, largest on top, too-small shares pooled into one grey band on top
function stack(letters, x, yBase, w, grey) {
  const st = buildStack(letters), out = []; let y = yBase;
  st.visible.forEach(({ letter, h }) => { out.push(glyph(letter, x, y, w, h, grey ? GREY_LETTER : aaColor(letter))); y -= h; });
  if (st.other > 0.3) out.push(sv('rect', { x: x + 2, y: y - st.other, width: w - 4, height: st.other, fill: OTHER_FILL, 'data-other': 1 },
    sv('title', {}, `too small to draw: ${st.otherLetters.join(' ')}`)));
  return out;
}
function viewToggle(classic, key, what) {
  const b = (c, t) => el('button', { type: 'button', 'aria-pressed': classic === c ? 'true' : 'false', 'data-view': c ? 'classic' : 'evidence',
    on: { click: () => { if (VIEW[key] !== c) setView({ [key]: c }); } } }, t);
  return el('div', { class: 'vt', role: 'group', 'aria-label': `${what}: logo view` }, b(true, 'Classic logo (bits)'), b(false, 'Evidence view'));
}
const classicEmptyNote = () => el('div', { class: 'classicNote', 'data-classic-empty': 1 }, 'Classic logo empty: each column comes from a single partner group — see the Evidence view.');
const axisNote = (maxIc, axisMax) => el('div', { class: 'axisNote', 'data-axis-max': axisMax, title: `tallest column ${maxIc.toFixed(2)} bits; axis 0–${axisMax.toFixed(1)} bits` });
function bitsTicks(top, h, max) {
  const step = max <= 2.5 ? 0.5 : 1, out = [];
  for (let b = 0; b <= max + 1e-9; b += step) { const v = Math.round(b * 10) / 10;
    out.push(el('div', { class: 'tick', style: { top: Math.max(top - 1, top + h - (v / max) * h - 5) + 'px' } }, step < 1 ? v.toFixed(1) : String(v))); }
  return out;
}
const gl = (y, h, t, title) => el('div', { class: 'gl', title: title || null, style: { top: y + 'px', height: h + 'px', lineHeight: h + 'px' } }, t);
// partner letters by residue number, from the interface's partner ranges ('?' never guessed: absent = undefined)
const SEQS = new Map();
function seqOf(fi) {
  if (!SEQS.has(fi)) { const out = new Map();
    (D.interfaces[fi].partners || []).forEach(p => { const m = new Map(); p.ranges.forEach(([a, b, , s]) => { for (let r = a; r <= b; r++) m.set(r, s[r - a] || '?'); }); out.set(p.id, m); });
    SEQS.set(fi, out); }
  return SEQS.get(fi);
}
const activeMotif = fi => { const e = D.interfaces[fi].motif; return S.iface === fi && S.motif && e ? e.motifs.find(m => m.id === S.motif) || null : null; };
// the bait-side columns: all of panel1, or -- with a motif picked -- only its positions, letters from ITS carriers
function baitCols(fi, active) {
  const e = D.interfaces[fi].motif;
  const norm = o => { const t = Object.values(o).reduce((a, b) => a + b, 0) || 1; return Object.entries(o).map(([l, v]) => [l, v / t]).sort((a, b) => b[1] - a[1]); };
  if (!active) return e.panel1.map(c => ({ pos: c.pos, nGroups: c.n_groups, q: c.q, agree: c.agree, shares: norm(c.letters),
    freq: norm(c.letters_grouped !== undefined ? c.letters_grouped : c.letters), n: c.n_groups }));
  const al = e.alignments[active.id], seqs = seqOf(fi), p1 = new Map(e.panel1.map(c => [c.pos, c]));
  return [...active.positions].sort((a, b) => a - b).map(pos => {
    const items = [];
    al.carriers.forEach(c => { const r = c.map[pos]; if (r !== undefined) { const s = seqs.get(c.id); items.push([(s && s.get(r)) || '?', c.group]); } });
    const counts = new Map(); items.forEach(([l]) => counts.set(l, (counts.get(l) || 0) + 1));
    const gw = groupWeighted(items), pc = p1.get(pos);
    return { pos, nGroups: gw.n, q: pc ? pc.q : null, agree: pc ? pc.agree : null,
      shares: [...counts].map(([l, k]) => [l, k / Math.max(1, items.length)]).sort((a, b) => b[1] - a[1]), freq: gw.freq, n: gw.n };
  });
}
function baitSide(fi) {
  const f = D.interfaces[fi], e = f.motif, classic = VIEW.baitClassic, active = activeMotif(fi), cols = baitCols(fi, active);
  const xs = []; let x = 0; cols.forEach((c, i) => { if (i > 0 && c.pos - cols[i - 1].pos > 1) x += JUMP_W; xs.push(x); x += COL_W; });
  const maxG = Math.max(1, ...cols.map(c => c.nGroups)), width = (xs.length ? xs[xs.length - 1] + COL_W : COL_W) + 2;
  const yLogo = ROW_GROUPS + ROW_SIG, yBait = yLogo + LOGO_H, yNum = yBait + ROW_BAIT, yLanes = yNum + ROW_NUM + 4;
  const height = yLanes + e.motifs.length * (LANE_H + LANE_GAP);
  const ics = cols.map(c => infoContent(c.freq.map(([, v]) => v), c.n)), axisMax = bitsAxisMax(ics);
  const classicEmpty = classic && cols.length > 0 && cols.every(c => c.n <= 1);
  const toggle = m => setSel(active && active.id === m.id ? { iface: fi, motif: null, carrier: null } : { iface: fi, motif: m.id, carrier: null });
  const gut = el('div', { class: 'gut', style: { height: height + 'px' } },
    gl(0, ROW_GROUPS, 'groups', 'independent partner groups touching the residue (related partners count once)'),
    gl(ROW_GROUPS, ROW_SIG, `✱ q < ${Q_SIG}`, 'significant against letter shuffles (BH q, all partners)'),
    ...(classic ? [el('div', { class: 'glAxis', style: { top: yLogo + 'px', height: LOGO_H + 'px', lineHeight: LOGO_H + 'px' } }, 'bits'), ...bitsTicks(yLogo, LOGO_H, axisMax)]
      : [gl(yLogo, LOGO_H, 'partner residues', 'column height ∝ √ groups (busiest column = full); letters = shares of the partner residues there')]),
    gl(yBait, ROW_BAIT, 'bait residue'), gl(yNum, ROW_NUM, 'position'));
  const svg = sv('svg', { width, height, 'data-bait-mode': classic ? 'classic' : 'evidence', role: 'img',
    'aria-label': `Bait-side ${classic ? 'classic logo' : 'evidence view'}${active ? ` of ${motifLabel(active.id)}` : ''}: ${cols.length} bait residues` });
  svg.append(sv('defs', {}, sv('pattern', { id: `hatch${fi}`, width: 5, height: 5, patternUnits: 'userSpaceOnUse', patternTransform: 'rotate(45)' },
    sv('line', { x1: 0, y1: 0, x2: 0, y2: 5, stroke: '#d6d9de', 'stroke-width': 2 }))));
  if (classic) svg.append(sv('line', { x1: 0, x2: width, y1: yLogo + LOGO_H, y2: yLogo + LOGO_H, stroke: '#d0d3d9', 'stroke-width': 1 }));
  cols.forEach((c, i) => {
    const X = xs[i], untested = c.q === null, sig = !untested && c.q < Q_SIG, hEvid = LOGO_H * Math.sqrt(c.nGroups / maxG);
    const ic = ics[i], colH = classic ? (ic / axisMax) * LOGO_H : hEvid, src = classic ? c.freq : c.shares;
    const bait = SEQ[c.pos - 1] || '?';
    const g = sv('g', { 'data-col': c.pos, 'data-ic': ic.toFixed(3), 'data-colh': colH.toFixed(2) },
      sv('title', {}, `bait ${bait}${c.pos} · ${c.nGroups} group${c.nGroups === 1 ? '' : 's'}${active ? ` among ${active.id}'s carriers` : ''} · ` +
        (untested ? 'untested (< 3 partner groups)' : `q ${fmtQ(c.q)}${sig ? ' (significant)' : ''}`) + (c.agree !== null ? ` · agreement ${c.agree.toFixed(2)}` : '') +
        (classic ? ` · ${ic.toFixed(2)} bits` : '') + `\n${classic ? 'group-weighted frequencies' : 'partner residues'}: ${src.map(([l, v]) => `${l} ${Math.round(v * 100)}%`).join(', ')}`));
    if (sig) g.append(sv('rect', { x: X, y: yLogo - 1, width: COL_W, height: LOGO_H + 2, fill: '#e6ebf2', rx: 2 }));
    if (untested) g.append(sv('rect', { x: X, y: yLogo, width: COL_W, height: LOGO_H, fill: `url(#hatch${fi})` }));
    if (sig) g.append(sv('rect', { x: X + 1, y: yLogo + LOGO_H, width: COL_W - 2, height: 2, fill: '#191c20' }));
    g.append(sv('text', { x: X + COL_W / 2, y: ROW_GROUPS - 3, 'text-anchor': 'middle', 'font-size': 9, 'font-weight': 600, fill: untested ? '#a3a9b1' : '#4a5058' }, c.nGroups));
    if (sig) g.append(sv('text', { x: X + COL_W / 2, y: ROW_GROUPS + ROW_SIG - 2, 'text-anchor': 'middle', 'font-size': 10, 'font-weight': 700, fill: '#191c20' }, '✱'));
    g.append(...stack(src.map(([l, v]) => [l, v * colH]), X, yLogo + LOGO_H, COL_W, untested));
    // the bait's own residue, Jalview-style: Zappo box, black letter
    g.append(sv('rect', { x: X + 1.5, y: yBait + 2, width: COL_W - 3, height: 12, rx: 2, fill: aaColor(bait) }));
    g.append(sv('text', { x: X + COL_W / 2, y: yBait + 12, 'text-anchor': 'middle', 'font-size': 11, 'font-weight': 700, fill: '#111418', style: 'font-family:var(--mono)' }, bait));
    g.append(sv('text', { transform: `translate(${X + COL_W / 2 + 3} ${yNum + 2}) rotate(90)`, 'font-size': 8.5, fill: '#7c828b', style: 'font-family:var(--mono)' }, c.pos));
    g.append(sv('rect', { x: X, y: 0, width: COL_W, height: yNum + ROW_NUM, fill: 'transparent' }));
    svg.append(g);
  });
  const colOf = new Map(cols.map((c, i) => [c.pos, i]));
  e.motifs.forEach((m, li) => {
    const y = yLanes + li * (LANE_H + LANE_GAP), idx = m.positions.map(p => colOf.get(p)).filter(v => v !== undefined), on = !!active && active.id === m.id;
    gut.append(el('button', { class: 'laneBtn', 'data-lane-btn': m.id, 'aria-pressed': on ? 'true' : 'false', title: on ? `${m.id}: show every residue again` : `${m.id}: show only this motif`,
      style: { top: y + 'px', height: LANE_H + 'px', lineHeight: (LANE_H - 2) + 'px' }, on: { click: () => toggle(m) } }, motifLabel(m.id)));
    let g;
    if (!idx.length) {      // filtered to another motif and no shared residue: an empty, still clickable lane
      g = sv('g', { style: 'cursor:pointer', 'data-lane': m.id }, sv('title', {}, `${motifLabel(m.id)}: not among the shown residues — click to show it`),
        sv('rect', { x: 0, y, width, height: LANE_H, fill: 'transparent' }),
        sv('line', { x1: 0, x2: width, y1: y + LANE_H / 2, y2: y + LANE_H / 2, stroke: '#e2e4e8', 'stroke-dasharray': '3 3' }));
    } else {
      const x0 = xs[Math.min(...idx)], x1 = xs[Math.max(...idx)] + COL_W;
      g = sv('g', { style: 'cursor:pointer', 'data-lane': m.id, 'data-on': on ? 1 : 0 }, sv('title', {}, `${motifLabel(m.id)}: bait ${compact(m.positions)} — click to ${on ? 'show every residue again' : 'show only it'}`),
        sv('rect', { x: x0, y: y - 1, width: x1 - x0, height: LANE_H + 2, fill: 'transparent' }),
        sv('line', { x1: x0, x2: x1, y1: y + LANE_H / 2, y2: y + LANE_H / 2, stroke: '#9aa1aa', 'stroke-width': 1 }),
        ...idx.map(ci => sv('rect', { x: xs[ci] + 0.5, y, width: COL_W - 1, height: LANE_H, rx: 2, fill: on ? HL : '#3b4552', stroke: on ? '#111418' : 'none', 'stroke-width': on ? 1.5 : 0 })));
    }
    g.addEventListener('click', () => toggle(m)); svg.append(g);
  });
  const nCar = active ? e.alignments[active.id].carriers.length : 0;
  return [active ? el('div', { class: 'filterBanner', 'data-filter-banner': 1 }, 'showing ', el('b', {}, motifLabel(active.id)),
      ` only · ${cols.length} residue${cols.length === 1 ? '' : 's'} · ${nCar} carrier${nCar === 1 ? '' : 's'} (letters from its carriers)`,
      el('button', { type: 'button', class: 'showAll', on: { click: () => setSel({ iface: fi, motif: null, carrier: null }) } }, 'show all')) : null,
    classicEmpty ? classicEmptyNote() : null, classic && !classicEmpty ? axisNote(Math.max(0, ...ics), axisMax) : null,
    el('div', { class: 'bs' }, gut, el('div', { class: 'scroll' }, svg))];
}
// one motif: the partner-sequence logo over the per-carrier table, in the same columns
function motifPanel(fi, m) {
  const f = D.interfaces[fi], e = f.motif, alBase = e.alignments[m.id], seqs = seqOf(fi);
  // always N->C (viz_int 2026-10-09): a binding sub-group whose own majority reads backwards is shown mirrored
  const mirrorAl = a => ({ ...a, flipped: !a.flipped,
    carriers: a.carriers.map(c => ({ ...c, dir: -c.dir, ref: c.ref ? [c.ref[0], -c.ref[1]] : null,
      cells: Object.fromEntries(Object.entries(c.cells || {}).map(([o, l]) => [String(-Number(o)), l])) })),
    columns: [...a.columns].reverse().map(c => c.kind === 'col' ? { ...c, offset: -c.offset } : c) });
  let al = alBase;
  const det = el('details', { class: 'motif', 'data-motif': m.id });
  const n1 = (k, w) => `${k} ${w}${k === 1 ? '' : 's'}`;
  det.append(el('summary', {}, el('b', {}, motifLabel(m.id)), el('span', { class: 'mmeta' }, `${m.positions.length} bait res · ${n1(m.n_groups, 'group')} · ${n1(m.carriers.length, 'carrier')}`), info('The partner side of this motif: one row per partner, aligned so equivalent residues share a column, always read N→C. "antiparallel" = its partners run against the bait; ← = a partner read the other way round. Classic logo = information in bits (every partner counts); evidence view = how many partners have a residue at each position. The bits here rate the motif: the information of its 5 most conserved positions, related partners counted once; “top N%” = its rank among 51 reference motifs (31 proteins); q = its significance against letter shuffles.'),
    alBase.quality !== undefined ? el('span', { class: 'mmeta', 'data-quality': alBase.quality, 'data-quality-top': alBase.quality_top,
      title: `motif quality: information of its 5 most conserved positions, related partners counted once${alBase.quality_top !== undefined ? `; top ${alBase.quality_top}% of the 51 reference motifs (31 proteins; median 4.3 bits)` : ''}` },
      `· ${alBase.quality.toFixed(1)} bits${alBase.quality_top !== undefined ? ` · top ${alBase.quality_top}%` : ''}`) : null,
    motifQ(e, m) !== null ? el('span', { class: 'mmeta', 'data-motif-q': motifQ(e, m), title: "the motif's significance: its strongest element against letter shuffles (Benjamini–Hochberg q)" }, `· q ${motifQ(e, m).toPrecision(1)}`) : null,
    alBase.flipped ? el('span', { class: 'badge', 'data-flipped': 1, title: 'its partners run against the bait (C→N along it); shown N→C' }, 'antiparallel') : null,
    m.bipartite ? el('span', { class: 'badge', title: 'two elements of this motif are carried together at a linker distance (bipartite)' }, 'bipartite') : null,
    ...m.linked_by.map(l => el('span', { class: 'badge', title: l === 'carriers' ? 'elements joined because the same partner groups carry them beyond chance'
      : 'elements joined because they share bait residues with similar preferred amino acids (BLOSUM62 ≥ 1)' }, `linked by ${l}`)), cxBtn(`view ${f.id}_${m.id}; select ${f.id}_${m.id}`)));
  const body = el('div', { class: 'mbody' }); det.append(body);
  const st = { table: false, lastCarrier: null, sub: null };
  // paralog groups numbered in order of first appearance; carriers ordered so a group sits together
  const groupNo = new Map(), groupSize = new Map();
  al.carriers.forEach(c => { if (!groupNo.has(c.group)) groupNo.set(c.group, groupNo.size + 1); groupSize.set(c.group, (groupSize.get(c.group) || 0) + 1); });
  const subs = alBase.subgroups || [];
  // per offset, over the shown carriers (a binding sub-group or all): evidence counts, independent groups, and the
  // classic frequencies -- in the classic view EVERY CARRIER COUNTS (viz_int 2026-10-08: paralogs counted once hid
  // VAPA FFAT's second F; the evidence view and each column's tooltip keep the independent-group count)
  const colsFor = shown => al.columns.map(c => { if (c.kind === 'gap') return null; const byCar = [], byGrp = [], cnt = new Map();
    shown.forEach(car => { const l = (car.cells || {})[c.offset]; if (!l) return; byCar.push([l, car.id]); byGrp.push([l, car.group]); cnt.set(l, (cnt.get(l) || 0) + 1); });
    const gw = groupWeighted(byCar);
    return { gw, nGroups: groupWeighted(byGrp).n, k: byCar.length, counts: [...cnt].sort((a, b) => b[1] - a[1]), ic: infoContent(gw.freq.map(([, v]) => v), gw.n) }; });
  const carrierOn = id => S.iface === fi && S.motif === m.id && S.carrier === id;
  const pick = id => setSel(carrierOn(id) ? { carrier: null, partner: null } : { iface: fi, motif: m.id, carrier: id, partner: id });
  const fill = () => {
    al = st.sub !== null && (alBase.subgroup_flip || [])[st.sub] ? mirrorAl(alBase) : alBase;
    const allRows = [...al.carriers].sort((a, b) => groupNo.get(a.group) - groupNo.get(b.group));
    const xs = []; let acc = 0; al.columns.forEach(c => { xs.push(acc); acc += c.kind === 'gap' ? GAP_W : PCOL_W; });
    const width = acc + 2;
    const subSet = st.sub !== null && subs[st.sub] ? new Set(subs[st.sub]) : null;
    const shown = subSet ? al.carriers.filter(c => subSet.has(c.id)) : al.carriers, rows = allRows.filter(c => !subSet || subSet.has(c.id));
    const nShown = shown.filter(c => c.aligned).length, colData = colsFor(shown);
    const pIcs = colData.map(d => d ? d.ic : 0), pAxisMax = bitsAxisMax(pIcs);
    const classic = VIEW.partnerClassic, PLOGO_H = classic ? 80 : Math.max(48, Math.min(96, nShown * 8));
    const yLogo = 12, yOff = yLogo + PLOGO_H + 12, height = yOff + 4, nFl = al.carriers.filter(c => c.dir === -1).length;
    body.replaceChildren();
    body.append(el('div', { class: 'mtop' }, viewToggle(classic, 'partnerClassic', `${motifLabel(m.id)} partner side`)));
    // consensus (A) over the carriers on show (viz_int 2026-10-09)
    const cons = st.sub !== null ? (alBase.subgroup_consensus || [])[st.sub] : alBase.consensus;
    // 2026-10-10: consensus as coloured residue boxes, best example as a chip (website-identical design)
    const consTok = tok => tok === 'x' ? el('span', { class: 'cDot' }, '·') : tok === '…' ? el('span', { class: 'cDot' }, '…')
      : el('span', { class: 'cBox' }, ...(tok.startsWith('[') ? tok.slice(1, -1).split('') : [tok]).map((l, i, a) =>
        el('span', { class: 'cPart', style: { background: aaColor(l), width: (20 / a.length) + 'px' } }, l)));
    if (subs.length) body.append(el('div', { class: 'subRow', role: 'group', 'aria-label': `${motifLabel(m.id)}: binding sub-groups`,
      title: 'carriers that touch the same bait residues with the same aligned residues (geometry only)' },
      el('span', { class: 'subLabel' }, 'binding sub-groups'),
      el('button', { type: 'button', 'aria-pressed': st.sub === null ? 'true' : 'false', 'data-subgroup': 'all', on: { click: () => { st.sub = null; fill(); } } }, `all ${al.n}`),
      ...subs.map((g, i) => el('button', { type: 'button', 'aria-pressed': st.sub === i ? 'true' : 'false', 'data-subgroup': i + 1,
        on: { click: () => { st.sub = st.sub === i ? null : i; fill(); } } }, `B${i + 1} · ${g.length}`))));
    const hasCol = colData.some(d => d !== null);
    if (classic && hasCol && colData.every(d => d === null || d.gw.n <= 1)) body.append(classicEmptyNote());
    else if (classic) body.append(axisNote(Math.max(0, ...pIcs), pAxisMax));
    const left = el('div', { class: 'pleft gut', style: { height: height + 'px', width: LEFT_W + 'px' } },
      ...(classic ? [gl(0, 12, 'carriers', 'carriers with a motif residue at this offset (in the classic view every carrier counts)'),
        el('div', { class: 'glAxis', style: { top: yLogo + 'px', height: PLOGO_H + 'px', lineHeight: PLOGO_H + 'px' } }, 'bits'), ...bitsTicks(yLogo, PLOGO_H, pAxisMax)]
        : [gl(0, 12, `k of ${nShown}`, 'carriers with a motif residue at this offset'), gl(yLogo, PLOGO_H, 'partner residues')]),
      gl(yLogo + PLOGO_H, 14, 'offset', `offset along the partner sequence; 0 = the consensus residue on bait ${al.anchor}`));
    const svg = sv('svg', { width, height, 'data-partner-mode': classic ? 'classic' : 'evidence', role: 'img',
      'aria-label': `${motifLabel(m.id)} partner-side ${classic ? 'classic logo' : 'evidence view'}: ${al.n} aligned carriers, anchored on bait residue ${al.anchor}` });
    if (classic) svg.append(sv('line', { x1: 0, x2: width, y1: yLogo + PLOGO_H, y2: yLogo + PLOGO_H, stroke: '#d0d3d9', 'stroke-width': 1 }));
    al.columns.forEach((c, i) => { const x = xs[i];
      if (c.kind === 'gap') { svg.append(sv('g', {}, sv('title', {}, `${c.skipped} partner residues where no carrier has a motif residue (collapsed)`),
        sv('text', { x: x + GAP_W / 2, y: yLogo + PLOGO_H - 4, 'text-anchor': 'middle', 'font-size': 13, 'font-weight': 700, fill: '#7c828b' }, '⋯'),
        sv('text', { x: x + GAP_W / 2, y: yOff, 'text-anchor': 'middle', 'font-size': 8.5, fill: '#7c828b', style: 'font-family:var(--mono)' }, `${c.skipped} aa`),
        sv('rect', { x, y: 0, width: GAP_W, height, fill: 'transparent' }))); return; }
      const d = colData[i], colH = classic ? (d.ic / pAxisMax) * PLOGO_H : (d.k / Math.max(1, nShown)) * PLOGO_H;
      const letters = classic ? d.gw.freq.map(([l, v]) => [l, v * colH]) : d.counts.map(([l, n]) => [l, (n / Math.max(1, nShown)) * PLOGO_H]);
      const g = sv('g', { 'data-pcol': c.offset, 'data-ic': d.ic.toFixed(3), 'data-colh': colH.toFixed(2) },
        sv('title', {}, `offset ${c.offset > 0 ? '+' : ''}${c.offset} · ${d.k}/${nShown} carriers · ${n1(d.nGroups, 'group')}` + (classic ? ` · ${d.ic.toFixed(2)} bits` : '') +
          (d.counts.length ? `\n${d.counts.map(([l, n]) => `${l} ${n}`).join(', ')}` : '')));
      g.append(sv('text', { x: x + PCOL_W / 2, y: 9, 'text-anchor': 'middle', 'font-size': 9, 'font-weight': 600, fill: d.k ? '#4a5058' : '#c3c7cd' }, d.k || '·'));
      g.append(...stack(letters, x, yLogo + PLOGO_H, PCOL_W, false));
      g.append(sv('text', { x: x + PCOL_W / 2, y: yOff, 'text-anchor': 'middle', 'font-size': c.offset === 0 ? 9 : 8.5, 'font-weight': c.offset === 0 ? 700 : 400, fill: '#7c828b', style: 'font-family:var(--mono)' }, c.offset > 0 ? `+${c.offset}` : c.offset));
      g.append(sv('rect', { x, y: 0, width: PCOL_W, height, fill: 'transparent' })); svg.append(g); });
    // known ELM class the partners match beyond chance, centred just above the logo (website-identical; 2026-10-12)
    const cands = m.candidates || [], elmHits = (cands[0] && cands[0].matches) || {};
    if (cands.length) {
      const c0 = cands[0], why = c => `${c.tier}: ${c.name.trim()} (${c.elm}) — ${c.groups_matched} of ${c.groups} partner groups match its ELM pattern, ≈${c.expected} expected by chance (p ${c.p.toPrecision(2)}); ` +
        (c.domain === 'supported' ? 'the bait has the domain ELM lists for this class'
          : c.domain === 'not found' ? 'the bait LACKS the domain ELM lists for this class — treat with caution'
          : 'no domain annotation to compare (--bait-domains)');
      // option C (2026-10-12, user): the badge links to the ELM class page, its pattern sits next to it, and the carrier
      // table marks the carriers whose sequence matches it
      const nHit = al.carriers.filter(c => elmHits[c.id]).length;
      body.append(el('div', { class: 'knownRow' },
        el('a', { href: elmHref(c0.elm), target: '_blank', rel: 'noopener noreferrer', 'data-elm-link': c0.elm,
          class: `kmBadge ${c0.tier === 'likely' ? 'kmLikely' : 'kmPossible'}`, 'data-candidate': c0.elm, 'data-tier': c0.tier,
          'data-domain': c0.domain, title: cands.map(why).join('\n') },
          c0.domain === 'not found' ? el('span', { class: 'kmWarn', 'aria-label': 'bait lacks its domain' }, '⚠') : null,
          c0.name.trim(), el('span', { class: 'kmTier' }, c0.tier), el('span', { class: 'kmExt', 'aria-hidden': 'true' }, '↗')),
        c0.regex ? el('span', { class: 'kmRegex', 'data-elm-regex': 1, title: `ELM pattern of ${c0.elm}; ${nHit} of ${al.carriers.length} carriers match it (shaded in the carrier table)` },
          el('code', {}, c0.regex), ' ', el('span', { class: 'kmCount' }, `${nHit}/${al.carriers.length} carriers`)) : null,
        cands.length > 1 && cands[1].name.trim() !== c0.name.trim()
          ? extA(elmHref(cands[1].elm), `also ${cands[1].name.trim()}`, { class: 'kmAlso', title: `${cands[1].tier}: ${cands[1].name.trim()} (${cands[1].elm}) — ${cands[1].groups_matched} of ${cands[1].groups} groups${cands[1].regex ? `; pattern ${cands[1].regex}` : ''}` }) : null,
        info('Known motif: an ELM class whose pattern the partners here match more often than shuffled versions of the same residues (likely = 5 or more independent groups and p ≤ 1e-6; possible = fewer). ⚠ = the bait lacks the domain ELM lists for it.')));
    }
    const tab = el('div', { class: 'ptab' }, el('div', { class: 'prow' }, left, svg));
    // consensus centred just below the logo (2026-10-10, user); best example removed (user)
    if (cons && cons.length) tab.append(el('div', { class: 'consRow', 'data-motif-summary-line': 1 },
      el('span', { class: 'consBoxes', 'data-consensus': cons.map(t => t[1]).join('-'), title: `consensus ${cons.map(t => t[1]).join('-')}` },
        ...cons.slice(0, 24).map(t => consTok(t[1])), cons.length > 24 ? el('span', { class: 'cDot' }, '…') : null), info('Consensus: a box = residue shared by most partners; a split box = two or three similar residues; a dot = variable.')));
    tab.append(el('button', { type: 'button', class: 'tableToggle', 'aria-expanded': st.table ? 'true' : 'false', 'data-carrier-table-toggle': 1,
      on: { click: () => { st.table = !st.table; fill(); } } }, `${st.table ? '▾' : '▸'} Carriers (${shown.length})`));
    if (st.table) {
      const grp = el('div', { role: 'group', 'aria-label': `${motifLabel(m.id)} carriers — select one to light it across the page` });
      rows.forEach(car => {
        const lo = car.residues[0], hi = car.residues[car.residues.length - 1], mres = new Set(car.residues), sq = seqs.get(car.id), seg = [];
        // the known class's matched stretch (ELM pattern) is shaded; the sequence widens to show all of it (option C, 2026-10-12)
        const hit = elmHits[car.id], inHit = r => !!hit && r >= hit[0] && r <= hit[1];
        for (let r = hit ? Math.min(lo, hit[0]) : lo; r <= (hit ? Math.max(hi, hit[1]) : hi); r++) { const l = sq && sq.get(r), x = inHit(r) ? ' segElm' : '';
          seg.push(l === undefined ? el('span', { class: 'segUnk' + x, title: `${r}: outside the partner's listed ranges` }, '·') : mres.has(r) ? el('b', { class: 'segM' + x }, l) : x ? el('span', { class: 'segElm' }, l) : l); }
        const gsz = groupSize.get(car.group) || 1, uni = (idn(car.id) || {}).uniprot, on = carrierOn(car.id);
        const cells = al.columns.map(c => { if (c.kind === 'gap') return el('span', { class: 'cgap', style: { width: GAP_W + 'px' } });
          const l = (car.cells || {})[c.offset];
          return l ? el('span', { class: 'cell', style: { width: PCOL_W + 'px' } }, el('span', { class: 'aabox', style: { background: aaColor(l) } }, l))
            : el('span', { class: 'cell cdot', style: { width: PCOL_W + 'px' } }, '·'); });
        grp.append(el('div', { class: 'crow', role: 'button', tabindex: 0, 'data-carrier': car.id, 'aria-pressed': on ? 'true' : 'false',
          title: `${nameOf(car.id)} · ${car.id} — click to select it (its ${motifLabel(m.id)} residues: ${compact(car.residues)})`,
          on: { click: () => pick(car.id), keydown: ev => { if (ev.key === 'Enter' || ev.key === ' ') { ev.preventDefault(); pick(car.id); } } } },
          el('div', { class: 'cleft' }, el('span', { class: 'cname' }, nameOf(car.id) === uni ? uniA(uni) : nameOf(car.id)), uni && nameOf(car.id) !== uni ? uniA(uni, 'cuni') : null,
            gsz > 1 ? el('span', { class: 'cgrp', title: `related partners: these ${gsz} carriers are one independent group and count once in the classic logo` }, `g${groupNo.get(car.group)}`) : null,
            el('span', { class: 'carrow', title: car.dir === 1 ? 'read N→C (parallel)' : 'read backwards, C→N (antiparallel, flipped)' }, car.dir === 1 ? '→' : '←'),
            car.aligned ? null : el('span', { class: 'cun' }, 'not in logo'),
            el('span', { class: 'crange' }, lo === hi ? lo : `${lo}–${hi}`)),
          el('div', { class: 'ccells' }, ...cells), el('div', { class: 'cseg' },
            hit ? el('span', { class: 'cElm', 'data-elm-hit': `${hit[0]}-${hit[1]}`, title: `matches the ${cands[0].name.trim()} pattern (${cands[0].elm}) at residues ${hit[0]}–${hit[1]} (shaded)` }, 'ELM') : null, ...seg)));
      });
      tab.append(grp);
    }
    body.append(tab);
  };
  det.addEventListener('toggle', () => { if (det.open) fill(); });
  onSel(() => {
    // a carrier picked elsewhere (partner-list badge) opens this panel and its table, and scrolls to the row
    const mine = S.iface === fi && S.motif === m.id && S.carrier && al.carriers.some(c => c.id === S.carrier) ? S.carrier : null;
    const fresh = mine && mine !== st.lastCarrier; st.lastCarrier = mine;
    if (fresh) { st.table = true; if (!det.open) { det.open = true; } }
    if (det.open) fill();
    if (fresh) setTimeout(() => { const r = body.querySelector(`[data-carrier="${CSS.escape(mine)}"]`); if (r) r.scrollIntoView({ block: 'nearest', inline: 'nearest' }); }, 60);
  });
  return det;
}
function motifView(fi) {
  const f = D.interfaces[fi], e = f.motif;
  if (!e) return el('div', { class: 'plate mplate' }, el('div', { class: 'summary' }, 'No motif-engine entry for this interface: none of its partners gave the engine a contact.'));
  const p1 = e.panel1, nT = p1.filter(c => c.q !== null).length, nS = p1.filter(c => c.q !== null && c.q < Q_SIG).length;
  const st = { open: true }, head = el('div', { class: 'mhead' }), bsBox = el('div');
  const draw = () => {
    head.replaceChildren(el('button', { type: 'button', class: 'headBtn', 'aria-expanded': st.open ? 'true' : 'false', on: { click: () => { st.open = !st.open; draw(); } } },
      `${st.open ? '▾' : '▸'} Interface – Bait-side view`), info('One column per bait residue touched by partners. Height grows with the number of independent partner groups (related partners count once); ✱ = significant against letter shuffles (q < 0.05); hatched = untested (fewer than 3 groups). The letters are the partner residues found there. The lanes below are the motifs: click one to show only its residues.'),
      el('span', { class: 'hmeta', title: `${p1.length} bait residues · ${nT} tested · ${nS} with q < ${Q_SIG} · partners used ${e.partners_used}/${e.partners_total} · shuffle control: ${e.control.false_elements} false element${e.control.false_elements === 1 ? '' : 's'}, columns ${e.control.false_columns.join('/')}` }, `${nS} of ${nT} tested residues significant`),
      st.open && p1.length ? viewToggle(VIEW.baitClassic, 'baitClassic', 'bait side') : null);
    bsBox.replaceChildren(...(!st.open ? [] : p1.length ? baitSide(fi).filter(Boolean) : [el('div', { class: 'summary' }, 'No partner residues here.')]));
  };
  draw(); onSel(draw);
  const rel = e.motifs.length === 0 ? el('div', { class: 'summary', 'data-motif-summary': 1 }, 'No motif found')
    : e.motifs.length > 1 ? el('div', { class: 'summary', 'data-motif-summary': 1 }, e.motif_relations.map(r => `${motifLabel(r.a)}–${motifLabel(r.b)} ${r.relation === 'same_pocket' ? 'same pocket' : 'separate sites'} · ${r.status === 'different' ? 'different (competing)' : 'unresolved, too few carriers to tell'}`).join(' · ')) : null;
  return el('div', { class: 'plate mplate', 'data-motif-view': f.id }, head, bsBox, rel, ...e.motifs.map(m => motifPanel(fi, m)));
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
        el('div', { class: 'pTop' }, el('span', { class: 'pName' }, nameOf(p.id)), subEl(p.id, 'pSub mono'),
          ...ms.map(x => { const on = S.iface === fi && S.motif === x.id && S.carrier === p.id;
            return el('button', { class: 'mBadge', 'aria-pressed': on ? 'true' : 'false', title: `carries ${motifLabel(x.id)} (read ${x.dir === 1 ? 'N→C' : 'backwards, C→N'})`,
              on: { click: ev => { ev.stopPropagation(); setSel(on ? { carrier: null } : { iface: fi, motif: x.id, carrier: p.id, partner: p.id }); } } }, `${motifLabel(x.id)} ${x.dir === 1 ? '→' : '←'}`); }),
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
  const sec = el('section', { class: 'card' }, el('div', { class: 'ctitle' }, `Interfaces (${D.interfaces.length})`, info('Bait sites merged in 3-D, one card each; open a card for its motifs and partners. Partners are listed by confidence: H = PAE ≤ 5 Å, M = 5–10 Å, L = above 10 Å (hidden unless asked); underlined residues are motif residues.')),
    pageKey());
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
      el('span', { class: 'imeta' }, `${n} · ${nm} motif${nm === 1 ? '' : 's'}`), cxBtn(`view ${f.id}; select ${f.id}`));
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
    // a div, not a <button>: the row holds a UniProt link, and a link cannot sit inside a button
    const selP = x => setSel({ partner: S.partner === x.id ? null : x.id, carrier: null });
    list.replaceChildren(...r.map(x => el('div', { class: 'brow', role: 'button', tabindex: 0, 'data-on': S.partner === x.id ? 1 : 0,
        on: { click: () => selP(x), keydown: ev => { if (ev.key === 'Enter' || ev.key === ' ') { ev.preventDefault(); selP(x); } } } },
      el('span', {}, el('b', {}, x.name), ' ', x.sub ? subEl(x.id, 'pSub mono') : el('span', { class: 'pSub mono' }, x.id)), el('span', { class: 'bc mono' }, `${x.count} res`)))); };
  const sb = (k, t) => { const b = el('button', { class: 'sortBtn', on: { click: () => { st.sort = k; draw(); sbs.forEach(([kk, bb]) => bb.setAttribute('data-on', st.sort === kk ? 1 : 0)); } } }, t);
    b.setAttribute('data-on', st.sort === k ? 1 : 0); return [k, b]; };
  const sbs = [sb('contacts', 'most contacts'), sb('name', 'name A–Z')];
  draw(); onSel(draw);
  // 2026-10-10 (boss): kept as the record of the proteins screened, named "Prey"; sort buttons labelled
  return el('section', { class: 'card' }, el('div', { class: 'ctitle' }, `Prey (${rows.length})`, info('Every prey that produced interface signal: the proteins this screen used. Search by gene or UniProt; click one to select it.')),
    el('div', { class: 'controls' }, el('input', { placeholder: 'Search gene or UniProt', on: { input: ev => { st.q = ev.target.value; draw(); } } }),
      el('span', { class: 'sortLab' }, 'Sort by'), ...sbs.map(x => x[1]), count), list);
}

document.getElementById('app').append(el('div', { class: 'wrap' }, header(),
  el('div', { class: 'grid' }, el('div', { class: 'side' }, side()), el('div', {}, heatmap(), interfaceCards(), allBinders())),
  el('div', { class: 'foot' }, `SCRIMS ${D.scrims_version} · schema ${D.schema} · ${D.screen.name} · the same data is in scrims.json`)));
"""


def html(doc):
    """The report for one `scrims/1` document, as a string."""
    data = json.dumps(doc, separators=(",", ":")).replace("</", "<\\/")
    idn = doc["bait"].get("identity") or {}
    bait = idn.get("gene") or idn.get("uniprot") or doc["screen"]["name"]
    return (f"<!doctype html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">"
            f"<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">"
            f"<title>{bait} — SCRIMS</title><style>{CSS}</style></head>"
            f"<body><div id=\"app\"></div>"
            f"<script type=\"application/json\" id=\"scrims-data\">{data}</script>"
            f"<script>{LOGO_JS}{JS}</script></body></html>\n")


def write_report(doc, path):
    path.write_text(html(doc), encoding="utf-8")
