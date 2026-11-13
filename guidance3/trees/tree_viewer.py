"""Self-contained interactive tree page (no external scripts or cookies): a phylogram whose branches are coloured by
alignment support (internal-node labels 0-100 in the Newick file), in rectangular or circular layout.
Opens in any browser, also offline. The page is ASCII-only (the web server re-encodes result pages)."""
import html
import json

CIRCULAR_FROM_LEAVES = 80  # larger trees open in circular layout, with support labels hidden (shown on hover)

_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
  :root { --ink:#1f2328; --muted:#59636e; --rule:#d0d7de; --bg:#ffffff; --link:#9198a1;
          --s95:#1a5fb4; --s75:#62a0ea; --s50:#e5a50a; --s0:#c01c28; --none:#8c959f; }
  @media (prefers-color-scheme: dark) {
    :root { --ink:#e6edf3; --muted:#9198a1; --rule:#3d444d; --bg:#0d1117; --link:#6e7681;
            --s95:#78aeed; --s75:#3584e4; --s50:#f6d32d; --s0:#ff7b72; --none:#6e7681; }
  }
  body { margin:0; padding:16px; font:14px/1.45 system-ui,-apple-system,"Segoe UI",sans-serif; color:var(--ink); background:var(--bg); }
  h1 { font-size:18px; margin:0 0 4px; }
  p.sub { margin:0 0 12px; color:var(--muted); max-width:75ch; }
  .bar { display:flex; flex-wrap:wrap; gap:10px 18px; align-items:center; margin:0 0 12px; }
  .legend { display:flex; flex-wrap:wrap; gap:6px 14px; }
  .legend span { display:inline-flex; align-items:center; gap:6px; }
  .legend i { display:inline-block; width:22px; height:4px; border-radius:2px; }
  .seg { display:inline-flex; border:1px solid var(--rule); border-radius:6px; overflow:hidden; }
  .seg button { border:0; border-radius:0; }
  .seg button.on { background:var(--rule); }
  label, button { font:inherit; }
  button { border:1px solid var(--rule); background:transparent; color:var(--ink); border-radius:6px; padding:3px 10px; cursor:pointer; }
  #wrap { overflow:auto; border:1px solid var(--rule); border-radius:8px; max-height:85vh; }
  svg text { fill:var(--ink); font-size:12px; }
  svg text.sup { font-size:11px; font-weight:600; fill:var(--ink);
                 paint-order:stroke; stroke:var(--bg); stroke-width:3.5px; stroke-linejoin:round; }
  .edge { stroke-linecap:round; }
  .hit { stroke:transparent; stroke-width:12; cursor:pointer; }
  .hit:hover + .edge, .edge.hl { stroke-width:6; }
  #tip { position:fixed; pointer-events:none; background:var(--bg); color:var(--ink); border:1px solid var(--rule);
         border-radius:6px; padding:6px 9px; font-size:12px; display:none; max-width:360px; box-shadow:0 4px 14px rgba(0,0,0,.15); }
  #tip b { font-weight:600; }
</style>
</head>
<body>
<h1>__TITLE__</h1>
<p class="sub">__SUBTITLE__</p>
<div class="bar">
  <div class="legend">
    <span><i style="background:var(--s95)"></i>&ge; 95%</span>
    <span><i style="background:var(--s75)"></i>75&ndash;94%</span>
    <span><i style="background:var(--s50)"></i>50&ndash;74%</span>
    <span><i style="background:var(--s0)"></i>&lt; 50%</span>
    <span><i style="background:var(--none)"></i>leaf / root edge</span>
  </div>
  <div class="seg" role="group" aria-label="Layout">
    <button id="lyRect" type="button">Rectangular</button><button id="lyCirc" type="button">Circular</button>
  </div>
  <label><input type="checkbox" id="showSup"> support values</label>
  <label title="Draw every branch exactly proportional to its length (short branches may then be too short for their support value)"><input type="checkbox" id="propLen"> proportional branch lengths</label>
  <span><button id="zoomOut" type="button" aria-label="Zoom out">&minus;</button>
        <button id="zoomIn" type="button" aria-label="Zoom in">+</button></span>
  <button id="dl" type="button">Download SVG</button>
</div>
<div id="wrap"></div>
<div id="tip"></div>
<script>
const NEWICK = __NEWICK__;
const INFO = __INFO__;          // per-branch aBayes / SH-aLRT / combined, keyed by branchKey()
const N_TREES = __NTREES__;
const CIRCULAR_FROM = __CIRCULAR_FROM__;

function parseNewick(s) {
  let i = 0;
  const root = {children: []};
  let node = root;
  const stack = [];
  function readName() {
    let out = "";
    if (s[i] === "'") { i++; while (i < s.length && s[i] !== "'") out += s[i++]; i++; return out; }
    while (i < s.length && !",():;[".includes(s[i])) out += s[i++];
    return out.trim();
  }
  while (i < s.length) {
    const c = s[i];
    if (c === "(") { const child = {children: []}; node.children.push(child); stack.push(node); node = child; i++; }
    else if (c === ",") { const sib = {children: []}; stack[stack.length - 1].children.push(sib); node = sib; i++; }
    else if (c === ")") { node = stack.pop(); i++; node.label = readName(); }
    else if (c === ":") { i++; let num = ""; while (i < s.length && /[0-9eE.+\-]/.test(s[i])) num += s[i++]; node.length = parseFloat(num); }
    else if (c === "[") { while (i < s.length && s[i] !== "]") i++; i++; }
    else if (c === ";" || /\s/.test(c)) { i++; }
    else { node.label = readName(); }
  }
  return root.children.length === 1 ? root.children[0] : root;
}

const root = parseNewick(NEWICK);
const leaves = [];
// branch lengths if the tree has them, otherwise unit lengths (cladogram)
const hasLen = (function has(n) { return n.length > 0 || n.children.some(has); })(root);
(function walk(n, x, parent) {
  n.x = x; n.parent = parent; n.label = n.label || "";
  if (!n.children.length) { n.idx = leaves.length; leaves.push(n); return; }
  n.children.forEach(c => walk(c, x + (hasLen ? (c.length || 0) : 1), n));
})(root, 0, null);
(function pos(n) {     // leaf order position (fractional for internal nodes)
  if (!n.children.length) { n.t = n.idx; return; }
  n.children.forEach(pos);
  n.t = (n.children[0].t + n.children[n.children.length - 1].t) / 2;
})(root);
const maxX = Math.max(...leaves.map(l => l.x)) || 1;
const allNames = leaves.map(l => l.label);
const firstName = allNames.slice().sort()[0];

function namesUnder(n) { return n.children.length ? n.children.flatMap(namesUnder) : [n.label]; }
function branchKey(n) {
  let side = namesUnder(n);
  if (side.includes(firstName)) { const s = new Set(side); side = allNames.filter(x => !s.has(x)); }
  return side.slice().sort().join("|");
}
function support(n) { const v = parseFloat(n.label); return n.children.length && !isNaN(v) ? v : null; }
function colorFor(n) {
  const v = support(n);
  if (v === null) return "var(--none)";
  return v >= 95 ? "var(--s95)" : v >= 75 ? "var(--s75)" : v >= 50 ? "var(--s50)" : "var(--s0)";
}
function esc(s) { return String(s).replace(/[&<>"']/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c])); }
function tipHtml(n) {
  if (!n.children.length) return "<b>" + esc(n.label) + "</b>";
  const names = namesUnder(n), v = support(n), extra = INFO[branchKey(n)] || {};
  const rows = [];
  if (v !== null) rows.push("Alignment support: <b>" + v + "%</b> of " + N_TREES + " trees");
  if (extra.abayes != null) rows.push("aBayes (base MSA): <b>" + extra.abayes.toFixed(3) + "</b>");
  if (extra.alrt != null) rows.push("SH-aLRT (base MSA): <b>" + extra.alrt.toFixed(1) + "</b>");
  if (extra.combined != null) rows.push("aBayes &times; alignment support: <b>" + extra.combined.toFixed(3) + "</b>");
  rows.push(names.length + " taxa: " + esc(names.slice(0, 12).join(", ")) + (names.length > 12 ? ", &hellip;" : ""));
  return rows.join("<br>");
}

let layout = leaves.length > CIRCULAR_FROM ? "circ" : "rect";
let zoom = 1;
document.getElementById("showSup").checked = leaves.length <= CIRCULAR_FROM;
const nodes = [];
(function collect(n) { nodes.push(n); n.children.forEach(collect); })(root);

function edgeMarkup(n, d, extraAttr) {
  const id = nodes.indexOf(n);
  return `<path class="hit" d="${d}" data-n="${id}"/>` +
         `<path class="edge" d="${d}" stroke="${colorFor(n)}" stroke-width="${extraAttr}" fill="none" data-n="${id}"/>`;
}

// Display positions: each branch is drawn at least minPx long (room for its support value); longer branches stay
// proportional to their lengths. The scale is searched so that the deepest tip lands at `avail` pixels.
function displayPositions(avail, minPx, proportional) {
  const len = c => (hasLen ? (c.length || 0) : 1);
  function compute(scale, mp) {
    let m = 0;
    (function w(n, px) { n.px = px; m = Math.max(m, px); n.children.forEach(c => w(c, px + Math.max(mp, len(c) * scale))); })(root, 0);
    return m;
  }
  const depth = (function d(n) { return n.children.length ? 1 + Math.max(...n.children.map(d)) : 0; })(root);
  let mp = proportional ? 0 : Math.min(minPx, avail / Math.max(1, depth) * 0.9);
  let lo = 0, hi = avail / maxX;
  if (compute(hi, mp) <= avail) { compute(hi, mp); return; }
  for (let it = 0; it < 40; it++) { const mid = (lo + hi) / 2; if (compute(mid, mp) > avail) hi = mid; else lo = mid; }
  compute(lo, mp);
}

function labelWidth(v) { return 7 * String(v).length + 6; }

function drawRect() {
  const rowH = 20, top = 22, left = 18;
  const longest = Math.max(...allNames.map(s => s.length));
  const labelW = Math.min(460, 7.4 * longest + 16);
  const avail = Math.max(560, (document.getElementById("wrap").clientWidth || 900) - 4);
  const width = Math.round(avail * zoom), plotW = Math.max(200, width - left - labelW);
  const height = Math.round(leaves.length * rowH * Math.max(1, zoom) + top * 2);
  const showSup = document.getElementById("showSup").checked;
  displayPositions(plotW, showSup ? 30 : 0, document.getElementById("propLen").checked);
  const X = n => left + n.px, Y = t => top + t * (height - 2 * top) / Math.max(1, leaves.length - 1);
  const p = [], sup = [];
  nodes.forEach(n => {
    if (n.children.length > 1)
      p.push(`<line x1="${X(n)}" y1="${Y(n.children[0].t)}" x2="${X(n)}" y2="${Y(n.children[n.children.length - 1].t)}" stroke="var(--link)" stroke-width="1.5"/>`);
    if (n.parent) {
      p.push(edgeMarkup(n, `M${X(n.parent)},${Y(n.t)}H${X(n)}`, 3));
      const v = support(n), blen = X(n) - X(n.parent);
      if (showSup && v !== null && blen >= labelWidth(v) - 2)   // centred above its own branch
        sup.push(`<text class="sup" x="${(X(n) + X(n.parent)) / 2}" y="${Y(n.t) - 6}" text-anchor="middle">${v}</text>`);
    }
  });
  leaves.forEach(l => p.push(`<text x="${X(l) + 6}" y="${Y(l.t) + 4}">${esc(l.label)}</text>`));
  return `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">${p.join("")}${sup.join("")}</svg>`;
}

function drawCirc() {
  const longest = Math.max(...allNames.map(s => s.length));
  const labelW = Math.min(300, 6.6 * longest + 14);
  const n = leaves.length;
  const avail = Math.max(560, Math.min((document.getElementById("wrap").clientWidth || 900) - 4, 1100));
  const minR = n * 11 / (2 * Math.PI);                 // room for one label line per leaf on the circle
  const R = Math.max(minR, avail / 2 - labelW) * zoom;
  const size = Math.round(2 * (R + labelW) + 20), c = size / 2;
  const step = 2 * Math.PI / Math.max(1, n);
  const A = t => -Math.PI / 2 + step * t;
  const P = (r, a) => [c + r * Math.cos(a), c + r * Math.sin(a)];
  const showSup = document.getElementById("showSup").checked;
  displayPositions(R, showSup ? 26 : 0, document.getElementById("propLen").checked);
  const leafCount = nd => nd.nLeaves || (nd.nLeaves = nd.children.length ? nd.children.reduce((s, ch) => s + leafCount(ch), 0) : 1);
  const p = [], sup = [];
  nodes.forEach(nd => {
    if (nd.children.length > 1) {
      const a1 = A(nd.children[0].t), a2 = A(nd.children[nd.children.length - 1].t), r = nd.px;
      const [x1, y1] = P(r, a1), [x2, y2] = P(r, a2);
      p.push(`<path d="M${x1},${y1}A${r},${r} 0 ${a2 - a1 > Math.PI ? 1 : 0} 1 ${x2},${y2}" stroke="var(--link)" stroke-width="1.2" fill="none"/>`);
    }
    if (nd.parent) {
      const a = A(nd.t), [x1, y1] = P(nd.parent.px, a), [x2, y2] = P(nd.px, a);
      p.push(edgeMarkup(nd, `M${x1},${y1}L${x2},${y2}`, 2.2));
      const v = support(nd);
      if (showSup && v !== null) {
        const rMid = (nd.px + nd.parent.px) / 2, radial = nd.px - nd.parent.px;
        const room = leafCount(nd) * step * rMid;          // arc space owned by this branch
        if (radial >= 16 && room >= labelWidth(v) + 4) {
          const [sx, sy] = P(rMid, a);
          sup.push(`<text class="sup" x="${sx}" y="${sy + 4}" text-anchor="middle">${v}</text>`);
        }
      }
    }
  });
  leaves.forEach(l => {
    const a = A(l.t), deg = a * 180 / Math.PI, flip = Math.cos(a) < 0;
    if (l.px < R - 2) {   // dotted leader from the tip to its label on the outer ring
      const [lx1, ly1] = P(l.px + 2, a), [lx2, ly2] = P(R + 2, a);
      p.push(`<line x1="${lx1}" y1="${ly1}" x2="${lx2}" y2="${ly2}" stroke="var(--rule)" stroke-width="0.8" stroke-dasharray="2 3"/>`);
    }
    const [x, y] = P(R + 6, a);
    p.push(`<text x="${x}" y="${y}" dy="4" font-size="11" text-anchor="${flip ? "end" : "start"}" transform="rotate(${flip ? deg + 180 : deg} ${x} ${y})">${esc(l.label)}</text>`);
  });
  return `<svg xmlns="http://www.w3.org/2000/svg" width="${size}" height="${size}" viewBox="0 0 ${size} ${size}">${p.join("")}${sup.join("")}</svg>`;
}

function draw() {
  document.getElementById("wrap").innerHTML = layout === "circ" ? drawCirc() : drawRect();
  document.getElementById("lyRect").classList.toggle("on", layout === "rect");
  document.getElementById("lyCirc").classList.toggle("on", layout === "circ");
}

const tip = document.getElementById("tip"), wrap = document.getElementById("wrap");
wrap.addEventListener("mousemove", e => {
  const id = e.target.getAttribute && e.target.getAttribute("data-n");
  wrap.querySelectorAll(".edge.hl").forEach(el => el.classList.remove("hl"));
  if (id !== null && id !== undefined) {
    const edge = wrap.querySelector(`.edge[data-n="${id}"]`); if (edge) edge.classList.add("hl");
    tip.innerHTML = tipHtml(nodes[+id]); tip.style.display = "block";
    const x = Math.min(e.clientX + 14, window.innerWidth - 370);
    tip.style.left = x + "px"; tip.style.top = (e.clientY + 14) + "px";
  } else tip.style.display = "none";
});
wrap.addEventListener("mouseleave", () => { tip.style.display = "none"; });
document.getElementById("lyRect").addEventListener("click", () => { layout = "rect"; zoom = 1; draw(); });
document.getElementById("lyCirc").addEventListener("click", () => { layout = "circ"; zoom = 1; draw(); });
document.getElementById("showSup").addEventListener("change", draw);
document.getElementById("propLen").addEventListener("change", draw);
document.getElementById("zoomIn").addEventListener("click", () => { zoom = Math.min(4, zoom * 1.25); draw(); });
document.getElementById("zoomOut").addEventListener("click", () => { zoom = Math.max(0.4, zoom / 1.25); draw(); });
document.getElementById("dl").addEventListener("click", () => {
  const svg = wrap.querySelector("svg").cloneNode(true);
  const css = getComputedStyle(document.documentElement);
  svg.querySelectorAll(".hit").forEach(el => el.remove());
  svg.querySelectorAll("[stroke^='var(']").forEach(el => el.setAttribute("stroke", css.getPropertyValue(el.getAttribute("stroke").slice(4, -1)).trim()));
  svg.querySelectorAll("text").forEach(el => {
    el.setAttribute("fill", "#1f2328");
    if (el.classList.contains("sup")) { el.setAttribute("stroke", "#ffffff"); el.setAttribute("stroke-width", "3.5");
      el.setAttribute("paint-order", "stroke"); el.setAttribute("font-size", "11"); el.setAttribute("font-weight", "600"); }
    else if (!el.getAttribute("font-size")) el.setAttribute("font-size", "12");
    el.setAttribute("font-family", "sans-serif");
  });
  const blob = new Blob([new XMLSerializer().serializeToString(svg)], {type: "image/svg+xml"});
  const a = document.createElement("a"); a.href = URL.createObjectURL(blob); a.download = "__DOWNLOAD__.svg"; a.click();
});
let rt; window.addEventListener("resize", () => { clearTimeout(rt); rt = setTimeout(draw, 150); });
draw();
</script>
</body>
</html>
"""


def write_tree_html(newick_path, out_html, title, n_trees, branch_info=None):
    """Write an interactive tree page for a Newick tree with alignment-support labels. branch_info: optional
    {branch key: {"abayes", "alrt", "combined"}} shown on hover (keys from alignment_trees._viewer_key)."""
    with open(newick_path) as fh:
        newick = fh.read().strip()
    subtitle = (f"Branch colours show alignment support: the percentage of the {n_trees} trees inferred from "
                f"alternative MSAs that contain the branch. Hover over a branch for its support values and taxa.")
    page = (_TEMPLATE
            .replace("__TITLE__", html.escape(title))
            .replace("__SUBTITLE__", html.escape(subtitle))
            .replace("__DOWNLOAD__", "".join(ch if ch.isalnum() else "_" for ch in title))
            .replace("__NTREES__", str(int(n_trees)))
            .replace("__CIRCULAR_FROM__", str(CIRCULAR_FROM_LEAVES))
            .replace("__INFO__", json.dumps(branch_info or {}).replace("</", "<\\/"))
            .replace("__NEWICK__", json.dumps(newick).replace("</", "<\\/")))
    # ASCII only (non-ASCII as HTML/JS escapes): the web server re-encodes result pages
    with open(out_html, "w", encoding="ascii", errors="xmlcharrefreplace") as fh:
        fh.write(page)
