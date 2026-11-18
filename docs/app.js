/* SAR Pass & Doppler Simulator — front end.
 * Reads data/passes.json (written by `python -m sarpass`) and renders it with Plotly.
 * Doppler is recomputed here from the stored radial velocity, so the carrier can be changed live:
 *   f_D = -2 f0 * rdot / c      (rdot < 0 = approaching = positive Doppler)
 */
const C = 299792458;
const css = (n) => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const $ = (id) => document.getElementById(id);

let DATA = null;
let selected = null; // pass object
let projection = "flat";

const MONTHS = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"];
function utcLabel(iso, withDate = true) {
  const d = new Date(iso);
  const p = (n) => String(n).padStart(2, "0");
  const t = `${p(d.getUTCHours())}:${p(d.getUTCMinutes())}:${p(d.getUTCSeconds())} UTC`;
  return withDate ? `${d.getUTCDate()} ${MONTHS[d.getUTCMonth()]} ${t}` : t;
}
function localLabel(iso) {
  return new Date(iso).toLocaleTimeString([], { hour: "2-digit", minute: "2-digit" });
}
const f0 = () => parseFloat($("band").value);
const dopplerKHz = (rdot) => (-2 * f0() * rdot) / C / 1e3;

function baseLayout(extra = {}) {
  return Object.assign(
    {
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(0,0,0,0)",
      font: { color: css("--text"), size: 12 },
      margin: { l: 62, r: 62, t: 10, b: 44 },
      xaxis: { gridcolor: css("--line"), zerolinecolor: css("--muted"), title: { text: "Time from closest approach [s]" } },
      legend: { orientation: "h", y: 1.12 },
      hovermode: "x unified",
    },
    extra
  );
}
const PLOT_CFG = { responsive: true, displaylogo: false };

async function init() {
  let res;
  try {
    res = await fetch("data/passes.json", { cache: "no-cache" });
    if (!res.ok) throw new Error(res.status);
    DATA = await res.json();
  } catch (e) {
    $("lede").textContent = "Could not load data/passes.json. Run `python -m sarpass --demo` first (see README).";
    return;
  }
  setupHeader();
  setupControls();
  renderList();
  if (DATA.passes.length) select(DATA.passes[0]);
  else $("chartDoppler").innerHTML = '<p class="empty">No SAR-eligible passes in this window.</p>';
}

function setupHeader() {
  const m = DATA.mission, t = DATA.target, p = DATA.params;
  $("lede").textContent =
    `${m.label} (${m.band}-band, ${(m.freq_hz / 1e9).toFixed(3)} GHz, λ = ${m.wavelength_cm} cm) over ${t.name} ` +
    `(${t.lat.toFixed(3)}°, ${t.lon.toFixed(3)}°). ${DATA.passes.length} SAR-eligible passes ` +
    `(incidence ${p.inc_min_deg}–${p.inc_max_deg}°) in the ${DATA.window_hours} h from ${utcLabel(DATA.window_start_utc)}.`;

  const notes = [];
  if (DATA.source === "demo")
    notes.push("<b>Demo data:</b> these orbits are synthetic (not real satellites). Run the tool with network access to use real CelesTrak TLEs.");
  const end = new Date(DATA.window_start_utc).getTime() + DATA.window_hours * 3600e3;
  if (Date.now() > end)
    notes.push(`<b>This window has ended</b> (generated ${utcLabel(DATA.generated_utc)}). Re-run the pipeline to refresh.`);
  if (DATA.skipped && DATA.skipped.length)
    notes.push(`${DATA.skipped.length} object(s) skipped (stale TLE or propagation error).`);
  if (notes.length) {
    $("banner").innerHTML = notes.join("<br>");
    $("banner").hidden = false;
  }

  // preselect the mission's own carrier
  const sel = $("band");
  const want = String(m.freq_hz);
  if (![...sel.options].some((o) => parseFloat(o.value) === m.freq_hz)) {
    const o = document.createElement("option");
    o.value = want; o.textContent = `${m.band} · ${(m.freq_hz / 1e9).toFixed(3)} GHz`;
    sel.appendChild(o);
  }
  sel.value = [...sel.options].find((o) => parseFloat(o.value) === m.freq_hz).value;

  const sats = [...new Set(DATA.passes.map((p) => p.sat))];
  for (const s of sats) {
    const o = document.createElement("option");
    o.value = s; o.textContent = s;
    $("satFilter").appendChild(o);
  }
}

function setupControls() {
  $("band").addEventListener("change", () => selected && renderDetail());
  $("satFilter").addEventListener("change", renderList);
  $("projFlat").addEventListener("click", () => setProjection("flat"));
  $("projGlobe").addEventListener("click", () => setProjection("globe"));
}
function setProjection(p) {
  projection = p;
  $("projFlat").classList.toggle("on", p === "flat");
  $("projGlobe").classList.toggle("on", p === "globe");
  if (selected) renderMap();
}

function renderList() {
  const filter = $("satFilter").value;
  const list = $("passList");
  list.innerHTML = "";
  const passes = DATA.passes.filter((p) => !filter || p.sat === filter);
  $("count").textContent = `(${passes.length})`;
  for (const p of passes) {
    const li = document.createElement("li");
    li.dataset.id = p.id;
    if (selected && selected.id === p.id) li.classList.add("sel");

    const t = document.createElement("div");
    t.className = "t";
    t.textContent = `${utcLabel(p.tca_utc)}`;
    const pill = document.createElement("span");
    pill.className = `pill ${p.side}`;
    pill.textContent = `${p.side}-looking`;
    t.appendChild(pill);

    const s = document.createElement("div");
    s.className = "s";
    s.textContent = `${p.sat} · local ${localLabel(p.tca_utc)} · max el ${p.max_elevation_deg.toFixed(0)}° · inc@TCA ${p.incidence_tca_deg.toFixed(0)}°`;

    li.append(t, s);
    li.addEventListener("click", () => select(p));
    list.appendChild(li);
  }
}

function select(p) {
  selected = p;
  for (const li of $("passList").children) li.classList.toggle("sel", Number(li.dataset.id) === p.id);
  renderDetail();
}

function renderDetail() {
  renderCards();
  renderDoppler();
  renderGeom();
  renderMap();
}

function renderCards() {
  const p = selected;
  const fd = p.range_rate_ms.map(dopplerKHz);
  const peak = Math.max(...fd.map(Math.abs));
  const ka = ((-2 * f0() * p.range_accel_ms2) / C) / 1e3; // kHz/s
  const win = p.imaging_windows_s.reduce((a, [s, e]) => a + (e - s), 0);
  const cards = [
    ["Peak |Doppler|", peak.toFixed(0), "kHz"],
    ["Doppler rate @ TCA", ka.toFixed(2), "kHz/s"],
    ["Min slant range", p.min_range_km.toFixed(0), "km"],
    ["Max elevation", p.max_elevation_deg.toFixed(1), "°"],
    ["Incidence @ TCA", p.incidence_tca_deg.toFixed(1), "°"],
    ["Off-nadir @ TCA", p.off_nadir_tca_deg.toFixed(1), "°"],
    ["Pass duration", Math.round(p.duration_s / 60 * 10) / 10, "min"],
    ["In SAR window", Math.round(win), "s"],
  ];
  $("cards").innerHTML = "";
  for (const [k, v, u] of cards) {
    const d = document.createElement("div");
    d.className = "card";
    d.innerHTML = `<div class="k"></div><div class="v"></div>`;
    d.querySelector(".k").textContent = k;
    d.querySelector(".v").textContent = v;
    const us = document.createElement("span");
    us.className = "u"; us.textContent = u;
    d.querySelector(".v").appendChild(us);
    $("cards").appendChild(d);
  }
}

const xRel = (p) => p.t_s.map((t) => t - p.t_tca_s);

function windowShapes(p) {
  return p.imaging_windows_s.map(([a, b]) => ({
    type: "rect", xref: "x", yref: "paper", x0: a - p.t_tca_s, x1: b - p.t_tca_s, y0: 0, y1: 1,
    fillcolor: "rgba(95,209,192,0.13)", line: { width: 0 }, layer: "below",
  }));
}

function renderDoppler() {
  const p = selected;
  const x = xRel(p);
  const y = p.range_rate_ms.map(dopplerKHz);
  const layout = baseLayout({
    yaxis: { gridcolor: css("--line"), zerolinecolor: css("--muted"), title: { text: "Doppler [kHz]" } },
    shapes: windowShapes(p),
    showlegend: false,
  });
  Plotly.react("chartDoppler", [
    { x, y, mode: "lines", line: { color: css("--accent"), width: 2.5 }, name: "Doppler",
      hovertemplate: "%{y:.1f} kHz<extra></extra>" },
  ], layout, PLOT_CFG);
}

function renderGeom() {
  const p = selected, x = xRel(p), prm = DATA.params;
  const layout = baseLayout({
    yaxis: { gridcolor: css("--line"), title: { text: "Incidence [°]" }, rangemode: "tozero" },
    yaxis2: { overlaying: "y", side: "right", title: { text: "Slant range [km]" }, showgrid: false },
    shapes: windowShapes(p).concat(
      [prm.inc_min_deg, prm.inc_max_deg].map((v) => ({
        type: "line", xref: "paper", x0: 0, x1: 1, yref: "y", y0: v, y1: v,
        line: { color: css("--muted"), width: 1, dash: "dot" },
      }))
    ),
  });
  Plotly.react("chartGeom", [
    { x, y: p.incidence_deg, mode: "lines", name: "Incidence", line: { color: css("--accent-2"), width: 2.5 },
      hovertemplate: "%{y:.1f}°<extra></extra>" },
    { x, y: p.range_km, mode: "lines", name: "Slant range", yaxis: "y2", line: { color: css("--warn"), width: 2 },
      hovertemplate: "%{y:.0f} km<extra></extra>" },
  ], layout, PLOT_CFG);
}

function renderMap() {
  const p = selected, tg = DATA.target;
  // sub-satellite points that lie inside an imaging window
  const lonIn = [], latIn = [];
  p.t_s.forEach((t, i) => {
    if (p.imaging_windows_s.some(([a, b]) => t >= a && t <= b)) { lonIn.push(p.sat_lon[i]); latIn.push(p.sat_lat[i]); }
    else if (lonIn.length && lonIn[lonIn.length - 1] !== null) { lonIn.push(null); latIn.push(null); }
  });
  const geo = {
    projection: projection === "globe" ? { type: "orthographic", rotation: { lon: tg.lon, lat: tg.lat } } : { type: "natural earth" },
    showland: true, landcolor: "#16233d", showocean: true, oceancolor: "#0a1222",
    showcountries: true, countrycolor: "#2a3c63", coastlinecolor: "#3a4f7d", showcoastlines: true,
    bgcolor: "rgba(0,0,0,0)", lakecolor: "#0a1222",
    lonaxis: { showgrid: true, gridcolor: "#1b2a4a" }, lataxis: { showgrid: true, gridcolor: "#1b2a4a" },
  };
  if (projection === "flat") geo.fitbounds = "locations";
  Plotly.react("chartMap", [
    { type: "scattergeo", lon: p.sat_lon, lat: p.sat_lat, mode: "lines", name: "Ground track (visible arc)",
      line: { width: 2, color: css("--muted") }, hoverinfo: "skip" },
    { type: "scattergeo", lon: lonIn, lat: latIn, mode: "lines", name: "In SAR incidence window",
      line: { width: 5, color: css("--accent-2") }, hoverinfo: "skip" },
    { type: "scattergeo", lon: [p.sat_lon[0], p.sat_lon[p.sat_lon.length - 1]], lat: [p.sat_lat[0], p.sat_lat[p.sat_lat.length - 1]],
      mode: "markers+text", text: ["AOS", "LOS"], textposition: "top center", name: "Rise / set",
      marker: { size: 7, color: css("--text") }, hoverinfo: "text" },
    { type: "scattergeo", lon: [tg.lon], lat: [tg.lat], mode: "markers+text", text: [tg.name], textposition: "bottom center",
      name: "Target", marker: { size: 11, symbol: "star", color: css("--warn") }, hoverinfo: "text" },
  ], {
    geo, paper_bgcolor: "rgba(0,0,0,0)", font: { color: css("--text"), size: 12 },
    margin: { l: 0, r: 0, t: 0, b: 0 }, legend: { orientation: "h", y: 0, bgcolor: "rgba(0,0,0,0)" },
  }, PLOT_CFG);
}

init();
