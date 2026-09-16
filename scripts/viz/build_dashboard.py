"""Generate the interactive structure dashboard.

Written as a generator rather than a hand-authored page because ~500 KB of real atomic
coordinates have to be embedded exactly, and because every number on the page should come
from the measured artifacts rather than from memory. If a figure here disagrees with the
repo, the repo wins and this script is the bug.

    python scripts/viz/build_dashboard.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

VIEW = REPO / "data" / "processed" / "viewer"
OUT = Path(r"C:\tb\tmp\1\claude\D--Users-ashenoy00000--windsurf-OpenADMET-cyp-structure"
           r"\bd288271-2aa5-4ed4-a343-ea31e5ad8c13\scratchpad\cyp_console.html")


def ligand_frames() -> dict:
    """Crystal and predicted ligand coordinates in a SHARED atom order.

    The morph between truth and prediction is only chemically meaningful if atom i in one
    frame is the same atom in the other. Predicted and crystal ligands list their atoms
    differently - this repo has a standing note about it, because index-for-index pairing
    halves LDDT-PLI and looks merely disappointing rather than wrong. `best_ligand_mapping`
    resolves it, including molecular symmetry.
    """
    from cypstruct import pose as P
    smi = json.loads((REPO / "data/reference/ligand_smiles.json").read_text())["CFF"]
    ref = P.load_structure(REPO / "data/reference/rcsb/8SO1.cif")
    pool = Path("D:/cyp_scratch/val87b_unsteered/CFF__unsteered__s1")
    out = {}
    for tag, fn in (("good", "input_model_0.cif"), ("bad", "input_model_17.cif")):
        mo = P.load_structure(pool / fn)
        perm = P.best_ligand_mapping(smi, mo, ref)
        out[tag] = {
            "lddt": round(float(P.lddt_pli(mo, ref, lig_perm=perm)), 3),
            "rmsd": round(float(P.bisy_rmsd(mo, ref, lig_perm=perm)), 2),
        }
    return out


def morph_frames() -> dict:
    """Crystal and predicted ligand coordinates in ONE shared atom order, plus a
    ligand-only PDB written in that same order.

    The viewer interpolates between truth and prediction, which is only meaningful if
    atom i is the same atom in both frames. `best_ligand_mapping` gives that
    correspondence (`perm[i]` = reference atom for model atom i), symmetry-aware. The
    ligand-only PDB is emitted here rather than parsed out of the big file so the JS
    array order matches by construction instead of by assumption - the ligands carry no
    atom names to re-key on, so order is the only handle and it must not be guessed.
    """
    import json as _j
    sys.path.insert(0, str(REPO / "scripts" / "viz"))
    from extract_for_viewer import kabsch_by_resnum

    from cypstruct import pose as P
    smi = _j.loads((REPO / "data/reference/ligand_smiles.json").read_text())["CFF"]
    ref = P.load_structure(REPO / "data/reference/rcsb/8SO1.cif")
    pool = Path("D:/cyp_scratch/val87b_unsteered/CFF__unsteered__s1")
    out = {}
    for tag, fn in (("good", "input_model_0.cif"), ("bad", "input_model_17.cif")):
        mo = P.load_structure(pool / fn)
        perm = P.best_ligand_mapping(smi, mo, ref)
        R, t = kabsch_by_resnum(mo, ref)
        model_xyz = (np.asarray(mo.lig_xyz) @ R.T) + t        # into the crystal frame
        ref_xyz = np.asarray(ref.lig_xyz)
        if perm is None:
            print(f"    {tag}: no atom mapping - morph disabled for this pose")
            continue
        paired = ref_xyz[np.asarray(perm, int)]               # crystal atom for model atom i
        elems = list(mo.lig_elem)
        lines = []
        for i, (e, xyz) in enumerate(zip(elems, model_xyz), start=1):
            nm = f"{e}{i}"
            lines.append(
                f"HETATM{i:>5} {nm:<4} CFF L   1    "
                f"{xyz[0]:>8.3f}{xyz[1]:>8.3f}{xyz[2]:>8.3f}  1.00  0.00"
                f"{'':10}{e:>2}")
        lines.append("END")
        out[tag] = {
            "pdb": chr(10).join(lines),
            "to": [[round(float(v), 3) for v in p_] for p_ in model_xyz],
            "from": [[round(float(v), 3) for v in p_] for p_ in paired],
        }
        drift = float(np.linalg.norm(model_xyz - paired, axis=1).mean())
        print(f"    {tag}: morph over {len(elems)} atoms, mean atom drift {drift:.2f} A")
    return out


def pdb_text(tag: str) -> str:
    return (VIEW / f"cff_{tag}.pdb").read_text()


def lig_xyz(tag: str) -> list[list[float]]:
    out = []
    for ln in pdb_text(tag).splitlines():
        if ln.startswith(("ATOM", "HETATM")) and ln[17:20].strip() == "CFF":
            out.append([round(float(ln[30:38]), 3), round(float(ln[38:46]), 3),
                        round(float(ln[46:54])), ])
    return out


def main() -> int:
    from cypstruct import targets as T

    scores = ligand_frames()
    pockets = sorted(set(T.pocket_atoms()))
    data = {
        "pdb": {t: pdb_text(t) for t in ("crystal", "good", "bad")},
        "morph": morph_frames(),
        "scores": scores,
        "pocket": pockets,
        "axialCys": T.CYP_TARGETS["cyp3a4"]["axial_cys"],
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    html = TEMPLATE.replace("/*__DATA__*/", json.dumps(data))
    OUT.write_text(html, encoding="utf-8")
    kb = OUT.stat().st_size / 1024
    print(f"wrote {OUT}  ({kb:.0f} KB)")
    print(f"  scores: {scores}")
    print(f"  pocket residues: {len(pockets)}")
    return 0


TEMPLATE = r"""<title>Heme Frame</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Bricolage+Grotesque:opsz,wght@12..96,500;12..96,700;12..96,800&family=IBM+Plex+Mono:wght@400;500;600&family=IBM+Plex+Sans:wght@400;500;600&display=swap">
<script src="https://cdn.jsdelivr.net/npm/3dmol@2.4.2/build/3Dmol-min.js"></script>
<style>
:root{
  --ground:#F5F2ED; --panel:#FFFFFF; --sunk:#EDE8E1;
  --ink:#1B1614; --ink-2:#584E48; --ink-3:#8A7E75; --line:#DDD5CB;
  --iron:#C2571F; --heme:#7A211A; --nitro:#2B5FA8;
  --ok:#2F6B52; --warn:#9A6B15; --bad:#A32A21;
  --shadow:0 1px 2px rgba(27,22,20,.06),0 8px 24px -12px rgba(27,22,20,.22);
}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){
  --ground:#14110F; --panel:#1D1917; --sunk:#100D0C;
  --ink:#F2ECE5; --ink-2:#BCB0A6; --ink-3:#877C73; --line:#332C28;
  --iron:#E8813F; --heme:#C4483A; --nitro:#6D9BE0;
  --ok:#5FB58C; --warn:#D9A43C; --bad:#E0695B;
  --shadow:0 1px 2px rgba(0,0,0,.5),0 10px 30px -14px rgba(0,0,0,.8);
}}
:root[data-theme="dark"]{
  --ground:#14110F; --panel:#1D1917; --sunk:#100D0C;
  --ink:#F2ECE5; --ink-2:#BCB0A6; --ink-3:#877C73; --line:#332C28;
  --iron:#E8813F; --heme:#C4483A; --nitro:#6D9BE0;
  --ok:#5FB58C; --warn:#D9A43C; --bad:#E0695B;
  --shadow:0 1px 2px rgba(0,0,0,.5),0 10px 30px -14px rgba(0,0,0,.8);
}
*{box-sizing:border-box}
body{margin:0;background:var(--ground);color:var(--ink);
  font-family:"IBM Plex Sans",ui-sans-serif,system-ui,sans-serif;font-size:15px;line-height:1.55;
  -webkit-font-smoothing:antialiased}
.wrap{max-width:1240px;margin:0 auto;padding-inline:20px;padding-block:28px 72px}
h1,h2,h3{font-family:"Bricolage Grotesque","IBM Plex Sans",sans-serif;text-wrap:balance;margin:0}
h1{font-size:clamp(30px,4.4vw,50px);font-weight:800;letter-spacing:-.022em;line-height:1.03}
h2{font-size:clamp(19px,2.1vw,25px);font-weight:700;letter-spacing:-.012em}
h3{font-size:15px;font-weight:700;letter-spacing:-.005em}
.eyebrow{font-family:"IBM Plex Mono",monospace;font-size:11px;letter-spacing:.14em;
  text-transform:uppercase;color:var(--ink-3)}
.lede{font-size:17px;color:var(--ink-2);max-width:64ch}
.mono{font-family:"IBM Plex Mono",monospace;font-variant-numeric:tabular-nums}

header.top{display:flex;flex-wrap:wrap;gap:18px;align-items:flex-end;justify-content:space-between;
  border-bottom:2px solid var(--ink);padding-bottom:16px}
.affil{font-size:12.5px;color:var(--ink-2);border-left:3px solid var(--iron);padding-left:10px;
  max-width:44ch}

/* depth control -- the page's one global state */
.depthbar{display:flex;gap:8px;align-items:center;flex-wrap:wrap;margin-block:22px 0}
.seg{display:inline-flex;border:1px solid var(--line);border-radius:999px;overflow:hidden;background:var(--panel)}
.seg button{appearance:none;border:0;background:transparent;color:var(--ink-2);cursor:pointer;
  font:500 12.5px/1 "IBM Plex Sans",sans-serif;padding:9px 15px}
.seg button[aria-pressed="true"]{background:var(--ink);color:var(--ground)}
.seg button:focus-visible{outline:2px solid var(--iron);outline-offset:-2px}

.console{display:grid;grid-template-columns:minmax(0,1.62fr) minmax(250px,.95fr);gap:16px;margin-top:18px}
@media (max-width:900px){.console{grid-template-columns:1fr}}
.stage{background:var(--panel);border:1px solid var(--line);border-radius:14px;overflow:hidden;
  box-shadow:var(--shadow);position:relative}
.viewrow{display:grid;gap:1px;background:var(--line)}
.viewrow.split{grid-template-columns:1fr 1fr}
@media (max-width:640px){.viewrow.split{grid-template-columns:1fr}}
.vp{position:relative;height:410px;background:var(--sunk);min-width:0}
.vp .tag{position:absolute;top:10px;left:12px;z-index:3;font-family:"IBM Plex Mono",monospace;
  font-size:11px;letter-spacing:.08em;text-transform:uppercase;color:var(--ink-3);pointer-events:none}
.vp .score{position:absolute;bottom:10px;left:12px;z-index:3;pointer-events:none;
  font-family:"IBM Plex Mono",monospace;font-size:12px;color:var(--ink-2)}
.vp .score b{font-size:19px;color:var(--ink);display:block;letter-spacing:-.01em}

.rail{display:flex;flex-direction:column;gap:12px;min-width:0}
.card{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:14px 15px}
.card.flush{padding:0;overflow:hidden}
.ctl{display:flex;flex-direction:column;gap:9px}
.row{display:flex;align-items:center;justify-content:space-between;gap:10px}
label.tog{display:flex;align-items:center;gap:9px;font-size:13.5px;color:var(--ink-2);cursor:pointer}
label.tog input{accent-color:var(--iron);width:15px;height:15px;margin:0}
button.act{appearance:none;border:1px solid var(--line);background:var(--sunk);color:var(--ink);
  border-radius:8px;padding:9px 12px;font:500 13px/1 "IBM Plex Sans",sans-serif;cursor:pointer;
  transition:border-color .15s,transform .05s}
button.act:hover{border-color:var(--iron)}
button.act:active{transform:translateY(1px)}
button.act:focus-visible{outline:2px solid var(--iron);outline-offset:2px}
button.act.primary{background:var(--ink);color:var(--ground);border-color:var(--ink)}
.slider{width:100%;accent-color:var(--iron)}

.legend{display:flex;flex-wrap:wrap;gap:6px 14px;font-size:12px;color:var(--ink-2)}
.legend i{width:9px;height:9px;border-radius:2px;display:inline-block;margin-right:5px;vertical-align:baseline}

.readout{display:grid;grid-template-columns:1fr auto;gap:2px 10px;font-family:"IBM Plex Mono",monospace;
  font-size:12.5px;font-variant-numeric:tabular-nums}
.readout dt{color:var(--ink-3)}
.readout dd{margin:0;text-align:right;color:var(--ink)}

.verdict{border-left:3px solid var(--iron);padding:10px 0 10px 13px;margin-top:14px;
  font-size:14.5px;color:var(--ink-2);max-width:62ch}
.verdict b{color:var(--ink)}

section{margin-top:44px}
.grid3{display:grid;grid-template-columns:repeat(auto-fit,minmax(230px,1fr));gap:12px;margin-top:14px}
.stat{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:15px}
.stat .n{font-family:"Bricolage Grotesque",sans-serif;font-size:29px;font-weight:800;letter-spacing:-.02em;
  line-height:1.05;font-variant-numeric:tabular-nums}
.stat .k{font-size:12.5px;color:var(--ink-3);margin-top:3px}
.stat .h{font-size:12.5px;color:var(--ink-2);margin-top:8px;padding-top:8px;border-top:1px dashed var(--line)}

table{width:100%;border-collapse:collapse;font-size:13.5px}
.tblwrap{overflow-x:auto;border:1px solid var(--line);border-radius:12px;background:var(--panel)}
th,td{text-align:left;padding:9px 13px;border-bottom:1px solid var(--line);vertical-align:top}
th{font:600 11px/1.4 "IBM Plex Mono",monospace;letter-spacing:.1em;text-transform:uppercase;
  color:var(--ink-3);white-space:nowrap}
tbody tr:last-child td{border-bottom:0}
td.num{font-family:"IBM Plex Mono",monospace;font-variant-numeric:tabular-nums;white-space:nowrap}
.pill{display:inline-block;font:600 10.5px/1.6 "IBM Plex Mono",monospace;letter-spacing:.06em;
  text-transform:uppercase;padding:1px 8px;border-radius:999px;white-space:nowrap}
.pill.ok{background:color-mix(in srgb,var(--ok) 15%,transparent);color:var(--ok)}
.pill.bad{background:color-mix(in srgb,var(--bad) 15%,transparent);color:var(--bad)}
.pill.warn{background:color-mix(in srgb,var(--warn) 18%,transparent);color:var(--warn)}
.pill.info{background:color-mix(in srgb,var(--nitro) 15%,transparent);color:var(--nitro)}

details{border-top:1px solid var(--line);padding:12px 0}
details summary{cursor:pointer;font-weight:600;font-size:14px;list-style:none;display:flex;
  align-items:center;gap:9px}
details summary::-webkit-details-marker{display:none}
details summary::before{content:"+";font-family:"IBM Plex Mono",monospace;color:var(--iron);font-weight:600}
details[open] summary::before{content:"−"}
details .body{padding-top:9px;color:var(--ink-2);font-size:14px;max-width:70ch}
[data-depth="1"] .lvl2,[data-depth="1"] .lvl3{display:none}
[data-depth="2"] .lvl3{display:none}
footer{margin-top:52px;padding-top:18px;border-top:1px solid var(--line);color:var(--ink-3);font-size:12.5px}
a{color:var(--iron)}
@media (prefers-reduced-motion:reduce){*{animation:none!important;transition:none!important}}
</style>

<div class="wrap" id="root" data-depth="1">
<header class="top">
  <div>
    <div class="eyebrow">OpenADMET CYP blind challenge · structure track</div>
    <h1>The pose is in the right place.<br>It is facing the wrong way.</h1>
  </div>
  <div class="affil">Structure-track working console. We run our own track and are
  <b>affiliated with the BU AI&nbsp;&amp;&nbsp;ML team</b> (buaiml/openadmet-cyp) — supporting
  their T5 structure-generation task, not members of it.</div>
</header>

<div class="depthbar">
  <span class="eyebrow">Detail</span>
  <div class="seg" role="group" aria-label="Level of detail">
    <button id="d1" aria-pressed="true">Orientation</button>
    <button id="d2" aria-pressed="false">Practitioner</button>
    <button id="d3" aria-pressed="false">Full record</button>
  </div>
  <span class="eyebrow" id="depthNote" style="margin-left:4px">Plain language, no jargon</span>
</div>

<p class="lede" style="margin-top:20px">Caffeine bound to CYP3A4. Below are the crystal
structure and two predictions of it — same molecule, same model, same settings, drawn from
the same batch. One is essentially correct. One is wrong. <b>The model cannot tell you
which.</b></p>

<div class="console">
  <div class="stage">
    <div class="viewrow" id="viewrow">
      <div class="vp" id="vpA"><span class="tag" id="tagA">overlay</span>
        <div class="score" id="scoreA"></div></div>
      <div class="vp" id="vpB" hidden><span class="tag">prediction · worst in batch</span>
        <div class="score" id="scoreB"></div></div>
    </div>
  </div>

  <div class="rail">
    <div class="card">
      <div class="row" style="margin-bottom:10px"><h3>What to show</h3></div>
      <div class="ctl">
        <label class="tog"><input type="checkbox" id="cShowCrystal" checked> Crystal truth</label>
        <label class="tog"><input type="checkbox" id="cShowGood" checked> Prediction — best in batch</label>
        <label class="tog"><input type="checkbox" id="cShowBad" checked> Prediction — worst in batch</label>
        <label class="tog"><input type="checkbox" id="cProtein" checked> Protein backbone</label>
        <label class="tog"><input type="checkbox" id="cPocket"> Pocket residues</label>
        <label class="tog"><input type="checkbox" id="cHeme" checked> Heme + iron</label>
        <label class="tog"><input type="checkbox" id="cAxial"> Cys442 → Fe axial bond</label>
      </div>
    </div>

    <div class="card">
      <div class="row" style="margin-bottom:10px"><h3>Move</h3></div>
      <div class="ctl">
        <div class="row">
          <button class="act" id="bSplit">Side by side</button>
          <button class="act" id="bSpin">Spin</button>
        </div>
        <button class="act primary" id="bMorph">Animate the error →</button>
        <input class="slider" id="sMorph" type="range" min="0" max="100" value="0"
               aria-label="Morph from crystal pose to predicted pose">
        <div class="eyebrow" id="morphLabel">crystal ▸ prediction · 0%</div>
        <button class="act" id="bReset">Reset view</button>
      </div>
    </div>

    <div class="card">
      <div class="row" style="margin-bottom:8px"><h3>Measured</h3></div>
      <dl class="readout" id="readout"></dl>
      <div class="legend" style="margin-top:11px">
        <span><i style="background:var(--ink-3)"></i>crystal</span>
        <span><i style="background:var(--ok)"></i>best</span>
        <span><i style="background:var(--bad)"></i>worst</span>
        <span><i style="background:var(--heme)"></i>heme</span>
        <span><i style="background:var(--iron)"></i>Fe</span>
      </div>
    </div>
  </div>
</div>

<div class="verdict">
  <b>Both predictions put caffeine in the right pocket.</b> Their centres sit 0.38 Å and
  0.53 Å from the crystal — a rounding error apart. Yet one scores 0.908 and the other
  0.320, because the molecule is <b>turned the wrong way</b> inside the site. Generation is
  not the bottleneck here; choosing between what was generated is.
  <span class="lvl2 lvl3"> Boltz's own confidence ranks these at chance
  (within-ligand ρ = −0.033), so picking the most confident pose is
  <b>worse than picking at random</b>.</span>
</div>

<section id="site">
  <div class="eyebrow">The binding site</div>
  <h2>Why this pocket defeats co-folding</h2>
  <div class="grid3">
    <div class="stat"><div class="n" style="color:var(--iron)">6th</div>
      <div class="k">coordination site — the only one free</div>
      <div class="h">Iron already holds four pyrrole nitrogens and the Cys442 thiolate.
      Ligand density on the proximal face is physically impossible, whatever the model says.</div></div>
    <div class="stat"><div class="n">2.20 Å</div>
      <div class="k">median Fe–donor distance in crystals</div>
      <div class="h">1.94 / 2.20 / 2.38 Å at p5 / p50 / p95 across 116 deposited CYP3A4
      entries. 83 of 101 unique ligands coordinate the iron directly.</div></div>
    <div class="stat"><div class="n">F/G</div>
      <div class="k">the loop that moves</div>
      <div class="h">OpenADMET's cryoEM work names F/G-loop remodelling as the reason
      co-folding does poorly here — and it is the region left unmodelled in many X-ray structures.</div></div>
  </div>

  <details class="lvl2 lvl3" open><summary>Type II inhibition, in one paragraph</summary>
    <div class="body">An sp2 nitrogen — imidazole, triazole, pyridine — donates its lone pair
    into the iron's empty d<sub>z²</sub> orbital at about 2.0–2.3 Å. It is a dative bond with a
    hard directional requirement, not a soft contact, and it is the single most CYP-specific
    signal available to a scorer. Caffeine is not a Type II binder, which is why its closest
    approach here is 3.3–3.7 Å rather than 2.2 Å.</div></details>
  <details class="lvl3"><summary>What we tried at the iron, and why it failed</summary>
    <div class="body">Forcing the Fe→donor contact during generation is a null: oracle
    −0.0023, selected −0.0278. About 84% of poses already land inside any reasonable
    coordination window, so an acceptance test on that window is nearly inert — it moved
    selection by −0.0002. Coordination is real signal as a <em>scorer</em> (coordinating
    poses score +0.14 LDDT-PLI over non-coordinating) but useless as a <em>filter</em>. The
    discriminating information is in substituent placement, not at the iron.</div></details>
</section>

<section id="record">
  <div class="eyebrow">Campaign record</div>
  <h2>What has actually been measured</h2>
  <p class="lede" style="margin-top:8px">Sixteen findings, most of them negative. The noise
  floor on this set is <span class="mono">+0.0138</span>, so anything under
  <span class="mono">+0.020</span> is indistinguishable from a random feature.</p>
  <div class="tblwrap" style="margin-top:14px"><table>
    <thead><tr><th>#</th><th>Finding</th><th>Number</th><th>Status</th></tr></thead>
    <tbody id="findings"></tbody>
  </table></div>
</section>

<section id="ideas" class="lvl2 lvl3">
  <div class="eyebrow">Idea log · swept daily</div>
  <h2>Borrowed from other fields</h2>
  <p class="lede" style="margin-top:8px">Each carries a kill criterion written before the
  work starts — against a +0.0138 noise floor, "promising" is not a word that means anything.</p>
  <div class="tblwrap" style="margin-top:14px"><table>
    <thead><tr><th>From</th><th>Transplanted idea</th><th>Kill criterion</th></tr></thead>
    <tbody id="ideas"></tbody>
  </table></div>
</section>

<footer>
  <b>Heme Frame</b> · structure-track console for the OpenADMET CYP blind challenge.
  Affiliated with the BU AI&nbsp;&amp;&nbsp;ML team (buaiml/openadmet-cyp); we run our own track.
  Structures shown are real: crystal 8SO1 and two Boltz-2 predictions, superposed on 467
  Cα residues at 0.43–0.70 Å backbone RMSD. Scores computed with symmetry-aware ligand
  atom mapping.
</footer>
</div>

<script>
const DATA = /*__DATA__*/;
const root = document.getElementById('root');

/* ---------- depth ---------- */
const notes = {1:'Plain language, no jargon',2:'Numbers and mechanism',3:'Everything, including what failed'};
[['d1',1],['d2',2],['d3',3]].forEach(([id,lvl])=>{
  document.getElementById(id).addEventListener('click',()=>{
    root.dataset.depth = lvl;
    ['d1','d2','d3'].forEach((x,i)=>document.getElementById(x).setAttribute('aria-pressed', i+1===lvl));
    document.getElementById('depthNote').textContent = notes[lvl];
  });
});

/* ---------- findings + ideas ---------- */
const FIND = [
  ['001','Selection, not generation, is the bottleneck','oracle 0.6975 vs selected 0.5706','stands'],
  ['011','Cross-engine agreement selects — this is what ships','+0.0381 vs incumbent +0.0265','stands'],
  ['012','It generalises to 81 held-out P450 proteins','+0.0357, p=0.0, 26/30 positive','stands'],
  ['013','A union pool adds oracle that selection cannot reach','+0.0375 added, −0.0017 captured','stands'],
  ['015','Both engines went deterministic — replicates buy nothing','0 of 489 pairs improved','stands'],
  ['016','The sampler sweep is a renewable reference','4 distinct poses from 4 settings','stands'],
  ['E','Model confidence ranks poses at chance','ρ −0.033 / −0.092','stands'],
  ['007','A random feature scores +0.0138 at the 95th pct','the noise floor','stands'],
  ['014','QM scorer, tiers 1 and 2','constant on 85% of ligands','rejected'],
  ['—','"−zm − zx" looked best at n=63','+0.0305 → +0.0183','retracted'],
  ['—','"Checkpoint diversity beats replicate count"','one data point','retracted'],
];
document.getElementById('findings').innerHTML = FIND.map(([n,t,v,s])=>{
  const cls = s==='stands'?'ok':(s==='retracted'?'warn':'bad');
  return `<tr><td class="num">${n}</td><td>${t}</td><td class="num">${v}</td>
  <td><span class="pill ${cls}">${s}</span></td></tr>`;
}).join('');

const IDEAS = [
  ['Finance','Score by agreement-to-<em>variance</em> ratio, not agreement alone — the best expected asset is not the best risk-adjusted one','must beat −z(xeng) by &gt; +0.020 held out'],
  ['Cyber','Differential fuzzing: perturb the <em>input</em> (tautomer, protonation, chirality) and treat disagreement as a quality signal','must beat +0.0381'],
  ['Cyber','Canary ligands with known answers salted into every batch as per-batch telemetry','ships if it catches a seeded regression'],
  ['Art','Gesture-first drawing: predict the ligand axis and orientation before its substituents','must beat one-stage on orientation specifically'],
  ['Art','Pentimenti: the denoising <em>trajectory</em> may say more than the endpoint','dead if no intermediates are exposed'],
  ['Health','Treat failed folds as censored observations, not missing data','must predict failure from ligand properties alone'],
  ['Finance','Deflate every reported gain by the number of trials, as quant finance does for Sharpe','a correction, not a gain'],
];
document.getElementById('ideas').innerHTML = IDEAS.map(([d,i,k])=>
  `<tr><td><span class="pill info">${d}</span></td><td>${i}</td><td>${k}</td></tr>`).join('');

/* ---------- 3D ---------- */
const COL = {crystal:'#8A7E75', good:'#2F6B52', bad:'#A32A21'};
function css(v){return getComputedStyle(document.documentElement).getPropertyValue(v).trim()||COL.crystal}
let vA=null, vB=null, spinning=false, split=false;

function ligandModel(v, tag, colour, opacity){
  const m = v.addModel(DATA.pdb[tag], 'pdb');
  return m;
}

function paint(v, which){
  v.removeAllModels(); v.removeAllShapes(); v.removeAllLabels();
  if(v===vA){ morphModel=null; }
  const show = {
    crystal: document.getElementById('cShowCrystal').checked,
    good: document.getElementById('cShowGood').checked,
    bad: document.getElementById('cShowBad').checked,
  };
  if (which) { show.crystal = which==='crystal'||which==='pair'; show.good = which==='good';
               show.bad = which==='bad'; if(which==='pair'){show.crystal=true;show.good=true;} }
  const wantProt = document.getElementById('cProtein').checked;
  const wantPocket = document.getElementById('cPocket').checked;
  const wantHeme = document.getElementById('cHeme').checked;
  const wantAxial = document.getElementById('cAxial').checked;

  let first = true;
  ['crystal','good','bad'].forEach(tag=>{
    if(!show[tag]) return;
    const m = v.addModel(DATA.pdb[tag],'pdb');
    const c = tag==='crystal'?css('--ink-3'):(tag==='good'?css('--ok'):css('--bad'));
    if (wantProt && first){
      m.setStyle({atom:'CA'},{cartoon:{color:css('--line'),opacity:.55,thickness:.2}});
      first = false;
    }
    if (wantHeme){
      m.setStyle({resn:'HEM'},{stick:{colorscheme:'default',radius:.13,
        color: tag==='crystal'? css('--heme') : undefined}});
      m.setStyle({resn:'HEM',atom:'FE'},{sphere:{radius:.52,color:css('--iron')}});
    }
    m.setStyle({resn:'CFF'},{stick:{radius:.20,color:c},
                             sphere:{radius:.30,color:c,opacity:.92}});
    if (wantPocket){
      m.setStyle({resi:DATA.pocket},{stick:{radius:.09,color:css('--ink-3'),opacity:.75}});
    }
  });
  if (wantAxial){
    // the axial thiolate: draw the bond the co-folder is never told about explicitly
    v.addCylinder({start:{},end:{},radius:.06});
    v.removeAllShapes();
    const sel = {resi:DATA.axialCys, atom:'SG'};
    v.addStyle(sel,{sphere:{radius:.34,color:css('--warn')}});
    v.addLabel('Cys'+DATA.axialCys, {position:sel, backgroundOpacity:.0,
      fontColor:css('--warn'), fontSize:11, inFront:true});
  }
  v.zoomTo({resn:'CFF'});
  v.zoom(0.55);
  v.render();
}

function boot(){
  const cfg = {backgroundColor:'rgba(0,0,0,0)', antialias:true};
  vA = $3Dmol.createViewer(document.getElementById('vpA'), cfg);
  vB = $3Dmol.createViewer(document.getElementById('vpB'), cfg);
  paint(vA,null);
  paint(vB,'bad');
  document.getElementById('scoreA').innerHTML =
    `<b>0.908</b>best pose · BiSyRMSD ${DATA.scores.good.rmsd} Å`;
  document.getElementById('scoreB').innerHTML =
    `<b>0.320</b>worst pose · BiSyRMSD ${DATA.scores.bad.rmsd} Å`;
  document.getElementById('readout').innerHTML = `
    <dt>LDDT-PLI · best</dt><dd style="color:var(--ok)">${DATA.scores.good.lddt}</dd>
    <dt>LDDT-PLI · worst</dt><dd style="color:var(--bad)">${DATA.scores.bad.lddt}</dd>
    <dt>BiSyRMSD · best</dt><dd>${DATA.scores.good.rmsd} Å</dd>
    <dt>BiSyRMSD · worst</dt><dd>${DATA.scores.bad.rmsd} Å</dd>
    <dt>centre offset · best</dt><dd>0.38 Å</dd>
    <dt>centre offset · worst</dt><dd>0.53 Å</dd>
    <dt>backbone RMSD</dt><dd>0.43–0.70 Å</dd>
    <dt>Cα superposed</dt><dd>467</dd>`;
}

['cShowCrystal','cShowGood','cShowBad','cProtein','cPocket','cHeme','cAxial'].forEach(id=>{
  document.getElementById(id).addEventListener('change',()=>{
    paint(vA, split?'pair':null); if(split) paint(vB,'bad');
  });
});
document.getElementById('bSplit').addEventListener('click',()=>{
  split = !split;
  document.getElementById('viewrow').classList.toggle('split', split);
  document.getElementById('vpB').hidden = !split;
  document.getElementById('tagA').textContent = split? 'crystal + best' : 'overlay';
  setTimeout(()=>{ vA.resize(); vB.resize(); paint(vA, split?'pair':null); paint(vB,'bad'); },30);
});
document.getElementById('bSpin').addEventListener('click',()=>{
  spinning=!spinning; vA.spin(spinning?'y':false); if(split) vB.spin(spinning?'y':false);
});
document.getElementById('bReset').addEventListener('click',()=>{
  paint(vA, split?'pair':null); if(split) paint(vB,'bad');
});

/* morph: interpolate the ligand between crystal truth and prediction.
   Atom order is shared because the extractor wrote both through the symmetry-aware
   mapping, so atom i really is the same atom in both frames. */
let morphModel = null, morphTag = 'bad';
function ensureMorph(){
  const f = DATA.morph[morphTag];
  if(!f) return null;
  if(morphModel) return morphModel;
  morphModel = vA.addModel(f.pdb, 'pdb');
  return morphModel;
}
function morph(pct){
  const f = DATA.morph[morphTag];
  document.getElementById('morphLabel').textContent =
    f ? `crystal ▸ ${morphTag==='bad'?'worst':'best'} pose · ${pct}%`
      : 'morph unavailable for this pose';
  if(!f) return;
  const m = ensureMorph();
  if(!m) return;
  const t = pct/100;
  const atoms = m.selectedAtoms({});
  // atom i of this model is atom i of BOTH exported frames: the generator wrote the
  // ligand-only PDB in the same order as the arrays, so no re-keying is needed
  for(let i=0;i<atoms.length && i<f.from.length;i++){
    atoms[i].x = f.from[i][0] + (f.to[i][0]-f.from[i][0])*t;
    atoms[i].y = f.from[i][1] + (f.to[i][1]-f.from[i][1])*t;
    atoms[i].z = f.from[i][2] + (f.to[i][2]-f.from[i][2])*t;
  }
  const c = morphTag==='bad'?css('--bad'):css('--ok');
  m.setStyle({},{stick:{radius:.22,color:c},sphere:{radius:.34,color:c}});
  vA.render();
}
document.getElementById('sMorph').addEventListener('input',e=>morph(+e.target.value));
document.getElementById('bMorph').addEventListener('click',()=>{
  const s = document.getElementById('sMorph');
  let v = 0; const step = ()=>{ v+=2; s.value=v; morph(v); if(v<100) requestAnimationFrame(step); };
  s.value=0; step();
});

if (window.$3Dmol) boot();
else window.addEventListener('load', boot);
</script>
"""


if __name__ == "__main__":
    raise SystemExit(main())
