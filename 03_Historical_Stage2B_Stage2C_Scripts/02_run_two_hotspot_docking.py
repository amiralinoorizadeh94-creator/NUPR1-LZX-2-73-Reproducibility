
from pathlib import Path
import subprocess, re, csv, statistics

vina=Path("vina.exe")
lig=Path("LZX-2-73.pdbqt")
if not vina.exists(): raise FileNotFoundError("vina.exe missing")
if not lig.exists(): raise FileNotFoundError("LZX-2-73.pdbqt missing")

def coords(path):
    pts=[]
    for line in Path(path).read_text(errors="ignore").splitlines():
        if line.startswith(("ATOM","HETATM")):
            try:
                pts.append((float(line[30:38]),float(line[38:46]),float(line[46:54])))
            except:
                pass
    return pts

outdir=Path("two_hotspot_docking_out")
outdir.mkdir(exist_ok=True)
rows=[]

recs=sorted(Path("combined_receptors").glob("COMBO_*.pdbqt"))
if not recs:
    raise RuntimeError("No combined receptor models found.")

for idx,rec in enumerate(recs,1):
    pts=coords(rec)
    xs=[p[0] for p in pts]; ys=[p[1] for p in pts]; zs=[p[2] for p in pts]
    center=(sum(xs)/len(xs),sum(ys)/len(ys),sum(zs)/len(zs))
    cfg=outdir/(rec.stem+".txt")
    cfg.write_text(f"""receptor = {rec.as_posix()}
ligand = {lig.as_posix()}
center_x = {center[0]:.3f}
center_y = {center[1]:.3f}
center_z = {center[2]:.3f}
size_x = 30
size_y = 30
size_z = 30
exhaustiveness = 8
num_modes = 10
energy_range = 6
""",encoding="utf-8")

    out=outdir/(rec.stem+"_out.pdbqt")
    cp=subprocess.run([str(vina.resolve()),"--config",str(cfg),"--out",str(out)],
                      capture_output=True,text=True)
    txt=(cp.stdout or "")+"\n"+(cp.stderr or "")
    (outdir/(rec.stem+".log")).write_text(txt,encoding="utf-8")
    if cp.returncode!=0:
        print(txt)
        raise SystemExit(cp.returncode)

    best=None
    for line in txt.splitlines():
        m=re.match(r"^\s*1\s+(-?\d+(?:\.\d+)?)\s+",line)
        if m:
            best=float(m.group(1)); break
    rows.append((rec.stem,best))
    print(f"[{idx}/{len(recs)}] {rec.stem} best = {best}")

with open(outdir/"summary.csv","w",newline="",encoding="utf-8") as f:
    w=csv.writer(f)
    w.writerow(["model","best_vina_kcal_mol"])
    w.writerows(rows)

vals=[v for _,v in rows if v is not None]
print("\nSUMMARY")
print("n =",len(vals))
if vals:
    print("best =",min(vals))
    print("median =",statistics.median(vals))
    print("mean =",round(statistics.mean(vals),3))
    for cutoff in [-5,-6,-7,-8]:
        n=sum(v<=cutoff for v in vals)
        print(f"fraction <= {cutoff:.1f} kcal/mol = {n}/{len(vals)} ({100*n/len(vals):.1f}%)")
    print("paper two-fragment reference ≈ -8.2 kcal/mol")
