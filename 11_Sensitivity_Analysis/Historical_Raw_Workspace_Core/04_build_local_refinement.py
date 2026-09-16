from pathlib import Path
import numpy as np, math, random, csv

random.seed(20260911)
np.random.seed(20260911)

SEEDS = ["COMBO_011","COMBO_031","COMBO_030","COMBO_044","COMBO_025"]
SRC = Path("combined_receptors")
OUT = Path("refined_receptors")
OUT.mkdir(exist_ok=True)

def parse_pdbqt(path):
    atoms=[]
    for line in Path(path).read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith(("ATOM","HETATM")):
            try:
                xyz=np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])],float)
            except Exception:
                continue
            chain=line[21:22] if len(line)>21 else " "
            atoms.append((line,xyz,chain))
    if not atoms:
        raise RuntimeError(f"No atoms in {path}")
    return atoms

def axis_angle(axis, angle_deg):
    axis=np.asarray(axis,float)
    axis=axis/np.linalg.norm(axis)
    a=math.radians(angle_deg)
    x,y,z=axis
    c=math.cos(a); s=math.sin(a); C=1-c
    return np.array([
        [c+x*x*C, x*y*C-z*s, x*z*C+y*s],
        [y*x*C+z*s, c+y*y*C, y*z*C-x*s],
        [z*x*C-y*s, z*y*C+x*s, c+z*z*C]
    ])

def rewrite(line, xyz):
    s=list(line)
    while len(s)<80: s.append(" ")
    x,y,z=xyz
    s[30:38]=list(f"{x:8.3f}")
    s[38:46]=list(f"{y:8.3f}")
    s[46:54]=list(f"{z:8.3f}")
    return "".join(s).rstrip()

def min_dist(A,B):
    d=A[:,None,:]-B[None,:,:]
    return np.sqrt((d*d).sum(axis=2)).min()

meta=[]
per_seed=30

for seed in SEEDS:
    path=SRC/f"{seed}.pdbqt"
    if not path.exists():
        raise FileNotFoundError(path)
    atoms=parse_pdbqt(path)
    A=[a for a in atoms if a[2]=="A"]
    B=[a for a in atoms if a[2]=="B"]
    if not A or not B:
        raise RuntimeError(f"{seed}: could not separate chains A/B")

    Axyz=np.array([x[1] for x in A])
    Bxyz=np.array([x[1] for x in B])
    Bcenter=Bxyz.mean(axis=0)

    made=0
    attempts=0
    while made<per_seed and attempts<5000:
        attempts += 1

        axis=np.random.normal(size=3)
        axis/=np.linalg.norm(axis)
        angle=np.random.uniform(-18,18)
        R=axis_angle(axis,angle)

        # small local translation
        trans=np.random.normal(scale=0.75,size=3)
        if np.linalg.norm(trans)>2.0:
            continue

        Bnew=(Bxyz-Bcenter)@R.T + Bcenter + trans

        dmin=min_dist(Axyz,Bnew)
        # reject hard clashes and very separated local environments
        if not (2.2 <= dmin <= 6.5):
            continue

        allxyz=np.vstack([Axyz,Bnew])
        span=allxyz.max(axis=0)-allxyz.min(axis=0)
        if span.max()>27:
            continue

        made += 1
        out=OUT/f"{seed}_R{made:02d}.pdbqt"
        lines=[]
        for line,xyz,_ in A:
            lines.append(rewrite(line,xyz))
        for (line,_,_),xyz in zip(B,Bnew):
            lines.append(rewrite(line,xyz))
        out.write_text("\n".join(lines)+"\n",encoding="utf-8")

        meta.append({
            "model":out.stem,
            "seed":seed,
            "angle_deg":angle,
            "translation_A":float(np.linalg.norm(trans)),
            "min_interfragment_A":float(dmin),
            "span_x":float(span[0]),
            "span_y":float(span[1]),
            "span_z":float(span[2]),
        })

    print(seed, "generated", made, "refined models")

with open(OUT/"refinement_metadata.csv","w",newline="",encoding="utf-8") as f:
    fields=["model","seed","angle_deg","translation_A","min_interfragment_A","span_x","span_y","span_z"]
    w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(meta)

print("Total refined models:",len(meta))
