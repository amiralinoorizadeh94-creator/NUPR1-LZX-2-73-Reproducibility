from pathlib import Path
import re, csv, numpy as np
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics import silhouette_score

ROOT=Path(__file__).resolve().parents[1]/'11_Sensitivity_Analysis'/'Historical_Raw_Workspace_Core'
REC=ROOT/'refined_receptors'; OUT=ROOT/'refinement_docking_out'
OUTDIR=Path(__file__).resolve().parent; OUTDIR.mkdir(exist_ok=True)

def atoms_pdbqt(path, model1=False):
    rows=[]; in_model=not model1
    for line in open(path, errors='ignore'):
        if model1 and line.startswith('MODEL'): in_model=line.split()[1]=='1'; continue
        if model1 and line.startswith('ENDMDL') and in_model: break
        if not in_model: continue
        if line.startswith(('ATOM  ','HETATM')):
            name=line[12:16].strip(); x=float(line[30:38]); y=float(line[38:46]); z=float(line[46:54])
            atype=line[77:].strip().split()[0] if len(line)>77 and line[77:].strip() else ''
            rows.append((name,atype,np.array([x,y,z],float)))
    return rows

def backbone(path):
    a=atoms_pdbqt(path)
    xyz=[r[2] for r in a if r[0] in {'N','CA','C'}]
    if len(xyz)!=30: raise ValueError((path,len(xyz)))
    return np.array(xyz)

def ligand_heavy(path):
    a=atoms_pdbqt(path, model1=True)
    xyz=[r[2] for r in a if r[0] != 'H' and not r[1].startswith('H')]
    if len(xyz)!=24: raise ValueError((path,len(xyz)))
    return np.array(xyz)

def kabsch_transform(P,Q): # map P -> Q
    pc=P.mean(0); qc=Q.mean(0); X=P-pc; Y=Q-qc
    U,S,Vt=np.linalg.svd(X.T@Y); R=U@Vt
    if np.linalg.det(R)<0:
        Vt[-1]*=-1; R=U@Vt
    t=qc-pc@R
    return R,t

names=sorted(p.stem for p in REC.glob('COMBO_*_R*.pdbqt'))
ref='COMBO_011_R28'; Q=backbone(REC/f'{ref}.pdbqt')
aligned=[]
for name in names:
    P=backbone(REC/f'{name}.pdbqt'); R,t=kabsch_transform(P,Q)
    L=ligand_heavy(OUT/f'{name}_out.pdbqt')
    aligned.append(L@R+t)
N=len(names); D=np.zeros((N,N))
for i in range(N):
    for j in range(i):
        d=np.sqrt(np.mean(np.sum((aligned[i]-aligned[j])**2,axis=1))); D[i,j]=D[j,i]=d
cond=squareform(D,checks=False); Z=linkage(cond,method='average')
scores=[]
for k in range(2,11):
    lab=fcluster(Z,k,criterion='maxclust')
    sc=silhouette_score(D,lab,metric='precomputed'); scores.append((k,sc))
bestk=max(scores,key=lambda x:x[1])[0]; labels=fcluster(Z,bestk,criterion='maxclust')
# relabel by descending population, tie by old label
old=sorted(set(labels), key=lambda c:(-np.sum(labels==c),c)); mp={c:i+1 for i,c in enumerate(old)}; labels=np.array([mp[x] for x in labels])
medoids={}
for c in sorted(set(labels)):
    idx=np.where(labels==c)[0]; sub=D[np.ix_(idx,idx)]; medoids[c]=idx[np.argmin(sub.mean(1))]
with open(OUTDIR/'refined_pose_clusters.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['model','cluster','cluster_size','is_medoid','is_COMBO_011_R28'])
    for i,n in enumerate(names):
        c=int(labels[i]); w.writerow([n,c,int(np.sum(labels==c)),int(i==medoids[c]),int(n==ref)])
with open(OUTDIR/'refined_pose_rmsd_matrix.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['model']+names)
    for n,row in zip(names,D): w.writerow([n]+[f'{x:.6f}' for x in row])
with open(OUTDIR/'refined_pose_clustering_silhouette.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['k','silhouette_score']); w.writerows([(k,f'{s:.8f}') for k,s in scores])
print('best',bestk,'scores',scores)
for c in sorted(set(labels)):
    print(c,np.sum(labels==c),names[medoids[c]], 'R28' if labels[names.index(ref)]==c else '')
