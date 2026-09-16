from pathlib import Path
import numpy as np, pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array, minimize_vectors
import sys
BASE=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent/'Raw_Rep1'
top=BASE/'COMBO_011_R28_OPC_015M_HMR.prmtop'
dcds=[BASE/'COMBO_011_R28_PROD_10ns.dcd',BASE/'COMBO_011_R28_PROD_10to50ns.dcd']
out=Path(__file__).resolve().parent; out.mkdir(exist_ok=True)
u=mda.Universe(str(top),*[str(x) for x in dcds])
a=u.select_atoms('resid 1:5'); b=u.select_atoms('resid 6:10')
rows=[]
for i,ts in enumerate(u.trajectory):
    box=ts.dimensions.copy()
    micmin=float(distance_array(a.positions,b.positions,box=box).min())
    rawmin=float(distance_array(a.positions,b.positions).min())
    ca=a.center_of_geometry(); cb=b.center_of_geometry(); v=cb-ca
    micv=minimize_vectors(v.reshape(1,3),box)[0]
    miccent=float(np.linalg.norm(micv)); rawcent=float(np.linalg.norm(v))
    rows.append((i,(i)*0.01,micmin,miccent,rawmin,rawcent,*box[:3]))
df=pd.DataFrame(rows,columns=['frame','time_ns','mic_min_atom_atom_A','mic_centroid_separation_A','raw_min_atom_atom_A','raw_centroid_separation_A','box_a_A','box_b_A','box_c_A'])
df.to_csv(out/'Figure_S2_Rep1_MIC_interfragment_source_v2.6.csv',index=False)
df['block_5ns']=(df['time_ns']//5).astype(int)+1
blk=df.groupby('block_5ns').agg(start_ns=('time_ns','min'),end_ns=('time_ns','max'),mean_mic_min_atom_atom_A=('mic_min_atom_atom_A','mean'),mean_mic_centroid_separation_A=('mic_centroid_separation_A','mean')).reset_index()
blk.to_csv(out/'Figure_S2_Rep1_5ns_block_means_v2.6.csv',index=False)
print('n',len(df),'min summary',df.mic_min_atom_atom_A.mean(),df.mic_min_atom_atom_A.median(),df.mic_min_atom_atom_A.std(ddof=1),df.mic_min_atom_atom_A.min(),df.mic_min_atom_atom_A.max())
print('cent',df.mic_centroid_separation_A.mean(),df.mic_centroid_separation_A.median(),df.mic_centroid_separation_A.std(ddof=1),df.mic_centroid_separation_A.min(),df.mic_centroid_separation_A.max())
print('diff counts',np.sum(np.abs(df.raw_min_atom_atom_A-df.mic_min_atom_atom_A)>1e-6), np.max(np.abs(df.raw_min_atom_atom_A-df.mic_min_atom_atom_A)), np.sum(np.abs(df.raw_centroid_separation_A-df.mic_centroid_separation_A)>1e-6),np.max(np.abs(df.raw_centroid_separation_A-df.mic_centroid_separation_A)))
print('blocks',blk.mean_mic_min_atom_atom_A.round(3).tolist())
