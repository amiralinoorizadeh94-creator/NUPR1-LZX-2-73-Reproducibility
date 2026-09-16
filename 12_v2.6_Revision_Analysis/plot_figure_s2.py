from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
src=Path(__file__).resolve().parent/'Figure_S2_Rep1_MIC_interfragment_source_v2.6.csv'
blk=Path(__file__).resolve().parent/'Figure_S2_Rep1_5ns_block_means_v2.6.csv'
df=pd.read_csv(src); b=pd.read_csv(blk)
fig,axs=plt.subplots(2,1,figsize=(9,6.2),sharex=True,constrained_layout=True)
axs[0].plot(df.time_ns,df.mic_min_atom_atom_A,lw=.65,label='MIC minimum atom–atom distance')
mid=(b.start_ns+b.end_ns)/2
axs[0].plot(mid,b.mean_mic_min_atom_atom_A,lw=2.0,marker='o',label='5-ns block mean')
axs[0].set_ylabel('Distance (Å)'); axs[0].legend(frameon=False)
axs[1].plot(df.time_ns,df.mic_centroid_separation_A,lw=.7,label='MIC centroid separation')
axs[1].set_ylabel('Separation (Å)'); axs[1].set_xlabel('Time (ns)'); axs[1].legend(frameon=False)
fig.savefig(Path(__file__).resolve().parent/'Figure_S2_REGENERATED.png',dpi=300)
