"""Plot exact certificate frequencies and measured phase profiles."""
from pathlib import Path
import argparse,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output-dir',type=Path,default=Path('plots'));args=parser.parse_args()
    output=args.output_dir.resolve();output.mkdir(parents=True,exist_ok=True)
    profiles=json.loads((Path(__file__).resolve().parent/'data/independent_validation_v1/artifacts_validation.json').read_text())['profiles']
    sizes=[p['r'] for p in profiles]
    phases=['reduction','real_roots','coordinates','propagate','rational_candidate','other']
    labels=['Reduction','Real roots','Coordinates','Propagation','Final checks','Other']
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(10.2,3.25),layout='constrained');x=np.arange(len(sizes));unique=np.array([p['unique'] for p in profiles])
    axes[0].bar(x,unique,color='#28618b',label='Certified singleton');axes[0].bar(x,16-unique,bottom=unique,color='#d4d9df',label='Unresolved')
    axes[0].set(xticks=x,xticklabels=sizes,xlabel='Six-cycle count r',ylabel='Inputs (16 at each size)',ylim=(0,18));axes[0].legend(frameon=False,loc='upper center',fontsize=8,ncol=2)
    for i,v in enumerate(unique):axes[0].text(i,16.4,str(v)+'/16',ha='center',fontsize=8)
    bottom=np.zeros(len(sizes));colors=['#235789','#4c9f70','#f0a202','#b86b77','#7460a8','#b5b9bd']
    for key,label,color in zip(phases,labels,colors):
        values=np.array([p['phases'][key] for p in profiles]);axes[1].bar(x,values,bottom=bottom,label=label,color=color);bottom+=values
    axes[1].set(xticks=x,xticklabels=sizes,xlabel='Six-cycle count r',ylabel='Mean median-run decode time (s)');axes[1].legend(frameon=False,fontsize=7,ncol=2,loc='upper left')
    axes[0].set_title('Exact certificates before a product join');axes[1].set_title('Cost by algorithmic stage')
    for ext in ['pdf','png']:fig.savefig(output/('coordinate_mechanism.'+ext),dpi=180)
    plt.close(fig)
if __name__=='__main__':main()
