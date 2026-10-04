"""Audit full-state entropy data; fits are descriptive until convergence passes."""
import glob,json,pathlib,csv,itertools
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

rows=[];data=[]
for file in sorted(glob.glob('results/*_L*_w*_chi*.json')):
    r=json.load(open(file))
    if 'records' not in r or 'family' not in r or 'length' not in r: continue
    data.append((file,r))
    for k,q in enumerate(r['records']):
        ee=q.get('sweep_energies',[]);tt=q.get('sweep_max_truncation_errors',[])
        rows.append({'file':file,'family':r['family'],'L':r['length'],'w':r['width'],
            'Ly':r['physical_circumference'],'spins':r['spins'],'ordering':r.get('ordering','axial'),
            'stage':k,'chi':q['cap'],'E':q['energy'],'S':q['entropy'],
            'last_sweep_dE':ee[-1]-ee[-2] if len(ee)>1 else '',
            'max_truncation_last_sweep':tt[-1] if tt else '',
            'refinement_dS':q.get('entropy_change_refinement','')})
if rows:
    with open('results/entropy_convergence.csv','w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n");writer.writeheader();writer.writerows(rows)

fits=[]
for family,ordering in itertools.product(['zigzag','armchair'],['axial','star']):
    for cap in [32,64,128,256,512]:
        selected=[]
        for w in [1,2,3]:
            candidates=[]
            for file,r in data:
                if r['family']==family and r['length']==4 and r['width']==w and r.get('ordering','axial')==ordering:
                    stages=[q for q in r['records'] if q['cap']==cap]
                    if stages:candidates.append((r['physical_circumference'],stages[-1]['entropy']))
            if candidates:selected.append(candidates[-1])
        if len(selected)==3:
            x,y=np.array(selected).T
            alpha,b=np.polyfit(x,y,1)
            fits.append({'family':family,'ordering':ordering,'length':4,'chi':cap,'alpha':float(alpha),
                'gamma_apparent':float(-b),'residual_rms':float(np.sqrt(np.mean((y-alpha*x-b)**2))),
                'status':'Descriptive only: narrow circumference, unconverged chi, untested length and sectors'})
pathlib.Path('results/apparent_gamma_fits.json').write_text(json.dumps(fits,indent=2))

fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
colors={32:'#aaaaaa',64:'#7b9acc',128:'#2266aa',256:'#cc7722',512:'#a52a2a'}
for ax,family in zip(axes,['zigzag','armchair']):
    for cap,ordering in itertools.product([32,64,128,256,512],['axial','star']):
        pts=[]
        for w in [1,2,3]:
            candidate=[q for q in rows if q['family']==family and q['L']==4 and q['w']==w and q['chi']==cap and q['ordering']==ordering]
            if candidate:pts.append((candidate[-1]['Ly'],candidate[-1]['S']))
        if pts:
            x,y=np.array(pts).T;ax.plot(x,y,'o'+('--' if ordering=='axial' else '-'),color=colors[cap],label=f'χ={cap}, {ordering}')
    ax.set(title=family+'; L=4 (exploratory)',xlabel='Circumference / nearest-neighbor distance',ylabel='Full-state entropy (nats)')
    ax.legend(fontsize=8);ax.grid(alpha=.2)
fig.suptitle('Exploratory finite-cylinder data; convergence required before topological inference',fontsize=10)
fig.savefig('results/entropy_circumference.png',dpi=180)
fig.savefig('results/entropy_circumference.pdf')
plt.close(fig)

fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
for ax,family in zip(axes,['zigzag','armchair']):
    for file,r in data:
        if r['family']!=family:continue
        q=r['records'];ax.plot([v['cap'] for v in q],[v['entropy'] for v in q],'.-',label=f"L={r['length']}, w={r['width']}, {r.get('ordering','axial')}")
    ax.set_xscale('log',base=2);ax.set(title=family,xlabel='Bond-dimension cap χ',ylabel='Full-state entropy (nats)')
    ax.legend(fontsize=6);ax.grid(alpha=.2)
fig.savefig('results/entropy_bond_convergence.png',dpi=180);fig.savefig('results/entropy_bond_convergence.pdf')
plt.close(fig)
print(json.dumps(fits,indent=2))
