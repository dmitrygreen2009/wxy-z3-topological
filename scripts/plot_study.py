"""Reproducible exploratory plots of raw microscopic convergence data."""
import json,pathlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from provenance import provenance
plt.rcParams.update({'font.size':10,'savefig.dpi':160})
records=json.load(open('results/entanglement_spectra.json'))['states']
selected={}
for p in records:
    if p['ordering']!='star':continue
    key=(p['family'],p['length'],p['width'],p['cap'])
    if key not in selected or p['energy']<selected[key]['energy']:selected[key]=p
points=list(selected.values())
fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
for ax,family in zip(axes,['zigzag','armchair']):
    for L in sorted(set(p['length'] for p in points if p['family']==family)):
        for cap in [128,256,512,1024]:
            group=sorted([p for p in points if p['family']==family and p['length']==L and p['cap']==cap],key=lambda p:p['width'])
            if len(group)<2:continue
            ax.plot([p['circumference'] for p in group],[p['spectrum']['entropy'] for p in group],'.--',label=f'L={L}, cap={cap}')
    ax.set(title=family,xlabel='Physical circumference (nearest-neighbor distance = 1)',ylabel='Full-state spatial entropy')
    ax.legend(fontsize=7)
fig.suptitle('Lowest recorded finite energies per width; number/flux branches may differ')
for suffix in ['png','pdf']:fig.savefig('results/full_state_circumference.'+suffix)
plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
for ax,family in zip(axes,['zigzag','armchair']):
    for L,w in sorted(set((p['length'],p['width']) for p in points if p['family']==family)):
        if L!=4 and w!=1:continue
        candidates=[p for p in records if p['ordering']=='star' and p['family']==family and p['length']==L and p['width']==w]
        for branch in sorted(set((p['nup'],p['initialization_branch']) for p in candidates)):
            raw=[p for p in candidates if (p['nup'],p['initialization_branch'])==branch]
            group=[min([p for p in raw if p['cap']==cap],key=lambda p:p['energy']) for cap in sorted(set(p['cap'] for p in raw))]
            if len(group)<2:continue
            label=f'L={L}, w={w}, Nup={branch[0]}, '+branch[1].split(f'_w{w}_',1)[-1]
            ax.plot([p['cap'] for p in group],[p['spectrum']['entropy'] for p in group],'.-',label=label)
    ax.set(title=family,xlabel='Maximum bond dimension',ylabel='Full-state spatial entropy');ax.set_xscale('log',base=2)
    ax.legend(fontsize=7)
fig.suptitle('Bond dependence within separate initialization and number branches')
for suffix in ['png','pdf']:fig.savefig('results/full_state_bond_convergence.'+suffix)
plt.close(fig)
infinite=json.load(open('results/infinite_analysis.json'))['selected_measurements']
fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
keys=set((p['family'],p['width'],p['ordering'],p['conserving_ansatz'],p['branch_id']) for p in infinite)
for key in sorted(keys,key=str):
    group=sorted([p for p in infinite if (p['family'],p['width'],p['ordering'],p['conserving_ansatz'],p['branch_id'])==key],key=lambda p:p['cap'])
    label=f'{key[0]} w={key[1]} {key[2]} U1={key[3]} {key[4]}'
    if any(p['energy_above_known_trial_state'] for p in group):label+=' *'
    axes[0].plot([p['actual_chi'] for p in group],[p['entropy'] for p in group],'.-',label=label)
    good=[p for p in group if isinstance(p['xi_cells'],(int,float)) and np.isfinite(p['xi_cells']) and p['xi_cells']>0]
    axes[1].plot([p['actual_chi'] for p in good],[2*p['xi_cells'] for p in good],'.-',label=label)
axes[0].set(xlabel='Achieved bond dimension',ylabel='Full-state spatial entropy')
axes[1].set(xlabel='Achieved bond dimension',ylabel='Transfer correlation length (axial slices)')
for ax in axes:ax.set_xscale('log',base=2);ax.legend(fontsize=6)
axes[1].set_yscale('log')
fig.suptitle('Separate infinite variational branches; * includes energy above a trial bound',fontsize=11)
for suffix in ['png','pdf']:fig.savefig('results/infinite_convergence.'+suffix)
plt.close(fig)
pathlib.Path('results/plot_study_audit.json').write_text(json.dumps(dict(audit=provenance(None,'Matplotlib saved-data visualization',{},'Recorded full-state spectra and infinite measurements',[]),
    sources=['results/entanglement_spectra.json','results/infinite_analysis.json'],interpretation='Exploratory plots; line connections do not certify convergence.'),indent=2))
print('Saved finite entropy and infinite entropy/correlation-length plots.')
