"""Fit recorded full-state entropies. Raw states and conclusions stay separate."""
import glob,json,pathlib,itertools
import numpy as np
from finite_results import finite_groups
LN3=float(np.log(3))

def fit(points):
    x=np.array([p['circumference'] for p in points]);y=np.array([p['entropy'] for p in points])
    design=np.column_stack([x,-np.ones(len(x))])
    alpha,gamma=np.linalg.lstsq(design,y,rcond=None)[0]
    residuals=y-design@np.array([alpha,gamma]);dof=len(x)-2
    covariance=np.linalg.inv(design.T@design)*sum(residuals**2)/dof if dof>0 else None
    fixed_alpha=float(x@(y+LN3)/(x@x));fixed_residuals=y-fixed_alpha*x+LN3
    return dict(circumferences=x.tolist(),widths=[p['width'] for p in points],
        alpha=float(alpha),gamma=float(gamma),gamma_apparent=float(gamma),
        gamma_minus_ln3=float(gamma-LN3),gamma_standard_error=float(np.sqrt(covariance[1,1])) if covariance is not None else None,
        alpha_standard_error=float(np.sqrt(covariance[0,0])) if covariance is not None else None,
        covariance=covariance.tolist() if covariance is not None else None,
        residuals=residuals.tolist(),residual_rms=float(np.sqrt(np.mean(residuals**2))),degrees_of_freedom=dof,
        fixed_ln3_alpha=fixed_alpha,fixed_ln3_residuals=fixed_residuals.tolist(),
        uncertainty_scope='OLS scatter only; unquantified length, chi, sectors, and circumference corrections prevent a thermodynamic confidence interval.')

def main():
    points=[]
    for path,r,stages in finite_groups():
        for cap in sorted(set(q['cap'] for k,q in stages)):
            index,q=[(k,q) for k,q in stages if q['cap']==cap][-1]
            points.append(dict(source_file=path,record_index=index,family=r['family'],length=r['length'],width=r['width'],
                ordering=r.get('ordering','axial'),nup=r.get('nup',(r['spins']+1)//2),number_sector_scope='Global cylinder minimum requires an independent sector search',chi=cap,actual_chi=q['maxlinkdim'],energy=q['energy'],entropy=q['entropy'],
                circumference=r['physical_circumference'],seed=r.get('seed'),truncation_error=q.get('sweep_max_truncation_errors',[None])[-1]))
    fits=[]
    for family,ordering,L,cap in sorted(set((p['family'],p['ordering'],p['length'],p['chi']) for p in points)):
        candidates=[p for p in points if (p['family'],p['ordering'],p['length'],p['chi'])==(family,ordering,L,cap)]
        selected=[min([p for p in candidates if p['width']==w],key=lambda p:p['energy']) for w in sorted(set(p['width'] for p in candidates))]
        if len(selected)<3:continue
        result=fit(selected)
        result.update(family=family,ordering=ordering,length=L,chi=cap,raw_points=selected,
            smallest_circumference_removed=fit(selected[1:]),
            eligible_for_topological_inference=False,selection_rule='Lowest recorded variational energy at each width; a common MES branch and global ground-sector selection are not certified.',status='Descriptive finite-size intercept only; current data do not supply a topological gamma confidence interval.')
        result['smallest_removed_gamma_change']=result['smallest_circumference_removed']['gamma']-result['gamma']
        fits.append(result)
    pathlib.Path('results/entropy_fits.json').write_text(json.dumps(dict(ln3=LN3,fits=fits,accepted_topological_fits=[],scientific_strategy_reaudit="results/scientific_strategy_reaudit.md"),indent=2))
    print('Wrote',len(fits),'fits with residuals, OLS errors and smallest-circumference sensitivity.')

if __name__=="__main__":main()
