"""Exact local number-preserving basis change for direct CGS winding QNs.
All eight matter states are retained; no clock truncation or effective theory.
"""
import json,pathlib,time
import numpy as np
from ed import W
from provenance import provenance
started=time.time();omega=np.exp(2j*np.pi/3);B=np.zeros((8,8),complex);B[0,0]=B[7,7]=1
for a in range(3):
 for b in range(3):
  B[1<<b,1<<a]=omega**((-a*b)%3)/np.sqrt(3)
  B[7^(1<<b),7^(1<<a)]=omega**((a*b)%3)/np.sqrt(3)
N=np.diag([b.bit_count() for b in range(8)]);Q=np.array([sum(a*((b>>a)&1) for a in range(3))%3 for b in range(8)])
assert np.linalg.norm(B.conj().T@B-np.eye(8))<1e-13
assert np.linalg.norm(B@N-N@B)<1e-13
records=[]
for c in range(3):
 for q in range(3):
  U=np.zeros((8,8),complex)
  for old in range(8):
   new=sum(((old>>a)&1)<<((a+q)%3) for a in range(3));U[new,old]=omega**((c*old.bit_count())%3)
  D=np.diag([omega**((c*b.bit_count()+q*Q[b])%3) for b in range(8)])
  err=np.linalg.norm(B.conj().T@U@B-D);assert err<1e-13
  records.append(dict(c=c,q=q,diagonalization_error=float(err)))
operators=[]
for i in range(3):
 A=np.zeros((8,8),complex)
 for a in range(3):
  for old in range(8):
   if (old>>a)&1:A[old^(1<<a),old]+=W[a,i]
 rotated=B.conj().T@A@B
 forbidden=[];terms=[]
 for out in range(8):
  for old in range(8):
   allowed=out.bit_count()==old.bit_count()-1 and (Q[out]-Q[old]+i)%3==0
   if not allowed:forbidden.append(abs(rotated[out,old]))
   if abs(rotated[out,old])>1e-14:
    assert allowed
    terms.append(dict(out=out,input=old,real=float(rotated[out,old].real),imag=float(rotated[out,old].imag)))
 err=max(forbidden,default=0);assert err<1e-13
 operators.append(dict(local_leg_i=i+1,nonzero_matrix_units=terms,forbidden_charge_amplitude=float(err)))
r=dict(basis_convention='ED bit 1 is up; B columns are transformed computational states in the physical matter basis',
 basis_real=B.real.tolist(),basis_imag=B.imag.tolist(),local_matter_states_retained=8,
 number_preserving_error=float(np.linalg.norm(B@N-N@B)),unitarity_error=float(np.linalg.norm(B.conj().T@B-np.eye(8))),
 transformed_number_charge_by_bitstring=Q.tolist(),local_symmetry_records=records,rotated_exchange_operators=operators,
 runtime_seconds=time.time()-started,
 audit=provenance(None,'Independent exact local eight-state basis algebra',dict(tolerance=1e-13,matrix_unit_roundoff_cutoff=1e-14),'All matter kets and all nine CGS actions',['physical matter N_up']),
 conclusion='After this exact local basis rotation, any specified microscopic CGS cycle is a product of onsite Z3 number rotations with matter weights c+q*a and gauge weights p_e. This permits direct MPS winding QNs without losing physical states; production implementation and ED/MPS validation remain pending.')
pathlib.Path('results/matter_charge_basis_audit.json').write_text(json.dumps(r,indent=2)+'\n')
print('All nine CGS actions diagonalized exactly; exchange selection rules passed without matter-space truncation')
