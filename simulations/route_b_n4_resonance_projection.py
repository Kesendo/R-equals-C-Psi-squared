"""Exact first-order resonant projection at N4 end weights sqrt(3)/2.

Gamma=1, XY, qCSharp; commutator frequencies measured in units q.
"""
import json,hashlib
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix as DM
from route_b_other_n_unfolding import parts


def main():
    field=s.QQ.algebraic_field(s.sqrt(3));a0=field.from_sympy(s.sqrt(3)/2)
    mat=lambda a:DM.from_Matrix(s.Matrix(a)).convert_to(field)
    d,c,left,right,reflection=parts(4)
    def integer(a):
        assert all(v.imag==0 and int(v.real)==v.real for row in a for v in row)
        return [[int(v.real) for v in row] for row in a]
    D=mat(integer(d));TL=mat(integer(left/1j));TR=mat(integer(right/1j))
    T=mat(integer((c-left-right)/1j))+(TL+TR).scalarmul(a0)
    I=DM.eye(24,field).to_dense()
    U=(T-I).nullspace().transpose().to_dense()
    assert U.shape==(24,5)
    W=(U.transpose()*U).inv()*U.transpose()
    R=mat([[int(i==reflection[j]) for j in range(24)] for i in range(24)])
    assert (R*U-U).is_zero_matrix
    DD=(W*D*U).to_Matrix()
    TT=(W*(TL+TR)*U).to_Matrix()
    u,v=s.symbols('u v',real=True)
    effective=DD+s.I*u/(4*s.sqrt(3))*TT
    polynomial=s.factor(((-4+s.I*v)*s.eye(5)-effective).det(method='domain-ge'))
    H=(u+4*v)*(u-4*v)**2+40*v-8*u
    assert s.expand(polynomial-s.I*((u-4*v)**2+64)*H/1024)==0
    ranks=[]
    for uu,vv in [(-2*s.sqrt(2),-s.sqrt(2)/4),(2*s.sqrt(2),s.sqrt(2)/4),(-5*s.sqrt(5)/4,-s.sqrt(5)/16),(5*s.sqrt(5)/4,s.sqrt(5)/16)]:
        numberfield=s.QQ.algebraic_field(s.I,uu,s.sqrt(3))
        M=((-4+s.I*vv)*s.eye(5)-effective.subs(u,uu)).applyfunc(s.simplify)
        rank=DM.from_Matrix(M).convert_to(numberfield).rank()
        assert rank==4
        ranks.append(dict(u=str(uu),v=str(vv),rank=rank))
    assert s.simplify(polynomial.subs({u:0,v:0}))==0
    # If the dissipator is replaced by a scalar, the same four points cease
    # to be roots: the Hamiltonian frequencies alone do not create these EPs.
    scalar_control=s.factor((s.I*v*s.eye(5)-s.I*u/(4*s.sqrt(3))*TT).det(method='domain-ge'))
    assert all(s.simplify(scalar_control.subs({u:s.sympify(row['u']),v:s.sympify(row['v'])}))!=0 for row in ranks)
    print('projected determinant:',polynomial,flush=True)
    print('quotient by H:',s.cancel(polynomial/H),flush=True)
    # Single-particle Hamiltonian in the same hopping convention.
    h=s.Matrix([[0,s.sqrt(3),0,0],[s.sqrt(3),0,2,0],[0,2,0,s.sqrt(3)],[0,0,s.sqrt(3),0]])
    eigenvalues=h.eigenvals();assert eigenvalues=={-3:1,-1:1,1:1,3:1}
    out=dict(scope=__doc__,resonant_dimension=5,reflection='+I',D_projection=str(DD),end_hop_projection=str(TT),
             effective_determinant=str(polynomial),residual_cubic=str(s.expand(H)),quotient=str(s.cancel(polynomial/H)),
             single_particle_energies=[-3,-1,1,3],script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    out.update(effective_ranks=ranks,scalar_dissipator_control_rejected=True)
    (Path(__file__).parent/'results/route_b_n4_resonance_projection.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
