"""Exact first-order resonant projection at N4 end weights sqrt(3)/2.

Gamma=1, XY, qCSharp; commutator frequencies measured in units q.

The projection also separates the two halves of the generator. Neither can
carry a Jordan block alone: the dephasing half is a diagonal matrix, the
absorption diagonal -2*gamma*hamming(bra, ket), and the coherent half is real
symmetric at EVERY end weight, so a scalar dephasing leaves the block normal at
every coupling. The dephasing does not create the coincidence either: flatten
it to the scalar -4 and the four isolated double roots are replaced by a
fourfold SEMISIMPLE line. What the spread does is resolve that line into
isolated points, four of which keep a Jordan block, and it can do so because
the resonant eigenspace MEETS both absorption strata and therefore carries the
extreme rates -2 and -6 themselves. Straddling the strata would not be enough;
every multiplet of the block straddles them evenly, which is why -4 is the
block's own centre rather than this eigenspace's property.
"""
import json
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix as DM
from route_b_artifact_provenance import artifact_provenance
from route_b_other_n_unfolding import parts
from framework.weight_coherence_block import weight_block_configs


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
        # Rank four is geometric multiplicity one. The DOUBLE root is the
        # factorization above; the two together are the Jordan block.
        control=(s.I*vv*s.eye(5)-s.I*uu/(4*s.sqrt(3))*TT).applyfunc(s.simplify)
        # Exact rank, not simplify(det) != 0: a structural inequality would
        # read an unreduced zero as nonzero.
        control_rank=DM.from_Matrix(control).convert_to(numberfield).rank()
        assert control_rank==5
        ranks.append(dict(u=str(uu),v=str(vv),rank=rank,scalar_dephasing_control_rank=control_rank))
    assert s.simplify(polynomial.subs({u:0,v:0}))==0
    print('projected determinant:',polynomial,flush=True)
    print('quotient by H:',s.cancel(polynomial/H),flush=True)
    # Single-particle Hamiltonian in the same hopping convention.
    h=s.Matrix([[0,s.sqrt(3),0,0],[s.sqrt(3),0,2,0],[0,2,0,s.sqrt(3)],[0,0,s.sqrt(3),0]])
    energies=h.eigenvals();assert energies=={-3:1,-1:1,1:1,3:1}

    # THE COHERENT HALF, at every end weight rather than at this one: the
    # frequency matrix is real symmetric for a SYMBOLIC weight, so -4*I + i*q*T
    # is normal and no coupling and no profile admits a Jordan block from it.
    weight=s.symbols('w',real=True)
    bulk=s.Matrix(integer((c-left-right)/1j));ends=s.Matrix(integer((left+right)/1j))
    assert (bulk-bulk.T).is_zero_matrix and (ends-ends.T).is_zero_matrix
    symbolic=bulk+weight*ends
    assert (symbolic-symbolic.T).is_zero_matrix,'coherent half is symmetric at every end weight'
    # Free-fermion additivity, not the equal spacing, is what makes the block
    # degenerate, and this is a derivation rather than a scan: over Q(w) the
    # characteristic polynomial carries two quadratic factors to the FOURTH
    # power, so a fourfold multiplet stands at every end weight.
    factors=s.factor_list(symbolic.charpoly().as_expr())[1]
    additivity=sorted(((s.Poly(f,s.Symbol('lambda')).degree(),m) for f,m in factors),reverse=True)
    assert [m for _,m in additivity]==[4,4,1,1,1,1] and {deg for deg,_ in additivity}=={2}

    # THE DEPHASING HALF is the absorption diagonal, hence also never defective
    # alone. Its two strata are what the resonant eigenspace can straddle.
    basis=[(a,b) for a in weight_block_configs(4,1) for b in weight_block_configs(4,2)]
    hamming=[bin(a^b).count('1') for a,b in basis]
    Dm=s.Matrix(D.to_Matrix())
    assert Dm==s.diag(*[-2*k for k in hamming]),'dissipator is the absorption diagonal -2*gamma*hamming'
    strata=sorted(set(hamming))
    assert strata==[1,3] and [hamming.count(k) for k in strata]==[12,12]
    rate={str(k):int(Dm[hamming.index(k),hamming.index(k)]) for k in strata}
    Um=U.to_Matrix()
    stratum_ranks=[Um[[i for i in range(24) if hamming[i]==k],:].rank() for k in strata]
    assert stratum_ranks==[4,4]
    # Straddling alone would NOT force a non-scalar restriction; the plane
    # spanned by two balanced sums straddles both and restricts to -4*I.
    balanced=s.zeros(24,2)
    for column,(one,three) in enumerate(zip([i for i in range(24) if hamming[i]==1][:2],
                                            [i for i in range(24) if hamming[i]==3][:2])):
        balanced[one,column]=1;balanced[three,column]=1
    plane=DM.from_Matrix(balanced).convert_to(field).to_dense()
    counterexample=s.simplify(((plane.transpose()*plane).inv()*plane.transpose()*D*plane).to_Matrix())
    assert counterexample==-4*s.eye(2),'straddling is not by itself non-scalarity'
    # What does force it: the resonant eigenspace MEETS each stratum, so the
    # restriction carries D's extreme rates themselves.
    meets={str(k):len(Um[[i for i in range(24) if hamming[i]!=k],:].nullspace()) for k in strata}
    assert meets=={'1':1,'3':1}
    restricted=s.simplify(DD).eigenvals()
    assert sum(restricted.values())==5 and len(restricted)==5
    assert set(restricted)>={s.Integer(-2),s.Integer(-4),s.Integer(-6)}
    mean=s.simplify(sum(k*m for k,m in restricted.items())/5)
    assert mean==-4

    # The -4 is not a property of this eigenspace. EVERY multiplet of the block
    # carries exactly half its weight in each stratum, so any Hamming-indexed
    # dephasing restricts to the mean of its two rung values.
    projector=DM.from_Matrix(s.diag(*[1 if k==1 else 0 for k in hamming])).convert_to(field).to_dense()
    means={};weights={}
    for k in range(-24,25):
        space=(T-I.scalarmul(field.convert(k))).nullspace().transpose().to_dense()
        if space.shape[1]==0:continue
        dual=(space.transpose()*space).inv()*space.transpose()
        means[k]=s.simplify((dual*D*space).to_Matrix().trace()/space.shape[1])
        weights[k]=s.simplify((dual*projector*space).to_Matrix().trace()/space.shape[1])
    assert set(weights.values())=={s.Rational(1,2)},'every multiplet splits evenly between the strata'
    comb={k:24-(T-I.scalarmul(field.convert(k))).rank() for k in range(-24,25)}
    comb={k:m for k,m in comb.items() if m}
    # Geometric multiplicities alone reaching the full dimension is both the
    # whole spectrum and the semisimplicity; the sum can only undercount.
    assert sum(comb.values())==24 and set(means)==set(comb)
    assert set(means.values())=={s.Integer(-4)}
    assert max(comb.values())==5 and comb[1]==5
    # The hand-written chain is not a separate object: its two-excitation minus
    # one-excitation differences ARE this comb, so the equal spacing -3,-1,1,3
    # and the integer comb are one statement.
    single=sorted(energies);pairs=[a+b for i,a in enumerate(single) for b in single[i+1:]]
    from collections import Counter
    assert Counter(int(x-y) for x in pairs for y in single)==Counter(comb),'the comb IS the spacing'
    assert len(comb)==8 and sum(comb.values())==24
    # And the fifth state is the equal spacing itself. The value +1 sits on a
    # multiplicity-four factor only at this end weight, and there, and only
    # there, a SIMPLE factor passes through +1 as well: four plus one.
    here=s.sqrt(3)/2
    at_resonance=[(m,s.solve(s.Eq(s.simplify(f.subs(s.Symbol('lambda'),1)),0),weight)) for f,m in factors]
    assert sum(1 for m,roots in at_resonance if m==4 and here in roots)==1
    assert sum(1 for m,roots in at_resonance if m==1 and here in roots)==1
    quadruple_at=sorted({str(root) for m,roots in at_resonance if m==4 for root in roots if root.is_positive})
    assert quadruple_at==['sqrt(3)/2'],'no other positive end weight gives +1 the fourfold factor'
    # Flattening the dephasing does not remove the coincidence, it un-resolves
    # it: a fourfold root along a whole line, and a semisimple one.
    flattened=s.factor((s.I*v*s.eye(5)-s.I*u/(4*s.sqrt(3))*TT).det(method='domain-ge'))
    assert s.expand(flattened-s.I*(u-4*v)**4*(u+4*v)/1024)==0
    shift=(TT-s.sqrt(3)*s.eye(5)).applyfunc(s.simplify)
    geometric=5-DM.from_Matrix(shift).convert_to(field).rank()
    algebraic=5-DM.from_Matrix((shift**5).applyfunc(s.simplify)).convert_to(field).rank()
    assert algebraic==geometric==4,'the flattened line is degenerate and semisimple'

    # The twin-scalar split of PROOF_CODIM1_BY_ADDITIVITY needs the two
    # coalescing directions to descend from ONE degenerate free-fermion
    # multiplet. They do not. T' sees the four-plus-one split as its own
    # spectrum, the two multiplets' differing slopes, and the coalescing plane
    # lies inside NEITHER of its eigenspaces: a plane contained in one would
    # join it at that eigenspace's own dimension, and every joint rank below
    # exceeds it. The premise fails, so the lemma reaches nothing here.
    slopes=s.simplify(TT).eigenvals()
    assert slopes=={s.sqrt(3):4,-s.sqrt(3):1}
    eigenspaces={str(k):s.Matrix.hstack(*(TT-k*s.eye(5)).applyfunc(s.simplify).nullspace())
                 for k in slopes}
    assert {k:m.cols for k,m in eigenspaces.items()}=={str(k):m for k,m in slopes.items()}
    assert len(ranks)==4
    twin=[]
    for row in ranks:
        uu,vv=s.sympify(row['u']),s.sympify(row['v'])
        numberfield=s.QQ.algebraic_field(s.I,uu,s.sqrt(3))
        shifted=((-4+s.I*vv)*s.eye(5)-effective.subs(u,uu)).applyfunc(s.simplify)
        coalescing=s.Matrix.hstack(*(shifted*shifted).applyfunc(s.simplify).nullspace()).applyfunc(s.simplify)
        assert coalescing.shape==(5,2),'the eigenvalue is not semisimple'
        joint=lambda A:DM.from_Matrix(s.Matrix.hstack(coalescing,(A*coalescing).applyfunc(s.simplify))).convert_to(numberfield).rank()
        entry=dict(u=row['u'],end_weight_term=joint((s.I*uu/(4*s.sqrt(3))*TT).applyfunc(s.simplify)),dephasing_term=joint(DD))
        assert entry['end_weight_term']==3 and entry['dephasing_term']==3
        # Intersection DIMENSIONS, not a joint rank: a two-plane meets a
        # one-dimensional eigenspace above its dimension for free, so only the
        # dimension of the meet is a measurement on both branches.
        meet={}
        for slope,space in eigenspaces.items():
            columns=s.Matrix.hstack(space,coalescing)
            meet[slope]=space.cols+2-DM.from_Matrix(columns).convert_to(numberfield).rank()
        assert meet=={'sqrt(3)':1,'-sqrt(3)':0}
        entry['multiplet_meet']=meet
        # The condition PROOF_CODIM1_BY_ADDITIVITY actually measures, which the
        # non-invariance above only implies: neither compression is scalar.
        inverse=(coalescing.T*coalescing).inv()*coalescing.T
        entry['compression_scalar']={name:s.simplify(s.simplify(inverse*A*coalescing)
                                                     -s.simplify(inverse*A*coalescing)[0,0]*s.eye(2)).is_zero_matrix
                                     for name,A in (('end_weight_term',s.I*uu/(4*s.sqrt(3))*TT),('dephasing_term',DD))}
        assert entry['compression_scalar']=={'end_weight_term':False,'dephasing_term':False}
        twin.append(entry)

    provenance=artifact_provenance(Path(__file__),Path(__file__))
    provenance.pop('source_sha256')  # This producer builds its object from code alone.
    out=dict(scope=__doc__,resonant_dimension=5,reflection='+I',D_projection=str(DD),end_hop_projection=str(TT),
             effective_determinant=str(polynomial),residual_cubic=str(s.expand(H)),quotient=str(s.cancel(polynomial/H)),
             single_particle_energies=[-3,-1,1,3],**provenance)
    out.update(effective_ranks=ranks)
    out.update(absorption_rate_per_stratum=rate,
               hamming_strata={str(k):hamming.count(k) for k in strata},
               resonance_rank_per_stratum={str(k):r for k,r in zip(strata,stratum_ranks)},
               resonance_stratum_intersections=meets,
               restricted_dissipator_rates=sorted(str(k) for k in restricted),
               restricted_dissipator_mean=int(mean),
               restricted_mean_per_multiplet={str(k):int(m) for k,m in sorted(means.items())},
               stratum_weight_per_multiplet={str(k):str(w) for k,w in sorted(weights.items())},
               coherent_half_symmetric_at_symbolic_weight=True,
               free_fermion_factor_multiplicities_over_Qw=[list(pair) for pair in additivity],
               end_weight_lifting_the_multiplet_to_five=quadruple_at,
               resonance_is_four_plus_one=True,
               coherent_half_frequency_comb={str(k):m for k,m in sorted(comb.items())},
               coherent_half_geometric_total=sum(comb.values()),
               scalar_dephasing_determinant=str(flattened),
               scalar_dephasing_line={'eigenvalue':'sqrt(3)','algebraic':algebraic,'geometric':geometric},
               end_weight_term_slopes={str(k):m for k,m in sorted(slopes.items(),key=str)},
               effective_term_plane_ranks=twin)
    (Path(__file__).parent/'results/route_b_n4_resonance_projection.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
