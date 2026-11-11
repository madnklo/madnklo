#!/usr/bin/env python
"""Exact Born/Born, Born/loop and loop/Born legacy-ID compatibility proof."""
from fractions import Fraction as F
import json
import os
import sys

sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from pyamplicol_legacy import PhysicalColorCache, Span, _physical_matrix, add_complex, immutable
from pyamplicol_conventions import outgoing_representations


def prove(record):
    born = [immutable(value) for value in record['born_color_basis']]
    loop = [immutable(value) for value in record['loop_color_basis']]
    if not born or not loop:
        raise ValueError('Virtual colour compatibility requires explicit Born and loop bases')
    colored = sorted(outgoing_representations(record))
    catalogue = [dict(id='born',bra=[],ket=[])]
    if record['required_correlations']['color_order']:
        catalogue += [dict(id='mnk-virtual/dipole/%d/%d'%(a,b),
            bra=[dict(kind='emit-gluon',emitter=a,emitted=-1)],
            ket=[dict(kind='emit-gluon',emitter=b,emitted=-1)]) for a in colored for b in colored]
    cache, span, matrices = PhysicalColorCache(), Span(), {}
    for request in catalogue:
        if not request['bra']:
            continue
        vector = (_physical_matrix(born,record,request,cache)+
                  _physical_matrix(born,record,request,cache,loop)+
                  _physical_matrix(loop,record,request,cache,born))
        matrices[request['id']] = vector
        add_complex(span,vector,request['id'])
    aliases = []
    for entry in record['legacy_colors']:
        target = [F(value) for block in ('bb','bl','lb') for value in entry[block]]
        direct = 'mnk-virtual/dipole/%d/%d' % tuple(entry['legs'])
        expression = {direct:(F(1),F(0))} if matrices.get(direct) == target else span.solve(target)
        if expression is None:
            raise ValueError('No Born-loop colour-space proof for legacy ID %s/%s' %
                             (entry['first'],entry['second']))
        aliases.append(dict(first=entry['first'],second=entry['second'],legs=entry['legs'],
            terms=[dict(id=key,real=str(value[0]),imag=str(value[1])) for key,value in sorted(expression.items())]))
    return dict(catalogue=catalogue,aliases=aliases,proof=dict(complete=True,
        spaces=['Born/Born','Born/loop','loop/Born'],born_dimension=len(born),
        loop_dimension=len(loop),matched_legacy_ids=len(aliases)))


if __name__ == '__main__':
    with open(sys.argv[1]) as source:
        result = prove(json.load(source))
    with open(sys.argv[2],'w') as output:
        json.dump(result,output,indent=2,sort_keys=True)
