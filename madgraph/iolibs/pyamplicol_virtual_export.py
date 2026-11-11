"""Python 2 export support for ordinary native pyAmpliCol virtuals."""
import json
import os
from fractions import Fraction

from madgraph.iolibs.pyamplicol import selected_provider


def color_data(matrix_element, process, order):
    """Retain the actual rectangular legacy matrices, not a tree-only proof."""
    from madgraph.core.color_amp import ColorMatrix
    born, loop = matrix_element.get('born_color_basis'), matrix_element.get('loop_color_basis')
    result = dict(born_color_basis=sorted(born), loop_color_basis=sorted(loop), legacy_colors=[])
    if not order:
        return result
    if order != 1:
        raise ValueError('NNLO subtraction requires at most one colour insertion on a virtual')
    args = (process.get('legs'),process.get('model'))
    bb, _ = ColorMatrix(born).build_color_correlated_matrices(*args, order='NLO',set_user_qqbar=True)
    lb, _ = ColorMatrix(loop,born).build_color_correlated_matrices(*args, order='NLO',set_user_qqbar=True)
    for pair in sorted(lb):
        def transposed(matrix, nbra, nket):
            return [str(value) for i in range(nbra) for j in range(nket)
                    for value in ((0,0) if matrix.col_matrix_fixed_Nc is None else
                                  matrix.col_matrix_fixed_Nc[(j,i)])]
        bl = transposed(lb[pair][1],len(born),len(loop))
        # The old virtual contracts 2 Re(B^dagger K_BL L). Keep its reverse
        # direction separately rather than assuming that tree matrices suffice.
        reverse = [bl[2*(i*len(loop)+j)+component] if not component else
                   str(-Fraction(bl[2*(i*len(loop)+j)+1]))
                   for j in range(len(loop)) for i in range(len(born)) for component in (0,1)]
        result['legacy_colors'].append(dict(first=pair[0],second=pair[1],
            legs=[int(x[0] if isinstance(x,(tuple,list)) else x) for x in lb[pair][2]],
            bb=transposed(bb[pair][1],len(born),len(born)), bl=bl, lb=reverse))
    return result


def virtual_code(matrix_element, directory, replacements, color_order=0, spin_order=0,
                 options=None):
    empty = dict(pyamp_virtual_helpers='', pyamp_virtual_threshold='',
                 pyamp_virtual_hel_threshold='')
    if selected_provider('virtual', options=options) != 'pyamplicol':
        return empty
    process = matrix_element.get('processes')[0]
    if process.get('n_loops') != 1:
        raise ValueError('pyAmpliCol supports one-loop virtuals, not two-loop VV; '
                         'regenerate with --ignore_contributions=VV')
    if not process.get('has_born'):
        raise ValueError('pyAmpliCol virtuals require Born-one-loop interference')
    if process.get('perturbation_couplings') != ['QCD'] or replacements['nSquaredSO'] != 1:
        raise ValueError('pyAmpliCol virtuals require QCD and one squared coupling order')
    legs = sorted(process.get('legs'), key=lambda leg: leg.get('number'))
    model = process.get('model')
    pdgs = [leg.get('id') for leg in legs]
    ninitial = process.get_ninitial()
    names = [model.get_particle(pdg).get_name() for pdg in pdgs]
    parameters = model.get('parameters')[('external',)]
    if any(parameter.type != 'real' for parameter in parameters):
        raise ValueError('pyAmpliCol virtual routing currently requires real external parameters')
    key = 'p%d_%s' % (ninitial, '_'.join(str(pdg) for pdg in pdgs))
    record = dict(key=key, pdgs=pdgs, ninitial=ninitial,
        process='%s > %s' % (' '.join(names[:ninitial]), ' '.join(names[ninitial:])),
        orders=dict((name, value) for name, value in process.get('orders').items() if name != 'WEIGHTED'),
        loop_coupling_orders=['QCD'],
        external_parameters=[p.name[4:] if p.name.startswith('mdl_') else p.name for p in parameters],
        required_correlations=dict(color_order=color_order, spin_order=spin_order))
    record['leg_colors'] = [model.get_particle(pdg).get_color() for pdg in pdgs]
    squared, amplitude_orders = matrix_element.get_split_orders_mapping()
    if len(amplitude_orders['born_amp_orders']) != 1 or len(squared) != 1:
        raise ValueError('pyAmpliCol virtuals require one Born and one virtual squared order')
    split_names = process.get('split_orders')
    record['born_orders'] = dict(zip(split_names,amplitude_orders['born_amp_orders'][0][0]))
    record['virtual_squared_order'] = dict(zip(split_names,squared[0][0]))
    record.update(color_data(matrix_element,process,color_order))
    if color_order and len(record['legacy_colors']) != replacements['n_color_correlators']:
        raise ValueError('Virtual colour manifest does not match the exported getter layout')
    root = os.path.join(directory, 'VirtualProvider')
    if not os.path.isdir(root):
        os.makedirs(root)
    with open(os.path.join(root, key+'.json'), 'w') as stream:
        json.dump(record, stream, indent=2, sort_keys=True)
    prefix = replacements['proc_prefix']
    values = dict(replacements)
    values.setdefault('n_color_correlators', 0)
    values.update(pyamp_n=len(pdgs), pyamp_ni=ninitial,
        pyamp_nc=len(record['legacy_colors']),
        pyamp_pdgs=','.join(str(pdg) for pdg in pdgs), pyamp_np=len(parameters),
        pyamp_parameters='\n'.join('      MNK_PARAMS(%d)=%s' % (i+1, p.name)
                                   for i, p in enumerate(parameters)))
    template = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            'template_files', 'pyamplicol_virtual_dispatch.inc')
    with open(template) as stream:
        helpers = stream.read() % values
    return dict(pyamp_virtual_helpers=helpers,
        pyamp_virtual_threshold='''
      CALL %sPYAMP_VIRTUAL(P,-1,ANS,PREC_ASKED,PREC_FOUND,RET_CODE)
''' % prefix,
        pyamp_virtual_hel_threshold='''
      CALL %sPYAMP_VIRTUAL(P,HEL,ANS,PREC_ASKED,PREC_FOUND,RET_CODE)
''' % prefix)
