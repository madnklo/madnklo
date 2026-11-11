"""Export-time support for the opt-in native pyAmpliCol Born provider."""
import json
import os
from fractions import Fraction


_DISPATCH_TEMPLATE = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                'template_files', 'pyamplicol_dispatch.inc')


def resolve_providers(provider=None, environ=None):
    """Resolve export options once, preserving the legacy tree-only switch."""
    environ = os.environ if environ is None else environ
    combined = environ.get('MADNKLO_ME_PROVIDER')
    legacy = environ.get('MADNKLO_BORN_PROVIDER')
    settings = (('me_provider', provider), ('MADNKLO_ME_PROVIDER', combined),
                ('MADNKLO_BORN_PROVIDER', legacy))
    for name, value in settings:
        if value is not None and value not in ('madgraph', 'pyamplicol'):
            raise ValueError('Unknown %s: %s' % (name, value))
    if len(set(value for name, value in settings if value is not None)) > 1:
        raise ValueError('Conflicting provider settings: '+', '.join(
            '%s=%s' % (name, value) for name, value in settings if value is not None))
    return dict(tree=provider or combined or legacy or 'madgraph',
                virtual=provider or combined or 'madgraph')


def selected_provider(kind='tree', environ=None, options=None):
    """Read a resolved export selection; retain standalone helper compatibility."""
    if kind not in ('tree', 'virtual'):
        raise ValueError('Unknown amplitude kind: '+kind)
    options = options or {}
    if 'me_providers' in options:
        return options['me_providers'][kind]
    return resolve_providers(options.get('me_provider'), environ)[kind]


def color_amplitudes(matrix_element):
    """Exact JAMP coefficients at Nc=3, before any kinematic evaluation."""
    rows = []
    for terms in matrix_element.get_color_amplitudes():
        row = []
        for (fermion, coefficient, imaginary, power), amplitude in terms:
            value = fermion * Fraction(coefficient) * Fraction(3)**power
            row.append(dict(amplitude=amplitude,
                            real=str(0 if imaginary else value),
                            imag=str(value if imaginary else 0)))
        rows.append(row)
    return rows


def provider_code(matrix_element, directory, prefix, n_squared_orders,
                  color_order=0, spin_order=0, options=None):
    result = dict((name, '') for name in (
        'pyamp_declarations', 'pyamp_evaluate', 'pyamp_color', 'pyamp_spin',
        'pyamp_table', 'pyamp_helpers', 'pyamp_hook_declarations', 'pyamp_hook_check',
        'pyamp_get_me_declarations', 'pyamp_get_me', 'pyamp_initialise',
        'pyamp_set_color_guard', 'pyamp_add_color_guard', 'pyamp_set_spin_guard'))
    if selected_provider(options=options) == 'madgraph':
        return result
    if n_squared_orders != 1:
        raise ValueError('pyAmpliCol currently requires one squared coupling order')
    if color_order not in (0, 1, 2) or spin_order not in (0, 1, 2):
        raise ValueError('pyAmpliCol supports correlation requirements through NNLO')
    result.update(pyamp_color_order=color_order, pyamp_spin_order=spin_order)
    process = matrix_element.get('processes')[0]
    for constraint in ('required_s_channels', 'forbidden_s_channels',
                       'forbidden_onsh_s_channels', 'forbidden_particles',
                       'is_decay_chain'):
        if process.get(constraint):
            raise ValueError('pyAmpliCol adapter does not translate %s' % constraint)
    legs = sorted(process.get('legs'), key=lambda leg: leg.get('number'))
    pdgs = [leg.get('id') for leg in legs]
    ninitial = process.get_ninitial()
    key = 'p%d_%s' % (ninitial, '_'.join(str(pdg) for pdg in pdgs))
    model = process.get('model')
    names = [model.get('particle_dict')[pdg].get_name() for pdg in pdgs]
    record = {'key': key, 'pdgs': pdgs, 'ninitial': ninitial,
              'process': '%s > %s' % (' '.join(names[:ninitial]),
                                      ' '.join(names[ninitial:])),
              'restriction_defaults': [],
              'required_correlations': dict(color_order=color_order, spin_order=spin_order),
              'orders': dict((name, value) for name, value in process.get('orders').items()
                             if name != 'WEIGHTED')}
    basis = matrix_element.get('color_basis')
    if basis is not None:
        record['color_basis'] = sorted(basis.keys())
        record['leg_colors'] = [model.get_particle(pdg).get_color() for pdg in pdgs]
        if basis:
            record['color_amplitudes'] = color_amplitudes(matrix_element)
    restriction = getattr(model,'restrict_card',None)
    if restriction:
        from models import check_param_card
        if isinstance(restriction,str):
            restriction = check_param_card.ParamCard(restriction)
        record['restriction_defaults'] = [
            [block,list(parameter.lhacode),parameter.value]
            for block,entries in sorted(restriction.items())
            if not block.startswith('decay_table')
            for parameter in entries
            if isinstance(parameter.value,(int,float))]
    manifest_dir = os.path.join(directory, 'BornProvider')
    if not os.path.isdir(manifest_dir):
        os.makedirs(manifest_dir)
    with open(os.path.join(manifest_dir, key + '.json'), 'w') as stream:
        json.dump(record, stream, indent=2, sort_keys=True)
    result['pyamp_declarations'] = '''
      INTEGER MNK_N,MNK_NINITIAL
      PARAMETER(MNK_N=%d,MNK_NINITIAL=%d)
      INTEGER MNK_PDGS(MNK_N)
      DATA MNK_PDGS/%s/
''' % (len(pdgs), ninitial, ','.join(str(pdg) for pdg in pdgs))
    result['pyamp_hook_declarations'] = '''
      INTEGER MNK_POLARIZATIONS(0:MNK_N,0:5)
      COMMON/%sBORN_BEAM_POL/MNK_POLARIZATIONS
''' % prefix
    result['pyamp_hook_check'] = '''
      IF(MNK_POLARIZATIONS(0,0).NE.-1) THEN
         WRITE(*,*) 'pyAmpliCol: selected beam polarizations are unsupported'
         STOP 1
      ENDIF
'''
    result['pyamp_get_me_declarations'] = '''
      REAL*8 MNK_ANS(0:1)
      INTEGER MNK_HEL
'''
    result['pyamp_get_me'] = '''
      MNK_HEL=NHEL
      IF(MNK_HEL.EQ.0) MNK_HEL=-1
      CALL %sME_ACCESSOR_HOOK(P,MNK_HEL,ALPHAS,MNK_ANS)
      ANS=MNK_ANS(0)
''' % prefix
    result['pyamp_evaluate'] = '''
      %sDEN=%d
      CALL MNK_PYAMP_EVAL(MNK_N,MNK_NINITIAL,MNK_PDGS,P,USER_ALPHAS,HEL,NSQAMPSO,ANS)
''' % (prefix, matrix_element.get_denominator_factor())
    result['pyamp_color'] = '''
      CALL MNK_PYAMP_COLOR(MNK_N,MNK_NINITIAL,MNK_PDGS,LEG_A,LEG_B,%sGET_CCBLO)
''' % prefix
    result['pyamp_spin'] = '''
      CALL MNK_PYAMP_SPIN(MNK_N,MNK_NINITIAL,MNK_PDGS,PB,PARENT_LEG,KT,%sGET_KKBLO)
''' % prefix
    return result


def correlation_code(replacements, denominator):
    """Redirect the selected-state standalone API without touching MadLoop."""
    if not replacements['pyamp_evaluate'] or not (
            replacements['pyamp_color_order'] or replacements['pyamp_spin_order']):
        return {}
    values = dict(replacements)
    values['pyamp_ncorr'] = max(1,values['n_color_correlators'])
    values['pyamp_nconnections'] = max(1,values['n_color_connections'])
    with open(_DISPATCH_TEMPLATE) as stream:
        helpers = stream.read() % values
    prefix = values['proc_prefix']
    return {'pyamp_evaluate': '''
      %sDEN=%s
      CALL %sPYAMP_EVALUATE(P,HEL,USER_ALPHAS,ANS)
''' % (prefix,denominator,prefix),
            'pyamp_table': '      CALL %sPYAMP_TABLE(CORRELATED_MES)\n'%prefix,
            'pyamp_color': '''
      CALL %sPYAMP_DIPOLE(LEG_A,LEG_B,%sGET_CCBLO)
''' % (prefix,prefix),
            'pyamp_spin': '''
      CALL %sPYAMP_SPIN(PARENT_LEG,PB,KT,%sGET_KKBLO)
''' % (prefix,prefix),
            'pyamp_initialise':'      CALL %sPYAMP_RESET()\n'%prefix,
            'pyamp_set_color_guard': '''
      IF(FIRST_CONNECTION.NE.0) THEN
         IF(FIRST_CONNECTION.LT.-1.OR.FIRST_CONNECTION.GT.NCOLORCONNECTIONS.OR.SECOND_CONNECTION.LT.-1.OR.SECOND_CONNECTION.EQ.0.OR.SECOND_CONNECTION.GT.NCOLORCONNECTIONS) THEN
            WRITE(*,*) 'pyAmpliCol: invalid colour connection selection'
            STOP 1
         ENDIF
         IF(FIRST_CONNECTION.GT.0.AND.SECOND_CONNECTION.GT.0) THEN
            IF(COLOR_CORRELATOR_TO_INDEX(FIRST_CONNECTION,SECOND_CONNECTION).LE.0) THEN
               WRITE(*,*) 'pyAmpliCol: unavailable colour connection pair'
               STOP 1
            ENDIF
         ENDIF
      ENDIF
''',
            'pyamp_add_color_guard': '''
      IF(FIRST_CONNECTION.LT.1.OR.FIRST_CONNECTION.GT.NCOLORCONNECTIONS.OR.SECOND_CONNECTION.LT.1.OR.SECOND_CONNECTION.GT.NCOLORCONNECTIONS) THEN
         WRITE(*,*) 'pyAmpliCol: invalid added colour connection'
         STOP 1
      ENDIF
      IF(COLOR_CORRELATOR_TO_INDEX(FIRST_CONNECTION,SECOND_CONNECTION).LE.0) THEN
         WRITE(*,*) 'pyAmpliCol: unavailable colour connection pair'
         STOP 1
      ENDIF
      IF(ANY(COLOR_CORRELATORS_TO_CONSIDER(1:COLOR_CORRELATORS_TO_CONSIDER(0)).EQ.COLOR_CORRELATOR_TO_INDEX(FIRST_CONNECTION,SECOND_CONNECTION))) RETURN
''',
            'pyamp_set_spin_guard': '''
      IF(LEG_INDEX.LT.1.OR.LEG_INDEX.GT.MNK_N) THEN
         WRITE(*,*) 'pyAmpliCol: invalid spin leg'
         STOP 1
      ENDIF
      IF(MNK_PDGS(LEG_INDEX).NE.21) THEN
         WRITE(*,*) 'pyAmpliCol: spin replacement requires a gluon'
         STOP 1
      ENDIF
      IF(N_VECTORS.LT.1.OR.N_VECTORS.GT.MAX_N_SPIN_CORR_VECTORS) THEN
         WRITE(*,*) 'pyAmpliCol: invalid spin vector count'
         STOP 1
      ENDIF
      IF(N_SPIN_CORR_VECTORS(LEG_INDEX).NE.0) THEN
         WRITE(*,*) 'pyAmpliCol: reset spin vectors before replacing a leg again'
         STOP 1
      ENDIF
      IF(COUNT(N_SPIN_CORR_VECTORS.GT.0).GE.MAX_LEGS_WITH_SPIN_CORR) THEN
         WRITE(*,*) 'pyAmpliCol: too many spin-correlated legs'
         STOP 1
      ENDIF
''', 'pyamp_helpers':helpers}
