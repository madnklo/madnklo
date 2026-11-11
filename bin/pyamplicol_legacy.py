#!/usr/bin/env python
"""Exact, preparation-only legacy colour compatibility (MadNkLO Python 2).

No legacy connection is identified with an adjointed connection by its name.
Both matrices are constructed in the exported Born basis and matched exactly.
When exported JAMP coefficients are available, equality is required on their
exact column space, for arbitrary complex diagram amplitudes (including spin
replacements), rather than on unreachable independent colour coefficients.
"""
import json
import itertools
import math
import os
import sys
import time
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from madgraph.core import base_objects, color_amp, color_algebra as ca
from madgraph.core.color_scalar import ClosedColorTraceCache
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from pyamplicol_conventions import (connection_signature, outgoing_representations,
                                   validate_correlator, dipole_product, relabel_correlator)


def immutable(value):
    return tuple((name, tuple(indices)) for name, indices in value)


class PhysicalColorCache(dict):
    def __init__(self):
        dict.__init__(self)
        self.closed_traces = ClosedColorTraceCache()
        self.trace_bases = {}
        self.trace_transitions = {}
        self.trace_words = {}
        self.trace_overlaps = {}

    def trace_permutations(self, basis, record):
        representations = outgoing_representations(record)
        key = (tuple(immutable(row) for row in basis), tuple(sorted(representations.items())),
               len(record['leg_colors']))
        if key not in self.trace_bases:
            self.trace_bases[key] = None
            n = len(record['leg_colors'])
            if (not 2 <= n <= 5 or sorted(representations) != list(range(1, n+1)) or
                    any(value != 8 for value in representations.values()) or
                    len(basis) != math.factorial(n-1)):
                return None
            words = []
            for row in key[0]:
                if (len(row) != 1 or row[0][0] != 'Tr' or len(row[0][1]) != n or
                        sorted(row[0][1]) != list(range(1, n+1))):
                    return None
                word = row[0][1]
                words.append(min(word[i:]+word[:i] for i in range(n)))
            lookup = dict((word, i) for i, word in enumerate(words))
            if len(lookup) != len(basis):
                return None
            permutations = []
            for order in itertools.permutations(range(1, n+1)):
                labels = dict((old, new) for new, old in enumerate(order, 1))
                indices = []
                for word in words:
                    image = tuple(labels[label] for label in word)
                    indices.append(lookup[min(image[i:]+image[:i] for i in range(n))])
                permutations.append((order, indices))
            self.trace_bases[key] = (permutations, {})
        return self.trace_bases[key]

    @staticmethod
    def normalise_traces(words):
        if any(len(word) == 1 for word in words):
            return (), F(0)
        factor = F(3**sum(not word for word in words))
        return tuple(sorted(min(word[i:]+word[:i] for i in range(len(word)))
                            for word in words if word)), factor

    def trace_dipole(self, state, pair):
        """Act on full external-index trace tensors, including generated multi-traces."""
        key = (state, tuple(sorted(pair)))
        if key in self.trace_transitions:
            return self.trace_transitions[key]
        if pair[0] == pair[1]:
            result = {state: F(3)}
        else:
            expanded = [(state, F(1))]
            for leg in pair:
                next_terms = []
                for words, coefficient in expanded:
                    i = next(i for i, word in enumerate(words) if leg in word)
                    word = words[i]
                    j = word.index(leg)
                    # -i f^{a b c} T^c = [T^b,T^a]. The temporary adjoint
                    # index 0 occurs twice after acting on the two endpoints.
                    for position, sign in ((j+1, 1), (j, -1)):
                        inserted = word[:position]+(0,)+word[position:]
                        next_terms.append((words[:i]+(inserted,)+words[i+1:], coefficient*sign))
                expanded = next_terms
            result = {}
            for words, coefficient in expanded:
                (i, p), (j, q) = [(i, p) for i, word in enumerate(words)
                                  for p, label in enumerate(word) if label == 0]
                if i == j:
                    first, second = words[i][p+1:q], words[i][q+1:]+words[i][:p]
                    rest = words[:i]+words[i+1:]
                    leading, subleading = (first, second)+rest, (first+second,)+rest
                else:
                    first = words[i][p+1:]+words[i][:p]
                    second = words[j][q+1:]+words[j][:q]
                    rest = words[:i]+words[i+1:j]+words[j+1:]
                    leading, subleading = (first+second,)+rest, (first, second)+rest
                for image, weight in ((leading, F(1,2)), (subleading, F(-1,6))):
                    normal, factor = self.normalise_traces(image)
                    if factor:
                        result[normal] = result.get(normal, F(0))+coefficient*weight*factor
            result = dict((state, value) for state, value in result.items() if value)
        self.trace_transitions[key] = result
        return result

    def trace_word(self, state, pairs):
        key = (state, pairs)
        if key not in self.trace_words:
            result = {}
            if not pairs:
                result[state] = F(1)
            else:
                # D_1 ... D_n acts on the ket from right to left. Intermediate
                # tensors are not projected onto the exported Born basis.
                for previous, coefficient in self.trace_word(state, pairs[1:]).items():
                    for current, value in self.trace_dipole(previous, pairs[0]).items():
                        result[current] = result.get(current, F(0))+coefficient*value
            self.trace_words[key] = dict((state, value) for state, value in result.items() if value)
        return self.trace_words[key]

    def dipole_matrix(self, basis, pairs):
        states = [self.normalise_traces((tuple(row[0][1]),))[0] for row in basis]
        columns = [self.trace_word(state, pairs) for state in states]
        result = []
        for bra in states:
            for column in columns:
                real, imag = F(0), F(0)
                for ket, coefficient in column.items():
                    key = (bra, ket)
                    if key not in self.trace_overlaps:
                        string = ca.ColorString([ca.Tr(*bra[0][::-1])]+[ca.Tr(*word) for word in ket])
                        self.trace_overlaps[key] = self.closed_traces.contract(string).set_Nc(3)
                    mr, mi = self.trace_overlaps[key]
                    real += coefficient*mr
                    imag += coefficient*mi
                result.extend((real, imag))
        return result


def connected(basis, record, operations, side):
    connection_signature(record, operations)
    cs = ca.ColorString()
    cs.from_immutable(basis)
    active = outgoing_representations(record)
    labels = set(active)
    for step in operations:
        labels.update(value for key, value in step.items() if key != 'kind')
    extra = sorted(labels-set(active))
    indices = dict((label, 1000+i) for i, label in enumerate(extra))
    indices.update((label,label) for label in active)
    dummy = -1000000*(side+1)
    # Born-basis internal indices must be disjoint on the two sides too.
    internal = sorted(set(index for obj in cs for index in obj if index < 0))
    cs.replace_indices(dict((index,dummy-i-1) for i,index in enumerate(internal)))
    dummy -= len(internal)+1
    for step in operations:
        parent = step.get('emitter',step.get('parent'))
        if parent not in active:
            raise ValueError('Connection uses a colour leg before its creation')
        representation = active[parent]
        outgoing = indices[parent]
        cs.replace_indices({outgoing:dummy})
        if step['kind'] == 'emit-gluon':
            emitted = step['emitted']
            if emitted in active:
                raise ValueError('Repeated emitted colour label')
            active[emitted] = 8
            if representation == 3:
                operator = ca.ColorString([ca.T(indices[emitted],outgoing,dummy)])
            elif representation == -3:
                operator = ca.ColorString([ca.T(indices[emitted],dummy,outgoing)],coeff=F(-1))
            elif representation == 8:
                operator = ca.ColorString([ca.f(indices[emitted],outgoing,dummy)],coeff=F(-1),is_imaginary=True)
            else:
                raise ValueError('Unsupported colour representation')
        elif step['kind'] == 'split-gluon':
            if representation != 8:
                raise ValueError('Only an adjoint can split into a pair')
            q,qb = step['quark'],step['antiquark']
            del active[parent]
            if q in active or qb in active or q == qb:
                raise ValueError('Repeated pair daughter label')
            active[q],active[qb] = 3,-3
            operator = ca.ColorString([ca.T(dummy,indices[q],indices[qb])])
        else:
            raise ValueError('Unknown colour operation')
        cs.product(operator)
        dummy -= 1
    return cs,active,indices


def physical_matrix(basis, record, request, cache):
    validate_correlator(record, request)
    if not basis:
        raise ValueError('Physical colour matrix needs an explicit Born basis')
    context = cache.trace_permutations(basis, record) if isinstance(cache, PhysicalColorCache) else None
    if context is None:
        return _physical_matrix(basis, record, request, cache)
    prefix = 'mnk-extra/dipole-product/'
    identifier = request.get('id')
    if isinstance(identifier, basestring) and identifier.startswith(prefix):
        try:
            suffix = identifier[len(prefix):]
            pairs = () if suffix == 'born' else tuple(tuple(int(leg) for leg in part.split('_'))
                                                      for part in suffix.split('/'))
            declared = dipole_product(record, pairs)
        except ValueError:
            declared = None
        if (declared is not None and declared['id'] == identifier and
                all(declared[side] == request[side] for side in ('bra', 'ket'))):
            return cache.dipole_matrix(basis, pairs)
    permutations, matrices = context
    best = None
    for order, indices in permutations:
        transformed = relabel_correlator(request, order)
        key = tuple(tuple(tuple(sorted(step.items())) for step in transformed[side])
                    for side in ('bra', 'ket'))
        if best is None or key < best[0]:
            best = (key, transformed, indices)
    key, transformed, indices = best
    if key not in matrices:
        matrices[key] = _physical_matrix(basis, record, transformed, cache)
    canonical = matrices[key]
    dimension = len(basis)
    # Relabel both the operator and basis, then restore the caller's row order.
    # This applies to physical matrices only, never to the legacy mirrored fill.
    return [canonical[2*(dimension*i+j)+component]
            for i in indices for j in indices for component in (0, 1)]


def _physical_matrix(basis, record, request, cache, ket_basis=None):
    left = [connected(key,record,request['bra'],0) for key in basis]
    right = [connected(key,record,request['ket'],1) for key in (basis if ket_basis is None else ket_basis)]
    if left[0][1] != right[0][1]:
        raise ValueError('Incompatible final colour representations')
    # Assign identical final indices even if the two histories have different retired labels.
    final = dict((label,10000+i) for i,label in enumerate(sorted(left[0][1])))
    result = []
    for lcs,lactive,lindices in left:
        for rcs,ractive,rindices in right:
            # The legacy base-class conjugation reverses f's indices as well.
            # Physical f and d tensors are real; only fundamental generators transpose.
            first = lcs.create_copy()
            first[:] = [obj if isinstance(obj,(ca.f,ca.d)) else obj.complex_conjugate() for obj in first]
            if first.is_imaginary:
                first.coeff = -first.coeff
            second = rcs.create_copy()
            first.replace_indices(dict((lindices[label],index) for label,index in final.items()))
            second.replace_indices(dict((rindices[label],index) for label,index in final.items()))
            first.product(second)
            # Every index is contracted here, so dummy-index renaming cannot change the value.
            key = (first.to_canonical()[0],first.coeff,first.is_imaginary,first.Nc_power)
            if key not in cache:
                scalar_cache = getattr(cache, 'closed_traces', None)
                contracted = scalar_cache.contract(first) if scalar_cache is not None else None
                if contracted is None:
                    contracted = ca.ColorFactor([first]).full_simplify()
                cache[key] = contracted.set_Nc(3)
            real,imag = cache[key]
            result.extend((F(real),F(imag)))
    return result


def add_expression(target, source, factor):
    for key,(real,imag) in source.items():
        previous = target.get(key,(F(0),F(0)))
        value = (previous[0]+factor*real,previous[1]+factor*imag)
        if value == (0,0):
            target.pop(key,None)
        else:
            target[key] = value


class Span(object):
    def __init__(self):
        self.rows = {}

    def add(self, vector, identifier, imaginary=False):
        vector = list(vector)
        expression = {identifier:(F(0),F(1)) if imaginary else (F(1),F(0))}
        for pivot,(row,terms) in sorted(self.rows.items()):
            factor = vector[pivot]
            if factor:
                for i, value in enumerate(row):
                    if value:
                        vector[i] -= factor*value
                add_expression(expression,terms,-factor)
        pivot = next((i for i,value in enumerate(vector) if value),None)
        if pivot is not None:
            factor = vector[pivot]
            self.rows[pivot] = ([value/factor for value in vector],
                dict((key,(value[0]/factor,value[1]/factor)) for key,value in expression.items()))
            return True
        return False

    def contains(self, vector):
        """Exact membership without constructing a usually unused alias."""
        vector = list(vector)
        for pivot, (row, _) in sorted(self.rows.items()):
            factor = vector[pivot]
            if factor:
                for i, value in enumerate(row):
                    if value:
                        vector[i] -= factor*value
        return not any(vector)

    def solve(self, vector):
        vector,expression = list(vector),{}
        for pivot,(row,terms) in sorted(self.rows.items()):
            factor = vector[pivot]
            if factor:
                for i, value in enumerate(row):
                    if value:
                        vector[i] -= factor*value
                add_expression(expression,terms,factor)
        return None if any(vector) else expression

    def separating_functional(self, vector):
        """Certify nonmembership over rational real/imaginary coordinates.

        Return w with w.dot(each spanning vector) == 0 and w.dot(vector) != 0.
        Elimination includes both M and i*M, so the certificate also excludes
        linear combinations with arbitrary complex coefficients.
        """
        residual = list(vector)
        rows = sorted(self.rows.items())
        for pivot, (row, _) in rows:
            factor = residual[pivot]
            if factor:
                residual = [a-factor*b for a, b in zip(residual, row)]
        candidates = []
        for coordinate, value in enumerate(residual):
            if not value:
                continue
            witness = [F(0)]*len(vector)
            witness[coordinate] = F(1)
            for pivot, (row, _) in reversed(rows):
                witness[pivot] -= sum(a*b for a, b in zip(witness, row))
            scale = next(v for v in witness if v)
            witness = [v/scale for v in witness]
            candidates.append(witness)
        if not candidates:
            return None
        return min(candidates, key=lambda w: (sum(v != 0 for v in w), w))


class LegacyMatchingError(ValueError):
    """An exact obstruction, retaining the original legacy matrix semantics."""

    def __init__(self, missing, targets, span):
        identifier = sorted(missing)[0]
        vector = targets[identifier]
        witness = span.separating_functional(vector)
        value = sum(a*b for a, b in zip(witness, vector))
        self.certificate = {
            'legacy_id': identifier,
            'missing_count': len(missing),
            'real_span_rank': len(span.rows),
            'real_matrix_dimension': len(vector),
            'coordinates': 'row-major real/imaginary entries after amplitude projection',
            'functional': [str(v) for v in witness],
            'target_value': str(value),
        }
        ValueError.__init__(self,
            'Physical catalogue cannot exactly represent legacy IDs: %s; '
            '%s is outside the exact complex span (real rank %d/%d). '
            'A rational separating functional annihilates every physical request '
            'but gives %s on this legacy target. Legacy values were not changed.' %
            (sorted(missing)[:8], identifier, len(span.rows), len(vector), value))


def amplitude_columns(record, dimension):
    """An exact complex basis for the exported JAMP coefficient column space."""
    rows = record.get('color_amplitudes')
    if rows is None:
        return None
    if len(rows) != dimension:
        raise ValueError('JAMP coefficients do not match the exported colour basis')
    columns = {}
    for i, terms in enumerate(rows):
        for term in terms:
            number = term['amplitude']
            if isinstance(number, bool) or not isinstance(number, int) or number < 1:
                raise ValueError('Invalid diagram amplitude number')
            column = columns.setdefault(number, [F(0)]*(2*dimension))
            column[2*i] += F(term['real'])
            column[2*i+1] += F(term['imag'])
    span, result = Span(), []
    for number, column in sorted(columns.items()):
        if span.solve(column) is not None:
            continue
        result.append(column)
        span.add(column, str(number))
        rotated = [value for real, imag in zip(column[::2], column[1::2])
                   for value in (-imag, real)]
        span.add(rotated, str(number), True)
    if not result:
        raise ValueError('Exported JAMP coefficient space is empty')
    return result


def project_matrix(vector, columns):
    """Compute C^dagger M C over complex rationals, with columns of C supplied."""
    if columns is None:
        return vector
    dimension = len(columns[0])//2
    if len(vector) != 2*dimension**2 or any(len(c) != 2*dimension for c in columns):
        raise ValueError('Invalid colour projection dimensions')
    nonzero = [[(i, column[2*i], column[2*i+1]) for i in range(dimension)
                if column[2*i] or column[2*i+1]] for column in columns]
    result = []
    for left in nonzero:
        for right in nonzero:
            real, imag = F(0), F(0)
            for i, ar, ai in left:
                for j, br, bi in right:
                    mr, mi = vector[2*(dimension*i+j):2*(dimension*i+j)+2]
                    cr, ci = ar*br+ai*bi, ar*bi-ai*br
                    real += cr*mr-ci*mi
                    imag += cr*mi+ci*mr
            result.extend((real, imag))
    return result


def extension_requests(record):
    """Stream short non-self dipole products without intermediate projection.

    Shared legs use the chronological cut proved by dipole_product.
    The search is deliberately incomplete and leaves an exact obstruction if
    it cannot span the original targets within the preparation budget.
    """
    legs = sorted(outgoing_representations(record))
    pairs = list(itertools.combinations(legs, 2))
    # Squared two-dipole words supply symmetric directions absent from the
    # NNLO physical catalogue and from all ordered triples for q g -> q g g.
    for first, second in itertools.product(pairs, repeat=2):
        yield dipole_product(record, (first, second, second, first))
    for choices in itertools.product(pairs, repeat=3):
        yield dipole_product(record, choices)


class Particle(object):
    def __init__(self,color):
        self.color = color

    def get_color(self):
        return self.color

    def get_name(self):
        return 'colour%d'%self.color


class Model(object):
    def __init__(self,record):
        self.colors = dict(zip(record['pdgs'],record['leg_colors']))

    def get_particle(self,pdg):
        return Particle(self.colors[pdg])


def legacy_matrices(record, basis, order):
    """Retain the literal matrices used by the original Fortran ID table."""
    cb = color_amp.ColorBasis()
    for key in basis:
        cb[key] = []
    matrix = color_amp.ColorMatrix(cb)
    legs = base_objects.LegList([base_objects.Leg({'id':pdg,'number':i+1,
        'state':i >= record['ninitial']}) for i,pdg in enumerate(record['pdgs'])])
    legacy,connections = matrix.build_color_correlated_matrices(legs,Model(record),
        order='N'*order+'LO',set_user_qqbar=True)
    offsets = {}
    offset = 1
    for level in range(1,order+1):
        for i,connection in enumerate(connections['N'*level+'LO']):
            offsets[offset+i] = connection['tuple_representation']
        offset += len(connections['N'*level+'LO'])
    targets = {}
    for pair,(_,current,_) in legacy.items():
        values = []
        for i in range(len(basis)):
            for j in range(len(basis)):
                # Original Fortran contracts A_i conj(A_j); physical requests use conj(A_i) A_j.
                values.extend(current.col_matrix_fixed_Nc[(j,i)] if current.col_matrix_fixed_Nc is not None else (F(0),F(0)))
        targets['legacy/%d/%d'%pair] = values
    return targets, offsets


def add_complex(span, vector, identifier):
    if span.contains(vector):
        return False
    changed = span.add(vector, identifier)
    rotated = [value for real, imag in zip(vector[::2], vector[1::2])
               for value in (-imag, real)]
    return span.add(rotated, identifier, True) or changed


def match_matrices(record, catalogue, targets, extension_catalogue=None, extension_limit=1024,
                   prefer_correlated=False):
    """Prove aliases over C^dagger M C, adding only necessary physical requests.

    None enables the bounded built-in search; an explicit empty iterable
    disables it. Supplied candidates are consumed lazily in caller order.
    """
    if isinstance(extension_limit, bool) or not isinstance(extension_limit, int) or extension_limit < 0:
        raise ValueError('Extension limit must be a nonnegative integer')
    basis = sorted(immutable(value) for value in record['color_basis'])
    columns = amplitude_columns(record, len(basis))
    projection_cache = {}

    def project(vector):
        key = tuple(vector)
        if key not in projection_cache:
            projection_cache[key] = project_matrix(vector, columns)
        return projection_cache[key]

    projected_targets = dict((key, project(value)) for key, value in targets.items())
    ncolored = len(outgoing_representations(record))

    def target_order(identifier):
        if not prefer_correlated:
            return 0
        parts = identifier.split('/')
        if len(parts) != 3 or parts[0] != 'legacy':
            return 0
        return 1 if max(map(int, parts[1:])) <= ncolored else 2

    order_spans = dict((level, Span()) for level in (1, 2))
    missing_by_order = dict((level, set(tuple(vector) for key, vector in projected_targets.items()
                                      if target_order(key) == level)) for level in (1, 2))
    span, cache = Span(), PhysicalColorCache()
    last_progress = time.time()
    # Different oriented IDs often export the same literal matrix.
    missing_vectors = set(tuple(vector) for vector in projected_targets.values())
    existing_ids = set()
    for request in catalogue:
        validate_correlator(record, request)
        identifier = request['id']
        if identifier in existing_ids:
            raise ValueError('Duplicate colour correlation ID: '+identifier)
        existing_ids.add(identifier)
    for completed, request in enumerate(sorted(catalogue,key=lambda entry:(len(entry['bra']),entry['id'])), 1):
        if prefer_correlated and not request['bra']:
            continue
        vector = project(physical_matrix(basis,record,request,cache))
        changed = add_complex(span, vector, request['id'])
        if changed:
            missing_vectors = set(vector for vector in missing_vectors if not span.contains(vector))
        level = len(request['bra'])
        if prefer_correlated and level in order_spans:
            if add_complex(order_spans[level], vector, request['id']):
                missing_by_order[level] = set(value for value in missing_by_order[level]
                                              if not order_spans[level].contains(value))
        if os.environ.get('MADNKLO_LEGACY_PROGRESS') == '1' and time.time()-last_progress >= 20:
            sys.stderr.write('Exact physical requests: %d/%d; complex rank: %d; unmatched matrices: %d\n' %
                             (completed, len(catalogue), len(span.rows)//2, len(missing_vectors)))
            sys.stderr.flush()
            last_progress = time.time()
        if not missing_vectors and not any(missing_by_order.values()):
            break
    base_rank = len(span.rows)//2
    base_rows = dict(span.rows)
    required = Span()
    required.rows = dict(span.rows)
    for index, vector in enumerate(sorted(missing_vectors)):
        add_complex(required, vector, 'target/%d' % index)
    minimum_additions = len(required.rows)//2-base_rank
    additional_requests = []
    additional_vectors = {}
    checked = 0
    if missing_vectors:
        if extension_catalogue is None:
            extension_catalogue = extension_requests(record)
        candidate_ids = set()
        for request in itertools.islice(extension_catalogue, extension_limit):
            checked += 1
            identifier = request['id']
            if not identifier.startswith('mnk-extra/'):
                raise ValueError('Extension request IDs must start with mnk-extra/')
            if identifier in existing_ids or identifier in candidate_ids:
                raise ValueError('Duplicate colour correlation ID: '+identifier)
            candidate_ids.add(identifier)
            validate_correlator(record,request)
            vector = project(physical_matrix(basis,record,request,cache))
            before = len(span.rows)
            add_complex(span, vector, identifier)
            if len(span.rows) != before:
                additional_requests.append(request)
                additional_vectors[identifier] = vector
                missing_vectors = set(value for value in missing_vectors if not span.contains(value))
                if not missing_vectors:
                    break
            if os.environ.get('MADNKLO_LEGACY_PROGRESS') == '1' and time.time()-last_progress >= 20:
                sys.stderr.write('Exact extension candidates: %d/%d; complex rank: %d; unmatched matrices: %d\n' %
                                 (checked, extension_limit, len(span.rows)//2, len(missing_vectors)))
                sys.stderr.flush()
                last_progress = time.time()
    if missing_vectors:
        missing = set(key for key, vector in projected_targets.items() if tuple(vector) in missing_vectors)
        error = LegacyMatchingError(missing, projected_targets, span)
        error.certificate['extension_candidates_checked'] = checked
        error.certificate['extension_limit'] = extension_limit
        augmented = Span()
        augmented.rows = dict(span.rows)
        for index, vector in enumerate(sorted(missing_vectors)):
            add_complex(augmented, vector, 'target/%d' % index)
        error.certificate['augmented_target_complex_rank'] = len(augmented.rows)//2
        error.certificate['base_target_complex_rank'] = base_rank+minimum_additions
        raise error
    aliases = {}
    expressions = {}
    for identifier,vector in sorted(targets.items()):
        level = target_order(identifier)
        projected = tuple(projected_targets[identifier])
        key = (level, projected)
        if key not in expressions:
            preferred = order_spans.get(level)
            expressions[key] = (preferred.solve(projected) if preferred is not None and
                                preferred.contains(projected) else span.solve(projected))
        expression = expressions[key]
        aliases[identifier] = {'terms':[{'id':key,'real':str(value[0]),'imag':str(value[1])}
                                       for key,value in sorted(expression.items())],
                               'matrix':[str(value) for value in vector]}
    used = set(term['id'] for alias in aliases.values() for term in alias['terms'])
    additional_requests = [request for request in additional_requests if request['id'] in used]
    span = Span()
    span.rows = base_rows
    for request in additional_requests:
        add_complex(span, additional_vectors[request['id']], request['id'])
    return {'process':record['key'],'aliases':aliases,'rank':len(span.rows),
            'additional_requests':additional_requests,
            'proof':{'complete':True,'legacy_target_semantics':
                     'current MadNkLO ColorMatrix including its symmetric fill',
                     'base_catalogue_count':len(catalogue),
                     'additional_request_count':len(additional_requests),
                     'additional_request_lower_bound':minimum_additions,
                     'minimum_additional_requests_proved':len(additional_requests) == minimum_additions,
                     'extension_candidates_checked':checked,
                     'extension_limit':extension_limit,
                     'complex_rank_before_extensions':base_rank,
                     'complex_rank_after_extensions':len(span.rows)//2,
                     'prefer_correlated_operators':prefer_correlated,
                     'matched_legacy_ids':len(aliases)},
            'matching_space': {'kind': 'full_color_basis' if columns is None else 'exported_amplitude_span',
                               'columns': None if columns is None else [[str(v) for v in c] for c in columns]},
            'reference':'current MadNkLO ColorMatrix with set_user_qqbar=True; no legacy bug corrections'}


def build(record,catalogue,order,extension_catalogue=None,extension_limit=1024):
    if 'color_basis' not in record or 'leg_colors' not in record:
        raise ValueError('Legacy compatibility needs a fresh provider-enabled export')
    outgoing_representations(record)
    for request in catalogue:
        validate_correlator(record, request)
    basis = sorted(immutable(value) for value in record['color_basis'])
    if not basis:
        return {'process':record['key'],'aliases':{},'connections':{},'rank':0,
                'additional_requests':[]}
    targets, offsets = legacy_matrices(record, basis, order)
    result = match_matrices(record, catalogue, targets, extension_catalogue, extension_limit,
                            prefer_correlated=True)
    result['connections'] = dict((str(key), value) for key,value in offsets.items())
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input')
    parser.add_argument('output')
    parser.add_argument('--physical-reference', help='Also save exact physical matrices for numerical audits')
    args = parser.parse_args()
    payload = json.load(open(args.input))
    result = build(payload['record'],payload['catalogue'],payload['order'],
                   payload.get('extension_catalogue'), payload.get('extension_limit',1024))
    with open(args.output,'w') as stream:
        json.dump(result,stream,indent=2,sort_keys=True)
    if args.physical_reference:
        basis = sorted(immutable(value) for value in payload['record']['color_basis'])
        cache = {}
        matrices = dict((entry['id'],[str(value) for value in physical_matrix(
            basis,payload['record'],entry,cache)]) for entry in
            payload['catalogue']+result['additional_requests'])
        with open(args.physical_reference,'w') as stream:
            json.dump(matrices,stream,indent=2,sort_keys=True)
