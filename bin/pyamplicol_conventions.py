"""Shared preparation conventions; usable by both Python 2 and Python 3.

Public leg labels are one-based and incoming-first. Only colour representations
are crossed here; physical momenta, spin sources and normalization are not.
"""
import itertools
import numbers


def _integer(value):
    return isinstance(value, numbers.Integral) and not isinstance(value, bool)


def outgoing_representations(record):
    """Return non-singlet charges, conjugating incoming 3/3bar exactly once."""
    colors = record.get('leg_colors')
    if not isinstance(colors, (list, tuple)):
        raise ValueError('Colour conventions require exported leg_colors')
    ninitial = record.get('ninitial')
    if not _integer(ninitial) or not 0 <= ninitial <= len(colors):
        raise ValueError('Invalid number of incoming colour legs')
    if 'pdgs' in record and len(record['pdgs']) != len(colors):
        raise ValueError('PDG and colour leg counts differ')
    active = {}
    for index, color in enumerate(colors):
        if not _integer(color) or color not in (1, 3, -3, 8):
            raise ValueError('Unsupported colour representation on leg %d: %s' %
                             (index + 1, color))
        if color != 1:
            active[index + 1] = -color if index < ninitial and abs(color) == 3 else color
    return active


def connection_signature(record, operations):
    """Validate a chronological history and return its final labelled charges.

    Reserve every Born label, including singlets, and every previously created
    label, including retired parents. A splitting represents one quark flavour.
    """
    active = outgoing_representations(record)
    seen = set(range(1, len(record['leg_colors']) + 1))
    if not isinstance(operations, (list, tuple)):
        raise ValueError('Colour connection must be a chronological list')
    for step in operations:
        if not isinstance(step, dict):
            raise ValueError('Colour operation must be a declaration')
        kind = step.get('kind')
        if kind == 'emit-gluon':
            fields = ('emitter', 'emitted')
            parent = step.get('emitter')
            daughters = (step.get('emitted'),)
        elif kind == 'split-gluon':
            fields = ('parent', 'quark', 'antiquark')
            parent = step.get('parent')
            daughters = (step.get('quark'), step.get('antiquark'))
        else:
            raise ValueError('Unknown colour operation: %s' % kind)
        if set(step) != set(('kind',) + fields):
            raise ValueError('Invalid fields for colour operation: %s' % kind)
        if any(not _integer(step[name]) or step[name] == 0 for name in fields):
            raise ValueError('Colour labels must be nonzero integers')
        if parent not in active:
            raise ValueError('Connection uses a colour leg before its creation or after retirement')
        if len(set(daughters)) != len(daughters) or any(label in seen for label in daughters):
            raise ValueError('New colour labels must be fresh, including retired and Born labels')
        if kind == 'emit-gluon':
            active[daughters[0]] = 8
        else:
            if active[parent] != 8:
                raise ValueError('Only an adjoint can split into a pair')
            del active[parent]
            active[daughters[0]], active[daughters[1]] = 3, -3
        seen.update(daughters)
    return active


def validate_correlator(record, request):
    """Validate physical bra-adjoint-times-ket, never an oriented legacy ID."""
    if not isinstance(request, dict) or 'bra' not in request or 'ket' not in request:
        raise ValueError('Colour correlator requires bra and ket histories')
    bra = connection_signature(record, request['bra'])
    ket = connection_signature(record, request['ket'])
    if len(request['bra']) != len(request['ket']):
        raise ValueError('Bra and ket must have the same colour order')
    if bra != ket:
        raise ValueError('Incompatible final colour representations')
    return bra


class UnorientableDipoleProduct(ValueError):
    """An ordered product cannot be represented by a single emission pair."""


def dipole_product(record, pairs, identifier=None):
    """Chronological histories for up to four full-colour dipoles in input order.

    [(a,b),(c,d)] means D_ab D_cd on the Born ket, without intermediate
    Born-space projection. Return a request with one emission per dipole on
    each side. Preserve an explicit identifier; otherwise derive a unique
    mnk-extra/dipole-product/ ID from the ordered endpoints.

    Orient non-self edges so each leg's operator word has a prefix of bra
    charges and a suffix of ket charges. At most sixteen orientations are tried:
    a graph with at most three non-self edges is bipartite or a triangle,
    and both admit this cut for every edge order. Four-edge words without a
    cut are rejected. Bra chronology follows the operator word; ket chronology
    reverses it.

    Self edges are scalar Casimirs. Appending their emission to both histories
    inserts sum_a T^a T^a between bra and ket, retaining the emission count.
    """
    pairs = list(pairs)
    if len(pairs) > 4:
        raise ValueError('Dipole product constructor supports at most four dipoles')
    active = outgoing_representations(record)
    edges, selves = [], []
    for index, pair in enumerate(pairs):
        if (not isinstance(pair, (list, tuple)) or len(pair) != 2 or
                any(not _integer(leg) or leg not in active for leg in pair)):
            raise ValueError('Dipole endpoints must be active Born colour labels')
        first, second = pair
        (selves if first == second else edges).append((index, first, second))
    oriented = None
    for flips in itertools.product((False, True), repeat=len(edges)):
        candidate, ket_started = [], set()
        for (index, first, second), flip in zip(edges, flips):
            bra, ket = (second, first) if flip else (first, second)
            if bra in ket_started:
                break
            candidate.append((index, bra, ket))
            ket_started.add(ket)
        else:
            oriented = candidate
            break
    if oriented is None:
        raise UnorientableDipoleProduct('Ordered dipoles have no chronological bra/ket cut')

    def emit(leg, index):
        return dict(kind='emit-gluon', emitter=leg, emitted=-index-1)

    bra = [emit(left, index) for index, left, right in oriented]
    ket = [emit(right, index) for index, left, right in reversed(oriented)]
    for index, leg, _ in selves:
        bra.append(emit(leg, index))
        ket.append(emit(leg, index))
    if identifier is None:
        identifier = 'mnk-extra/dipole-product/' + ('/'.join('%d_%d' % tuple(pair) for pair in pairs) or 'born')
    request = dict(id=identifier, bra=bra, ket=ket)
    validate_correlator(record, request)
    return request


def relabel_correlator(request, order):
    """Relabel external legs: order[new-1] is the old one-based leg label.

    This changes labels, not incoming/outgoing states. Auxiliary labels and
    chronological order are preserved. A physical crossing needs its own record
    and separately generated artifact, not just this permutation.
    """
    if any(not _integer(leg) for leg in order) or sorted(order) != list(range(1, len(order) + 1)):
        raise ValueError('External leg order must be a one-based permutation')
    labels = dict((old, new) for new, old in enumerate(order, 1))
    result = dict(request)
    for side in ('bra', 'ket'):
        result[side] = [dict((key, labels.get(value, value) if key != 'kind' else value)
                             for key, value in step.items()) for step in request[side]]
    return result
