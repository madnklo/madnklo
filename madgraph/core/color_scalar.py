"""Exact closed SU(N) contractions, with bounded memoization of Fierz steps.

This specializes scalar products, not open colour tensors. Unsupported or open
networks return None so callers can retain the general colour algebra.
"""
from collections import Counter
from fractions import Fraction
import itertools

from madgraph.core import color_algebra as ca


class ClosedColorTraceCache(object):
    def __init__(self, max_entries=32768):
        if max_entries < 1:
            raise ValueError('The trace cache must have positive capacity')
        self.max_entries = max_entries
        self.cache = {}

    @staticmethod
    def _canonical(words):
        # Only cyclic rotations, factor ordering and dummy renaming are used.
        # This need not identify every isomorphism to be a valid cache key.
        ordered = sorted((min(word[i:]+word[:i] for i in range(len(word)))
                          if word else () for word in words), key=lambda w: (len(w), w))
        labels = {}
        return tuple(tuple(labels.setdefault(label, len(labels)) for label in word)
                     for word in ordered)

    def _traces(self, words):
        """Integer Laurent coefficients with one factor 1/2 per adjoint pair omitted."""
        key = self._canonical(words)
        if key in self.cache:
            return self.cache[key]
        if any(len(word) == 1 for word in key):
            result = ()
        else:
            empty = sum(not word for word in key)
            nonempty = tuple(word for word in key if word)
            if not nonempty:
                result = ((empty, 1),)
            elif empty:
                result = tuple((power+empty, value) for power, value in self._traces(nonempty))
            else:
                # Contract the first generator with its unique partner using
                # sum_a T^a_ij T^a_kl = (delta_il delta_kj-delta_ij delta_kl/N)/2.
                first = nonempty[0]
                label = first[0]
                partner = next((i, word.index(label, 1 if i == 0 else 0))
                               for i, word in enumerate(nonempty)
                               if label in word[1 if i == 0 else 0:])
                i, j = partner
                if i == 0:
                    leading = (first[1:j], first[j+1:])+nonempty[1:]
                    subleading = (first[1:j]+first[j+1:],)+nonempty[1:]
                else:
                    second = nonempty[i]
                    rest = nonempty[1:i]+nonempty[i+1:]
                    leading = (first[1:]+second[j+1:]+second[:j],)+rest
                    subleading = (first[1:], second[j+1:]+second[:j])+rest
                terms = {}
                for power, value in self._traces(leading):
                    terms[power] = terms.get(power, 0)+value
                for power, value in self._traces(subleading):
                    terms[power-1] = terms.get(power-1, 0)-value
                result = tuple(sorted((power, value) for power, value in terms.items() if value))
        if len(self.cache) >= self.max_entries:
            self.cache.clear()
        self.cache[key] = result
        return result

    def contract(self, string):
        """Return a scalar ColorFactor, preserving the full Laurent polynomial in Nc."""
        traces, vertices, lines = [], [], {}
        heads, tails = [], []
        adjoint = []
        for obj in string:
            if type(obj) is ca.Tr:
                traces.append(tuple(obj))
                adjoint.extend(obj)
            elif type(obj) in (ca.f, ca.d):
                vertices.append((type(obj) is ca.f, tuple(obj)))
                adjoint.extend(obj)
            elif type(obj) is ca.T:
                head, tail = obj[-2:]
                heads.append(head)
                tails.append(tail)
                lines[head] = (tail, tuple(obj[:-2]))
                adjoint.extend(obj[:-2])
            elif type(obj) is not ca.ColorOne:
                return None
        # Each fundamental index must connect one head to one tail; adjoint
        # indices occur twice, in a separate namespace. Otherwise use the general path.
        if (len(set(heads)) != len(heads) or len(set(tails)) != len(tails) or
                set(heads) != set(tails) or set(heads).intersection(adjoint) or
                any(count != 2 for count in Counter(adjoint).values())):
            return None
        while lines:
            start = next(iter(lines))
            label, word = start, ()
            while True:
                label, segment = lines.pop(label)
                word += segment
                if label == start:
                    break
            traces.append(word)

        terms = {}
        imaginary_power = int(string.is_imaginary)+sum(is_f for is_f, _ in vertices)
        # Every Fierz branch removes exactly one adjoint pair. Keeping this
        # common denominator outside the recursion avoids repeated Fraction work.
        coefficient = string.coeff * Fraction(2**len(vertices), 2**(len(adjoint)//2))
        coefficient *= -1 if imaginary_power % 4 >= 2 else 1
        # f(a,b,c)=-2i[Tr(a,b,c)-Tr(c,b,a)], d(a,b,c)=2[Tr(a,b,c)+Tr(c,b,a)].
        for choices in itertools.product((False, True), repeat=len(vertices)):
            words = list(traces)
            sign = 1
            for (is_f, word), reverse in zip(vertices, choices):
                words.append(word[::-1] if reverse else word)
                if is_f and not reverse:
                    sign = -sign
            for power, value in self._traces(words):
                terms[power] = terms.get(power, 0)+sign*value
        return ca.ColorFactor([ca.ColorString([], coeff=coefficient*value,
            is_imaginary=bool(imaginary_power % 2), Nc_power=power+string.Nc_power,
            loop_Nc_power=string.loop_Nc_power)
            for power, value in sorted(terms.items()) if coefficient*value])
