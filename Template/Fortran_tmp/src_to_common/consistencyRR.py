#!/usr/bin/env python3
"""Consistency check of the counterterms of the double real (testRR).
"""
import argparse
import os
import re
import sys
from collections import OrderedDict

def sector_type(tag):
    """ijjk, ijkj or ijkl from the tag isec_jsec_c3p_d3p."""
    isec, jsec, c3p, d3p = [int(x) for x in tag.split("_")]
    if jsec == c3p and jsec != d3p:
        return "ijjk"
    if jsec != c3p and jsec == d3p:
        return "ijkj"
    return "ijkl"

def operators(label):
    """K<op>_<op>_...: the operators of a piece, e.g. KSS_SC -> {SS, SC}."""
    return frozenset(label[1:].split("_"))

def by_set(labels):
    d = {}
    for lab in labels:
        d.setdefault(operators(lab), []).append(lab)
    return d

def set_rule(target_of):
    def relations(labels):
        sets = by_set(labels)
        pairs = []
        for lab in labels:
            t = operators(lab)
            target = frozenset(target_of(t))
            if target == t:
                continue
            for partner in sets.get(target, [None]):
                pairs.append((lab, partner))
        return pairs
    return relations

def pair_rule(table):
    """Explicit relations: table = [(A, B), ...], those with both pieces
    in the log."""
    def relations(labels):
        have = set(labels)
        return [(a, b) for a, b in table if a in have and b in have]
    return relations


# ijjk (Sij, Cijk, SCijk), sector 3_4_4_5 of e+e- > g g c cx
IJJK = {
    # K1 soft limit of i
    "Si": set_rule(lambda t: t | {"S"}),
    # soft limit of the pair
    "Sij": set_rule(lambda t: t | {"SS"}),
    # collinear limit of ij: the SC piece becomes CC
    "Cij": set_rule(lambda t: t | ({"CC"} if "SC" in t else {"C"})),
    # soft-collinear limit of ijk: S where there is already a C
    "SCijk": set_rule(lambda t: t | ({"S"} if "C" in t else {"SC"})),
    # triple collinear limit
    "Cijk": set_rule(lambda t: t | {"CC"}),
    # spurious limits: pieces that coincide
    "Cir": pair_rule([
        ("KC", "KS_C"), ("KCC", "KS_CC"), ("KC_CC", "KS_C_CC"),
        ("KC_SS", "KS_C_SS"), ("KC_SS_CC", "KS_C_SS_CC"),
        ("KSC", "KS_SC"), ("KSC_CC", "KS_SC_CC"),
        ("KSS_CC", "KS_SS_CC"), ("KSS_SC", "KS_SS_SC"),
        ("KSS_SC_CC", "KS_SS_CC_SC")]),
    "Cjr": pair_rule([
        ("KC", "KS_C"), ("KC_CC", "KS_C_CC"), ("KC_SS", "KS_C_SS"),
        ("KC_SS_CC", "KS_C_SS_CC"), ("KSC", "KSS_SC"),
        ("KSC_CC", "KSS_SC_CC"), ("KS_CC", "KS_SS_CC")]),
    "Cijr": pair_rule([
        ("KCC", "KSS_CC"), ("KC_CC", "KC_SS_CC"), ("KSC", "KSS_SC"),
        ("KSC_CC", "KSS_SC_CC"), ("KS_CC", "KS_SS_CC"),
        ("KS_C_CC", "KS_C_SS_CC"), ("KS_SC", "KS_SS_SC"),
        ("KS_SC_CC", "KS_SS_CC_SC")]),
}

RELATIONS = {"ijjk": IJJK}



def read_testrr(path):
    """{event number: {limit: {lambda: {label: value}}}} from
    testRR_<sector>.log: for each lambda, the pieces logged by ct_log
    followed by the line 'lambda RR LIM |RR-LIM|/|LIM|'."""
    events = OrderedDict()
    event = limit = None
    cur = OrderedDict()
    with open(path) as f:
        for line in f:
            t = line.split()
            if not t or t[0].startswith("---") or t[0] == "lambda":
                continue
            if t[0] == "EVENT":
                event = events.setdefault(int(t[2]), OrderedDict())
                limit = None
                continue
            if t[0] == "LIM":
                limit, cur = t[2], OrderedDict()
                if event is not None:
                    event[limit] = OrderedDict()
                continue
            if limit is None or event is None:
                continue
            if t[0].startswith("K") and len(t) == 2:
                cur[t[0]] = float(t[1].replace("D", "E"))
                continue
            try:
                lam = float(t[0])
            except ValueError:
                continue
            if len(t) >= 4 and cur:
                event[limit][lam] = cur
                cur = OrderedDict()
    return events



TINY = 1e-8      # converged: best |A-B|/|B| below this
DROP = 0.5       # else: it must be this factor below the one at 1e3 times larger lambda
LOOSE = 1e-2     # and below this
PLATEAU = 1e-5   # or: below this and 1e3 times below the difference at the largest lambda
                 # (the difference can stop dropping at a precision floor of the pieces)


def reldiff(a, b):
    return abs(a - b) / abs(b) if b != 0 else (0.0 if a == 0 else float("inf"))


def judge(lams, diffs):
    if not diffs:
        return None, None
    k = min(range(len(diffs)), key=lambda i: diffs[i])
    best = diffs[k]
    if best < TINY:
        return True, best
    if best < PLATEAU and diffs[0] > 1e3 * best:
        return True, best
    ref = [i for i in range(len(lams)) if lams[i] >= 300 * lams[k]]
    if not ref or best >= LOOSE:
        return False, best
    r = min(ref, key=lambda i: abs(lams[i] / lams[k] / 1e3 - 1))
    return best <= DROP * diffs[r], best


def full_lambdas(points):
    labels = OrderedDict()
    for vals in points.values():
        for lab in vals:
            labels[lab] = None
    lams = sorted([lam for lam, vals in points.items() if len(vals) == len(labels)],
                  reverse=True)
    return list(labels), lams


def write_summary(f, limit, rels):
    f.write("\n\n SUMMARY %s (events ok / events; best |A-B|/|B| over lambda, worst event;\n"
            " ok: below %g; or below %g and 1e3 times below the one at the largest lambda;\n"
            " or below %g and %g times smaller than at 1e3 times larger lambda)\n"
            % (limit, TINY, PLATEAU, LOOSE, DROP))
    nbad = 0
    for (a, b), res in rels.items():
        nok = sum(1 for _, ok, _ in res if ok)
        ds = [d for _, _, d in res if d is not None]
        bad = [str(i) for i, ok, _ in res if not ok]
        flag = "ok" if nok == len(res) else "FAIL"
        nbad += flag == "FAIL"
        f.write(" %-14s -> %-14s %-5s %3d/%-3d  worst %-10s %s\n"
                % (a, b or "-", flag, nok, len(res), "%.2e" % max(ds) if ds else "-",
                   ("events " + ",".join(bad)) if bad else ""))
    return nbad


def check(path):
    name = os.path.basename(path)
    m = re.match(r"testRR_(\d+_\d+_\d+_\d+)\.log$", name)
    if not m:
        sys.exit("%s: expected testRR_<isec>_<jsec>_<c3p>_<d3p>.log" % name)
    tag = m.group(1)
    stype = sector_type(tag)
    events = read_testrr(path)
    out = os.path.join(os.path.dirname(path), "consistencyRR_%s.log" % tag)
    summary = OrderedDict()      # limit -> {(A, B): [(event, ok, best)]}
    skipped = OrderedDict()
    dash = "-" * 40
    with open(out, "w") as f:
        f.write(" consistency of the counterterms, sector %s (type %s), from %s\n" % (tag, stype, name))
        f.write(" relation REL = A -> B: under the limit, piece A tends to piece B\n")
        for evno, limits in events.items():
            f.write("\n %s\n %s\n  EVENT NUMBER %12d\n %s\n %s\n" % (dash, dash, evno, dash, dash))
            for limit, points in limits.items():
                if limit not in RELATIONS.get(stype, {}):
                    skipped[limit] = None
                    continue
                labels, lams = full_lambdas(points)
                if not lams:
                    continue
                f.write("\n\n LIM = %s\n" % limit)
                for a, b in RELATIONS[stype][limit](labels):
                    res = summary.setdefault(limit, OrderedDict()).setdefault((a, b), [])
                    if b is None:
                        f.write(" REL = %s -> ?   *** no piece to compare with\n" % a)
                        res.append((evno, False, None))
                        continue
                    f.write(" REL = %s -> %s\n" % (a, b))
                    f.write("%24s%24s%24s%24s\n" % ("lambda", a, b, "|A-B|/|B|"))
                    diffs = []
                    for lam in lams:
                        x, y = points[lam][a], points[lam][b]
                        diffs.append(reldiff(x, y))
                        f.write("%24.15E%24.15E%24.15E%24.15E\n" % (lam, x, y, diffs[-1]))
                    ok, best = judge(lams, diffs)
                    res.append((evno, bool(ok), best))
                    f.write(" %s (best |A-B|/|B| = %.3e)\n\n" % ("ok" if ok else "*** FAIL", best))
        nbad = 0
        for limit, rels in summary.items():
            nbad += write_summary(f, limit, rels)
        for limit in skipped:
            f.write("\n\n LIM = %s: no relations defined for sectors of type %s, not checked\n"
                    % (limit, stype))
    nrel = sum(len(rels) for rels in summary.values())
    if nbad:
        print("%s -> %s: %d of %d relations FAILED" % (name, os.path.basename(out), nbad, nrel))
    else:
        print("%s -> %s: all good (%d relations)" % (name, os.path.basename(out), nrel))
    return nbad


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("logs", nargs="+", help="testRR_<sector>.log")
    args = parser.parse_args()
    sys.exit(1 if sum(check(p) for p in args.logs) else 0)
