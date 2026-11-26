#!/usr/bin/env python3
"""consistency check of the counterterms of the double real
"""
import argparse
import os
import re
import sys
from collections import ordereddict

def sector_type(tag):
    """ijjk, ijkj or ijkl from the tag isec_jsec_c3p_d3p."""
    isec, jsec, c3p, d3p = [int(x) for x in tag.split("_")]
    if jsec == c3p and jsec != d3p:
        return "ijjk"
    if jsec != c3p and jsec == d3p:
        return "ijkj"
    return "ijkl"

def operators(label):
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
            for partner in sets.get(target, [none]):
                pairs.append((lab, partner))
        return pairs
    return relations


def pair_rule(table):
    def relations(labels):
        have = set(labels)
        return [(a, b) for a, b in table if a in have and b in have]
    return relations


# ijjk (sij, cijk, scijk), sector 3_4_4_5 of e+e- > g g c cx
ijjk = {
    # k1 soft limit of i
    "si": set_rule(lambda t: t | {"s"}),
    # soft limit of the pair
    "sij": set_rule(lambda t: t | {"ss"}),
    # collinear limit of ij: the sc piece becomes cc
    "cij": set_rule(lambda t: t | ({"cc"} if "sc" in t else {"c"})),
    # soft-collinear limit of ijk: s where there is already a c
    "scijk": set_rule(lambda t: t | ({"s"} if "c" in t else {"sc"})),
    # triple collinear limit
    "cijk": set_rule(lambda t: t | {"cc"}),
    # spurious limits: pieces that coincide
    "cir": pair_rule([
        ("kc", "ks_c"), ("kcc", "ks_cc"), ("kc_cc", "ks_c_cc"),
        ("kc_ss", "ks_c_ss"), ("kc_ss_cc", "ks_c_ss_cc"),
        ("ksc", "ks_sc"), ("ksc_cc", "ks_sc_cc"),
        ("kss_cc", "ks_ss_cc"), ("kss_sc", "ks_ss_sc"),
        ("kss_sc_cc", "ks_ss_cc_sc")]),
    "cjr": pair_rule([
        ("kc", "ks_c"), ("kc_cc", "ks_c_cc"), ("kc_ss", "ks_c_ss"),
        ("kc_ss_cc", "ks_c_ss_cc"), ("ksc", "kss_sc"),
        ("ksc_cc", "kss_sc_cc"), ("ks_cc", "ks_ss_cc")]),
    "cijr": pair_rule([
        ("kcc", "kss_cc"), ("kc_cc", "kc_ss_cc"), ("ksc", "kss_sc"),
        ("ksc_cc", "kss_sc_cc"), ("ks_cc", "ks_ss_cc"),
        ("ks_c_cc", "ks_c_ss_cc"), ("ks_sc", "ks_ss_sc"),
        ("ks_sc_cc", "ks_ss_cc_sc")]),
}

relations = {"ijjk": ijjk}

def read_log(path):
    """{limit: [event, ...]}, event = ordereddict {lambda: {label: value}}."""
    events = ordereddict()
    current = none
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith("#"):
                current = none          # new event: its limit is on the next line
                continue
            lim, lam, label, val = line.split()
            if current is none:
                current = ordereddict()
                events.setdefault(lim, []).append(current)
            current.setdefault(float(lam), ordereddict())[label] = float(val)
    return events


def read_testrr(path):
    """same as read_log, from a testrr log of the former format."""
    events = ordereddict()
    limit = none
    cur = ordereddict()
    event = none
    with open(path) as f:
        for line in f:
            t = line.split()
            if not t or t[0].startswith("---") or t[0] == "lambda":
                continue
            if t[0] == "event":
                limit, event = none, none
                continue
            if t[0] == "lim":
                limit, cur, event = t[2], ordereddict(), none
                continue
            if limit is none:
                continue
            if t[0].startswith("k") and len(t) == 2:
                cur[t[0]] = float(t[1].replace("d", "e"))
                continue
            try:
                lam = float(t[0])
            except valueerror:
                continue
            if len(t) >= 4 and cur:
                if event is none:
                    event = ordereddict()
                    events.setdefault(limit, []).append(event)
                event[lam] = cur
                cur = ordereddict()
    return events



tiny = 1e-8      # converged: rel. difference below this at the smallest lambda
drop = 0.5       # else: it must drop by this factor over 3 decades of lambda
loose = 1e-2     # and be below this


def reldiff(a, b):
    return abs(a - b) / abs(b) if b != 0 else (0.0 if a == 0 else float("inf"))


def judge(lams, diffs):
    if not diffs:
        return None, None
    k = min(range(len(diffs)), key=lambda i: diffs[i])
    best = diffs[k]
    if best < TINY:
        return True, best
    ref = [i for i in range(len(lams)) if lams[i] >= 300 * lams[k]]
    if not ref or best >= LOOSE:
        return False, best
    # nearest to 1e3 times larger
    r = min(ref, key=lambda i: abs(lams[i] / lams[k] / 1e3 - 1))
    return best <= DROP * diffs[r], best


def full_lambdas(event):
    labels = OrderedDict()
    for vals in event.values():
        for lab in vals:
            labels[lab] = None
    lams = sorted([lam for lam, vals in event.items() if len(vals) == len(labels)], reverse=True)
    return list(labels), lams


def check_limit(f, stype, limit, events):
    summary = OrderedDict()
    f.write("\n\n LIM = %s\n" % limit)
    for iev, event in enumerate(events, 1):
        labels, lams = full_lambdas(event)
        if not lams:
            continue
        for a, b in RELATIONS[stype][limit](labels):
            if b is None:
                f.write("\n EVENT %d  RELATION  %s -> ?   *** no piece to compare with\n" % (iev, a))
                summary.setdefault((a, b), []).append((iev, False, None))
                continue
            f.write("\n EVENT %d  RELATION  %s -> %s\n" % (iev, a, b))
            f.write("%24s%24s%24s%24s\n" % ("lambda", a, b, "|A-B|/|B|"))
            diffs = []
            for lam in lams:
                x, y = event[lam][a], event[lam][b]
                diffs.append(reldiff(x, y))
                f.write("%24.15E%24.15E%24.15E%24.15E\n" % (lam, x, y, diffs[-1]))
            ok, d = judge(lams, diffs)
            summary.setdefault((a, b), []).append((iev, bool(ok), d))
    return summary


def write_summary(f, limit, summary):
    f.write("\n\n SUMMARY %s (events ok / events; best rel. diff over lambda, worst event;\n"
            " ok: below %g, or below %g and %g times smaller than at 1e3 times larger lambda)\n"
            % (limit, TINY, LOOSE, DROP))
    nbad = 0
    for (a, b), res in summary.items():
        nok = sum(1 for _, ok, _ in res if ok)
        ds = [d for _, _, d in res if d is not None]
        bad = [str(i) for i, ok, _ in res if not ok]
        flag = "ok" if nok == len(res) else "FAIL"
        nbad += flag == "FAIL"
        f.write(" %-14s -> %-14s %-5s %3d/%-3d  worst %-10s %s\n"
                % (a, b or "-", flag, nok, len(res), "%.2e" % max(ds) if ds else "-",
                   ("events " + ",".join(bad)) if bad else ""))
    return nbad


def check(path, testrr):
    name = os.path.basename(path)
    m = re.match(r"(?:testRR_|consistency_(?:[A-Za-z0-9]+)_)(\d+_\d+_\d+_\d+)", name)
    if not m:
        sys.exit("%s: cannot read the sector from the name" % name)
    tag = m.group(1)
    stype = sector_type(tag)
    events = (read_testrr if testrr else read_log)(path)
    out = os.path.join(os.path.dirname(path), "consistencyRR_" + name.split("_", 1)[1]
                       if not testrr else "consistencyRR_" + tag + ".log")
    nbad = 0
    with open(out, "w") as f:
        f.write(" sector %s, type %s, source %s\n" % (tag, stype, name))
        summaries = []
        for limit, evs in events.items():
            if limit not in RELATIONS.get(stype, {}):
                f.write("\n\n LIM = %s: no relations defined for sectors of type %s, skipped\n"
                        % (limit, stype))
                print("%s: no relations for the limit %s (%s)" % (name, limit, stype))
                continue
            summaries.append((limit, check_limit(f, stype, limit, evs)))
        for limit, summary in summaries:
            nbad += write_summary(f, limit, summary)
    print("%s -> %s: %d relations not ok" % (name, os.path.basename(out), nbad))
    return nbad


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("logs", nargs="+")
    parser.add_argument("--testrr", action="store_true",
                        help="the logs are testRR_<sector>.log of the former format")
    args = parser.parse_args()
    sys.exit(1 if sum(check(p, args.testrr) for p in args.logs) else 0)
