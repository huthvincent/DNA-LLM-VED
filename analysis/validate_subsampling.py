#!/usr/bin/env python3
"""Validate our Python re-implementation of the July-12 subsampling procedure against
the reference outputs shipped with the fix package.

2a (genome-wide, pruned HM3): stratify by AF_bin, draw ceil(0.8*n) [>=1000, capped]
without replacement from every bin, B=2000, Spearman(delta, AgeMedian_Jnt).
The reference table used AgeMedian_Jnt, so we validate with AgeMedian_Jnt here even
though the database itself standardises on AgeMean_Jnt.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
import config

import sqlite3, numpy as np, time
from scipy import stats

DB = str(config.DB_PATH)
PRUNE = str(config.PRUNE_IN)
AF_LEVELS = ["AF_0-0.2", "AF_0.2-0.4", "AF_0.4-0.6", "AF_0.6-0.8", "AF_0.8-1.0"]

REF = {  # from Bootstrap2000_..._AF_bin.CI.tsv  (GEVA, Evo2_40B_AvgRC, AgeMedian_Jnt)
    "All":    dict(N=40863, rho=0.0162817093697937, lo=0.0118527414760614, hi=0.0206857966638352),
    "LCR":    dict(N=13244, rho=0.0186045085144159, lo=0.0110375548938476, hi=0.0259614635157215),
    "nonLCR": dict(N=27806, rho=0.0133042215728781, lo=0.00805056219806831, hi=0.0186524408272904),
}

t0 = time.time()
db = sqlite3.connect(DB)
rows = db.execute("""SELECT v.RSID, v.AF_bin, v.Repeat_Family, v."Evo2_40B_AvgRC_Delta", a.age_median
                     FROM variants v JOIN allele_age a ON a.rsid=v.RSID
                     WHERE a.age_median IS NOT NULL""").fetchall()
print(f"loaded {len(rows):,} rows with age_median ({time.time()-t0:.0f}s)")

prune = set(l.strip() for l in open(PRUNE) if l.strip())
rsid = np.array([r[0] for r in rows], dtype=object)
afb = np.array([r[1] or "" for r in rows], dtype=object)
lcr = np.array([r[2] is not None for r in rows], bool)
delta = np.array([np.nan if r[3] is None else float(r[3]) for r in rows])
agem = np.array([np.nan if r[4] is None else float(r[4]) for r in rows])
inp = np.array([r in prune for r in rsid], bool)

base = inp & np.isfinite(delta) & np.isfinite(agem) & (afb != "")
print(f"pruned & complete: {int(base.sum()):,}")


def subsample(mask, B=2000, frac=0.8, min_bin=1000, seed=123):
    idxs = np.where(mask)[0]
    tiers = [idxs[afb[idxs] == lv] for lv in AF_LEVELS]
    tiers = [t for t in tiers if t.size]
    rng = np.random.default_rng(seed)
    rhos = []
    n_each = 0
    for _ in range(B):
        picks = []
        for sub in tiers:
            na = sub.size
            k = min(max(int(np.ceil(frac * na)), min_bin), na)
            picks.append(sub if k >= na else rng.choice(sub, size=k, replace=False))
        idx = np.concatenate(picks)
        n_each = idx.size
        r = stats.spearmanr(delta[idx], agem[idx]).statistic
        if r == r:
            rhos.append(float(r))
    r = np.array(rhos)
    return dict(N=n_each, rho=float(r.mean()), lo=float(np.quantile(r, .025)),
                hi=float(np.quantile(r, .975)))


print(f"\n{'region':<8} {'field':<7} {'ours':>12} {'reference':>12} {'diff':>10}")
ok = True
for region, mask in (("All", base), ("LCR", base & lcr), ("nonLCR", base & ~lcr)):
    got = subsample(mask)
    ref = REF[region]
    for f in ("N", "rho", "lo", "hi"):
        d = got[f] - ref[f]
        flag = ""
        if f == "N":
            flag = "  <-- MISMATCH" if got[f] != ref[f] else ""
            ok &= got[f] == ref[f]
        else:
            # subsampling is stochastic: agree well within the Monte-Carlo error
            flag = "  <-- LARGE" if abs(d) > 3e-4 else ""
            ok &= abs(d) <= 3e-4
        print(f"{region:<8} {f:<7} {got[f]:>12.6f} {ref[f]:>12.6f} {d:>+10.2e}{flag}")
print(f"\n{'VALIDATION PASSED' if ok else 'VALIDATION FAILED'}  ({time.time()-t0:.0f}s)")
