#!/usr/bin/env python3
"""Validate the region-set (2b) procedure against the shipped HAR reference.

2b differs from 2a: strata are CHR x AF_bin (not AF_bin alone) and the per-stratum
draw is max(1, floor(0.8*n)) rather than ceil(0.8*n) with a 1000 floor. The reported
rho is the FULL-data point estimate; the CI comes from the subsampling distribution;
p_boot = mean(rho_b <= 0) floored at 1/(B+1). 2b uses AgeMean_Jnt.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
import config

import sqlite3, numpy as np, time
from scipy import stats

DB = str(config.DB_PATH)
AF_LEVELS = ["AF_0-0.2", "AF_0.2-0.4", "AF_0.4-0.6", "AF_0.6-0.8", "AF_0.8-1.0"]
REF = {  # HAR_raw_SiteSubsample_80pctWithinChrxAF_2000fold_OverAllSummary.tsv
    "All":         dict(N=1519, rho=0.0390473658235253, lo=0.00789858428012878, hi=0.0624606131221699),
    "AF_0-0.2":    dict(N=571,  rho=-0.00764773561907364, lo=-0.0551331566929922, hi=0.0347385549357049),
    "AF_0.2-0.4":  dict(N=378,  rho=-0.0179317543524622, lo=-0.0751809202838259, hi=0.0401005683995998),
    "AF_0.4-0.6":  dict(N=260,  rho=0.0809786671855637, lo=0.0209404832274272, hi=0.145562509235998),
    "AF_0.6-0.8":  dict(N=188,  rho=0.0895754743602859, lo=0.00744708606636309, hi=0.17060840545732),
    "AF_0.8-1.0":  dict(N=122,  rho=-0.0242580653688938, lo=-0.142597206887291, hi=0.0609079854376175),
}

t0 = time.time()
db = sqlite3.connect(DB)
regions = db.execute("SELECT DISTINCT chrom,start,end FROM region_sets WHERE set_name='HAR'").fetchall()
print(f"HAR regions: {len(regions)}")

seen, rows = set(), []
for chrom, s, e in regions:
    for r in db.execute("""SELECT v.RSID, v.Chromosome, v.AF_bin, v."Evo2_40B_AvgRC_Delta", a.geva
                           FROM variants v JOIN allele_age a ON a.rsid=v.RSID
                           WHERE v.Chromosome=? AND v.Position BETWEEN ? AND ?""", (chrom, s, e)):
        if r[0] not in seen:
            seen.add(r[0]); rows.append(r)
print(f"variants in HAR (deduped, with allele-age row): {len(rows):,} ({time.time()-t0:.0f}s)")

chrom = np.array([str(r[1]) for r in rows], dtype=object)
afb = np.array([r[2] or "" for r in rows], dtype=object)
delta = np.array([np.nan if r[3] is None else float(r[3]) for r in rows])
age = np.array([np.nan if r[4] is None else float(r[4]) for r in rows])
ok = np.isfinite(delta) & np.isfinite(age) & (afb != "")
print(f"complete cases: {int(ok.sum()):,}   (reference All N = {REF['All']['N']})")

chrom, afb, delta, age = chrom[ok], afb[ok], delta[ok], age[ok]


def point_rho(i):
    if i.size < 3:
        return np.nan
    return stats.spearmanr(age[i], delta[i]).statistic


def run(B=2000, prop=0.8, seed=999):
    idx_all = np.arange(len(delta))
    strata = {}
    for k in idx_all:
        strata.setdefault((chrom[k], afb[k]), []).append(k)
    strata = {k: np.array(v) for k, v in strata.items()}
    sizes = {k: max(1, int(np.floor(prop * v.size))) for k, v in strata.items()}
    rng = np.random.default_rng(seed)
    rall, rbin = [], {lv: [] for lv in AF_LEVELS}
    n_each = 0
    for _ in range(B):
        picks = [(v if sizes[k] >= v.size else rng.choice(v, size=sizes[k], replace=False))
                 for k, v in strata.items()]
        sub = np.concatenate(picks); n_each = sub.size
        r = point_rho(sub)
        if r == r:
            rall.append(float(r))
        for lv in AF_LEVELS:
            s2 = sub[afb[sub] == lv]
            r2 = point_rho(s2)
            if r2 == r2:
                rbin[lv].append(float(r2))
    return rall, rbin, n_each


rall, rbin, n_each = run()
print(f"subsample size per draw (N_each): {n_each}  (reference summary.tsv: 1180)\n")
print(f"{'stratum':<12} {'field':<5} {'ours':>12} {'reference':>12} {'diff':>10}")
for name, rhos, sel in [("All", rall, np.ones(len(delta), bool))] + \
                       [(lv, rbin[lv], afb == lv) for lv in AF_LEVELS]:
    ref = REF[name]
    got = dict(N=int(sel.sum()), rho=float(point_rho(np.where(sel)[0])),
               lo=float(np.quantile(rhos, .025)), hi=float(np.quantile(rhos, .975)))
    for f in ("N", "rho", "lo", "hi"):
        print(f"{name:<12} {f:<5} {got[f]:>12.6f} {ref[f]:>12.6f} {got[f]-ref[f]:>+10.2e}")
print(f"\n({time.time()-t0:.0f}s)")
