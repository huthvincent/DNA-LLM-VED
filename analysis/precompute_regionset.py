#!/usr/bin/env python3
"""Pre-compute the genome-wide Region-Set analyses for ALL 6 Evo2 models and archive them
to regionset_genomewide.json (baked into the Space image) so /regionset can show and export
genome-wide results for any model instantly, with no per-request computation.

Methods follow the maintainer's reference R code:

  Task 4A phastCons : Spearman( |delta| , phastCons100way )        -- overall + per MAF_tier  (3A.R)
  Task 4B phyloP    : Spearman( delta * sign , phyloP100way )      -- RAW delta, sign=-1 if NFE_AF>0.5,
                                                                      overall + per MAF_tier  (3B.R)
  Task 5  alleleage : Spearman( delta , {GEVA age, SDS} ) on the pruned-HM3 set,
                      region = All / LCR / non-LCR   (LCR = RepeatMasker repFamily present)
  Task 6  by AF_bin : the same, split by AF_bin
  Task 7G cCRE      : |delta| by ENCODE cCRE category, Mann-Whitney vs Intergenic (BH)   (3G.R)
  Task 7H coding    : |delta| by coding consequence,  Mann-Whitney vs the rest (BH)

Task 5/6 use the July-12 **subsampling** procedure (correlation_analysis_hm3_prune_boostrap.R),
NOT bootstrap resampling: each of B=2000 draws takes ceil(80%) of the sites -- without
replacement -- from every AF_bin stratum, then Spearman is computed on the pooled draw.
Reported rho is the FULL-data point estimate; the interval is the 2.5/97.5 percentile of the
subsampling distribution; significance is judged by the interval excluding zero. The directional
rate p_boot = mean(rho_b <= 0) is floored at 1/(B+1) and is NOT a conventional p-value.

Allele-age fields: geva = AgeMean_Jnt (the database standardises on the mean), sds = SDS_Final
(SDS re-signed to this database's allele orientation).
Validated against the shipped reference tables by validate_subsampling.py / validate_har.py.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
import config

import sqlite3, numpy as np, json, time, os
from scipy import stats

DB = str(config.DB_PATH)
PRUNE = str(config.PRUNE_IN)
OUT = str(config.REGIONSET_GW_JSON)

MODELS = ["Evo2_7B_NoRC_Delta", "Evo2_7B_AvgRC_Delta", "Evo2_7B_WeightedRC_Delta",
          "Evo2_40B_NoRC_Delta", "Evo2_40B_AvgRC_Delta", "Evo2_40B_WeightedRC_Delta"]
LABEL = {"Evo2_7B_NoRC_Delta": "7B · noRC", "Evo2_7B_AvgRC_Delta": "7B · avgRC", "Evo2_7B_WeightedRC_Delta": "7B · wtRC",
         "Evo2_40B_NoRC_Delta": "40B · noRC", "Evo2_40B_AvgRC_Delta": "40B · avgRC", "Evo2_40B_WeightedRC_Delta": "40B · wtRC"}
MAF_TIERS = ["Q1_0-20%", "Q2_20-40%", "Q3_40-60%", "Q4_60-80%", "Q5_80-100%"]
AF_LEVELS = ["AF_0-0.2", "AF_0.2-0.4", "AF_0.4-0.6", "AF_0.6-0.8", "AF_0.8-1.0"]
CODING = ["synonymous SNV", "nonsynonymous SNV", "stopgain", "stoploss", "startloss",
          "nonframeshift substitution", "frameshift substitution"]
B_SUB = 2000


def spear(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 10:
        return {"n": n, "rho": None, "p": None}
    r = stats.spearmanr(x[m], y[m])
    return {"n": n, "rho": round(float(r.statistic), 4), "p": float(r.pvalue)}


def bh(ps):
    ps = [p if p is not None else 1.0 for p in ps]
    n = len(ps); order = np.argsort(ps); adj = np.empty(n); prev = 1.0
    for rank, i in enumerate(reversed(order)):
        v = ps[i] * n / (n - rank); prev = min(prev, v); adj[i] = min(prev, 1.0)
    return [round(float(a), 6) for a in adj]


def _empty(n=0):
    return {"n": int(n), "n_each": 0, "B": 0, "rho": None, "rho_mean": None,
            "ci_lo": None, "ci_hi": None, "p_boot": None, "excludes_zero": False,
            "direction": "Interval includes zero"}


def subsample(delta, metric, strata, B=B_SUB, frac=0.8, min_bin=1000, seed=123):
    """July-12 stratified site-subsampling (2a rule: ceil(frac*n), floored at min_bin, capped
    at n, drawn WITHOUT replacement from each stratum). `strata` is a list of index arrays.
    Returns the full-data point rho plus the subsampling interval."""
    strata = [s for s in strata if s.size]
    if not strata:
        return _empty()
    full = np.concatenate(strata)
    n = full.size
    if n < 10:
        return _empty(n)
    rho_hat = stats.spearmanr(delta[full], metric[full]).statistic
    rng = np.random.default_rng(seed)
    rhos, n_each = [], 0
    for _ in range(B):
        picks = []
        for s in strata:
            k = min(max(int(np.ceil(frac * s.size)), min_bin), s.size)
            picks.append(s if k >= s.size else rng.choice(s, size=k, replace=False))
        idx = np.concatenate(picks); n_each = idx.size
        if idx.size <= 2:
            continue
        r = stats.spearmanr(delta[idx], metric[idx]).statistic
        if r == r:
            rhos.append(float(r))
    if not rhos:
        return _empty(n)
    r = np.array(rhos)
    lo, hi = float(np.quantile(r, .025)), float(np.quantile(r, .975))
    ex = lo > 0 or hi < 0
    return {"n": int(n), "n_each": int(n_each), "B": len(r),
            "rho": round(float(rho_hat), 5) if rho_hat == rho_hat else None,
            "rho_mean": round(float(r.mean()), 5),
            "ci_lo": round(lo, 5), "ci_hi": round(hi, 5),
            "p_boot": round(max(float((r <= 0).mean()), 1.0 / (len(r) + 1)), 5),
            "excludes_zero": bool(ex),
            "direction": "Robust positive" if lo > 0 else ("Robust negative" if hi < 0 else "Interval includes zero")}


def strata_by(valid, key):
    """Index arrays for `valid` rows grouped by the given per-row label array."""
    idxs = np.where(valid)[0]
    return [idxs[key[idxs] == lv] for lv in AF_LEVELS]


t0 = time.time()
db = sqlite3.connect(DB)
print("loading columns ...")
cols = (["Non_Finnish_European_AF", "MAF_tier", "AF_bin", "PhastCons_100way", "PhyloP_100way",
         "Repeat_Family", "ENCODE_cCRE_Annotation", "Exonic_Function", "RSID"] + MODELS)
colsel = ",".join(f'v."{c}"' for c in cols)
rows = db.execute(f'SELECT {colsel}, a.geva, a.sds '
                  f'FROM variants v LEFT JOIN allele_age a ON a.rsid=v.RSID').fetchall()
arr = list(zip(*rows)); ci = {c: i for i, c in enumerate(cols)}
n_total = len(rows)
print(f"  {n_total:,} rows ({time.time()-t0:.0f}s)")

f = lambda name: np.array([np.nan if v is None else float(v) for v in arr[ci[name]]])
af = f("Non_Finnish_European_AF"); phast = f("PhastCons_100way"); phylo = f("PhyloP_100way")
geva = np.array([np.nan if v is None else float(v) for v in arr[-2]])
sds = np.array([np.nan if v is None else float(v) for v in arr[-1]])
maf_tier = np.array([x or "" for x in arr[ci["MAF_tier"]]], dtype=object)
af_bin = np.array([x or "" for x in arr[ci["AF_bin"]]], dtype=object)
lcr = np.array([x is not None for x in arr[ci["Repeat_Family"]]], bool)
ccre = np.array([x or "NA" for x in arr[ci["ENCODE_cCRE_Annotation"]]], dtype=object)
exonic = np.array([x or "NA" for x in arr[ci["Exonic_Function"]]], dtype=object)
sign = np.where(af > 0.5, -1.0, 1.0)
prune = set(l.strip() for l in open(PRUNE) if l.strip())
in_prune = np.array([r in prune for r in arr[ci["RSID"]]], bool)
n_pruned_geva = int((in_prune & np.isfinite(geva)).sum())
n_pruned_sds = int((in_prune & np.isfinite(sds)).sum())
print(f"  pruned-HM3 {int(in_prune.sum()):,} (GEVA {n_pruned_geva:,}, SDS {n_pruned_sds:,}); "
      f"LCR(repFamily) {int(lcr.sum()):,}")

out = {"n_total": n_total, "models": {},
       "n_pruned": n_pruned_geva, "n_pruned_sds": n_pruned_sds,
       "method": ("Task 5/6: stratified site subsampling (B=2000, 80% of sites per AF_bin, without "
                  "replacement); rho is the full-data point estimate, the interval is the 2.5/97.5 "
                  "percentile of the subsampling distribution; significance = interval excludes zero. "
                  "Task 4A/4B remain stratified by MAF_tier."),
       "age_field": "AgeMean_Jnt", "sds_field": "SDS_Final",
       "maf_cutpoints_note": "MAF_tier = NFE-AF minor-allele quintiles; AF_bin = equal-width NFE-AF bins"}

for model in MODELS:
    delta = f(model)
    adelta = np.abs(delta)

    # ---- Task 4A / 4B : conservation, stratified by MAF_tier (unchanged, per 3A/3B.R)
    t4a = [dict(stratum="All", **spear(adelta, phast))]
    for tier in MAF_TIERS:
        sel = maf_tier == tier
        t4a.append(dict(stratum=tier, **spear(adelta[sel], phast[sel])))
    signed = delta * sign
    t4b = [dict(stratum="All", **spear(signed, phylo))]
    for tier in MAF_TIERS:
        sel = maf_tier == tier
        t4b.append(dict(stratum=tier, **spear(signed[sel], phylo[sel])))

    # ---- Task 5 / 6 : allele age & SDS on the pruned-HM3 set, AF_bin-stratified subsampling
    t5, t6 = {}, {}
    for mname, metric, mseed in (("GEVA", geva, 100), ("SDS", sds, 200)):
        base = in_prune & np.isfinite(metric) & np.isfinite(delta) & (af_bin != "")
        t5[mname] = [{"region": reg, **subsample(delta, metric, strata_by(m, af_bin), seed=mseed + j)}
                     for j, (reg, m) in enumerate((("All", base), ("LCR", base & lcr), ("non-LCR", base & ~lcr)))]
        # per-AF_bin: subsample within the single bin
        t6[mname] = []
        for j, lv in enumerate(AF_LEVELS):
            m = base & (af_bin == lv)
            t6[mname].append({"af_bin": lv, **subsample(delta, metric, [np.where(m)[0]], seed=mseed + 50 + j)})

    # ---- Task 7G / 7H : functional stratification
    ref = adelta[(ccre == "Intergenic") & np.isfinite(adelta)]
    cats = [c for c in ["Coding", "Splice_Region", "Promoter/Flanking", "UTR", "ncRNA", "Intronic", "Intergenic"]
            if (ccre == c).sum() > 0]
    t7g, praw = [], []
    for c in cats:
        g = adelta[(ccre == c) & np.isfinite(adelta)]
        p = None
        if c != "Intergenic" and len(g) >= 3 and len(ref) >= 3:
            p = float(stats.mannwhitneyu(g, ref, alternative="two-sided").pvalue)
        praw.append(p if p is not None else 1.0)
        t7g.append({"category": c, "n": int(len(g)), "mean": round(float(g.mean()), 3),
                    "median": round(float(np.median(g)), 3), "p": p})
    adj = bh(praw)
    for i, r in enumerate(t7g):
        r["p_adj"] = adj[i] if r["p"] is not None else None

    t7h, praw = [], []
    for c in [c for c in CODING if (exonic == c).sum() >= 3]:
        g = adelta[(exonic == c) & np.isfinite(adelta)]
        rest = adelta[(exonic != c) & (exonic != "NA") & np.isfinite(adelta)]
        p = float(stats.mannwhitneyu(g, rest, alternative="two-sided").pvalue) if len(rest) >= 3 else None
        praw.append(p if p is not None else 1.0)
        t7h.append({"category": c, "n": int(len(g)), "mean": round(float(g.mean()), 3),
                    "median": round(float(np.median(g)), 3), "p": p})
    adj = bh(praw)
    for i, r in enumerate(t7h):
        r["p_adj"] = adj[i] if r["p"] is not None else None

    out["models"][model] = {"label": LABEL[model], "task4A_phastCons": t4a, "task4B_phyloP": t4b,
                            "task5_GEVA": t5["GEVA"], "task5_SDS": t5["SDS"],
                            "task6_GEVA": t6["GEVA"], "task6_SDS": t6["SDS"],
                            "task7G_cCRE": t7g, "task7H_coding": t7h}
    print(f"  {LABEL[model]}: 4A/4B/5/6/7G/7H done ({time.time()-t0:.0f}s)", flush=True)

json.dump(out, open(OUT, "w"))
print(f"[done] wrote {OUT} ({os.path.getsize(OUT)/1e3:.0f} KB) in {time.time()-t0:.0f}s")
db.close()
