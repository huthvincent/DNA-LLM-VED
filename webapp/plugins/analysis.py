"""Datasette plugin: server-side compute for the analysis modules.

Adds JSON API routes that the analysis pages call. Keeping this inside Datasette
(via the register_routes hook) means we reuse the baked SQLite db, the existing
templates/static/canned-queries, and the same deploy pipeline -- we only add a
compute layer (numpy/scipy) on top.

Module 1 (this file): /api/region/analyze  -- Regional-Level Analysis (Panels A-D)
"""
import os, json, math
import numpy as np
from scipy import stats
from datasette import hookimpl
from datasette.utils.asgi import Response

# ---- model configurations (6 = 2 model sizes x 3 strand strategies) -----------
CONFIGS = [
    ("Evo2_7B_NoRC_Delta",        "7B · noRC"),
    ("Evo2_7B_AvgRC_Delta",       "7B · avgRC"),
    ("Evo2_7B_WeightedRC_Delta",  "7B · wtRC"),
    ("Evo2_40B_NoRC_Delta",       "40B · noRC"),
    ("Evo2_40B_AvgRC_Delta",      "40B · avgRC"),
    ("Evo2_40B_WeightedRC_Delta", "40B · wtRC"),
]
CONFIG_COLS = [c for c, _ in CONFIGS]
CONFIG_LABEL = dict(CONFIGS)
DEFAULT_CONFIG = "Evo2_40B_NoRC_Delta"

FUNC_ORDER = ["exonic", "UTR5", "UTR3", "ncRNA_exonic", "upstream", "downstream",
              "intronic", "ncRNA_intronic", "intergenic"]
EXONIC_ORDER = ["synonymous SNV", "nonsynonymous SNV", "stopgain", "stoploss", "startloss",
                "nonframeshift substitution", "frameshift substitution", "unknown"]
EXONIC_LABEL = {"synonymous SNV": "synonymous", "nonsynonymous SNV": "nonsynonymous",
                "nonframeshift substitution": "nonframeshift", "frameshift substitution": "frameshift"}

MAX_VARIANTS = 200000      # hard cap pulled into memory
LANDSCAPE_MAX = 8000       # downsample scatter (keep the most extreme |score|)
MIN_GROUP_N = 8            # below this we skip the statistical test for a group

_GENOME_STATS = None


def genome_stats():
    global _GENOME_STATS
    if _GENOME_STATS is None:
        p = os.path.join(os.path.dirname(__file__), "..", "genome_stats.json")
        try:
            _GENOME_STATS = json.load(open(p))
        except Exception:
            _GENOME_STATS = {"configs": {}, "dataset_version": "unknown"}
    return _GENOME_STATS


# ---- small stats helpers ------------------------------------------------------
def bh_adjust(pvals):
    """Benjamini-Hochberg FDR. None entries are passed through as None."""
    idx = [i for i, p in enumerate(pvals) if p is not None and not (isinstance(p, float) and math.isnan(p))]
    out = [None] * len(pvals)
    if not idx:
        return out
    p = np.array([pvals[i] for i in idx], float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.clip(ranked, 0, 1)
    for k, i in enumerate(idx):
        out[i] = float(adj[k])
    return out


def one_sample_test(vals, ref):
    """One-sample location test of |score| vs a reference value.
    t-test by default; Wilcoxon signed-rank if the data look non-normal."""
    v = np.asarray(vals, float)
    v = v[np.isfinite(v)]
    n = len(v)
    if n < MIN_GROUP_N or np.allclose(v, v[0]):
        return {"test": "n/a", "stat": None, "p": None, "n": int(n)}
    normal = True
    if n >= 20:
        try:
            normal = stats.normaltest(v).pvalue >= 0.05
        except Exception:
            normal = True
    try:
        if normal:
            r = stats.ttest_1samp(v, popmean=ref)
            return {"test": "t", "stat": float(r.statistic), "p": float(r.pvalue), "n": int(n)}
        else:
            d = v - ref
            d = d[d != 0]
            if len(d) < MIN_GROUP_N:
                return {"test": "n/a", "stat": None, "p": None, "n": int(n)}
            r = stats.wilcoxon(d)
            return {"test": "Wilcoxon", "stat": float(r.statistic), "p": float(r.pvalue), "n": int(n)}
    except Exception:
        return {"test": "n/a", "stat": None, "p": None, "n": int(n)}


def box_stats(vals):
    v = np.asarray(vals, float)
    v = v[np.isfinite(v)]
    if len(v) == 0:
        return None
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    iqr = q3 - q1
    lo_candidates = v[v >= q1 - 1.5 * iqr]
    hi_candidates = v[v <= q3 + 1.5 * iqr]
    lowerfence = float(lo_candidates.min()) if len(lo_candidates) else float(v.min())
    upperfence = float(hi_candidates.max()) if len(hi_candidates) else float(v.max())
    return {
        "n": int(len(v)),
        "q1": float(q1), "median": float(med), "q3": float(q3),
        "lowerfence": lowerfence, "upperfence": upperfence,
        "mean": float(v.mean()), "sem": float(v.std(ddof=1) / math.sqrt(len(v))) if len(v) > 1 else 0.0,
    }


def stratify(absvals, groups, order, ref_sel, ref_gw, label_map=None):
    """Build per-group box stats + one-sample tests vs selected-set mean and genome-wide mean."""
    out = []
    p_sel, p_gw = [], []
    present = [g for g in order if np.any(groups == g)]
    for g in present:
        v = absvals[groups == g]
        if len(v) < 3:
            continue
        bs = box_stats(v)
        t_sel = one_sample_test(v, ref_sel)
        t_gw = one_sample_test(v, ref_gw)
        p_sel.append(t_sel["p"]); p_gw.append(t_gw["p"])
        out.append({"name": (label_map or {}).get(g, g), "raw": g, "box": bs,
                    "test_vs_selected": t_sel, "test_vs_genomewide": t_gw})
    adj_sel = bh_adjust(p_sel); adj_gw = bh_adjust(p_gw)
    for i, o in enumerate(out):
        o["test_vs_selected"]["p_adj"] = adj_sel[i]
        o["test_vs_genomewide"]["p_adj"] = adj_gw[i]
    return out


def spearman_matrix(cols_data):
    """Pairwise-complete Spearman rho among the 6 config columns."""
    k = len(CONFIG_COLS)
    M = np.full((k, k), np.nan)
    for i in range(k):
        M[i, i] = 1.0
        for j in range(i + 1, k):
            a = cols_data[:, i]; b = cols_data[:, j]
            m = np.isfinite(a) & np.isfinite(b)
            if m.sum() >= 5:
                rho = stats.spearmanr(a[m], b[m]).statistic
                M[i, j] = M[j, i] = float(rho) if rho == rho else None
    return [[None if (x != x) else float(x) for x in row] for row in M]


# ---- data pull ----------------------------------------------------------------
def _pull(conn, where, params, model_col):
    cols = ["Variant_ID", "RSID", "Position", "Functional_Annotation", "Exonic_Function"] + CONFIG_COLS
    sql = (f'SELECT {",".join(chr(34)+c+chr(34) for c in cols)} FROM variants '
           f'WHERE {where} AND "{model_col}" IS NOT NULL LIMIT {MAX_VARIANTS + 1}')
    cur = conn.execute(sql, params)
    return cur.fetchall(), cols


async def region_analyze(request, datasette):
    args = request.args
    model = args.get("model") or DEFAULT_CONFIG
    if model not in CONFIG_COLS:
        return Response.json({"error": f"unknown model '{model}'"}, status=400)
    gene = (args.get("gene") or "").strip()
    chrom = (args.get("chr") or "").strip()
    start = args.get("start"); end = args.get("end")

    if gene:
        where = 'Gene LIKE ?'
        params = [f"%{gene}%"]
        label = f"gene {gene}"
    elif chrom and start and end:
        try:
            start = int(start); end = int(end)
        except ValueError:
            return Response.json({"error": "start/end must be integers"}, status=400)
        where = 'Chromosome = ? AND Position BETWEEN ? AND ?'
        params = [chrom, start, end]
        label = f"chr{chrom}:{start:,}-{end:,}"
    else:
        return Response.json({"error": "provide gene, or chr+start+end"}, status=400)

    db = datasette.get_database("evo2")
    rows, cols = await db.execute_fn(lambda conn: _pull(conn, where, params, model))
    if not rows:
        return Response.json({"error": f"No variants found for {label}."}, status=404)
    truncated = len(rows) > MAX_VARIANTS
    rows = rows[:MAX_VARIANTS]

    arr = list(zip(*rows))  # column-major
    ci = {c: i for i, c in enumerate(cols)}
    pos = np.array(arr[ci["Position"]], float)
    func = np.array([x if x is not None else "NA" for x in arr[ci["Functional_Annotation"]]], dtype=object)
    exonic = np.array([x if x is not None else "NA" for x in arr[ci["Exonic_Function"]]], dtype=object)
    rsid = arr[ci["RSID"]]
    config_data = np.array([[None if v is None else float(v) for v in arr[ci[c]]] for c in CONFIG_COLS],
                           dtype=float).T  # shape (n, 6)
    sel = config_data[:, CONFIG_COLS.index(model)]
    absvals = np.abs(sel)
    finite = np.isfinite(absvals)

    gw = genome_stats()["configs"].get(model, {})
    ref_gw = gw.get("mean_abs")
    ref_sel = float(np.nanmean(absvals[finite])) if finite.any() else None

    panelA = stratify(absvals[finite], func[finite], FUNC_ORDER, ref_sel, ref_gw)
    panelB = stratify(absvals[finite], exonic[finite], EXONIC_ORDER, ref_sel, ref_gw, EXONIC_LABEL)
    panelC = {"labels": [CONFIG_LABEL[c] for c in CONFIG_COLS], "matrix": spearman_matrix(config_data)}

    # Panel D landscape: |score| vs position, downsample keeping the most extreme
    n = int(finite.sum())
    idx = np.where(finite)[0]
    if n > LANDSCAPE_MAX:
        order = idx[np.argsort(-absvals[idx])]
        keep = set(order[:LANDSCAPE_MAX // 2].tolist())
        rest = [i for i in idx.tolist() if i not in keep]
        step = max(1, len(rest) // (LANDSCAPE_MAX // 2))
        keep.update(rest[::step])
        idx = np.array(sorted(keep))
        d_truncated = True
    else:
        d_truncated = False
    panelD = {
        "downsampled": d_truncated,
        "x": [pos[i] / 1e6 for i in idx],
        "y": [float(absvals[i]) for i in idx],
        "signed": [float(sel[i]) for i in idx],
        "func": [func[i] for i in idx],
        "rsid": [rsid[i] for i in idx],
    }

    return Response.json({
        "query": {"label": label, "gene": gene or None, "chr": chrom or None,
                  "start": start if not gene else None, "end": end if not gene else None,
                  "model": model, "model_label": CONFIG_LABEL[model], "n": n, "truncated": truncated},
        "configs": [{"col": c, "label": CONFIG_LABEL[c]} for c in CONFIG_COLS],
        "selected_mean_abs": ref_sel,
        "genomewide_mean_abs": ref_gw,
        "dataset_version": genome_stats().get("dataset_version"),
        "panelA": panelA, "panelB": panelB, "panelC": panelC, "panelD": panelD,
    })


# =============================================================================
# Module 2 -- Region Set-Based Analysis (Tasks 4-7)
# Triggered by a built-in region set (e.g. HAR) or an uploaded BED file.
# =============================================================================
LCR_TOKENS = ("Low_complexity", "Simple_repeat", "Satellite")
AF_BINS = [("All", None), ("0.05-0.2", (0.05, 0.2)), ("0.2-0.4", (0.2, 0.4)),
           ("0.4-0.6", (0.4, 0.6)), ("0.6-0.8", (0.6, 0.8)), ("0.8-0.95", (0.8, 0.95))]
FUNC_ORDER7 = ["splicing", "exonic", "UTR5", "UTR3", "ncRNA_exonic", "ncRNA_splicing",
               "upstream", "downstream", "intronic", "ncRNA_intronic", "intergenic"]
EXONIC_ORDER7 = ["stopgain", "nonsynonymous SNV", "startloss", "stoploss", "synonymous SNV",
                 "nonframeshift substitution", "frameshift substitution", "unknown"]
MAX_SET_VARIANTS = 300000
MAX_REGIONS = 60000
GW_SAMPLE = 80000
BOOT_B = 200
BOOT_SUBSAMPLE = 20000

PULL_COLS = (["RSID", "Repeat_Class", "Repeat_Family", "Chromosome", "AF_bin",
              "Functional_Annotation", "Exonic_Function",
              "PhastCons_100way", "PhyloP_100way", "Joint_AF"] + CONFIG_COLS + ["__geva", "__sds"])

AF_LEVELS = ["AF_0-0.2", "AF_0.2-0.4", "AF_0.4-0.6", "AF_0.6-0.8", "AF_0.8-1.0"]
SUB_B = 2000          # subsampling draws (reference R uses B = 2000)
SUB_FRAC = 0.8        # 80% of the sites in each stratum, drawn WITHOUT replacement


def is_lcr(repclass):
    if not repclass:
        return False
    return any(tok in repclass for tok in LCR_TOKENS)


def boot_spearman(x, y, B=BOOT_B, seed=0, sub=BOOT_SUBSAMPLE, samples=False):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]; n = len(x)
    if n < 10:
        return {"rho": None, "mean": None, "sem": None, "n": int(n), "p": None, "samples": []}
    if n > sub:
        rng = np.random.default_rng(seed); idx = rng.choice(n, sub, replace=False)
        x, y = x[idx], y[idx]; n = sub
    rho_obs = stats.spearmanr(x, y).statistic
    rng = np.random.default_rng(seed + 1)
    boots = []
    for _ in range(B):
        bi = rng.integers(0, n, n)
        r = stats.spearmanr(x[bi], y[bi]).statistic
        if r == r:
            boots.append(float(r))
    arr = np.array(boots) if boots else np.array([0.0])
    p = float(2 * min((arr >= 0).mean(), (arr <= 0).mean()))
    out = {"rho": float(rho_obs) if rho_obs == rho_obs else None,
           "mean": float(arr.mean()), "sem": float(arr.std(ddof=1)) if len(arr) > 1 else 0.0,
           "n": int(n), "p": p}
    if samples:
        out["samples"] = [round(v, 4) for v in boots[:300]]
    return out


def conservation_panel(absmat, cons, seed):
    """Spearman(|delta|, conservation) per config with paired bootstrap + key pairwise tests."""
    cons = np.asarray(cons, float)
    valid = np.isfinite(cons)
    idx_all = np.where(valid)[0]
    n = len(idx_all)
    if n < 10:
        return {"configs": [], "pairwise": [], "n": int(n)}
    rng = np.random.default_rng(seed)
    if n > BOOT_SUBSAMPLE:
        idx_all = rng.choice(idx_all, BOOT_SUBSAMPLE, replace=False); n = BOOT_SUBSAMPLE
    c = cons[idx_all]
    boot_idx = [rng.integers(0, n, n) for _ in range(BOOT_B)]   # shared across configs (paired)
    per = []
    boot_rhos = {}
    for k, col in enumerate(CONFIG_COLS):
        a = np.abs(absmat[idx_all, k])
        rho_obs = stats.spearmanr(a, c, nan_policy="omit").statistic
        br = np.array([stats.spearmanr(a[bi], c[bi], nan_policy="omit").statistic for bi in boot_idx])
        br = br[np.isfinite(br)]
        boot_rhos[col] = br
        per.append({"config": col, "label": CONFIG_LABEL[col],
                    "rho": float(rho_obs) if rho_obs == rho_obs else None,
                    "mean": float(br.mean()) if len(br) else None,
                    "sem": float(br.std(ddof=1)) if len(br) > 1 else 0.0})
    pairs = [("Evo2_7B_NoRC_Delta", "Evo2_40B_NoRC_Delta"),
             ("Evo2_7B_AvgRC_Delta", "Evo2_40B_AvgRC_Delta"),
             ("Evo2_40B_NoRC_Delta", "Evo2_40B_AvgRC_Delta"),
             ("Evo2_7B_NoRC_Delta", "Evo2_7B_AvgRC_Delta")]
    pairwise = []
    praw = []
    for i, j in pairs:
        d = boot_rhos[i][:min(len(boot_rhos[i]), len(boot_rhos[j]))] - boot_rhos[j][:min(len(boot_rhos[i]), len(boot_rhos[j]))]
        p = float(2 * min((d >= 0).mean(), (d <= 0).mean())) if len(d) else None
        praw.append(p)
        pairwise.append({"a": CONFIG_LABEL[i], "b": CONFIG_LABEL[j], "delta_rho": float(d.mean()) if len(d) else None, "p": p})
    adj = bh_adjust(praw)
    for k, pw in enumerate(pairwise):
        pw["p_adj"] = adj[k]
    return {"configs": per, "pairwise": pairwise, "n": int(n)}


def site_subsample(score, metric, strata, B=SUB_B, frac=SUB_FRAC, seed=1, samples=False):
    """July-12 stratified site-subsampling used for every region-set allele-age analysis
    (HAR_hg38_V2.R / RBA_ASD_test.R): in each of B draws take max(1, floor(frac*n)) sites
    -- WITHOUT replacement -- from every CHR x AF_bin stratum, pool them, and take
    Spearman(metric, score).

    This is subsampling, NOT bootstrap resampling: rho is the FULL-data point estimate and
    the interval is the 2.5/97.5 percentile of the subsampling distribution. Significance is
    read off the interval (excludes zero -> robust). `p_boot` is the directional subsampling
    rate mean(rho_b <= 0), floored at 1/(B+1); it is not a conventional p-value.
    """
    score = np.asarray(score, float); metric = np.asarray(metric, float)
    strata = [np.asarray(s) for s in strata if len(s)]
    empty = {"n": 0, "n_each": 0, "B": 0, "rho": None, "rho_mean": None, "ci_lo": None,
             "ci_hi": None, "p_boot": None, "excludes_zero": False,
             "direction": "Interval includes zero", "samples": []}
    if not strata:
        return empty
    full = np.concatenate(strata)
    if full.size < 10:
        return {**empty, "n": int(full.size)}
    rho_hat = stats.spearmanr(metric[full], score[full]).statistic
    rng = np.random.default_rng(seed)
    sizes = [max(1, int(np.floor(frac * s.size))) for s in strata]
    rhos, n_each = [], 0
    for _ in range(B):
        picks = [s if k >= s.size else rng.choice(s, size=k, replace=False)
                 for s, k in zip(strata, sizes)]
        idx = np.concatenate(picks); n_each = idx.size
        if idx.size <= 2:
            continue
        r = stats.spearmanr(metric[idx], score[idx]).statistic
        if r == r:
            rhos.append(float(r))
    if not rhos:
        return {**empty, "n": int(full.size)}
    r = np.array(rhos)
    lo, hi = float(np.quantile(r, .025)), float(np.quantile(r, .975))
    out = {"n": int(full.size), "n_each": int(n_each), "B": len(r),
           "rho": float(rho_hat) if rho_hat == rho_hat else None,
           "rho_mean": float(r.mean()),
           "ci_lo": lo, "ci_hi": hi,
           "p_boot": max(float((r <= 0).mean()), 1.0 / (len(r) + 1)),
           "excludes_zero": bool(lo > 0 or hi < 0),
           "direction": "Robust positive" if lo > 0 else ("Robust negative" if hi < 0 else "Interval includes zero")}
    if samples:
        out["samples"] = [round(v, 5) for v in rhos[:300]]
    return out


def _chr_af_strata(mask, chrom, af_bin):
    """Index arrays for the CHR x AF_bin strata among the rows selected by `mask`."""
    idxs = np.where(mask)[0]
    groups = {}
    for i in idxs:
        groups.setdefault((chrom[i], af_bin[i]), []).append(i)
    return [np.array(v) for v in groups.values()]


def allele_age_panel(score, metric, lcr, chrom, af_bin, seed):
    """Allele-age / SDS correlation for All / LCR / non-LCR, via CHR x AF_bin site subsampling."""
    score = np.asarray(score, float); metric = np.asarray(metric, float); lcr = np.asarray(lcr, bool)
    ok = np.isfinite(score) & np.isfinite(metric) & (af_bin != "")
    out = []
    for j, (name, mask) in enumerate((("All", ok), ("LCR", ok & lcr), ("non-LCR", ok & ~lcr))):
        r = site_subsample(score, metric, _chr_af_strata(mask, chrom, af_bin), seed=seed + j, samples=True)
        r["subset"] = name
        out.append(r)
    by = {r["subset"]: r for r in out}
    cmp = None
    a = by["LCR"].get("samples") or []; b = by["non-LCR"].get("samples") or []
    if len(a) >= 20 and len(b) >= 20:
        aa = np.array(a); bb = np.array(b); m = min(len(aa), len(bb))
        diff = aa[:m] - bb[:m]
        lo, hi = float(np.quantile(diff, .025)), float(np.quantile(diff, .975))
        cmp = {"delta_rho": float(aa.mean() - bb.mean()), "ci_lo": lo, "ci_hi": hi,
               "excludes_zero": bool(lo > 0 or hi < 0)}
    n_with = int(np.isfinite(metric).sum())
    return {"subsets": out, "compare_lcr_vs_nonlcr": cmp,
            "n_with": n_with, "n_excluded": int(len(metric) - n_with)}


def af_binned_panel(score, metric, chrom, af_bin, seed):
    """The same subsampling, split by AF_bin (strata are the chromosomes within each bin)."""
    score = np.asarray(score, float); metric = np.asarray(metric, float)
    ok = np.isfinite(score) & np.isfinite(metric) & (af_bin != "")
    out = []
    for j, lv in enumerate(["All"] + AF_LEVELS):
        mask = ok if lv == "All" else (ok & (af_bin == lv))
        r = site_subsample(score, metric, _chr_af_strata(mask, chrom, af_bin), seed=seed + j)
        r["bin"] = lv
        out.append(r)
    return {"bins": out, "trend": None}


def _intersect(conn, regions, model):
    """Per-region indexed range query (uses idx_chr_pos), deduped by RSID.
    Deterministic and fast; avoids the planner choosing a full-scan range join."""
    colsql = ", ".join(f'v."{c}"' for c in PULL_COLS[:-2])
    stmt = (f'SELECT {colsql}, a.geva, a.sds FROM variants v '
            f'LEFT JOIN allele_age a ON a.rsid=v.RSID '
            f'WHERE v.Chromosome=? AND v.Position BETWEEN ? AND ? AND v."{model}" IS NOT NULL')
    seen = {}
    for chrom, start, end in regions[:MAX_REGIONS]:
        for row in conn.execute(stmt, (chrom, start, end)):
            rid = row[0]  # RSID is PULL_COLS[0]
            if rid not in seen:
                seen[rid] = row
        if len(seen) > MAX_SET_VARIANTS:
            break
    return list(seen.values())


def _sample_genomewide(conn, model, n):
    """Genome-wide example: a uniform sample (via rowid modulo) PLUS every coding
    variant, so each functional-region and coding-consequence class has enough n."""
    colsql = ", ".join(f'v."{c}"' for c in PULL_COLS[:-2])
    total = conn.execute("SELECT max(rowid) FROM variants").fetchone()[0] or n
    step = max(1, total // n)
    stmt = (f'SELECT {colsql}, a.geva, a.sds FROM variants v LEFT JOIN allele_age a ON a.rsid=v.RSID '
            f'WHERE (v.rowid % {step} = 0 OR v.Exonic_Function IS NOT NULL) AND v."{model}" IS NOT NULL '
            f'LIMIT {MAX_SET_VARIANTS}')
    return [tuple(r) for r in conn.execute(stmt)]


async def regionset_analyze(request, datasette):
    args = request.args
    model = args.get("model") or DEFAULT_CONFIG
    if model not in CONFIG_COLS:
        return Response.json({"error": f"unknown model '{model}'"}, status=400)
    set_name = (args.get("set") or "").strip()
    db = datasette.get_database("evo2")
    n_regions = 0
    if set_name == "genomewide":
        label = "Genome-wide (random sample + all coding)"
        rows = await db.execute_fn(lambda conn: _sample_genomewide(conn, model, GW_SAMPLE))
    else:
        if set_name:
            label = set_name
            rs = await db.execute("SELECT chrom,start,end FROM region_sets WHERE set_name=?", [set_name])
            regions = [(r["chrom"], r["start"], r["end"]) for r in rs]
            if not regions:
                return Response.json({"error": f"region set '{set_name}' not found"}, status=404)
        else:
            label = args.get("label") or "uploaded BED"
            body = (await request.post_body()).decode("utf-8", "replace") if request.method == "POST" else ""
            regions = []
            for line in body.splitlines():
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                p = line.split()
                if len(p) < 3:
                    continue
                try:
                    regions.append((p[0].replace("chr", ""), int(p[1]), int(p[2])))
                except ValueError:
                    continue
            if not regions:
                return Response.json({"error": "no valid BED regions (need chrom start end; tab/space separated)"}, status=400)
        n_regions = len(regions)
        rows = await db.execute_fn(lambda conn: _intersect(conn, regions, model))
    if not rows:
        return Response.json({"error": f"No variants in {label}."}, status=404)
    truncated = len(rows) > MAX_SET_VARIANTS
    rows = rows[:MAX_SET_VARIANTS]

    arr = list(zip(*rows)); ci = {c: i for i, c in enumerate(PULL_COLS)}
    chrom = np.array([str(x) if x is not None else "" for x in arr[ci["Chromosome"]]], dtype=object)
    af_bin = np.array([x if x is not None else "" for x in arr[ci["AF_bin"]]], dtype=object)
    func = np.array([x if x is not None else "NA" for x in arr[ci["Functional_Annotation"]]], dtype=object)
    exonic = np.array([x if x is not None else "NA" for x in arr[ci["Exonic_Function"]]], dtype=object)
    phast = np.array([np.nan if v is None else float(v) for v in arr[ci["PhastCons_100way"]]])
    phylo = np.array([np.nan if v is None else float(v) for v in arr[ci["PhyloP_100way"]]])
    af = np.array([np.nan if v is None else float(v) for v in arr[ci["Joint_AF"]]])
    geva = np.array([np.nan if v is None else float(v) for v in arr[ci["__geva"]]])
    sds = np.array([np.nan if v is None else float(v) for v in arr[ci["__sds"]]])
    # LCR follows the reference R code: RepeatMasker repFamily present (not repClass tokens)
    lcr = np.array([x is not None for x in arr[ci["Repeat_Family"]]], bool)
    config_data = np.array([[np.nan if v is None else float(v) for v in arr[ci[c]]] for c in CONFIG_COLS]).T
    sel = config_data[:, CONFIG_COLS.index(model)]

    n = len(rows)
    n_geva = int(np.isfinite(geva).sum()); n_sds = int(np.isfinite(sds).sum())

    # Task 4 -- conservation (Panel A phastCons, Panel B phyloP)
    task4 = {"phastCons": conservation_panel(config_data, phast, seed=11),
             "phyloP": conservation_panel(config_data, phylo, seed=12)}
    # Task 5 -- allele age (Panel C GEVA, Panel D SDS) by All/LCR/non-LCR
    task5 = {"GEVA": allele_age_panel(sel, geva, lcr, chrom, af_bin, seed=21),
             "SDS": allele_age_panel(sel, sds, lcr, chrom, af_bin, seed=24),
             "n_lcr": int(lcr.sum()), "n_nonlcr": int((~lcr).sum())}
    # Task 6 -- the same subsampling, split by AF_bin (equal-width Non-Finnish-European AF bins)
    task6 = {"GEVA": af_binned_panel(sel, geva, chrom, af_bin, seed=31),
             "SDS": af_binned_panel(sel, sds, chrom, af_bin, seed=38)}
    # Task 7 -- functional property: log2|delta| by functional region (G) and consequence (H)
    absL = np.log2(np.abs(sel) + 1e-6)
    fin = np.isfinite(absL)
    ref_sel = float(np.nanmean(absL[fin])) if fin.any() else None
    gw_abs = genome_stats()["configs"].get(model, {}).get("mean_abs")
    ref_gw = float(np.log2(gw_abs)) if gw_abs else None
    task7 = {"unit": "log2|Evo2 Δ|", "ref_selected": ref_sel, "ref_genomewide": ref_gw,
             "panelG": stratify(absL[fin], func[fin], FUNC_ORDER7, ref_sel, ref_gw),
             "panelH": stratify(absL[fin], exonic[fin], EXONIC_ORDER7, ref_sel, ref_gw, EXONIC_LABEL)}

    return Response.json({
        "query": {"label": label, "set": set_name or None, "model": model, "model_label": CONFIG_LABEL[model],
                  "n_regions": n_regions, "n": n, "n_geva": n_geva, "n_sds": n_sds, "truncated": truncated},
        "configs": [{"col": c, "label": CONFIG_LABEL[c]} for c in CONFIG_COLS],
        "dataset_version": genome_stats().get("dataset_version"),
        "task4": task4, "task5": task5, "task6": task6, "task7": task7,
    })


# =============================================================================
# Module 3 -- GWAS / Heritability (Tasks 8-9): S-LDSC coefficient z-score
# UKB / FinnGen / MVP, never mixed. Data baked as gwas_ldsc.json (no db needed).
# =============================================================================
_GWAS = None


def gwas_data():
    global _GWAS
    if _GWAS is None:
        p = os.path.join(os.path.dirname(__file__), "..", "gwas_ldsc.json")
        try:
            _GWAS = json.load(open(p))
        except Exception:
            _GWAS = {"sources": {}, "metric": "S-LDSC coefficient z-score"}
    return _GWAS


def _rank_quantile(zs, z):
    """Rank (1 = highest z) and percentile of z within the list zs."""
    total = len(zs)
    if not total:
        return None, 0, None
    greater = sum(1 for v in zs if v > z)
    less_eq = sum(1 for v in zs if v <= z)
    return greater + 1, total, round(100.0 * less_eq / total, 1)


async def gwas_meta(request, datasette):
    g = gwas_data()["sources"]
    meta = {}
    for s, recs in g.items():
        cats = {}
        for r in recs:
            cats[r["category"]] = cats.get(r["category"], 0) + 1
        meta[s] = {"n": len(recs), "categories": sorted(cats.items())}
    return Response.json({"sources": meta, "metric": gwas_data().get("metric")})


async def gwas_search(request, datasette):
    src = request.args.get("source") or "UKB"
    q = (request.args.get("q") or "").lower().strip()
    want_all = request.args.get("all")
    recs = gwas_data()["sources"].get(src, [])
    rank = {r["id"]: i + 1 for i, r in enumerate(recs)}   # recs are sorted z-descending
    res = [r for r in recs if (not q) or q in r["name"].lower() or q in r["id"].lower() or q in r["category"].lower()]
    lim = len(res) if want_all else 60
    return Response.json({"source": src, "n_total": len(recs),
                          "results": [{"id": r["id"], "name": r["name"], "category": r["category"],
                                       "z": r["z"], "p": r.get("p"), "n": r.get("n"), "rank": rank[r["id"]]} for r in res[:lim]]})


async def gwas_trait(request, datasette):
    src = request.args.get("source") or "UKB"
    tid = request.args.get("id")
    recs = gwas_data()["sources"].get(src, [])
    if not recs:
        return Response.json({"error": f"unknown GWAS database '{src}'"}, status=400)
    rec = next((r for r in recs if r["id"] == tid), None)
    if rec is None:
        return Response.json({"error": "The selected trait was not found in the selected GWAS database."}, status=404)
    z_all = [r["z"] for r in recs]
    cat = rec["category"]
    z_cat = [r["z"] for r in recs if r["category"] == cat]
    ro, to, po = _rank_quantile(z_all, rec["z"])
    rc, tc, pc = _rank_quantile(z_cat, rec["z"])
    return Response.json({
        "source": src, "metric": gwas_data().get("metric"), "trait": rec,
        "overall": {"rank": ro, "total": to, "percentile": po, "z": z_all},
        "category": {"name": cat, "rank": rc, "total": tc, "percentile": pc, "z": z_cat},
    })


async def gwas_overview(request, datasette):
    src = request.args.get("source") or "UKB"
    recs = gwas_data()["sources"].get(src, [])
    if not recs:
        return Response.json({"error": f"unknown GWAS database '{src}'"}, status=400)
    bycat = {}
    for r in recs:
        bycat.setdefault(r["category"], []).append(r)
    cats = []
    for c, rs in bycat.items():
        a = np.array([r["z"] for r in rs])
        cats.append({"category": c, "n": len(rs), "mean": float(a.mean()), "median": float(np.median(a)),
                     "frac_pos": float((a > 0).mean()),
                     "traits": [{"id": r["id"], "name": r["name"], "z": r["z"], "p": r.get("p"), "n": r.get("n")} for r in rs]})
    cats.sort(key=lambda d: -d["median"])
    return Response.json({"source": src, "metric": gwas_data().get("metric"), "n": len(recs), "categories": cats})


# =============================================================================
# Module 4 -- Evo2 score x GWAS P-value co-localization (v1.2)
# Regional view: Evo2 track + GWAS -log10(P) track(s) + annotation track.
# GWAS associations (P<=1e-4) live in gwas_assoc, linked to variants by rsID.
# =============================================================================
COLOC_MAX_VARIANTS = 25000


async def coloc_traits(request, datasette):
    """Search registered GWAS studies (for the multi-trait selector)."""
    q = (request.args.get("q") or "").lower().strip()
    cohort = request.args.get("cohort")
    db = datasette.get_database("evo2")
    rows = await db.execute("SELECT study_id, cohort, trait_name, phenotype_category, sample_size, ancestry FROM gwas_studies")
    res = []
    for r in rows:
        co = r["cohort"] or ""
        if cohort and co != cohort:
            continue
        name = (r["trait_name"] or "").lower(); cat = (r["phenotype_category"] or "").lower()
        if q and q not in name and q not in cat and q not in co.lower():
            continue
        res.append({"study_id": r["study_id"], "cohort": co, "name": r["trait_name"],
                    "category": r["phenotype_category"], "n": r["sample_size"], "ancestry": r["ancestry"]})
    return Response.json({"n_total": len(res), "results": res[:80]})


async def coloc_region(request, datasette):
    """Region co-localization data: region variants (position, Evo2, annotation) +
    per selected trait the GWAS -log10(P) points that fall in the region."""
    import math
    args = request.args
    model = args.get("model") or DEFAULT_CONFIG
    if model not in CONFIG_COLS:
        return Response.json({"error": f"unknown model '{model}'"}, status=400)
    studies = [int(s) for s in (args.get("studies") or "").split(",") if s.strip().lstrip("-").isdigit()][:8]
    db = datasette.get_database("evo2")

    gene = (args.get("gene") or "").strip()
    if gene:
        gr = (await db.execute("SELECT Chromosome c, min(Position) s, max(Position) e FROM variants WHERE Gene LIKE ?", ["%" + gene + "%"])).first()
        if not gr or gr["s"] is None:
            return Response.json({"error": f"gene '{gene}' not found"}, status=404)
        chrom, start, end = gr["c"], gr["s"], gr["e"]; label = f"gene {gene}"
    else:
        try:
            chrom = args.get("chr"); start = int(args.get("start")); end = int(args.get("end"))
        except (TypeError, ValueError):
            return Response.json({"error": "need gene, or chr + start + end"}, status=400)
        label = f"chr{chrom}:{start:,}-{end:,}"

    rows = await db.execute(
        f'SELECT Position p, RSID r, Reference_Allele ref, Alternate_Allele alt, "{model}" e, Functional_Annotation f '
        f'FROM variants WHERE Chromosome=? AND Position BETWEEN ? AND ? ORDER BY Position LIMIT {COLOC_MAX_VARIANTS + 1}', [chrom, start, end])
    variants = [{"pos": x["p"], "rsid": x["r"], "ref": x["ref"], "alt": x["alt"], "evo2": x["e"], "func": x["f"]} for x in rows]
    truncated = len(variants) > COLOC_MAX_VARIANTS
    variants = variants[:COLOC_MAX_VARIANTS]
    rsid_pos = {v["rsid"]: v["pos"] for v in variants}

    traits = []
    p_rsids = set()
    for sid in studies:
        meta = (await db.execute("SELECT trait_name, cohort, phenotype_category, sample_size FROM gwas_studies WHERE study_id=?", [sid])).first()
        if not meta:
            continue
        pts = []
        for a in await db.execute("SELECT rsid, p_value FROM gwas_assoc WHERE study_id=?", [sid]):
            pos = rsid_pos.get(a["rsid"])
            if pos is not None and a["p_value"] is not None:
                pv = a["p_value"]
                pts.append({"pos": pos, "rsid": a["rsid"], "p": pv,
                            "nlp": round(-math.log10(pv), 2) if pv > 0 else 320.0})
                p_rsids.add(a["rsid"])
        pts.sort(key=lambda d: d["pos"])
        traits.append({"study_id": sid, "name": meta["trait_name"], "cohort": meta["cohort"],
                       "category": meta["phenotype_category"], "n": meta["sample_size"], "points": pts})

    # only display sites that have BOTH an Evo2 score AND a GWAS P-value (so every track
    # shares the same set of variants). Applies only when >=1 trait is selected.
    if studies:
        variants = [v for v in variants if v["rsid"] in p_rsids]

    return Response.json({
        "query": {"label": label, "chr": chrom, "start": start, "end": end,
                  "model": model, "model_label": CONFIG_LABEL[model], "n": len(variants), "truncated": truncated,
                  "both_only": bool(studies)},
        "variants": variants, "traits": traits,
        "note": "Only variants with both an Evo2 score and a stored GWAS P-value (P<=1e-4) are shown.",
    })


# --- pre-computed genome-wide Region-Set results (R-aligned, all 6 models) ---
_RSGW = None


def rsgw_data():
    global _RSGW
    if _RSGW is None:
        p = os.path.join(os.path.dirname(__file__), "..", "regionset_genomewide.json")
        try:
            _RSGW = json.load(open(p))
        except Exception:
            _RSGW = {"models": {}}
    return _RSGW


async def regionset_genomewide(request, datasette):
    model = request.args.get("model") or DEFAULT_CONFIG
    d = rsgw_data()
    md = d.get("models", {}).get(model)
    if not md:
        return Response.json({"error": f"no pre-computed genome-wide results for '{model}'"}, status=404)
    return Response.json({"n_total": d.get("n_total"), "n_pruned": d.get("n_pruned"),
                          "n_pruned_sds": d.get("n_pruned_sds"),
                          "method": d.get("method"), "age_field": d.get("age_field"),
                          "sds_field": d.get("sds_field"),
                          "model": model, "model_label": md.get("label"),
                          "models_available": list(d.get("models", {}).keys()), **md})


@hookimpl
def register_routes():
    return [
        (r"^/api/region/analyze$", region_analyze),
        (r"^/api/regionset/analyze$", regionset_analyze),
        (r"^/api/regionset/genomewide$", regionset_genomewide),
        (r"^/api/gwas/meta$", gwas_meta),
        (r"^/api/gwas/search$", gwas_search),
        (r"^/api/gwas/trait$", gwas_trait),
        (r"^/api/gwas/overview$", gwas_overview),
        (r"^/api/coloc/traits$", coloc_traits),
        (r"^/api/coloc/region$", coloc_region),
    ]
