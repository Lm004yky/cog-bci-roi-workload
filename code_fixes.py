# =============================================================================
# IJP REVISION — corrected / new code blocks
# Paste block-by-block into EEG_v2_cleanup.ipynb, replacing the originals.
# Covers: Fig 2A fix + checklist #3 (ICC CI), #5 (effect sizes),
#         #6 (multiple-comparison correction), #7 (mixed-effects sensitivity).
#
# Requires (already in your setup): numpy, pandas, scipy.stats, matplotlib, seaborn
# New deps:
#   from statsmodels.stats.multitest import multipletests   # ships with statsmodels
#   import statsmodels.formula.api as smf                    # ships with statsmodels
#   (optional cross-check) pip install pingouin
# =============================================================================


# =============================================================================
# BLOCK A — FIG 2 (replaces PART 1 plotting).
# Change vs original: panel A no longer draws the pairwise "**" bracket
# (which implied an unreported 0-vs-2 pairwise test). It now shows the omnibus
# one-way ANOVA as a text annotation, matching panel B. The pairwise p-values
# are still COMPUTED and printed so you can decide whether to report them.
# =============================================================================

order = ['Low', 'Medium', 'High']
palette = {'Low': '#4A90D9', 'Medium': '#50B86C', 'High': '#E05D5D'}
FIG_DPI = 300

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# ---- A: N-Back ----
nback = df[df['task'] == 'N-Back']
sns.boxplot(data=nback, x='filename', y='roi_ab_ratio',
            order=['zeroBACK', 'oneBACK', 'twoBACK'],
            palette=['#4A90D9', '#50B86C', '#E05D5D'], ax=axes[0],
            fliersize=2, linewidth=1)
axes[0].set_title('(A) N-Back: ROI \u03b1/\u03b2 by difficulty', fontsize=12, fontweight='bold')
axes[0].set_xticklabels(['0-back\n(Low)', '1-back\n(Medium)', '2-back\n(High)'])
axes[0].set_ylabel('ROI \u03b1/\u03b2 Ratio', fontsize=11)
axes[0].set_xlabel('')

nback_groups = [nback[nback['filename'] == f]['roi_ab_ratio']
                for f in ['zeroBACK', 'oneBACK', 'twoBACK']]

# Omnibus ANOVA -> text annotation (same style as MATB panel B)
f_nback, p_nback = stats.f_oneway(*nback_groups)
axes[0].text(0.95, 0.95, f'ANOVA: p={p_nback:.3f}',
             transform=axes[0].transAxes, ha='right', va='top', fontsize=9,
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# Pairwise tests still computed + printed (NOT drawn as a bracket).
# Use these if you decide to report 0-vs-2 explicitly; check it against the
# Bonferroni/Holm threshold from BLOCK B before claiming significance.
_, p01 = stats.ttest_ind(nback_groups[0], nback_groups[1])
_, p12 = stats.ttest_ind(nback_groups[1], nback_groups[2])
t02, p02 = stats.ttest_ind(nback_groups[0], nback_groups[2])
n0, n2 = len(nback_groups[0]), len(nback_groups[2])
pooled02 = np.sqrt(((n0 - 1) * nback_groups[0].var() + (n2 - 1) * nback_groups[2].var())
                   / (n0 + n2 - 2))
d02 = (nback_groups[0].mean() - nback_groups[2].mean()) / pooled02
print(f'N-Back omnibus: F={f_nback:.3f}, p={p_nback:.4f}')
print(f'N-Back pairwise (raw, uncorrected):')
print(f'   0 vs 1: p={p01:.4f}')
print(f'   1 vs 2: p={p12:.4f}')
print(f'   0 vs 2: t={t02:.3f}, p={p02:.4f}, d={d02:+.3f}, n0={n0}, n2={n2}')

# ---- B: MATB (unchanged; this is the target pattern) ----
matb = df[df['task'] == 'MATB']
sns.boxplot(data=matb, x='filename', y='roi_ab_ratio',
            order=['MATBeasy', 'MATBmed', 'MATBdiff'],
            palette=['#4A90D9', '#50B86C', '#E05D5D'], ax=axes[1],
            fliersize=2, linewidth=1)
axes[1].set_title('(B) MATB: ROI \u03b1/\u03b2 by difficulty', fontsize=12, fontweight='bold')
axes[1].set_xticklabels(['Easy\n(Low)', 'Medium', 'Difficult\n(High)'])
axes[1].set_ylabel('ROI \u03b1/\u03b2 Ratio', fontsize=11)
axes[1].set_xlabel('')

matb_groups = [matb[matb['filename'] == f]['roi_ab_ratio']
               for f in ['MATBeasy', 'MATBmed', 'MATBdiff']]
_, p_matb = stats.f_oneway(*matb_groups)
axes[1].text(0.95, 0.95, f'ANOVA: p={p_matb:.3f} (ns)',
             transform=axes[1].transAxes, ha='right', va='top', fontsize=9,
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.suptitle(f'Task-Specific Workload Analysis (N={n_subjects})',
             fontsize=14, fontweight='bold', y=1.03)
plt.tight_layout()
plt.savefig('fig2_nback_matb.png', dpi=FIG_DPI, bbox_inches='tight')
plt.savefig('fig2_nback_matb.pdf', dpi=FIG_DPI, bbox_inches='tight')
plt.show()


# =============================================================================
# BLOCK B — full_stats v2  (replaces PART 2a full_stats).
# Adds checklist #5 (partial eta^2, omega^2) and #6 (Holm + FDR correction of
# the 3 pairwise post-hocs within each metric). Prints raw AND adjusted p.
# For a one-way design partial eta^2 == eta^2 (single factor) — reported as such.
# =============================================================================

from statsmodels.stats.multitest import multipletests


def full_stats(data, metric, label, correct='holm'):
    print(f'\n{"=" * 70}')
    print(label)
    print(f'{"=" * 70}')

    groups = {wl: data[data['workload'] == wl][metric].dropna()
              for wl in ['Low', 'Medium', 'High']}
    vals = list(groups.values())

    for wl in ['Low', 'Medium', 'High']:
        g = groups[wl]
        print(f'  {wl:8s}: {g.mean():.4g} \u00b1 {g.std():.4g} '
              f'(median={g.median():.4g}, n={len(g)})')

    # --- omnibus ANOVA ---
    f_stat, p_val = stats.f_oneway(*vals)

    # --- effect sizes ---
    allv = np.concatenate([g.values for g in vals])
    N, k = len(allv), len(vals)
    gm = allv.mean()
    ss_b = sum(len(g) * (g.mean() - gm) ** 2 for g in vals)
    ss_t = ((allv - gm) ** 2).sum()
    ss_w = ss_t - ss_b
    df_b, df_w = k - 1, N - k
    ms_w = ss_w / df_w
    eta2 = ss_b / ss_t                                   # == partial eta^2 (one-way)
    omega2 = (ss_b - df_b * ms_w) / (ss_t + ms_w)        # less biased

    omni_mark = '\u2713' if p_val < 0.05 else '\u2717'
    print(f'\n  ANOVA: F({df_b},{df_w})={f_stat:.3f}, p={p_val:.2e} {omni_mark}')
    print(f'  partial \u03b7\u00b2 = {eta2:.4f}   (\u03b7\u00b2 == partial \u03b7\u00b2 for one-way), '
          f'\u03c9\u00b2 = {omega2:.4f}')

    # --- pairwise post-hocs + multiple-comparison correction ---
    pairs = [('Low', 'Medium'), ('Medium', 'High'), ('Low', 'High')]
    raw_p, tvals, dvals = [], [], []
    for a, b in pairs:
        t, p = stats.ttest_ind(groups[a], groups[b])
        na, nb = len(groups[a]), len(groups[b])
        pooled = np.sqrt(((na - 1) * groups[a].var() + (nb - 1) * groups[b].var())
                         / (na + nb - 2))
        d = (groups[a].mean() - groups[b].mean()) / pooled if pooled > 0 else 0.0
        raw_p.append(p); tvals.append(t); dvals.append(d)

    rej_h, p_holm, _, _ = multipletests(raw_p, alpha=0.05, method=correct)
    rej_f, p_fdr, _, _ = multipletests(raw_p, alpha=0.05, method='fdr_bh')

    print(f'\n  Post-hoc (3 comparisons, corrected within metric):')
    print(f'  {"comparison":16s} {"t":>8s} {"p_raw":>10s} '
          f'{"p_"+correct:>10s} {"p_FDR":>10s} {"d":>8s}')
    for i, (a, b) in enumerate(pairs):
        size = 'small' if abs(dvals[i]) < 0.5 else 'medium' if abs(dvals[i]) < 0.8 else 'large'
        star = '\u2713' if p_holm[i] < 0.05 else '\u2717'
        print(f'  {a+" vs "+b:16s} {tvals[i]:8.3f} {raw_p[i]:10.2e} '
              f'{p_holm[i]:10.2e} {p_fdr[i]:10.2e} {dvals[i]:+8.3f} {size} {star}')

    return {'F': f_stat, 'df_between': df_b, 'df_within': df_w, 'p': p_val,
            'partial_eta2': eta2, 'omega2': omega2,
            'pairs': pairs, 'p_raw': raw_p,
            f'p_{correct}': list(p_holm), 'p_fdr': list(p_fdr),
            't': tvals, 'd': dvals}


full_stats(df_cog, 'roi_ab_ratio',    'ROI \u03b1/\u03b2 RATIO (parietal \u03b1 / frontal \u03b2)')
full_stats(df_cog, 'global_ab_ratio', 'GLOBAL \u03b1/\u03b2 RATIO (old method)')
full_stats(df_cog, 'theta_beta_ratio','\u03b8/\u03b2 RATIO (frontal \u03b8 / frontal \u03b2)')
full_stats(df_cog, 'frontal_theta',   'FRONTAL MIDLINE THETA POWER')

# Reporting note for the paper (Table 1 / Abstract):
#   report as F(df_b, df_w) = ..., p = ..., partial eta^2 = ...
#   e.g. "ROI alpha/beta: F(2, N) = 37.7, p < 0.001, partial eta^2 = <value>"
#   For the 2.8x-gain claim, cite partial eta^2 for BOTH ROI and global side by side.


# =============================================================================
# BLOCK C — icc_31 with 95% CI  (replaces PART 2c icc_31), checklist #3.
# Adds F-test based CI (McGraw & Wong 1996, consistency, single measures) and
# a complete-case guard. Prints paper-ready string. Optional pingouin cross-check.
# =============================================================================

def icc_31(df_nri_sessions, alpha=0.05):
    """ICC(3,1): two-way mixed, consistency, single measurement (Shrout&Fleiss).
    Complete-case only: subjects lacking any session are dropped.
    Returns point estimate + (1-alpha) CI via the F-distribution."""
    pivot = df_nri_sessions.pivot(index='subject', columns='session',
                                  values='nri').dropna()
    arr = pivot.values                      # (n_subjects, k_sessions)
    n, k = arr.shape
    if n < 2 or k < 2:
        raise ValueError(f'Need >=2 subjects and >=2 sessions; got n={n}, k={k}')

    grand_mean = arr.mean()
    ss_between_subj = k * np.sum((arr.mean(axis=1) - grand_mean) ** 2)
    ss_between_ses  = n * np.sum((arr.mean(axis=0) - grand_mean) ** 2)
    ss_total        = np.sum((arr - grand_mean) ** 2)
    ss_error        = ss_total - ss_between_subj - ss_between_ses

    df_subj  = n - 1
    df_error = df_subj * (k - 1)
    ms_subj  = ss_between_subj / df_subj
    ms_error = ss_error / df_error

    icc31 = (ms_subj - ms_error) / (ms_subj + (k - 1) * ms_error)
    icc3k = (ms_subj - ms_error) / ms_subj

    # F-test CI for consistency single-measures ICC
    F_obs = ms_subj / ms_error
    F_L = F_obs / stats.f.ppf(1 - alpha / 2, df_subj, df_error)
    F_U = F_obs * stats.f.ppf(1 - alpha / 2, df_error, df_subj)
    ci_low  = (F_L - 1) / (F_L + k - 1)
    ci_high = (F_U - 1) / (F_U + k - 1)

    return {'icc_31': icc31, 'icc_3k': icc3k,
            'ci_low': ci_low, 'ci_high': ci_high,
            'n_subjects': n, 'k_sessions': k,
            'ms_subj': ms_subj, 'ms_error': ms_error}


icc_result = icc_31(df_nri[['subject', 'session', 'nri']])
print(f"ICC(3,1) = {icc_result['icc_31']:.3f}, "
      f"95% CI [{icc_result['ci_low']:.3f}, {icc_result['ci_high']:.3f}]  "
      f"(n={icc_result['n_subjects']} complete-case subjects, "
      f"k={icc_result['k_sessions']} sessions)")
print(f"ICC(3,k) = {icc_result['icc_3k']:.3f}")
# -> paste the printed CI into the paper:
#    "ICC(3,1) = 0.075, 95% CI [X, Y] (classification: poor; Koo & Li 2016)"
#    A wide CI spanning ~0 reinforces the state-level (not trait) interpretation.

# Optional cross-check (pip install pingouin):
# import pingouin as pg
# icc_pg = pg.intraclass_corr(data=df_nri, targets='subject', raters='session',
#                             ratings='nri').set_index('Type')
# print(icc_pg.loc['ICC3', ['ICC', 'CI95%']])   # should match within rounding


# =============================================================================
# BLOCK D — non-independence sensitivity check (NEW), checklist #7.
# 4,736 epochs from 20 subjects violate the independence assumption of a plain
# epoch-level ANOVA. Each workload effect is re-tested THREE ways:
#   (1) SUBJECT-LEVEL repeated-measures ANOVA on epoch means (N=20)  <- headline
#       the textbook fix: collapse epochs to one value per subject x condition,
#       then a within-subject ANOVA. Always converges, most conservative.
#   (2) OLS with by-subject CLUSTER-ROBUST standard errors (keeps epoch-level)
#   (3) linear mixed model, random intercept per subject  [optional confirmation;
#       MixedLM's convergence flag is conservative, so report it honestly]
# (1) and (2) always work and each independently answers the reviewer. y is
# z-scored inside for numerical stability; this does not affect the F/Wald tests.
# =============================================================================

import statsmodels.formula.api as smf
from statsmodels.stats.anova import AnovaRM


def _workload_term_idx(names):
    return [i for i, n in enumerate(names) if n.startswith('C(workload)')]


def workload_sensitivity(data, metric):
    d = data[['subject', 'workload', metric]].dropna().copy()
    d = d.rename(columns={metric: 'y'})
    d['y'] = (d['y'] - d['y'].mean()) / d['y'].std()
    d['workload'] = pd.Categorical(d['workload'],
                                   categories=['Low', 'Medium', 'High'])

    # ---- (1) subject-level repeated-measures ANOVA on epoch means ----
    # keep only subjects present in all three levels (AnovaRM needs balance)
    agg = d.groupby(['subject', 'workload'], observed=True)['y'].mean().reset_index()
    complete = agg.groupby('subject')['workload'].nunique()
    keep = complete[complete == 3].index
    agg = agg[agg['subject'].isin(keep)]
    rm = AnovaRM(agg, depvar='y', subject='subject', within=['workload']).fit()
    t = rm.anova_table
    rm_F = float(t['F Value'].iloc[0])
    rm_p = float(t['Pr > F'].iloc[0])
    rm_dfn = float(t['Num DF'].iloc[0])
    rm_dfd = float(t['Den DF'].iloc[0])
    rm_peta = (rm_dfn * rm_F) / (rm_dfn * rm_F + rm_dfd)     # partial eta^2
    rm_nsub = agg['subject'].nunique()

    # ---- (2) cluster-robust OLS (epoch-level) ----
    ols = smf.ols('y ~ C(workload)', d).fit(
        cov_type='cluster', cov_kwds={'groups': d['subject']})
    names = ols.params.index.tolist()
    idx = _workload_term_idx(names)
    L = np.zeros((len(idx), len(names)))
    for r, i in enumerate(idx):
        L[r, i] = 1.0
    w = ols.wald_test(L, scalar=True)
    ols_chi2, ols_p = float(np.squeeze(w.statistic)), float(w.pvalue)

    # ---- (3) mixed model (optional confirmation, guarded) ----
    lmm_chi2 = lmm_p = np.nan
    lmm_ok = False
    for method in ('lbfgs', 'cg', 'powell', 'bfgs'):
        try:
            mf = smf.mixedlm('y ~ C(workload)', d, groups=d['subject']).fit(
                reml=False, method=method, maxiter=200)
            lidx = _workload_term_idx(mf.fe_params.index.tolist())
            Ll = np.zeros((len(lidx), len(mf.fe_params)))
            for r, i in enumerate(lidx):
                Ll[r, i] = 1.0
            stat = float(np.squeeze(mf.wald_test(Ll, scalar=True).statistic))
            if mf.converged and np.isfinite(mf.llf) and np.isfinite(stat):
                lmm_chi2 = stat
                lmm_p = float(mf.wald_test(Ll, scalar=True).pvalue)
                lmm_ok = True
                break
        except Exception:
            continue

    return {'metric': metric, 'n_epochs': int(len(d)),
            'rm_F': rm_F, 'rm_dfn': rm_dfn, 'rm_dfd': rm_dfd,
            'rm_p': rm_p, 'rm_partial_eta2': rm_peta, 'rm_n_subjects': rm_nsub,
            'ols_chi2': ols_chi2, 'ols_p': ols_p,
            'lmm_chi2': lmm_chi2, 'lmm_p': lmm_p, 'lmm_converged': lmm_ok}


print(f'\n{"=" * 92}')
print('NON-INDEPENDENCE SENSITIVITY CHECK  (workload effect re-tested per metric)')
print(f'{"=" * 92}')
print(f'  {"metric":18s} | {"RM-ANOVA F(subj)":>16s} {"p":>10s} {"pEta2":>7s} '
      f'| {"cluster chi2":>12s} {"p":>10s} | {"LMM":>10s}')
for metric in ['roi_ab_ratio', 'global_ab_ratio', 'theta_beta_ratio', 'frontal_theta']:
    r = workload_sensitivity(df_cog, metric)
    lmm = (f'{r["lmm_chi2"]:.1f}/{r["lmm_p"]:.1e}' if r['lmm_converged'] else 'n/c')
    print(f'  {metric:18s} | F({r["rm_dfn"]:.0f},{r["rm_dfd"]:.0f})={r["rm_F"]:7.2f} '
          f'{r["rm_p"]:10.2e} {r["rm_partial_eta2"]:7.3f} '
          f'| {r["ols_chi2"]:12.2f} {r["ols_p"]:10.2e} | {lmm:>10s}')

# Reporting note (Methods/Results, 1-2 sentences):
#   "To confirm the epoch-level results are not an artifact of treating epochs as
#    independent, each workload effect was re-tested at the subject level with a
#    repeated-measures ANOVA on per-subject condition means (N=20) and, at the
#    epoch level, with by-subject cluster-robust standard errors. The workload
#    factor remained significant for <metrics> (RM-ANOVA F(2,38)=..., p=...,
#    partial eta^2=...; cluster-robust Wald chi2(2)=..., p=...)."
# The subject-level partial eta^2 here is the honest between-SUBJECT effect size
# and will be much larger than the per-epoch Cohen's d (0.09-0.44) — mention both:
# the small per-epoch d reflects inter-subject noise, the subject-level effect is
# strong. 'n/c' = mixed model did not converge for that metric; report RM-ANOVA
# and cluster-robust results for it. Report any effect that weakens honestly.

# =============================================================================
# END
# =============================================================================
