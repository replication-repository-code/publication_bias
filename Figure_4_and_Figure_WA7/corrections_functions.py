"""Precision-based corrections and trim-and-fill, as in Table_WA10/1_corrections.R, without R packages.

- WLS with weights 1/SE^2 (UWLS, WAAP, PET, PEESE): CR2 cluster-robust standard errors with Satterthwaite degrees of freedom,
  following clubSandwich 0.7.0 for weighted lm fits (working variance = identity, the default when a lm fit has weights).
- Two-level random-effects model (REML) and trim-and-fill (L0 estimator, side chosen from the sign of the slope of a REML
  meta-regression on the standard error), following metafor 4.4.0's rma.uni and trimfill.
"""
import numpy as np
from scipy import optimize, stats


def _sym_power(a, p):
    val, vec = np.linalg.eigh((a + a.T) / 2)
    val_p = np.where(val > 1e-12, np.abs(val) ** p, 0.0)
    return (vec * val_p) @ vec.T


def wls_cr2(y, X, w, cl):
    """WLS fit with CR2 cluster-robust covariance and Satterthwaite df per coefficient (clubSandwich, target = identity)."""
    y = np.asarray(y, float); X = np.asarray(X, float); w = np.asarray(w, float); cl = np.asarray(cl)
    N, p = X.shape
    XtW = X.T * w
    M = np.linalg.inv(XtW @ X)
    beta = M @ (XtW @ y)
    e = y - X @ beta
    groups = [np.where(cl == g)[0] for g in dict.fromkeys(cl)]
    MUWTWUM = M @ ((X.T * w ** 2) @ X) @ M
    meat = np.zeros((p, p)); E = []
    for idx in groups:
        Xj = X[idx]; wj = w[idx]
        Hjj = Xj @ M @ (Xj.T * wj)
        G = np.eye(len(idx)) - Hjj - Hjj.T + Xj @ MUWTWUM @ Xj.T
        A = _sym_power(G, -0.5)
        Ej = M @ (Xj.T * wj) @ A                       # p x n_j
        E.append(Ej)
        u = Ej @ e[idx]; meat += np.outer(u, u)
    V = meat
    se = np.sqrt(np.diag(V))
    WXM = (X * w[:, None]) @ M                          # N x p
    df = np.empty(p)
    for c in range(p):
        Gm = np.zeros((N, len(groups)))
        for j, idx in enumerate(groups):
            pj = E[j][c]                                # n_j
            col = -WXM @ (X[idx].T @ pj)
            col[idx] += pj
            Gm[:, j] = col
        gram = Gm.T @ Gm
        df[c] = np.trace(gram) ** 2 / np.sum(gram ** 2)
    t = beta / se
    pval = 2 * stats.t.sf(np.abs(t), df)
    crit = stats.t.ppf(0.975, df)
    return dict(beta=beta, se=se, df=df, p=pval, lo=beta - crit * se, hi=beta + crit * se)


def row(fit, c=0, **extra):
    return dict(g=fit['beta'][c], se=fit['se'][c], lo=fit['lo'][c], hi=fit['hi'][c], df=fit['df'][c], p=fit['p'][c], **extra)


def precision_methods(y, s, cl):
    y = np.asarray(y, float); s = np.asarray(s, float); cl = np.asarray(cl).astype(str); w = 1 / s ** 2
    one = np.ones((len(y), 1))
    out = {}
    u = wls_cr2(y, one, w, cl); out['UWLS'] = row(u)
    ap = s <= abs(u['beta'][0]) / 2.8
    if ap.sum() >= 2 and len(set(cl[ap])) >= 2:
        wa = wls_cr2(y[ap], one[ap], w[ap], cl[ap])
        out['WAAP'] = row(wa, n_estimates=int(ap.sum()), n_trials=len(set(cl[ap])))
    else:
        out['WAAP'] = dict(g=np.nan, se=np.nan, lo=np.nan, hi=np.nan, df=np.nan, p=np.nan, n_estimates=int(ap.sum()),
                           n_trials=len(set(cl[ap])))
    pet = wls_cr2(y, np.column_stack([one, s]), w, cl); out['PET'] = row(pet, slope=pet['beta'][1], slope_p=pet['p'][1])
    pe = wls_cr2(y, np.column_stack([one, s ** 2]), w, cl); out['PEESE'] = row(pe)
    pos = pet['beta'][0] > 0 and pet['p'][0] / 2 < .05
    out['PET-PEESE'] = dict(out['PEESE'] if pos else out['PET'], rule='PEESE' if pos else 'PET')
    return out


# ------------------------------------------------------------------ random effects (REML) and trim-and-fill
def _reml_fit(y, v, X, tau2):
    wi = 1 / (v + tau2); XtW = X.T * wi
    cov = np.linalg.inv(XtW @ X); b = cov @ (XtW @ y); r = y - X @ b
    ll = -0.5 * (np.sum(np.log(v + tau2)) + np.linalg.slogdet(XtW @ X)[1] + np.sum(wi * r ** 2))
    return b, cov, ll


def rma_reml(y, v, X=None):
    y = np.asarray(y, float); v = np.asarray(v, float)
    X = np.ones((len(y), 1)) if X is None else np.asarray(X, float)
    f = lambda t2: -_reml_fit(y, v, X, t2)[2]
    hi = max(10 * np.var(y), 1.0)
    grid = np.concatenate([[0.0], np.geomspace(1e-8, hi, 200)])
    vals = [f(t) for t in grid]; i = int(np.argmin(vals))
    if i == 0:
        res = optimize.minimize_scalar(f, bounds=(0.0, grid[1]), method='bounded', options={'xatol': 1e-14})
        t2 = res.x if res.fun < vals[0] else 0.0
    else:
        a = grid[i - 1]; b = grid[min(i + 1, len(grid) - 1)]
        t2 = optimize.minimize_scalar(f, bounds=(a, b), method='bounded', options={'xatol': 1e-14}).x
    b, cov, ll = _reml_fit(y, v, X, t2)
    se = np.sqrt(np.diag(cov))
    return dict(beta=b, se=se, tau2=t2, lo=b - 1.959963984540054 * se, hi=b + 1.959963984540054 * se, ll=ll)


def _rank_first(x):
    order = np.argsort(x, kind='stable'); r = np.empty(len(x)); r[order] = np.arange(1, len(x) + 1); return r


def trimfill(y, v, maxiter=100):
    y0 = np.asarray(y, float); v0 = np.asarray(v, float)
    reg = rma_reml(y0, v0, np.column_stack([np.ones(len(y0)), np.sqrt(v0)]))
    side = 'right' if reg['beta'][1] < 0 else 'left'
    yi = -y0 if side == 'right' else y0.copy()
    idx = np.argsort(yi, kind='stable'); yi = yi[idx]; vi = v0[idx]; k = len(yi)
    k0_sav, k0, it = -1, 0, 0
    while abs(k0 - k0_sav) > 0:
        k0_sav = k0; it += 1
        if it > maxiter: raise RuntimeError('trim-and-fill did not converge')
        res = rma_reml(yi[:k - k0], vi[:k - k0]); beta = res['beta'][0]
        yc = yi - beta
        rs = np.sign(yc) * _rank_first(np.abs(yc))
        Sr = rs[rs > 0].sum()
        k0 = max(0, int(round((4 * Sr - k * (k + 1)) / (2 * k - 1))))
    if k0 > 0:
        yc = -1 * (yc - beta) if side == 'right' else yc - beta
        yfill = np.concatenate([y0, -1 * yc[k - k0:]]); vfill = np.concatenate([v0, vi[k - k0:]])
        fit = rma_reml(yfill, vfill)
    else:
        fit = rma_reml(y0, v0)
    before = rma_reml(y0, v0)
    return dict(g=fit['beta'][0], se=fit['se'][0], lo=fit['lo'][0], hi=fit['hi'][0], k0=k0, side=side, tau2=fit['tau2'],
                g_before=before['beta'][0], se_before=before['se'][0], lo_before=before['lo'][0], hi_before=before['hi'][0],
                tau2_before=before['tau2'])
