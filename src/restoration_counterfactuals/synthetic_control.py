"""Univariate port of R/utils_scm_fit.R (quadprog, raw outcome units)."""
import numpy as np
import pandas as pd
import quadprog


def prepare_timeseries(treated, donors, intervention_year, min_pre_years=3):
    """Align one indicator; retain missing years explicitly, never impute outcomes.

    treated: year/value table; donors: donor_id/year/value table.
    Pre-fit years must be finite for the treatment and every supplied donor.
    """
    if treated.duplicated('year').any() or donors.duplicated(['year', 'donor_id']).any():
        raise ValueError('Expected one indicator and one observation per unit/year.')
    if donors.empty or donors.donor_id.isna().any():
        raise ValueError('No eligible donors or missing donor IDs.')
    y = treated.set_index('year').value.astype(float)
    x = donors.pivot(index='year', columns='donor_id', values='value').astype(float)
    years = y.index.union(x.index).sort_values()
    y, x = y.reindex(years), x.reindex(years)
    complete = np.isfinite(y) & np.isfinite(x).all(axis=1)
    pre = complete & (years < intervention_year)
    if pre.sum() < min_pre_years:
        raise ValueError(f'Insufficient complete pre-intervention years: {pre.sum()} < {min_pre_years}.')
    return y, x, pre


def solve_scm_ridge(y, x, penalty=1e-4, ridge_pd=1e-8):
    """Minimize squared error + penalty * sum(w²), w >= 0, sum(w) = 1."""
    y, x = np.asarray(y, dtype=float), np.asarray(x, dtype=float)
    if x.ndim != 2 or x.shape[0] != y.size or x.shape[1] == 0:
        raise ValueError('Expected an observations × donors matrix aligned with y.')
    if not np.isfinite(x).all() or not np.isfinite(y).all() or penalty < 0:
        raise ValueError('Fit inputs must be finite and penalty non-negative.')
    n = x.shape[1]
    dmat = 2 * (x.T @ x + penalty * np.eye(n)) + ridge_pd * np.eye(n)
    w = quadprog.solve_qp(dmat, 2 * x.T @ y,
                         np.column_stack([np.ones(n), np.eye(n)]),
                         np.r_[1., np.zeros(n)], meq=1)[0]
    if not np.isfinite(w).all() or w.min() < -1e-7 or abs(w.sum() - 1) > 1e-7:
        raise RuntimeError('Ridge solver returned infeasible weights.')
    w[w < 1e-10] = 0
    w /= w.sum()
    return dict(weights=w, penalty=float(penalty),
                rmspe=float(np.sqrt(np.mean((y - x @ w) ** 2))),
                effective_n=float(1 / (w @ w)), max_weight=float(w.max()))


def fit_synthetic_control(treated, donors, intervention_year, lambda_grid=None):
    """Return trajectories, weights, diagnostics and the full penalty search table."""
    y, x, pre = prepare_timeseries(treated, donors, intervention_year)
    grid = np.logspace(-6, 3, 30) if lambda_grid is None else np.asarray(lambda_grid)
    if len(grid) == 0 or not np.isfinite(grid).all() or (grid < 0).any():
        raise ValueError('Provide a nonempty finite non-negative penalty grid.')
    fits = [solve_scm_ridge(y[pre], x.loc[pre], p) for p in grid]
    search = pd.DataFrame([{k: v for k, v in fit.items() if k != 'weights'} for fit in fits])
    best = search.rmspe.min()
    # Multiplication also handles a zero best RMSPE without a 0/0 ratio.
    selected = search.loc[search.rmspe <= best * 1.02].penalty.max()
    fit = solve_scm_ridge(y[pre], x.loc[pre], selected)
    synthetic = x.to_numpy() @ fit['weights']  # NaNs propagate; no partial weighted sums.
    trajectories = pd.DataFrame({'year': y.index, 'observed': y.to_numpy(),
                                 'synthetic': synthetic, 'effect': y.to_numpy() - synthetic,
                                 'used_for_fit': pre.to_numpy()})
    trajectories['period'] = np.where(trajectories.year < intervention_year, 'pre', 'post')
    weights = pd.DataFrame({'donor_id': x.columns, 'weight': fit['weights']})
    diagnostics = pd.DataFrame([{k: v for k, v in fit.items() if k != 'weights'}])
    diagnostics['n_pre_years'] = int(pre.sum())
    diagnostics['n_donors'] = x.shape[1]
    diagnostics['n_incomplete_years'] = int((~(np.isfinite(y) & np.isfinite(x).all(axis=1))).sum())
    return dict(trajectories=trajectories, weights=weights, diagnostics=diagnostics, penalty_search=search)
