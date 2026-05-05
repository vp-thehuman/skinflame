import numpy as np

from skinflame.utils import bca_ci, fdr_bh, ols, percentile_ci


def test_ols_recovers_known_coefficients():
    rng = np.random.default_rng(0)
    n = 500
    X = rng.normal(size=(n, 3))
    Xd = np.column_stack([np.ones(n), X])
    true = np.array([1.0, 2.0, -1.0, 0.5])
    y = Xd @ true + rng.normal(size=n) * 0.5
    beta, cov, sigma2 = ols(y, Xd)
    np.testing.assert_allclose(beta, true, atol=0.1)
    assert sigma2 == np.float64(sigma2)
    assert cov.shape == (4, 4)


def test_percentile_and_bca_ci_contain_mean():
    rng = np.random.default_rng(0)
    samples = rng.normal(loc=2.0, scale=1.0, size=2000)
    lo, hi = percentile_ci(samples)
    assert lo < 2.0 < hi
    lo, hi = bca_ci(samples, theta_hat=2.0)
    assert lo < 2.0 < hi


def test_fdr_bh_identity_when_all_significant():
    p = np.array([0.001, 0.002, 0.003])
    rejected, q = fdr_bh(p, q=0.05)
    assert rejected.all()
    assert (q <= 0.05).all()


def test_fdr_bh_controls_false_discoveries():
    # With 10 strong signals (p=1e-6) among 90 nulls (p=0.5):
    # BH adjusted p for the 10 signals is 1e-6 * 100 / k <= 0.05 for k >= 1, so all reject.
    p = np.concatenate([np.full(10, 1e-6), np.full(90, 0.5)])
    rejected, q = fdr_bh(p, q=0.05)
    assert rejected[:10].all()
    assert rejected[10:].sum() == 0
