import numpy as np
import pytest

from earsim import layouts
from earsim.baselines import spherical_spline
from earsim.head import load_template
from earsim.transfer import bipolar, build_transfer

SCALP = ["Fp1", "Fp2", "F7", "F3", "Fz", "F4", "F8", "T7", "C3", "Cz", "C4", "T8",
         "P7", "P3", "Pz", "P4", "P8", "O1", "Oz", "O2", "FT9", "FT10", "TP9", "TP10"]


@pytest.fixture(scope="module")
def tr():
    src = layouts.standard(SCALP)
    tgt = layouts.in_ear_like("L") | layouts.in_ear_like("R")
    return build_transfer(src, tgt, method="MNE", snr=3.0)


def test_shape_and_rank(tr):
    assert tr.T.shape == (4, len(SCALP))
    assert np.linalg.matrix_rank(tr.T) <= len(SCALP)


def test_linear(tr):
    rng = np.random.default_rng(0)
    a, b = rng.standard_normal((2, len(SCALP), 50))
    np.testing.assert_allclose(tr.apply(2 * a + b), 2 * tr.apply(a) + tr.apply(b), atol=1e-10)


def test_reference_invariant(tr):
    """Adding a common signal to all scalp channels (a reference change) is a no-op."""
    rng = np.random.default_rng(1)
    x = rng.standard_normal((len(SCALP), 50))
    common = rng.standard_normal((1, 50))
    np.testing.assert_allclose(tr.apply(x + common), tr.apply(x), atol=1e-8)


def test_bipolar_rows_sum_zero():
    R = bipolar(["a", "b", "c"], [("a", "b"), ("c", "a")])
    np.testing.assert_allclose(R.sum(1), 0)


def test_layouts_on_scalp():
    head = load_template()
    for lay in layouts.wearable_layouts().values():
        P = np.array(list(lay.values()))
        d = np.linalg.norm(P[:, None] - head.scalp_rr[None], axis=2).min(1)
        # snapped layouts sit exactly on vertices; 10-05 positions within a few mm
        assert d.max() < 0.02


def test_spline_reproduces_source_channels():
    src = layouts.standard(SCALP)
    W = spherical_spline(src, src, alpha=None)
    np.testing.assert_allclose(W, np.eye(len(SCALP)), atol=1e-6)
