import numpy as np

from earsim import aad, earsaad, sources
from earsim.head import load_template


def test_ceegrid_runs_top_to_bottom_and_sits_on_scalp():
    pos = earsaad.ceegrid_positions()
    head = load_template()
    for names in (earsaad.CEEGRID_L, earsaad.CEEGRID_R):
        z = np.array([pos[n][2] for n in names])
        assert z[0] > z[len(z) // 2] > z[-1], "electrode 1 must be the top, last the bottom (manual)"
        for n in names:
            assert np.linalg.norm(head.scalp_rr - pos[n], axis=1).min() < 1e-9
    # left/right on the correct sides of the head
    assert all(pos[n][0] < 0 for n in earsaad.CEEGRID_L) and all(pos[n][0] > 0 for n in earsaad.CEEGRID_R)


def test_channel_inventory_matches_dataset():
    assert len(earsaad.SCALP) == 29
    assert len(earsaad.EAR_CHANNELS) == 31 and len(set(earsaad.EAR_CHANNELS)) == 31
    assert len(earsaad.CEEGRID_L) == 9 and len(earsaad.CEEGRID_R) == 10


def test_cortical_bases_are_psd_unit_trace_and_cover_parcels():
    pos = earsaad.scalp_positions()
    names, S = sources.cortical_bases(pos)
    assert len(names) >= 60 and S.shape[1:] == (29, 29)
    n = S.shape[1]
    P = np.eye(n) - 1.0 / n
    np.testing.assert_allclose([np.trace(P @ s @ P) for s in S], 1.0, atol=1e-8)
    # MNE forwards are single precision: allow ~1e-6 relative negative eigenvalues
    assert min(np.linalg.eigvalsh(s).min() for s in S) > -1e-6 * S.max()


def test_raw_and_normalised_bases_are_consistent():
    pos = earsaad.scalp_positions()
    _, raw = sources.cortical_raw(pos)
    _, unit = sources.cortical_bases(pos)
    np.testing.assert_allclose(raw / sources.avgref_trace(raw)[:, None, None], unit)


def test_references_are_reference_free_and_handle_nan():
    rng = np.random.default_rng(0)
    Y = rng.standard_normal((50, 10))
    groups = [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9]]    # a group needs >= 3 good channels to serve as reference
    common = rng.standard_normal((50, 1))
    for fn in (aad.local_reference, aad.cross_reference):
        np.testing.assert_allclose(fn(Y + common, groups), fn(Y, groups), atol=1e-10)  # common ref cancels
    Y2 = Y.copy(); Y2[:, 0] = np.nan
    out = aad.cross_reference(Y2, groups)
    assert np.isnan(out[:, 0]).all() and np.isfinite(out[:, 1:]).all()
