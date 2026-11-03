import numpy as np

from earsim import aad, placement

NAMES = ["a", "b", "c", "d"]


def make_subject(seed, n_trials=4, T=1500, informative=("a", "b")):
    """Channels carry noise; the bipolar pair (a-b) additionally carries the attended envelope."""
    rng = np.random.default_rng(seed)
    trials, envs, att = [], [], []
    for _ in range(n_trials):
        env = np.cumsum(rng.standard_normal((T, 2)), 0)
        env = env - np.convolve(env[:, 0], np.ones(20) / 20, "same")[:, None]
        env = np.stack([np.convolve(rng.standard_normal(T), np.ones(5) / 5, "same") for _ in range(2)], 1)
        x = rng.standard_normal((T, len(NAMES)))
        x[:, NAMES.index(informative[0])] += 2.0 * env[:, 0]
        trials.append(x)
        envs.append(env)
        att.append(0)
    return trials, envs, att


def test_lagged_shape_and_content():
    X = np.arange(12.0).reshape(6, 2)
    F = aad.lagged(X)
    assert F.shape == (6, 2 * aad.LAGS)
    np.testing.assert_array_equal(F[0, :2], X[0])
    np.testing.assert_array_equal(F[0, 2:4], X[1])      # lag 1 = next sample
    assert np.all(F[-1, 2:4] == 0)                      # zero padded at the end


def test_pair_columns_counts():
    pairs, cols = placement.pair_columns(NAMES, NAMES)
    assert len(pairs) == 6 and cols.shape == (6, 2)


def test_gram_score_equals_direct_ridge():
    trials, envs, att = make_subject(0)
    pairs, cols = placement.pair_columns(NAMES, NAMES)
    st = placement.subject_stats(trials, envs, att, cols)
    subset = (0, 3)
    got = placement.loto_score(st, subset, c=0.1)
    feats = [aad.lagged(np.nan_to_num(placement.pair_features(t, cols)[:, list(subset)] / st.scale[list(subset)]))
             for t in trials]
    ys = [aad.zs(e[:, a]) for e, a in zip(envs, att)]
    exp = []
    for k in range(len(trials)):
        tr = [i for i in range(len(trials)) if i != k]
        w = aad.fit_ridge(np.vstack([feats[i] for i in tr]), np.concatenate([ys[i] for i in tr]), 0.1)
        pred = feats[k] @ w
        exp.append(pred @ ys[k] / (np.linalg.norm(pred) * np.linalg.norm(ys[k])))  # uncentred, as scored
    np.testing.assert_allclose(got, np.mean(exp), atol=1e-6)


def test_greedy_picks_informative_pair_and_is_leave_one_out():
    pairs, cols = placement.pair_columns(NAMES, NAMES)
    stats = {s: placement.subject_stats(*make_subject(s), cols) for s in range(4)}
    score = placement.ScoreCache(stats)
    chosen = placement.greedy_select(score, len(pairs), 1, exclude=0)
    # the informative channel 'a' is in every pair that has it; first pick must involve 'a'
    assert "a" in pairs[chosen[0]]
    # excluding a subject must change the training mean (no leakage of the excluded subject)
    assert score.mean_excluding((chosen[0],), exclude=0) != score.mean_excluding((chosen[0],), exclude=1)


def test_window_accuracy_perfect_when_signal_is_clean():
    rng = np.random.default_rng(1)
    T = 2000
    env = np.stack([np.convolve(rng.standard_normal(T), np.ones(5) / 5, "same") for _ in range(2)], 1)
    recon = env[:, 0] + 0.01 * rng.standard_normal(T)
    ok, tot = aad.window_accuracy(recon, env, 0, 10)
    assert ok == tot and tot == T // (10 * aad.FS)


def test_dead_channel_pair_scores_zero_and_is_never_selected():
    pairs, cols = placement.pair_columns(NAMES, NAMES)
    subjects = {}
    for s in range(3):
        trials, envs, att = make_subject(s)
        for t in trials:
            t[:, NAMES.index("d")] = np.nan                  # channel d is dead for everyone
        subjects[s] = placement.subject_stats(trials, envs, att, cols)
    score = placement.ScoreCache(subjects)
    dead = [j for j, p in enumerate(pairs) if "d" in p]
    assert all(score((j,))[0] == 0.0 for j in dead)
    chosen = placement.greedy_select(score, len(pairs), 2, exclude=0)
    assert not any(j in dead for j in chosen)
