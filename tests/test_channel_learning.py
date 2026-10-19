"""Tomography in the algebra as training (scripts/channel_learning.py, Section sec:learn-tomo of the paper): a tiny
run of the pipeline. The true process of the noisy XY register (n = 3), simulated counts, one model fitted by maximum
likelihood: the loss decreases and the learned channel is closer to the truth than the starting point."""
import importlib.util
import sys
from pathlib import Path

import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))                     # the script imports scripts/common.py
_spec = importlib.util.spec_from_file_location('channel_learning', ROOT / 'scripts' / 'channel_learning.py')
cl = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cl)


def test_problem_is_consistent():
    prob = cl.Problem(3)
    # the outputs computed in the algebra agree with an independent dense Lindblad evolution in M_8
    assert prob.dense_check() < 1e-12
    # the fast training paths agree with the layers' forward passes and with POVM.probabilities
    assert cl.check_fast_paths(prob) < 1e-12
    # the true Choi matrices are states, and the counts of the exact data are the probabilities
    for J in prob.true_choi:
        assert abs(torch.diagonal(J).sum().real.item() - 1) < 1e-12
        assert torch.linalg.eigvalsh(J).min() > -1e-10
    assert torch.allclose(prob.p.sum(-1), torch.ones(prob.n_train, dtype=torch.float64))


def test_maximum_likelihood_fit_learns_the_noise_process():
    torch.manual_seed(0)
    prob = cl.Problem(3)
    counts = prob.counts(10000, dense=False, seed=1)
    start = cl.evaluate('lindblad', cl.build('lindblad', prob, torch.Generator().manual_seed(4)), prob)
    net, loss, history, _ = cl.fit('lindblad', prob, counts, (60, 40), seed=4)       # same start (seed 4)
    end = cl.evaluate('lindblad', net, prob)
    assert loss < 0.1 * history[0]
    assert end['choi'] < 0.3 * start['choi']
    assert end['trace_distance'] < 0.3 * start['trace_distance']
    rates = cl.project_rates(prob, net)
    assert abs(rates['down'] - cl.RATES['down']) < 0.3 * cl.RATES['down']
