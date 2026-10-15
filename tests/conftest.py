"""Every test starts from the same random state, whatever ran before it: torch, numpy and random are seeded from
the test's id, so tests that draw from the global generators are deterministic and independent of the order."""
import random
import zlib

import numpy as np
import pytest
import torch


@pytest.fixture(autouse=True)
def _seed_global_generators(request):
    seed = zlib.crc32(request.node.nodeid.encode())
    torch.manual_seed(seed)
    np.random.seed(seed % 2 ** 32)
    random.seed(seed)
    yield
