"""torch_vn_algebra: finite-dimensional Type I von Neumann algebras in PyTorch."""
from .algebra import TypeIAlgebra
from .hilbert_space import HilbertSpace, tensor_product_hilbert, direct_sum_hilbert
from .channels import Channel, InterSectorChannel, SectorJump
from .states import DensityMatrix
from .chains import SpinChain
from .composite import tensor_product, fused_tensor_product, kron, partial_trace
from . import chains, channels, composite, krylov, cost, dynamics, optimize, states, functional, inclusions, decompose, su2, divergences, jordan, interop
from .inclusions import Inclusion, ConditionalExpectation
from .decompose import Decomposition, RealStructure, wedderburn, commutant_basis
from .su2 import SU2Chain
from .jordan import SpinFactor, jordan_algebra_generated
from .functional import caching

Operator = TypeIAlgebra.Operator

__all__ = ["TypeIAlgebra", "Operator", "HilbertSpace", "tensor_product_hilbert", "direct_sum_hilbert",
           "Channel", "InterSectorChannel", "SectorJump", "DensityMatrix", "tensor_product", "fused_tensor_product", "kron", "partial_trace",
           "SpinChain", "SU2Chain", "su2", "Inclusion", "ConditionalExpectation", "inclusions", "Decomposition", "RealStructure", "wedderburn", "commutant_basis", "decompose", "divergences", "SpinFactor", "jordan_algebra_generated", "jordan", "chains", "channels", "composite", "krylov", "cost", "dynamics", "optimize", "states", "functional", "caching", "interop"]
__version__ = "1.1.0"
from . import permutational  # noqa: E402  (PermutationalSystem, PermutationalLindbladian)
