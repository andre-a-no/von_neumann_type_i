"""torch_vn_algebra: finite-dimensional Type I von Neumann algebras in PyTorch."""
from .algebra import TypeIAlgebra
from .hilbert_space import HilbertSpace, tensor_product_hilbert, direct_sum_hilbert
from .channels import Channel, InterSectorChannel
from .states import DensityMatrix
from .composite import tensor_product, kron, partial_trace
from . import channels, composite, cost, dynamics, optimize, states

Operator = TypeIAlgebra.Operator

__all__ = ["TypeIAlgebra", "Operator", "HilbertSpace", "tensor_product_hilbert", "direct_sum_hilbert",
           "Channel", "InterSectorChannel", "DensityMatrix", "tensor_product", "kron", "partial_trace",
           "channels", "composite", "cost", "dynamics", "optimize", "states"]
__version__ = "0.5.0"
