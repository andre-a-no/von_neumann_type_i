"""torch_vn_algebra: finite-dimensional Type I von Neumann algebras in PyTorch."""
from .algebra import TypeIAlgebra
from .hilbert_space import HilbertSpace, tensor_product_hilbert, direct_sum_hilbert
from .channels import Channel
from .states import DensityMatrix
from . import channels, dynamics, states

Operator = TypeIAlgebra.Operator

__all__ = ["TypeIAlgebra", "Operator", "HilbertSpace", "tensor_product_hilbert", "direct_sum_hilbert",
           "Channel", "DensityMatrix", "channels", "dynamics", "states"]
__version__ = "0.3.0"
