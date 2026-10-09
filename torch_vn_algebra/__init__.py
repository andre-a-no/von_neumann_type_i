"""torch_vn_algebra: finite-dimensional Type I von Neumann algebras in PyTorch."""
from .algebra import TypeIAlgebra
from .hilbert_space import HilbertSpace, tensor_product_hilbert, direct_sum_hilbert

Operator = TypeIAlgebra.Operator

__all__ = ["TypeIAlgebra", "Operator", "HilbertSpace", "tensor_product_hilbert", "direct_sum_hilbert"]
__version__ = "0.2.0"
