---
name: torch-vn-algebra
description: Write and check code that uses the Python library torch_vn_algebra (block-diagonal operators in symmetry sectors, finite-dimensional von Neumann algebras ⊕ M_k(C) in PyTorch). Covers quantum states, channels, Lindblad dynamics, spin and anyon chains, Wedderburn decompositions, inclusions and Jones indices, measurements, and channels as torch.nn layers. Use when the user's code imports torch_vn_algebra, or when the task is batched computation with block-diagonal matrices, superselection or charge sectors, or quantum channels between sectors on CPU or GPU.
---

# torch_vn_algebra

1. Make sure the library is importable:

   ```bash
   python -c "import torch_vn_algebra as t; print(t.__version__)"
   ```

   If it is not, run `pip install torch_vn_algebra` (PyTorch must be installed first).
2. Read the guide before writing code:

   ```bash
   python -m torch_vn_algebra.guide
   ```

   It explains:
   - the layout `(batch, C, k_max, k_max)` with zero padding;
   - the three traces;
   - which attributes are properties and which are methods;
   - the usual mistakes;
   - a table from tasks to functions;
   - tested recipes.
3. Look up every function you use in the API index of the installed version. Do not guess names or arguments:

   ```bash
   python -m torch_vn_algebra.guide --api <module>
   ```

   For example, `--api states`. Without a module the command lists everything.
4. Run the code and check an invariant of the result, for example:
   - `rho.Tr_blunt()` = 1;
   - `channel.is_trace_preserving()`;
   - conserved `sector_probabilities()`;
   - `Decomposition.verify()`.

   Compare with a dense `torch.linalg` computation on a small case when in doubt.
5. Use `precision='double'` unless the quantity is a Monte Carlo statistic. Batch the samples instead of looping in
   Python.

Install this skill by copying this directory to `~/.claude/skills/` (personal) or `.claude/skills/` (one project).
