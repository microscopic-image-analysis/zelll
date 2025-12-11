[![Crates.io](https://img.shields.io/crates/v/zelll.svg)](https://crates.io/crates/zelll)
[![PyPI](https://img.shields.io/pypi/v/zelll.svg)](https://pypi.python.org/pypi/zelll)
[![Documentation](https://img.shields.io/badge/docs-Python_API-green)](https://microscopic-image-analysis.github.io/zelll)

# Python bindings for zelll

Bindings for [`zelll`](https://github.com/microscopic-image-analysis/zelll) with the aim to allow for idiomatic use in Python.
The latest release is documented [here](https://microscopic-image-analysis.github.io/zelll).
Pre-built wheels can be found on [PyPI](https://pypi.python.org/pypi/zelll).

## Usage

```python
from zelll import CellGrid
import numpy as np

# CellGrid accepts any iterable object and converts its elements if possible.
# Values that can't be interpreted as something like [<float>, <float>, <float>]
# will be silently omitted.
points = np.random.random_sample((10, 3))
cg = CellGrid(points, 0.5)

# rebuild() accepts an optional cutoff parameter
cg.rebuild(points, 1.0)

# CellGrid objects are iterable, so you can use them like any other Iterable in Python:
pairs = list(cg)
pairs = [(p, q) for p, q in cg]
pairs = []
for p, q in cg:
    pairs.append((p, q))

# Note that while CellGrid produces unique ordered index pairs, it visits its cells in arbitrary order.
# So if you want to check whether the index pairs after `rebuild()` changed,
# prefer `set(cg)` over `list(cg)` 

# you can keep a CellGridIter object:
it = iter(cg)
# however, CellGrid can't be mutated while there are iterators of it alive
# i.e. `cg.rebuild(...)` throws a RuntimeError as long as `it` is alive
# either use `del it` or just use iterators implicitly/in local scopes
# (see above for examples)
# Additionally, CellGridIter is not thread-safe 
# but CellGrid is and can be sent between threads instead.

# The index pairs produced by CellGridIter also contain pairs
# with distance > cutoff.
# Here's an example dropping pairs with distance > cutoff.
# Note that there are faster ways to compute the (squared) euclidean distance.
pairs = [((i, p), (j, q)) for (i, p), (j, q) in cg 
    if np.linalg.norm(np.array(p) - np.array(q)) <= 0.5]
```

### Case Study

`examples/psssh.py` illustrates how the bindings can be used for prototyping purposes
by replicating the core design implemented in 
[`surface-sampling/`](https://github.com/microscopic-image-analysis/zelll/tree/main/surface-sampling):

```sh
maturin develop --release
uv venv examples/.venv
source examples/.venv/bin/activate
uv pip install -r examples/requirements.txt
# download some protein structures to test
# e.g. from here:
# https://dockground.compbio.ku.edu/unbound/unbound-docking-benchmarks.php
python psssh.py <PDB> -o psssh.pdb
# you can visualize the output file using e.g. PyMol
```

## Building

The following steps assume a working Rust toolchain.

1. clone this repository and `cd ./zelll/python/`
2. install [`maturin`](https://www.maturin.rs/installation)
3. create and activate a virtual environment, eg. `python3 -m venv .venv && source .venv/bin/activate`
4. (optionally install `numpy` in your environment for testing purposes)
5. run `maturin develop --release` to build and install an optimized `.whl` into the current virtual environment

