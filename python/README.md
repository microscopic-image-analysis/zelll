# `python-zelll`

proof-of-concept bindings for `zelll` with the aim to allow for idiomatic use in Python.

## Try it Yourself

The following steps assume a working Rust toolchain.

1. clone this repository and `cd ./zelll/python/`
2. install [`maturin`](https://www.maturin.rs/installation)
3. create and activate a virtual environment, eg. `python3 -m venv .venv && source .venv/bin/activate`
4. (optionally install `numpy` in your environment for testing purposes)
5. run `maturin develop --release` to build and install an optimized `.whl` into the current virtual environment
6. open a Python REPL and start playing around (currently w/o any documentation):
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
pairs = [(i, j) for i, j in cg]
pairs = []
for i, j in cg:
    pairs.append((i, j))

# Note that while CellGrid produces unique ordered index pairs, it visits its cells in arbitrary order.
# So if you want to check whether the index pairs after `rebuild()` changed,
# prefer `set(cg)` over `list(cg)` 

# you can keep a CellGridIterator object:
it = iter(cg)
# however, CellGrid can't be mutated while there are iterators of it alive
# i.e. `cg.rebuild(...)` throws a RuntimeError as long as `it` is alive
# either use `del it` or just use iterators implicitly/in local scopes
# (see above for examples)
# Additionally, CellGridIterator is not thread-safe 
# but CellGrid is and can be sent between threads instead.

# The index pairs produced by CellGridIterator also contain pairs
# with distance > cutoff.
# Here's an example dropping the index pairs with distance > cutoff.
# Note that there are faster ways to compute the (squared) euclidean distance.
pairs = [(i,j) for i, j in cg if np.linalg.norm(points[i] - points[j]) <= 0.5]
```

## TODO

- [ ] measure bindings performance overhead
- [ ] fully convince myself that usage of std::mem::transmute() is sound
- [ ] complete API
- [ ] documentation
    * [ ] integrate with readthedocs.io
- [ ] distribution
    * [ ] build `*.whl`s on CI
    * [ ] publish to PyPI on release

