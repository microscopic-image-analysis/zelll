# `python-zelll`

## Try it Yourself

The following steps assume a working Rust toolchain.

1. clone this repository and `cd ./zelll/python/`
2. install [`maturin`](https://www.maturin.rs/tutorial#install-and-configure-maturin-in-a-virtual-environment) in a virtual environment
3. (optionally install `numpy` in the same environment for testing purposes)
4. run `maturin develop` to build and install a `.whl` into the current virtual environment
5. start a Python REPL and start playing around (currently w/o any documentation):
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
```

