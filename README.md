# `zelll`: a Rust implementation of the cell lists algorithm.

[![Crates.io](https://img.shields.io/crates/v/zelll.svg)](https://crates.io/crates/zelll)
[![Documentation](https://docs.rs/zelll/badge.svg)](https://docs.rs/zelll)

Particle simulations usually require to compute interactions between those particles.
Considering all _pairwise_ interactions of _`n`_ particles would be of time complexity _`O(n²)`_.\
Cell lists facilitate _linear-time_ enumeration of particle pairs closer than a certain
cutoff distance by dividing the enclosing bounding box into (cuboid) grid cells.

## Caveats

`zelll`[^etymology] is motivated by _coarse-grained_ (bio-)molecular simulations but is not restricted to that.\
This is reflected by a few points:

- internally, the simulation box is represented by a (sparse) hash map only storing non-empty grid cells,
  which gives an upper bound for memory usage of _`n`_
- bounding boxes are assumed to change and are computed from particle data\
  (future APIs may be added to set a fixed bounding box)
- instead of cell _lists_, slices into a contiguous storage buffer are used
- periodic boundary conditions are currently not supported
- parts of this implementation are more cache-aware than others, which becomes noticeable with
  larger data sets\
  (at `10⁶` -- `10⁷` particles, mostly depending on last-level cache size)
  but is less pronounced with structured data[^structureddata]

## Usage

The general pattern in which this crate is intended to be used is roughly:

1. construct `CellGrid` from particle positions
2. enumerate pairs in order to compute particle interactions
3. simulate particle motion
4. rebuild `CellGrid` from updated particle positions

This crate only provides iteration over particle pairs.
It is left to the user to filter (e.g. by distance) and compute interaction potentials.
The `rayon` feature enables parallel iteration. Performance gains depend on data size and
computational cost per pair though. Benchmarks are encouraged.
The `serde` feature flag enables serialization.

This crate is intended for simulations where performance is often paramount.
The rust compiler offers [codegen options](https://doc.rust-lang.org/rustc/codegen-options/index.html#target-cpu) 
that can be useful in these settings, e.g. like this:

```sh
RUSTFLAGS="-C target-cpu=native" cargo bench --features rayon
```

Experimental Python bindings can be found in the
[`python/`](https://github.com/microscopic-image-analysis/zelll/tree/main/python)
directory.
These bindings are not intended for productive use because they are
incomplete, likely contain unsound code and carry significant overhead.
They are, however, suitable for exploratory purposes.

### Examples

```rust
use zelll::CellGrid;

let data = vec![[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
let mut cg = CellGrid::new(data.iter().copied(), 1.0);

for ((i, p), (j, q)) in cg.particle_pairs() {
    /* do some work */
}

cg.rebuild_mut(data.iter().copied(), Some(0.5));
```


### Benchmarks

In addition to the `rayon` feature flag, benchmarks also read `quick_bench`
for reduced sample sizes as full benchmarks may take quite some time.

```sh
# only runs the "Iteration" benchmark (the other valid choice is "CellGrid") 
RUSTFLAGS="-C target-cpu=native" cargo bench --features quick_bench,rayon -- Iteration
```

Cache misses are measured via `scripts/cachemisses.sh`:

```sh
# this requires a Valgrind installation
# presorted data: false, f32: false
./scripts/cachemisses.sh false false > cachemisses.csv
```

## Case Study

Information for a self-contained example can be found in the 
[`surface-sampling/`](https://github.com/microscopic-image-analysis/zelll/tree/main/surface-sampling)
directory.

## Roadmap

These are improvements we want to make eventually:

- [ ] parallel `CellGrid` construction
    * might help a bit with cache awareness
    * possible approach: merging 2 `CellGrid`s into one
        - cell indices maximum bounding box might help here
    * explore [`cubecl`](https://crates.io/crates/cubecl)
- [ ] periodic boundaries
- [ ] revisit flat cell indices 
    * maximum bounding box 
    * other hashing approaches
- [ ] redo `CellStorage`, this is rather hacky at the moment

[^etymology]: abbrv. from German _Zelllisten_ /ˈʦɛlɪstən/, for cell lists.
[^structureddata]: Usually, (bio-)molecular data files are not completely unordered
    even though they could be.
    In practice, it may be a reasonable assumption that sequentially proximate
    particles often have spatially clustered coordinates as well.
