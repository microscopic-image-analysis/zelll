# `zelll`: a Rust implementation of the cell lists algorithm.

<!--
[![Crates.io](https://img.shields.io/crates/v/zelll.svg)](https://crates.io/crates/zelll)
[![Documentation](https://docs.rs/zelll/badge.svg)](https://docs.rs/zelll)
-->

Particle simulations usually require to compute interactions between those particles.
Considering all _pairwise_ interactions of _`n`_ particles would be of time complexity _`O(n²)`_.\
Cell lists facilitate _linear-time_ enumeration of particle pairs closer than a certain
cutoff distance by dividing the enclosing bounding box into (cuboid) grid cells.

## Caveats

`zelll`[^etymology] is motivated by _coarse-grained_ (bio-)molecular simulations but is not restricted to that.\
This is reflected by a few things:

- internally, the simulation box is represented by a (sparse) hash map only storing non-empty grid cells,
  which gives an upper bound for memory usage given by _`n`_
- bounding boxes are assumed to change and are computed from particle data\
  (future APIs may be added to set a fixed bounding box)
- instead of cell _lists_, slices into a contiguous storage buffer are used
- periodic boundary conditions are currently not supported
- parts of this implementation are more cache-aware than others, which becomes noticeable with
  larger data sets\
  (at `10⁶` -- `10⁷` particles, mostly depending on L2 cache size)
  but is less pronounced with structured data[^structureddata]

## Usage

The general pattern in which this crate is intended to be used is roughly:

1. construct `CellGrid` from particle positions
2. enumerate pairs in order to compute particle interactions
3. simulate particle motion
4. rebuild `CellGrid` from updated particle positions

This crate only provides iteration over particle pairs.
It is left to the user to filter (eg. by distance) and compute interaction potentials.
The `rayon` feature enables parallel iteration. Performance gains depend on data size and
computational cost per pair though. Benchmarks are encouraged.

This crate is intended for simulations where performance is often paramount.
The rust compiler offers [codegen options](https://doc.rust-lang.org/rustc/codegen-options/index.html#target-cpu) 
that can be useful in these settings, eg. like this:

```sh
RUSTFLAGS="-C target-cpu=native" cargo bench --features rayon
```

## Examples
```rust
use zelll::CellGrid;

let data = vec![[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
let mut cg = CellGrid::new(data.iter().copied(), 1.0);

for ((i, p), (j, q)) in cg.particle_pairs() {
    /* do some work */
}

cg.rebuild_mut(data.iter().copied(), Some(0.5));
```

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
    * completely different "hashing" method
    * hashing/flat indices are cheap, we don't really need to allocate memory for that
- [ ] redo `CellStorage`, this is rather hacky at the moment


## ToDo (temporary)

- [ ] finish python bindings + docs
    * only expose subset of API
- [ ] clean up examples
    * [ ] remove rendering dependencies
- [ ] showcase application: docking decoy sets
    * put into subcrate
- [ ] benchmarks
    * [x] cache misses (cf. `minimal_new.rs`)
    * [ ] clean up
    * [ ] figures for the above

[^etymology]: abbrv. from German _Zelllisten_ /ˈʦɛlɪstən/, for cell lists.
[^structureddata]: Usually, (bio-)molecular data files are not completely unordered
    even though they could be.
    In practice, it may be a reasonable assumption that sequentially proximate
    particles often have spatially clustered coordinates as well.
