# `zelll`

## Outline & Goals
- implement "naive" cell list algorithm following Allen and Tildesley
- clean up implementation and make it ergonomic to use
    * work on docs
    * fill and improve README.md
- explore various optimization of memory locality:
    * manual prefetch instructions
    * sort point cloud by point location/multi index (points in the same cell should be close to each other in the cell list)
        - might make iteration simpler, naive algorithm relies heavily on individual indexing
        - sort either on updating the data structure or after a certain (entropic?) criterion (Christian said sth. about a cheap way...  )
        - also see "presortedness"
    * replace cell "list" and "head" array by other backing storage (octrees, z-trees, etc.; quite common structures in this context)
        - requires changes to multi index approach?
        - probably could just use individual "cell list" Vecs that are allocated on a (known) fixed-size arena instead of head+list arrays?
- maintain `O(n)` time complexity
- python bindings
- benchmarks (run time, cache hit/miss rate; see e.g. cachegrind)
- allow "minimal" bounding boxes to define the grid/lattice (currently using axis-aligned bounding box)
    * more general lattices
    * lattice boundary conditions?
    
## Examples

Caution, this is not very polished yet.
Parameters are hard-coded atm and memory usage isn't profiled yet.
```
cargo run --example borbs --release
```
