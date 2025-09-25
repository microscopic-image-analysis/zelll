# `psssh` 🤫: Protein Structure Surface Sampling using HMC

Example application of `zelll`.
This contains the case study accompanying `zelll`.
Here, `CellGrid` neighborhood queries are used to approximate a smooth
distance function which is then used to define a prior distribution to
sample from on protein structure surfaces using the No-U-Turn-Sampler (NUTS).

See also
[`../python/examples/psssh.py`](https://github.com/microscopic-image-analysis/zelll/tree/main/python/examples/psssh.py)
for a condensed and illustrative implementation.

## Building

```sh
RUSTFLAGS="-C target-cpu=native" cargo build --release --example psssh
```

## Usage

```sh
RUSTFLAGS="-C target-cpu=native" cargo run --release --example psssh -- help sample
```

```
Sample points on the surface of a single protein structure

Usage: psssh sample [OPTIONS] <PDB> [OUT]

Arguments:
  <PDB>  Protein structure file to sample on
  [OUT]  File path to save sampled surface to. defaults to input path + ".psssh.pdb"

Options:
  -c, --cutoff <CUTOFF>
          Neighborhood cutoff treshold used for sampling [default: 10]
  -n, --samples <N>
          Number of samples to produce [default: 2000]
  -b, --burn-in <B>
          Number of samples to discard before sampling 'n' samples [default: 1000]
  -l, --surface-level <SURFACE_LEVEL>
          Distance to the protein structure at which the surface will be sampled [default: 1.05]
  -f, --force-constant <FORCE_CONSTANT>
          Force constant used for sampling on the protein surface. Smaller values might work better with smaller cutoff radii but might also require adjusting the surface level [default: 10]
  -d, --nuts-depth <NUTS_DEPTH>
          Maximum tree depth for NUTS. The default value is robust enough for this application. Lower values are cheaper and may suffice if it's not required to cover the complete surface with the sampled points or the sample size is large enough [default: 7]
  -h, --help
          Print help
  -V, --version
          Print version
```

### Benchmark

`scripts/sdf_queries.sh` contains an ad-hoc benchmark
of the internally used approximate smooth distance function:

```sh
scripts/sdf_queries.sh > sdf_queries.csv
```

