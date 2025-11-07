#!/usr/bin/env julia
using ArgTools
using CSV
using BenchmarkTools
using CellListMap
using DataFrames
using Printf


#ARGS[1] for simple file name input
# input "parsing"
df = arg_read(ARGS[1]) do data_io
    DataFrame(CSV.File(data_io; header=false, skipto=11, select=[3, 4, 5]))
end

n = length(df[!, 1])

function lj(dsq)
    tmp = (1 / dsq)^3
    4.0 * tmp * (tmp - 1.0)
end

cutoff = 10.0
concentration = 10 / cutoff^3
a = 3.0 * cutoff
b = 3.0 * cutoff
c = max((n / concentration) / a / b, 3.0 * cutoff)
sides = [a, b, c]

particles = df |> Tables.matrix |> transpose


function compute(data)
    box = Box(sides, cutoff)
    cl = CellList(data, box)

    map_pairwise!(
        (x, y, i, j, dsq, acc) -> lj(dsq) + acc,
        0.0,
        box,
        cl,
        parallel=false,
    ) / n # total energy per atom (as done by LAMMPS)
end

b = @benchmarkable compute(x) setup = (x = copy(particles))
results = run(b)


total_energy = compute(particles)
@printf(
    "%d 1 %f %f %f \"CellListMap.jl\"\n",
    n,
    total_energy,
    mean(results.times) * 1e-9, # seconds
    results.memory / 1024^2, # MB
)
