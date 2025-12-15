import numpy as np
import pymc as pm
import nutpie
import pytensor
import pytensor.tensor as pt
from pytensor.graph import Apply, Op
from zelll import CellGrid
from num_dual import gradient, first_derivative
from Bio.PDB import *
from argparse import ArgumentParser
from enum import Enum


class Element(Enum):
    C = 1.70
    H = 1.09
    O = 1.52
    N = 1.55
    S = 1.80
    Se = 1.90


class SdfOp(Op):
    # By default only show the first output, and "hide" the other ones
    default_output = 0

    def __init__(self, points, cutoff=1.0, atom_radii=None):
        self.inner = CellGrid(points, cutoff)
        self.radii = atom_radii

    # we use a parameter `neighbors` to make sure it's treated as constant during autodiff
    def _sdf(self, pos, neighbors):
        scaled_exp_dists = 0.0
        atom_radii = 0.0
        total_exp_dists = 0.0

        for radius, coords in neighbors:
            dist = np.linalg.norm(np.asarray(pos) - np.asarray(coords))
            if dist.value != 0.0:
                scaled_exp_dists += np.exp(-dist / radius)
                atom_radii += np.exp(-dist) * radius
                total_exp_dists += np.exp(-dist)
            else:
                scaled_exp_dists += 1.0
                atom_radii += radius
                total_exp_dists += 1.0

        sigma = atom_radii / total_exp_dists
        return -sigma * np.log(scaled_exp_dists)

    def eval_sdf(self, pos):
        # NOTE: in practice, this is not well suited to tensor libraries for reasons described e.g. here:
        # NOTE: https://pytensor.readthedocs.io/en/latest/extending/creating_a_c_op.html#methods-the-c-op-needs-to-define
        neighbors = self.inner.neighbors(pos)
        if not neighbors:
            return -np.inf, [-np.inf, -np.inf, -np.inf]

        if self.radii:
            neighbors = [(self.radii[i].value, coords) for i, coords in neighbors]
        else:
            neighbors = [(1.0, coords) for i, coords in neighbors]

        return gradient(lambda x: self._sdf(x, neighbors), pos)

    def make_node(self, x):
        x = pt.as_tensor(x)
        inputs = [x]
        outputs = [pt.dscalar()] + [x.type()]
        return Apply(self, inputs, outputs)

    def perform(self, node, inputs, outputs):
        result, grad_results = self.eval_sdf(inputs[0])
        outputs[0][0] = np.asarray(result, dtype=node.outputs[0].dtype)
        outputs[1][0] = np.asarray(grad_results, dtype=node.outputs[1].dtype)

    def grad(self, inputs, output_gradients):
        value = self(inputs[0])
        gradient = value.owner.outputs[1]
        return [output_gradients[0] * gradient]


# essentially (a special case of) a normal logp
def harmonic_potential(x, r, k):
    return -k * (x - r) ** 2


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-n", nargs="?", default=2000, type=int)
    parser.add_argument("-b", "--burn-in", nargs="?", default=1000, type=int)
    parser.add_argument("-d", "--nuts-depth", nargs="?", default=7, type=int)
    parser.add_argument("-l", "--surface-level", nargs="?", default=1.05, type=float)
    parser.add_argument("-f", "--force-constant", nargs="?", default=10.0, type=float)
    parser.add_argument("-c", "--cutoff", nargs="?", default=10.0, type=float)
    parser.add_argument("PDB")
    parser.add_argument("-o", "--out", nargs="?", default="psssh.pdb", type=str)

    args = parser.parse_args()

    pdbparser = PDBParser()
    structure = pdbparser.get_structure("input", args.PDB)
    points = [atom.coord for atom in structure.get_atoms()]
    atom_radii = [Element[atom.element] for atom in structure.get_atoms()]

    sdf_op = SdfOp(points, args.cutoff, atom_radii)

    with pm.Model() as model:
        # using padded bounding box for limits
        padding = np.array([5.0, 5.0, 5.0]) + args.surface_level
        l, u = sdf_op.inner.aabb()
        l, u = np.array(l) - padding, np.array(u) + padding
        x = pm.Uniform("x", shape=(3,), initval=points[0], lower=l, upper=u)

        # we're interested in uniformly sampling x under the constraints given by the potentials below:
        sdf = pm.Potential("sdf", sdf_op(x))

        k = args.force_constant
        r = args.surface_level
        harmonic = pm.Potential("harmonic", harmonic_potential(sdf, r, k))

        compiled = nutpie.compile_pymc_model(model)
        # zelll's Python bindings don't benefit from nutpie's multithreading
        # unless a free-threaded Python build is used
        # however, num_dual currently has no pre-built free-threaded wheels
        trace = nutpie.sample(
            compiled,
            save_warmup=False,
            chains=1,
            draws=args.n,
            maxdepth=args.nuts_depth,
            use_grad_based_mass_matrix=False,
            tune=args.burn_in,
            cores=1,
        )

        # alternatively use PyMC's internal NUTS implementation
        # trace = pm.sample(draws=100)

    chain = Chain.Chain("A")
    model = Model.Model(0)
    surface = Structure.Structure("PSSSH")
    model.add(chain)
    surface.add(model)

    atoms = np.vstack(trace.posterior["x"].data)

    for i, atom in enumerate(atoms):
        residue = Residue.Residue(("H", i, "A"), "R", "R")
        residue.add(Atom.Atom("H", atom, 1.0, 1.0, "H", "H", i, "H"))
        chain.add(residue)

    io = PDBIO()
    io.set_structure(surface)
    io.save(args.out)

