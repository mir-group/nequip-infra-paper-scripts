import numpy as np
from openff.toolkit.topology import Molecule
from openmm.unit import *
from collections import defaultdict
import h5py

typeDict = {
    ("B", -1): 0,
    ("B", 0): 1,
    ("B", 1): 2,
    ("Br", -1): 3,
    ("Br", 0): 4,
    ("C", -1): 5,
    ("C", 0): 6,
    ("C", 1): 7,
    ("Ca", 2): 8,
    ("Cl", -1): 9,
    ("Cl", 0): 10,
    ("F", -1): 11,
    ("F", 0): 12,
    ("H", 0): 13,
    ("I", -1): 14,
    ("I", 0): 15,
    ("K", 1): 16,
    ("Li", 1): 17,
    ("Mg", 2): 18,
    ("N", -1): 19,
    ("N", 0): 20,
    ("N", 1): 21,
    ("Na", 1): 22,
    ("O", -1): 23,
    ("O", 0): 24,
    ("O", 1): 25,
    ("P", -1): 26,
    ("P", 0): 27,
    ("P", 1): 28,
    ("S", -1): 29,
    ("S", 0): 30,
    ("S", 1): 31,
    ("Si", 0): 32,
}
appearing_types = {}
main_infile = h5py.File(
    "/lustre/orion/mat281/proj-shared/24_11_04_BigData/SPICE-2/data/SPICE-test-2.hdf5"
)
# First pass: group the samples by total number of atoms.

groupsSplit = {}
groupsSplit["dimers"] = defaultdict(list)


for name in main_infile:
    group = main_infile[name]
    count = len(group["atomic_numbers"])

    groupsSplit["dimers"][count].append(group)

# Create the output file.


posScale = 1 * bohr / angstrom
energyScale = 1 * hartree / item / (kilojoules_per_mole)
forceScale = energyScale / posScale

for split, groupsByAtomCount in groupsSplit.items():
    filename = f"SPICE-test-2-processed_charges_fullprecision_totalenergy_{split}.hdf5"
    outfile = h5py.File(filename, "w")

    # One pass for each number of atoms, creating a group for it.

    print(sorted(list(groupsByAtomCount.keys())))
    appearingtypes = set()
    for count in sorted(groupsByAtomCount.keys()):
        print(count)
        smiles = []
        subset = []
        pos = []
        types = []
        energy = []
        forces = []
        charges = []
        for g in groupsByAtomCount[count]:
            molSmiles = g["smiles"][0]
            molSubset = g["subset"][0]

            mol = Molecule.from_mapped_smiles(molSmiles, allow_undefined_stereo=True)

            molTypes = [
                typeDict[(atom.symbol, atom.formal_charge.magnitude)]
                for atom in mol.atoms
            ]
            print("molTypes", molTypes)
            assert len(molTypes) == count

            for i, atom in enumerate(mol.atoms):
                assert atom.atomic_number == g["atomic_numbers"][i]

            appearingtypes.update(set(molTypes))
            assert len(molTypes) == count
            numConfs = g["conformations"].shape[0]
            for i in range(numConfs):
                smiles.append(molSmiles)
                subset.append(molSubset)
                pos.append(g["conformations"][i])
                types.append(molTypes)
                energy.append(g["dft_total_energy"][i])
                forces.append(g["dft_total_gradient"][i])
                if "mbis_charges" in g.keys():
                    charges.append(g["mbis_charges"][i])

                else:
                    charges.append(np.ones((count, 1)) * -1000)
                    print(
                        f"At count: {count}, molecule {molSmiles} has no mbis charges"
                    )
        print("appearing_types", appearing_types)
        group = outfile.create_group(f"samples{count}")
        group.create_dataset("smiles", data=smiles, dtype=h5py.string_dtype())
        group.create_dataset("subset", data=subset, dtype=h5py.string_dtype())
        group.create_dataset("types", data=np.array(types), dtype=np.int8)
        group.create_dataset("pos", data=np.array(pos) * posScale, dtype=np.float64)
        group.create_dataset(
            "energy", data=np.array(energy) * energyScale, dtype=np.float64
        )
        group.create_dataset(
            "forces", data=-np.array(forces) * forceScale, dtype=np.float64
        )
        group.create_dataset("mbis_charges", data=np.array(charges), dtype=np.float64)
