import numpy as np

from ase import Atoms
from ase.io import read, write
from ase.calculators.singlepoint import SinglePointCalculator
from ase.geometry import get_distances

# This script is designed to take the ion pair frames from SPICE and remove any frames where the ions are
# farther apart than the allegro model's cutoff is going to be.
# This is needed so that we avoid the error about having no neighbor atoms for a
# system when building the nequip dataset

infile = "../01_CreateDataset/SPICE-2.0.1-processed_charges_element-types_all_all-totalcharge_use-dfttotalenergy_units-eV-A.extxyz"

# Create the output file.
cutoff = 6.0
filename = f"SPICE-2.0.1-processed_element-types_all_within-{cutoff}_all-totalcharge_use-dfttotalenergy_units-eV-A.extxyz"

frames = read(infile, index=":", format="extxyz")

for frame in frames:
    curr_atoms = Atoms(
        positions=frame.get_positions(), symbols=frame.symbols, pbc=False
    )
    calculator = SinglePointCalculator(
        curr_atoms, energy=frame.info["energy"], forces=frame.get_forces()
    )
    curr_atoms.info["subset"] = frame.info["subset"]
    curr_atoms.calc = calculator

    _, distance_matrix = get_distances(frame.get_positions())

    mask = distance_matrix > 0
    distance_matrix_masked = np.where(mask, distance_matrix, np.nan)
    # Minimum neighbor distance for each atom
    row_min_distances = np.nanmin(distance_matrix_masked, axis=1)
    # Furthest of the nearest neighbors in the frame.
    minimum_neighbor_distance = np.nanmax(row_min_distances)

    if minimum_neighbor_distance < cutoff:
        write(filename, curr_atoms, format="extxyz", append=True)
    else:
        print("skipping", curr_atoms.positions)
