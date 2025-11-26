# Simon, modified for new work by Marc:
import numpy as np
import h5py

from ase import Atoms
from ase.io import write
from ase import units
from ase.calculators.singlepoint import SinglePointCalculator

from collections import defaultdict


for split in ["small", "large", "pentapeptides"]:
    infile = h5py.File(
        f"/lustre/orion/mat281/proj-shared/24_11_04_BigData/SPICE-2/Revise-TotalEnergy/01_CreateDataset/SegmentedTest/SPICE-test-processed_charges_fullprecision_totalenergy_{split}.hdf5"
    )

    # Create the output file.
    filename = f"/lustre/orion/mat281/proj-shared/24_11_04_BigData/SPICE-2/Revise-TotalEnergy/01_CreateDataset/SegmentedTest/SPICE-test-processed_element-types_all-totalcharge_use-dfttotalenergy_units-eV-A_{split}.extxyz"

    # From dataset creation
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
    # For easily determining formal charge
    reverseTypeDict = dict((value, key) for key, value in typeDict.items())

    # Define the filters
    filter_force = False
    max_force = 51.42208619083232  # Filter any structure out with a force component over 1 Hartree/Bohr = 51.422 eV/Ang

    filter_molecular_charge = False
    allowed_charge = [0]

    filter_formal_charge = False
    allowed_formal_charge = [0]

    filter_num_atoms = False
    allowed_num_atoms = [2]

    save_charges = False

    # First pass: group the samples by total number of atoms.

    groupsByAtomCount = defaultdict(list)
    for name in infile:
        group = infile[name]
        count = len(group["types"])
        groupsByAtomCount[count].append(group)

    # One pass for each number of atoms, creating a group for it.
    print(sorted(list(groupsByAtomCount.keys())))

    KJMOL2EV = units.kJ / units.mol
    skip_counter = 0
    for count in sorted(groupsByAtomCount.keys()):
        print(count)
        for g in groupsByAtomCount[count]:
            numConfs = g["pos"].shape[0]
            for i in range(numConfs):
                if filter_num_atoms:
                    num_atoms = len(g["types"][i])
                    if num_atoms not in allowed_num_atoms:
                        skip_counter += 1
                        continue
                sum_charge = None
                if filter_force:
                    if np.max(g["forces"][i] * units.kJ / units.mol) > max_force:
                        skip_counter += 1
                        continue

                if filter_formal_charge:
                    sum_charge = 0
                    forbidden = False
                    for datasettype in g["types"][i]:
                        formal_charge = reverseTypeDict[datasettype][1]
                        sum_charge += formal_charge
                        if formal_charge not in allowed_formal_charge:
                            forbidden = True

                    if forbidden:
                        skip_counter += 1
                        continue

                if filter_molecular_charge:
                    if sum_charge is None:
                        sum_charge = sum(
                            [
                                reverseTypeDict[datasettype][1]
                                for datasettype in g["types"][i]
                            ]
                        )
                    if sum_charge not in allowed_charge:
                        skip_counter += 1
                        continue

                symbols = [
                    reverseTypeDict[datasettype][0] for datasettype in g["types"][i]
                ]
                curr_atoms = Atoms(positions=g["pos"][i], symbols=symbols, pbc=False)

                if save_charges:
                    # This is the case that I added if the dataset was missing an mbis_charges field.
                    if (g["mbis_charges"][i] == -1000).all():
                        skip_counter += 1
                        continue
                    else:
                        curr_atoms.new_array("mbis_charges", g["mbis_charges"][i])

                calculator = SinglePointCalculator(
                    curr_atoms,
                    energy=g["energy"][i] * units.kJ / units.mol,
                    forces=g["forces"][i] * units.kJ / units.mol,
                )
                curr_atoms.calc = calculator
                print(g["subset"][i])
                curr_atoms.info["subset"] = g["subset"][i].decode("utf-8")
                write(filename, curr_atoms, format="extxyz", append=True)

    print("Skipped", skip_counter, "frames")
