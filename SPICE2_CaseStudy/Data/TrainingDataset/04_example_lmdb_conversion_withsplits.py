"""
Example of how to use NequIPLMDBDataset to convert an xyz file to LMDB format, using ase and some nequip ase utilities. This can be adapted to convert data from other formats, as long as there one writes code to convert data from the custom format to nequip's AtomicDataDict format.
"""

import torch
import ase
from nequip.data.dataset import NequIPLMDBDataset
from nequip.data.ase import from_ase
import numpy as np
import json
from tqdm import tqdm


train_lmdb_file_path = "SPICE-2.0.1-processed_element-types_all_within-6.0_0-totalcharge_use-dfttotalenergy_units-eV-A_train_v2.lmdb"
val_lmdb_file_path = "SPICE-2.0.1-processed_element-types_all_within-6.0_0-totalcharge_use-dfttotalenergy_units-eV-A_val_v2.lmdb"
test_lmdb_file_path = "SPICE-2.0.1-processed_element-types_all_within-6.0_0-totalcharge_use-dfttotalenergy_units-eV-A_test_v2.lmdb"
json_file_path = "SPICE-2.0.1-processed_element-types_all_within-6.0_0-totalcharge_use-dfttotalenergy_units-eV-A_subsets_v2.json"

xyz_file_path = "SPICE-2.0.1-processed_element-types_all_within-6.0_0-totalcharge_use-dfttotalenergy_units-eV-A.extxyz"


# make sure the data is saved in float64!
torch.set_default_dtype(torch.float64)

# === ase.Atoms -> AtomicDataDict ===
atoms_list = list(
    tqdm(
        ase.io.iread(filename=xyz_file_path, parallel=False),
        desc="Reading dataset with ASE...",
    )
)

train_split = 0.90
val_split = 0.05
test_split = 0.05
seed = 20244202

np.random.seed(seed)
assert train_split + val_split + test_split == 1


indices = np.arange(len(atoms_list))
np.random.shuffle(indices)
train_end = int(train_split * len(atoms_list))
val_end = train_end + int(val_split * len(atoms_list))

train_indices = indices[:train_end]
val_indices = indices[train_end:val_end]
test_indices = indices[val_end:]

comparison1 = set(train_indices)
comparison2 = set(val_indices)
for index in test_indices:
    assert index not in comparison1
    assert index not in comparison2
comparison2 = set(test_indices)
for index in val_indices:
    assert index not in comparison1
    assert index not in comparison2


print(f"Train end is: {train_end}, val end is {val_end}, length is {len(atoms_list)}")
print(f"This should be a shuffled order: {indices[:10]}")

np.save("split_train.npy", train_indices)
np.save("split_val.npy", val_indices)
np.save("split_test.npy", test_indices)

# Label subset:
subset_map = {}
subset_index = 0
elements = set()
for frame in atoms_list:
    subset = frame.info["subset"]
    if subset not in subset_map.keys():
        subset_map[subset] = subset_index
        subset_index += 1
    frame.info["subset"] = torch.tensor(subset_map[subset], dtype=torch.long)
    elements.update(set(np.unique(frame.symbols)))
with open(json_file_path, "w") as file:
    json.dump(subset_map, file)

print(
    f"Included elements in {xyz_file_path},{train_lmdb_file_path},{val_lmdb_file_path},{test_lmdb_file_path} are {elements}"
)

# Test:
test_atomic_data_dicts = (
    from_ase(atoms_list[index], include_keys=["subset"])
    for index in tqdm(test_indices, desc="Saving Test to LMDB...")
)

# === convert to LMDB ===
NequIPLMDBDataset.save_from_iterator(
    file_path=test_lmdb_file_path,
    iterator=test_atomic_data_dicts,
    write_frequency=5000,  # increase this from default 1000 to speed up writing of very large datasets
)

# Validation:
val_atomic_data_dicts = (
    from_ase(atoms_list[index], include_keys=["subset"])
    for index in tqdm(val_indices, desc="Saving Validation to LMDB...")
)

# === convert to LMDB ===
NequIPLMDBDataset.save_from_iterator(
    file_path=val_lmdb_file_path,
    iterator=val_atomic_data_dicts,
    write_frequency=5000,  # increase this from default 1000 to speed up writing of very large datasets
)


# Train:
train_atomic_data_dicts = (
    from_ase(atoms_list[index], include_keys=["subset"])
    for index in tqdm(train_indices, desc="Saving Train to LMDB...")
)

# === convert to LMDB ===
NequIPLMDBDataset.save_from_iterator(
    file_path=train_lmdb_file_path,
    iterator=train_atomic_data_dicts,
    write_frequency=5000,  # increase this from default 1000 to speed up writing of very large datasets
)
