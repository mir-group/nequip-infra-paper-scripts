# High-performance training and inference for deep equivariant interatomic potentials

This repository contains publicly released scripts and examples to reproduce the results presented in the paper:

> Chuin Wei Tan, Marc L. Descoteaux, Mit Kotak, Gabriel de Miranda Nascimento, Seán R. Kavanagh, Laura Zichi, Menghang Wang, Aadit Saluja, Yizhong R. Hu, Tess Smidt, Anders Johansson, William C. Witt, Boris Kozinsky, Albert Musaelian. <br/>
> "High-performance training and inference for deep equivariant interatomic potentials." <br/>
> https://doi.org/10.48550/arXiv.2504.16068

---

## Repository Organization

### Benchmarking

#### Training
- `Fig2_TrainScaling/` - Training scalability benchmarks

#### Inference
- `Fig3_SingleRank/` - Single-rank inference benchmarks
- `Fig4_SingleRankWater/` - Single-rank water system benchmarks
- `Fig5_MultiNodeBioBench/` - Multi-node biological system benchmarks

### SPICE2 Case Study

#### Data
- `TrainingDataset/` - Training dataset files
- `TestingDatasets/` - Testing dataset files

#### Models
- `Configs/` - Model configuration files
- `Jobscripts/` - Job submission scripts
