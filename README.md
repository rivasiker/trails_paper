# Code for reproducing the analyses in [Rivas-González et al. 2024](https://doi.org/10.1371/journal.pgen.1010836)

TRAILS is a coalescent hidden Markov model (HMM) that reconstructs demographic parameters (ancestral Ne and split times) for three species and an outgroup. After model fitting, TRAILS can also perform posterior decoding to identify incomplete lineage sorting fragments. 

In this repository, you can find the code for reproducing the analyses in the original TRAILS publication ([Rivas-González et al. 2024](https://doi.org/10.1371/journal.pgen.1010836)):

- `trails_paper/parameter_estimation` contains the scripts for reproducing Fig 2.
- `trails_paper/simulated_posterior` contains the scripts for reproducing Fig 3 and Fig 4C.
- `trails_paper/selection_posterior` contains the scripts for reproducing Fig 4A and B.
- `trails_paper/chr1` contains the scripts for reproducing Fig 5. 

The scripts in each of the directories are unified with a [gwf workflow](https://gwf.app). The packages neccessary to run the workflows can be installed by running the following:

```{bash}
conda create -n trails -c conda-forge python=3.8 "ray-default" pandas numba numpy scipy biopython gwf msprime
```

The TRAILS python package can be installed from pip:

```{bash}
pip install trails-rivasiker
```

After running the scripts, the analyses can be performed using the following packages:

```{bash}
conda create -n trails_plot python=3.8 numpy pandas scipy r-essentials r-ggrastr r-patchwork pip tqdm jupyterlab rpy2 ray
```
