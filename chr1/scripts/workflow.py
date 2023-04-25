from gwf import Workflow

gwf = Workflow()


gwf.target(f"chr1_maffilter",
    inputs=['../data/chr1.maf.gz'],
    outputs=['../data/chr1.filtered.maf'],
    cores=1,
    memory='4g',
    walltime= '12:00:00',
    account='Primategenomes') << f"""
conda activate maffilter
maffilter param=options_2.txt file=../data/chr1 sp1=hg38 sp2=panTro5 sp3=gorGor5 sp4=ponAbe2 ref=hg38
"""

start = 25000000
end   = 75000000

gwf.target(f"chr1_extract_region",
    inputs=['../data/chr1.filtered.maf'],
    outputs=['../data/chr1.filtered.region.maf'],
    cores=1,
    memory='4g',
    walltime= '1:00:00',
    account='Primategenomes') << f"""
python maf_extract_region.py ../data/chr1.filtered hg38.chr1 {start} {end}
"""


n_int_AB = 3
n_int_ABC = 3

N_AB = 200000
N_ABC = 200000
t_1 = 300000
t_2 = 100000
t_3 = 500000
r = 1e-08
mu = 2e-08

t_A = t_1
t_B = t_1
t_C = t_1+t_2

for method in ['L-BFGS-B', 'Nelder-Mead']:
    # print(f"../results/chr21_{method}.csv")
    gwf.target(f"chr1_{method.replace('-', '')}_newparam",
        inputs=['../data/chr1.filtered.region.maf'],
        outputs=[f"../results/chr1_{method}.csv"],
        cores=4,
        memory='12g',
        walltime= '{}:00:00'.format(72),
        account='Primategenomes') << f"""
    python optimize_4.py {t_A} {t_B} {t_C} {t_2} {t_3} {N_AB} {N_ABC} {r} {mu} {n_int_AB} {n_int_ABC} {method}
    """

n_int_AB = 1
n_int_ABC = 5

for method in ['L-BFGS-B', 'Nelder-Mead']:
    # print(f"../results/chr21_{method}.csv")
    gwf.target(f"chr1_{method.replace('-', '')}_newparam_second_run",
        inputs=['../data/chr1.filtered.region.maf', "../results/chr1_L-BFGS-B.csv"],
        outputs=[f"../results/chr1_{method}_second_run.csv"],
        cores=6,
        memory='12g',
        walltime= '{}:00:00'.format(120),
        account='Primategenomes') << f"""
    python optimize_5.py {n_int_AB} {n_int_ABC} {method}
    """

