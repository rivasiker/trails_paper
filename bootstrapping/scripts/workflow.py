from gwf import Workflow

gwf = Workflow()

n_int_AB = 1
n_int_ABC = 5
nsim = 30

for method in ['Nelder-Mead']:
    for seed in range(nsim):
        # print(f"../results/chr21_{method}.csv")
        gwf.target(f"bootstrap_{method.replace('-', '')}_{seed}",
            inputs=["../results/chr1_Nelder-Mead_third_run.csv"],
            outputs=[f"../results/sim_{n_int_AB}_{n_int_ABC}_{seed}_{method}.csv"],
            cores=6,
            memory='12g',
            walltime= '{}:00:00'.format(120),
            account='Primategenomes') << f"""
        python optimize_bootstrap.py {n_int_AB} {n_int_ABC} {method} {seed}
        """