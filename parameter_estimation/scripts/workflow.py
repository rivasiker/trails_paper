from gwf import Workflow
from os.path import exists
import numpy as np

gwf = Workflow()

ILS = 32

t_1 = 200000
N_AB = 80000
t_2 = -N_AB*np.log(3/2*ILS/100)
N_ABC = 70000
t_3 = t_1*5
r = 0.5e-8
mu = 1.5e-8

dct = {1:4, 2:30,  3:57, 4:90, 5:121, 7:251}
# dct = {1:4, 2:15,  3:60, 4:45, 5:60, 7:125, 9:250}
tot_lst = []

for n_int_AB in [1, 3, 5]:
    for n_int_ABC in [1, 3, 5]:
        for seed in range(1, 6):
            if exists(f'../results/sim_{n_int_AB}_{n_int_ABC}_{seed}_fixed_model.csv'):
                continue
            tot_lst.append('simulate_{}_{}_{}_fixed_model'.format(n_int_AB, n_int_ABC, seed))
            gwf.target('simulate_{}_{}_{}_fixed_model'.format(n_int_AB, n_int_ABC, seed), 
                inputs=['optimize_2.py'], 
                outputs=['../results/{}_{}_{}_{}_fixed_model.csv'.format(x, n_int_AB, n_int_ABC, seed) for x in ['trees', 'sim']],  
                cores=n_int_ABC,
                memory='{}g'.format(n_int_ABC*4),
                walltime= '{}:00:00'.format(dct[n_int_ABC]),
                account='Primategenomes') << """
            python optimize_2.py {} {} {} {} {} {} {} {} {} {} fixed_model
            """.format(seed, t_1, t_2, t_3, N_AB, N_ABC, r, mu, n_int_AB, n_int_ABC)

dct = {1:15, 2:30,  3:57, 4:90, 5:180, 7:251}

t_A = t_1
t_B = t_1
t_C = t_1+t_2

for n_int_AB in [1, 5]:
    for n_int_ABC in [1]:
        n_int_ABC = n_int_AB
        for seed in range(1, 21):
            if exists(f'../results/sim_{n_int_AB}_{n_int_ABC}_{seed}_error_model.csv'):
                continue
            tot_lst.append('simulate_{}_{}_{}_error_model'.format(n_int_AB, n_int_ABC, seed))
            gwf.target('simulate_{}_{}_{}_error_model'.format(n_int_AB, n_int_ABC, seed),
                inputs=['optimize_4.py'],
                outputs=['../results/{}_{}_{}_{}_error_model.csv'.format(x, n_int_AB, n_int_ABC, seed) for x in ['trees', 'sim']],
                cores=n_int_ABC,
                memory='{}g'.format(n_int_ABC*4),
                walltime= '{}:00:00'.format(dct[n_int_ABC]),
                account='Primategenomes') << f"""
            python optimize_4.py {seed} {t_A} {t_B} {t_C} {t_2} {t_3} {N_AB} {N_ABC} {r} {mu} {n_int_AB} {n_int_ABC} error_model
            """

#[print(i, end = ' ') for i in tot_lst]
