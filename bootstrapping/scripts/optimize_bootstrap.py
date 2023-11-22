import sys
from trails.optimizer import trans_emiss_calc
from trails.cutpoints import cutpoints_ABC
import numpy as np
from trails.optimizer import loglik_wrapper, write_list, optimizer
import pandas as pd

####################### Loading parameters #######################

print('Loading params')

n_int_AB = int(sys.argv[1])
n_int_ABC = int(sys.argv[2])
model = sys.argv[3]
seed = int(sys.argv[4])

# Read the output from TRAILS
df = pd.read_csv('../../chr1/results/chr1_Nelder-Mead_third_run.csv')
# Find iteration with largest likelihood
df = df[df['loglik'] == df['loglik'].max()]
# Convert parameter estimates into dictionary
dct = dict(zip(list(df.columns), df.iloc[0].to_list()))

n_sites = 50_000_000

cut_ABC = cutpoints_ABC(n_int_ABC, 1)
t_out = (((dct["t_A"]+dct["t_B"])/2+dct["t_2"])+dct["t_C"])/2 + cut_ABC[n_int_ABC-1]*dct["N_ABC"] + dct["t_upper"] + 2*dct["N_ABC"]

####################### Simulation #######################

print('Simulating')

t_A = dct["t_A"]
t_B = dct["t_B"]
t_C = dct["t_C"]
t_2 = dct["t_2"]
t_upper = dct["t_upper"]
t_out = t_out
N_AB = dct["N_AB"]
N_ABC = dct["N_ABC"]
r = dct["r"]

transitions, emissions, starting, hidden_states, observed_states = trans_emiss_calc(
    t_A, t_B, t_C, t_2, t_upper, t_out,
    N_AB, N_ABC,
    r, n_int_AB, n_int_ABC)

np.random.seed(seed)
idx_lst = list(range(0, len(starting)))
idx = np.random.choice(idx_lst, size = 1, p = starting)[0]
# hid = np.zeros(n_sites, dtype = np.int16)
obs = np.zeros(n_sites, dtype = np.int16)
# hid[0] = idx
obs[0] = np.random.choice(list(range(256)), size = 1, p = emissions[idx])[0]
for i in range(1, n_sites):
    idx = np.random.choice(idx_lst, size = 1, p = transitions[idx])[0]
    # hid[i] = idx
    obs[i] = np.random.choice(list(range(256)), size = 1, p = emissions[idx])[0]

####################### Optimization #######################

print('Optimizing')

loglik = loglik_wrapper(transitions, emissions, starting, [obs])

write_list([-1, t_A, t_B, t_C, t_2, t_upper, N_AB, N_ABC, r, loglik, 0], '../results/sim_{}_{}_{}_{}.csv'.format(n_int_AB, n_int_ABC, seed, model))

t_init_A = t_A
t_init_B = t_B
t_init_C = t_C
t_init_2 = t_2
t_init_upper = t_upper
N_init_AB = N_AB
N_init_ABC = N_ABC
r_init = r

dct = {
    't_A':     [t_init_A,     t_init_A/10, t_init_A*10], 
    't_B':     [t_init_B,     t_init_B/10, t_init_B*10], 
    't_C':     [t_init_C,     t_init_C/10, t_init_C*10], 
    't_2':     [t_init_2,     t_init_2/10, t_init_2*10], 
    't_upper': [t_init_upper, t_init_upper/10, t_init_upper*10], 
    'N_AB':    [N_init_AB,    N_init_AB/10,  N_init_AB*10], 
    'N_ABC':   [N_init_ABC,   N_init_ABC/10,  N_init_ABC*10], 
    'r':       [r_init,       r_init/10,  r_init*10]
    }

dct2 = {'n_int_AB':n_int_AB, 'n_int_ABC':n_int_ABC}
res = optimizer(
    optim_params = dct, 
    fixed_params = dct2, 
    V_lst = [obs], 
    res_name = f'../results/sim_{n_int_AB}_{n_int_ABC}_{seed}_{model}.csv', 
    header = False,
    method = "Nelder-Mead"
    )
