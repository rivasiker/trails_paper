import sys
from trails.cutpoints import cutpoints_ABC, cutpoints_AB
from trails.optimizer import optimizer
from trails.read_data import maf_parser
import pandas as pd



####################### Model parameters #######################

n_int_AB = int(sys.argv[1])
n_int_ABC = int(sys.argv[2])
method = sys.argv[3]
interval = sys.argv[4]

df = pd.read_csv('../results/chr1_L-BFGS-B.csv')
val_dct = df[df['loglik'] == df['loglik'].max()].iloc[0].to_dict()

t_A = val_dct['t_A']
t_B = val_dct['t_B']
t_C = val_dct['t_C']
t_2 = val_dct['t_2']
t_upper = val_dct['t_upper']
N_AB = val_dct['N_AB']
N_ABC = val_dct['N_ABC']
r = val_dct['r']

t_3 = t_upper + cutpoints_ABC(3, 1/N_ABC)[-2]
t_upper = t_3-cutpoints_ABC(n_int_ABC, 1/N_ABC)[-2]

t_init_A = t_A
t_init_B = t_B
t_init_C = t_C
t_init_2 = t_2
t_init_upper = t_upper
N_init_AB = N_AB
N_init_ABC = N_ABC
r_init = r

dct = {
    't_A':     [t_init_A,     t_A/10, t_A*10], 
    't_B':     [t_init_B,     t_B/10, t_B*10], 
    't_C':     [t_init_C,     t_C/10, t_C*10], 
    't_2':     [t_init_2,     t_2/10, t_2*10], 
    't_upper': [t_init_upper, t_upper/10, t_upper*10], 
    'N_AB':    [N_init_AB,    N_AB/10,  N_AB*10], 
    'N_ABC':   [N_init_ABC,   N_ABC/10,  N_ABC*10], 
    'r':       [r_init,       r/10,  r*10]
    }

print(dct)

dct2 = {'n_int_AB':n_int_AB, 'n_int_ABC':n_int_ABC}

sp_lst = ['hg38','panTro5','gorGor5','ponAbe2']
E = maf_parser(f'../data/chr1.filtered.{interval}.region.maf', sp_lst)

res = optimizer(
    optim_params = dct, 
    fixed_params = dct2, 
    V_lst = E, 
    res_name = f'../results/chr1_{method}_second_run_{interval}.csv',
    method = method, 
    header = True
    )
