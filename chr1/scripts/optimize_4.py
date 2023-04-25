import sys
from trails.cutpoints import cutpoints_ABC, cutpoints_AB
from trails.optimizer import optimizer
from trails.read_data import maf_parser



####################### Model parameters #######################

t_A = float(sys.argv[1])
t_B = float(sys.argv[2])
t_C = float(sys.argv[3])
t_2 = float(sys.argv[4])
t_3 = float(sys.argv[5])
N_AB = float(sys.argv[6])
N_ABC = float(sys.argv[7])
r = float(sys.argv[8])
mu = float(sys.argv[9])
n_int_AB = int(sys.argv[10])
n_int_ABC = int(sys.argv[11])
method = sys.argv[12]

t_upper = t_3-cutpoints_ABC(n_int_ABC, 1/N_ABC)[-2]

t_A = t_A*mu
t_B = t_B*mu
t_C = t_C*mu
t_2 = t_2*mu
t_upper = t_upper*mu
N_AB = N_AB*mu
N_ABC = N_ABC*mu
r = r/mu

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
E = maf_parser('../data/chr1.filtered.region.maf', sp_lst)

res = optimizer(
    optim_params = dct, 
    fixed_params = dct2, 
    V_lst = E, 
    res_name = f'../results/chr1_{method}.csv',
    method = method, 
    header = True
    )
