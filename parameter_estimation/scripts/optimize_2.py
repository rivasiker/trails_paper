import sys
from trails.optimizer import trans_emiss_calc
from trails.cutpoints import cutpoints_ABC, cutpoints_AB
import numpy as np
from trails.optimizer import loglik_wrapper, write_list, optimizer
import pandas as pd
import time
import re
import msprime


####################### Model parameters #######################

n_sites = 10_000_000

seed = int(sys.argv[1])
t_1 = float(sys.argv[2])
t_2 = float(sys.argv[3])
t_3 = float(sys.argv[4])
N_AB = float(sys.argv[5])
N_ABC = float(sys.argv[6])
r = float(sys.argv[7])
mu = float(sys.argv[8])
n_int_AB = int(sys.argv[9])
n_int_ABC = int(sys.argv[10])
model = str(sys.argv[11])

t_out = t_1+t_2+t_3+2*N_ABC

N_ref = N_ABC

coal_ABC = N_ref/N_ABC
coal_AB = N_ref/N_AB
t_upper = t_3-cutpoints_ABC(n_int_ABC, coal_ABC)[-2]*N_ref
t_AB = t_2/N_ref

cut_AB = t_1+cutpoints_AB(n_int_AB, t_AB, coal_AB)*N_ref
cut_ABC = t_1+t_2+cutpoints_ABC(n_int_ABC, coal_ABC)*N_ref

transitions, emissions, starting, hidden_states, observed_states = trans_emiss_calc(
    t_1, t_1, t_1, t_2, t_upper, t_out,
    N_AB, N_ABC, 
    r, n_int_AB, n_int_ABC)

dct_hid = {v: k for k, v in hidden_states.items()}
dct = {v: k for k, v in observed_states.items()}

####################### Add demography #######################

demography = msprime.Demography()
demography.add_population(name="A", initial_size=N_AB)
demography.add_population(name="B", initial_size=N_AB)
demography.add_population(name="C", initial_size=N_AB)
demography.add_population(name="D", initial_size=N_AB)
demography.add_population(name="AB", initial_size=N_AB)
demography.add_population(name="ABC", initial_size=N_ABC)
demography.add_population(name="ABCD", initial_size=N_ABC)
demography.add_population_split(time=t_1, derived=["A", "B"], ancestral="AB")
demography.add_population_split(time=t_1+t_2, derived=["AB", "C"], ancestral="ABC")
demography.add_population_split(time=t_1+t_2+t_3, derived=["ABC", "D"], ancestral="ABCD")

ts = msprime.sim_ancestry(
    {"A": 1, "B": 1, "C": 1, 
     "D": 1
    }, 
    demography=demography, 
    recombination_rate=r,
    sequence_length=n_sites,
    ploidy=1, 
    random_seed=seed
)


####################### Feature extraction #######################

left_lst = []
right_lst = []
tree_lst = []
for t in ts.trees():
    # Append start coordinate
    left_lst.append(t.interval.left)
    # Append end coordinate
    right_lst.append(t.interval.right-1)
    # Find topology
    tree_lst.append(t.as_newick(include_branch_lengths=True))
        
len_lst = (np.array(right_lst)-np.array(left_lst)+1)

dat_sim = pd.DataFrame({'tree':tree_lst, 'len':len_lst})
dat_sim.to_csv('../results/trees_{}_{}_{}_{}.csv'.format(n_int_AB, n_int_ABC, seed, model), index = False)

#### Add mutations

mutated_ts = msprime.sim_mutations(ts, rate=mu, random_seed=seed)

nochange_lst = [dct['AAAA'], dct['CCCC'], dct['TTTT'], dct['GGGG']]
sim_genome = np.random.choice(nochange_lst, n_sites)

mut_lst = []
mut_loc = []
for variant in mutated_ts.variants():
    mut_loc.append(variant.site.position)
    mut_lst.append(''.join([variant.alleles[i] for i in variant.genotypes]))

for i in range(len(mut_loc)):
    sim_genome[int(mut_loc[i])] = dct[mut_lst[i]]

E = sim_genome

pd.DataFrame({'E':E}).to_csv('../results/obs_{}_{}_{}_{}.csv'.format(n_int_AB, n_int_ABC, seed, model), index = False)

####################### Optimization #######################

np.random.seed(seed)

# t_init_1 = np.random.normal(t_1, t_1/5)
# t_init_2 = np.random.normal(t_2, t_2/5)
# t_init_upper = np.random.normal(t_upper, t_upper/5)
# N_init = np.random.normal(np.mean([N_AB, N_ABC]), np.mean([N_AB, N_ABC])/5)
# r_init = np.random.normal(r, r/5)
# mu_init = mu

t_1 = t_1*mu
t_2 = t_2*mu
t_upper = t_upper*mu
t_out = t_out*mu
N_AB = N_AB*mu
N_ABC = N_ABC*mu
r = r/mu

transitions, emissions, starting, hidden_states, observed_states = trans_emiss_calc(
    t_1, t_1, t_1, t_2, t_upper, t_out,
    N_AB, N_ABC,
    r, n_int_AB, n_int_ABC)

loglik = loglik_wrapper(transitions, emissions, starting, [E])

write_list([-1, t_1, t_2, t_upper, N_AB, N_ABC, r, loglik, 0], '../results/sim_{}_{}_{}_{}.csv'.format(n_int_AB, n_int_ABC, seed, model))

np.random.seed(seed)

t_init_1 =  np.random.normal(t_1, t_1/5)
t_init_2 =  np.random.normal(t_2, t_2/5)
t_init_upper =  np.random.normal(t_upper, t_upper/5)
N_init_AB =  np.random.normal(N_AB, N_AB/5)
N_init_ABC =  np.random.normal(N_ABC, N_ABC/5)
r_init =  np.random.normal(r, r/5)


dct = {
    't_1':     [t_init_1,     t_init_1/10, t_init_1*10], 
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
    V_lst = [E], 
    res_name = f'../results/sim_{n_int_AB}_{n_int_ABC}_{seed}_{model}.csv', 
    header = False,
    method = 'Nelder-Mead'
    )
