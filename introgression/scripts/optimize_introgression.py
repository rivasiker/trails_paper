import sys
from trails.optimizer import trans_emiss_calc_introgression
from trails.cutpoints import cutpoints_ABC, cutpoints_AB
import numpy as np
from trails.optimizer import loglik_wrapper, write_list, optimizer_introgression
import pandas as pd
import time
import re
import msprime
import sys

####################### Model parameters #######################

print('Model parameters')

n_sites = 10_000_000

seed = int(sys.argv[1])
t_A = float(sys.argv[2])
t_B = float(sys.argv[3])
t_C = float(sys.argv[4])
t_2 = float(sys.argv[5])
t_3 = float(sys.argv[6])
N_AB = float(sys.argv[7])
N_ABC = float(sys.argv[8])
r = float(sys.argv[9])
mu = float(sys.argv[10])
n_int_AB = int(sys.argv[11])
n_int_ABC = int(sys.argv[12])
model = str(sys.argv[13])

t_m = float(sys.argv[14])
m = float(sys.argv[15])

t_C_prime = t_C-t_2
t_1 = max([t_A, t_B, t_C_prime])
t_out = t_1+t_2+t_3+2*N_ABC

N_ref = N_ABC

coal_ABC = N_ref/N_ABC
coal_AB = N_ref/N_AB
t_upper = t_3-cutpoints_ABC(n_int_ABC, 1/N_ABC)[-2]
t_AB = t_2/N_ref

cut_AB = t_1+cutpoints_AB(n_int_AB, t_AB, coal_AB)*N_ref
cut_ABC = t_1+t_2+cutpoints_ABC(n_int_ABC, coal_ABC)*N_ref

transitions, emissions, starting, hidden_states, observed_states = trans_emiss_calc_introgression(
    t_A, t_B, t_C, t_2, t_upper, t_out, t_m,
    N_AB, N_ABC, 
    r, m, n_int_AB, n_int_ABC)

dct_hid = {v: k for k, v in hidden_states.items()}
dct = {v: k for k, v in observed_states.items()}

####################### Add demography #######################

print('Simulating demography')

demography = msprime.Demography()
demography.add_population(name="A", initial_size=N_AB, default_sampling_time=t_1-t_A)
demography.add_population(name="B", initial_size=N_AB, default_sampling_time=t_1-t_B)
demography.add_population(name="B_anc", initial_size=N_AB, initially_active=False)
demography.add_population(name="C", initial_size=N_AB, default_sampling_time=t_1+t_2-t_C)
demography.add_population(name="D", initial_size=N_AB, default_sampling_time=t_1-t_1)
demography.add_population(name="AB", initial_size=N_AB)
demography.add_population(name="ABC", initial_size=N_ABC)
demography.add_population(name="ABCD", initial_size=N_ABC)
demography.add_admixture(time = t_1-t_m, derived="B", ancestral=["B_anc", "C"], proportions=(1-m, m))
demography.add_population_split(time=t_1, derived=["A", "B_anc"], ancestral="AB")
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

print('Extracting features')

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
# dat_sim.to_csv('../results/trees_{}_{}_{}_{}.csv'.format(n_int_AB, n_int_ABC, seed, model), index = False)

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

# pd.DataFrame({'E':E}).to_csv('../results/obs_{}_{}_{}_{}.csv'.format(n_int_AB, n_int_ABC, seed, model), index = False)

####################### Optimization #######################

print('Optimizing')

# t_init_1 = np.random.normal(t_1, t_1/5)
# t_init_2 = np.random.normal(t_2, t_2/5)
# t_init_upper = np.random.normal(t_upper, t_upper/5)
# t_init_out = np.random.normal(t_out, t_out/5)
# N_init = np.random.normal(np.mean([N_AB, N_ABC]), np.mean([N_AB, N_ABC])/5)
# r_init = np.random.normal(r, r/5)
# mu_init = mu

t_A = t_A*mu
t_B = t_B*mu
t_C = t_C*mu
t_2 = t_2*mu
t_upper = t_upper*mu
t_out = t_out*mu
t_m = t_m*mu
N_AB = N_AB*mu
N_ABC = N_ABC*mu
r = r/mu

print(f't_m = {t_m}')

transitions, emissions, starting, hidden_states, observed_states = trans_emiss_calc_introgression(
    t_A, t_B, t_C, t_2, t_upper, t_out, t_m,
    N_AB, N_ABC,
    r, m, n_int_AB, n_int_ABC)

loglik = loglik_wrapper(transitions, emissions, starting, [E])

write_list([-1, t_A, t_B, t_C, t_2, t_upper, t_m, N_AB, N_ABC, r, m, loglik, 0], '../results/sim_{}_{}_{}_{}.csv'.format(n_int_AB, n_int_ABC, seed, model))

np.random.seed(seed)
t_init_A = np.random.normal(t_A, t_A/5)
t_init_B = np.random.normal(t_B, t_B/5)
t_init_C = np.random.normal(t_C, t_C/5)
t_init_2 = np.random.normal(t_2, t_2/5)
t_init_upper = np.random.normal(t_upper, t_upper/5)
t_init_m = np.random.normal(t_m, t_m/5)
N_init_AB = np.random.normal(N_AB, N_AB/5)
N_init_ABC = np.random.normal(N_ABC, N_ABC/5)
r_init = np.random.normal(r, r/5)
m_init = m

acc = 1
while (t_init_B/2 <= t_init_m) or ((t_init_C/2-t_init_2/2) <= t_init_m):
    print(acc)
    np.random.seed(seed*10000+acc)
    t_init_A = np.random.normal(t_A, t_A/5)
    t_init_B = np.random.normal(t_B, t_B/5)
    t_init_C = np.random.normal(t_C, t_C/5)
    t_init_2 = np.random.normal(t_2, t_2/5)
    t_init_upper = np.random.normal(t_upper, t_upper/5)
    t_init_m = np.random.normal(t_m, t_m/5)
    N_init_AB = np.random.normal(N_AB, N_AB/5)
    N_init_ABC = np.random.normal(N_ABC, N_ABC/5)
    r_init = np.random.normal(r, r/5)
    m_init = m
    acc += 1

dct = {
    't_A':     [t_init_A,     t_init_A/2, t_init_A*2], 
    't_B':     [t_init_B,     t_init_B/2, t_init_B*2], 
    't_C':     [t_init_C,     t_init_C/2, t_init_C*2], 
    't_2':     [t_init_2,     t_init_2/2, t_init_2*2], 
    't_upper': [t_init_upper, t_init_upper/2, t_init_upper*2], 
    't_m':     [t_init_m,     0, min([t_init_B/2, t_init_2/2+t_init_C/2])], 
    'N_AB':    [N_init_AB,    N_init_AB/2,  N_init_AB*2], 
    'N_ABC':   [N_init_ABC,   N_init_ABC/2,  N_init_ABC*2], 
    'r':       [r_init,       r_init/2,  r_init*2],
    'm':       [m_init,       0.0001,  0.9999]
    }

print(dct)

dct2 = {'n_int_AB':n_int_AB, 'n_int_ABC':n_int_ABC}
res = optimizer_introgression(
    optim_params = dct, 
    fixed_params = dct2, 
    V_lst = [E], 
    res_name = f'../results/sim_{n_int_AB}_{n_int_ABC}_{seed}_{model}.csv', 
    header = False,
    method = "Nelder-Mead"
    )
