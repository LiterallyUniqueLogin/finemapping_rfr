import argparse

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('original_betas')
parser.add_argument('ld')
parser.add_argument('sigma', type=float)
parser.add_argument('num_samples', type=int)
parser.add_argument('n_replicates', type=int)
args = parser.parse_args()

original_betas = np.genfromtxt(args.original_betas).squeeze()
ld = np.genfromtxt(args.ld)

n_betas = original_betas.shape[0]
assert ld.shape == (n_betas, n_betas)

rng = np.random.default_rng(13)

new_betas = rng.multivariate_normal(
    original_betas,
    ld*(args.sigma**2)/args.num_samples,
    size=args.n_replicates
)

assert new_betas.shape == (n_betas, args.n_replicates)
np.savetxt('replicate_betas.tab', new_betas, delimiter='\t')


