## Phase estimation with GHZ blocks

This repository contains code for numerical simulations of Bayesian phase estimation with blocks of GHZ states.
The preprint can be found here: arxiv.org/abs/2407.06006. 

## Usage

The main code scripts are the following:

- GHZBlocks_adaptive_analytical.m: For computing the Bayesian mean squared error (BMSE) of phase estimation using a partition of GHZ blocks.
- optimize_over_N.m: Finds the BMSE of the optimal quantum interferometer (OQI) for various number of qubits.
- optimize_over_delta_phi.m: Finds the BMSE of the OQI for a fixed number of qubits, for different prior widths.
- wiseman_semianalytic.m: Simulates the protocol outlined in B. L. Higgins et al., New Journal of Physics 11, 073023 (2009).
- sampling_lukin_analytical.m: Simulates the protocol outlined in E. M. Kessler et al., Phys. Rev. Lett. 112, 190403 (2014).
- rosenband.m: Simulates the protocol outlined in T. Rosenband and D. R. Leibrandt, arXiv:1303.6357 (2013).

The rest of the files contain helper methods.
