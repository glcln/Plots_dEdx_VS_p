# Description of the project

These are the codes used for the computation of the mean energy deposit per unit of length ($dE/dx$) inside strip modules in the CMS tracker as a function of the particle momentum $p$. The mean energy deposit is computed by a harmonic sum of rank 2 on the charge clusters along the track, it is labelled as the $dE/dx$ *estimator*.

Three plots are computed here:
- The $dE/dx$ estimator as a function of the charge sign times the track momentum.
- The $dE/dx$ estimator as a function of the track momentum.
- The $dE/dx$ estimator as a function of the track momentum fitted to identify the mass lines of the pions, kaons, protons and deuterons.

# Environment

> CMSSW_15_0_X/

# How to run the codes

1. Compute the `.root` files on the data of 2025:

> ./Script_dEdx_vs_p.sh

2. Fit the *dE/dx VS p* distribution to compute later the mass lines of the pions, kaons, protons and deuterons:

> root doFitOndEdx.C

3. Display and store the resulting plots:

> root doDisplay.C
