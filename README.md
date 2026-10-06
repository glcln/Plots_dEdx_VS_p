# dE/dx versus momentum in the CMS strip tracker

This project computes the mean energy deposit per unit of length (d$E$/d$x$) in the strip modules of the CMS tracker as a function of the particle momentum $p$.

For each track, the d$E$/d$x$ *estimator* is the harmonic mean of power $-2$ of the charge per unit path length $c_i$ (in MeV/cm) of its $N$ strip clusters:

$$I_h = \left( \frac{1}{N} \sum_{i=1}^{N} c_i^{-2} \right)^{-1/2}$$

The charge of the clusters with saturated strips is corrected beforehand. Two corrections are compared, and every output exists once for each:
- `NewCorr`: the shape saturation correction, included here as the git submodule `SaturationCorrection/`;
- `OldCorr`: the former correction, kept in `code_dedx.C` for the comparison.

The plots produced for each correction are:
- the d$E$/d$x$ estimator as a function of the charge sign times the track momentum;
- the d$E$/d$x$ estimator as a function of the track momentum;
- the same plot with the mass lines of the pions, kaons, protons and deuterons, for each of the two parametrisations described [below](#parametrisations-of-the-mass-lines).


# Disclaimer

I wrote all the code in this project myself. However, this README and the comments in the scripts were generated using Claude.


# Environment and installation

The macros are plain ROOT macros. The environment used is `CMSSW_15_0_X`.

The saturation correction is a git submodule, so the repository must be cloned with its submodules:

```bash
git clone --recurse-submodules https://github.com/glcln/EPR-CorrectionAlgorithm
```

In a clone made without this option, run `git submodule update --init`.

The input ntuples (`TTree` named `stage/ttree`) are listed in `Script_dEdx_vs_p.sh`. They are stored locally at Strasbourg: the paths are hard-coded and must be edited to run elsewhere.


# Repository structure

| File | Role |
| --- | --- |
| `Script_dEdx_vs_p.sh` | Runs `code_dedx.C` on the ntuples. |
| `code_dedx.C` | Loop on the events: applies the two saturation corrections, computes the d$E$/d$x$ estimator of each selected track and fills the histograms. |
| `run2analysis.h` | `TTree::MakeClass()` header describing the tree of the ntuples. |
| `doFitOndEdx.C` | Fits the d$E$/d$x$ versus $p$ distribution and extracts the parameters of the mass lines. |
| `doDisplay.C` | Draws and saves the plots. |
| `SaturationCorrection/` | Git submodule with the saturation correction (`CorrFunctions.h` and its templates). |
| `Template_correction` | Symbolic link to `SaturationCorrection/Template_correction`, where `CorrFunctions.h` looks for its templates. |
| `ROOT_histograms/` | Histograms written by step 1. |
| `Results/` | Fit results and plots written by steps 2 and 3. |


# How to run the codes

The three steps are run from the root of the repository, in this order. The directories `ROOT_histograms/` and `Results/` must exist (`mkdir -p ROOT_histograms Results`).

Each macro starts with a `SETTINGS` section holding the file names and the values that can be changed.

## 1. Fill the histograms

```bash
./Script_dEdx_vs_p.sh
```

This runs `code_dedx.C` on the data and writes `ROOT_histograms/dEdx_output.root`, which holds for each correction (`<corr>` = `NewCorr` or `OldCorr`):

| Histogram | Content |
| --- | --- |
| `dEdX0stripVsP_lowp_<corr>` | Estimator versus momentum, below 5 GeV. |
| `dEdX0stripVsP_charge_<corr>` | Estimator versus charge sign times momentum, below 5 GeV. |
| `dEdX0stripVsP_<corr>` | Estimator versus momentum, up to 50 GeV. |
| `dEdX0stripVsP_NewCorr__large` | Estimator versus momentum, up to 4 TeV (new correction only). |

For a quick test on a limited number of events, set `kMaxEntries` in `code_dedx.C`.

## 2. Fit the mass lines

```bash
root -l -b -q doFitOndEdx.C
```

The d$E$/d$x$ versus $p$ distribution is cut in momentum slices. In each slice, the peak of each visible species is fitted with a Landau convoluted with a Gaussian, and its most probable value is taken as the d$E$/d$x$ of the species at that momentum. These values are then used to fit the parameters of the two parametrisations.

The outputs are written in `Results/`:

| File | Content |
| --- | --- |
| `K_C_fit_<corr>.txt` | Parameters $K$ and $C$. |
| `Atlas_fit_<corr>.txt` | Parameters $p_1$ to $p_5$. |
| `dEdx_fit_<corr>.root` | Most probable values versus momentum, fitted function and control canvases. |
| `Atlas_fit_<corr>.pdf`, `Atlas_fit2_<corr>.pdf` | "Atlas" fit: view of the proton points and zoom on the pion points. |

The momentum slices, the starting values and the ranges of the fits are tuned for the data and for the binning of `code_dedx.C`. To check the fit of every slice, set `kSaveSliceFits` to `true`: the fits are then saved in `dEdx_fit_<corr>.root`.

## 3. Draw the plots

```bash
root -l -b -q doDisplay.C
```

The plots are written in `Results/`, each one as `.pdf`, `.C`, `.root` and `.png`:

| Plot | Content |
| --- | --- |
| `dEdX0stripVsP_charge_<corr>` | Estimator versus charge sign times momentum. |
| `dEdX0stripVsP_lowp_<corr>` | Estimator versus momentum. |
| `dEdX0stripVsP_lowp_FIT_<corr>` | Estimator versus momentum, with the "K and C" mass lines. |
| `dEdX0stripVsP_lowp_FIT_A_<corr>` | Estimator versus momentum, with the "Atlas" mass lines. |


# Parametrisations of the mass lines

The mass lines give the d$E$/d$x$ expected for a particle of mass $m$ as a function of its momentum $p$. Two parametrisations are fitted in step 2 and drawn in step 3.

**"K and C"**, with two parameters:

$$\frac{dE}{dx} = K \frac{m^2}{p^2} + C$$

$C$ is fitted on the pions at high momentum, then $K$ on the protons at low momentum.

**"Atlas"**, with five parameters, as a function of $\beta\gamma = p/m$:

$$\frac{dE}{dx} = p_1 X^{p_2/2} \ln\left(1 + (p_3 \beta\gamma)^{p_4}\right) - p_5$$

with

$$X = \frac{\sqrt{(\beta\gamma)^4 + 4 (\beta\gamma)^2} - (\beta\gamma)^2}{2}$$

The five parameters are fitted on the protons at low momentum and the pions at high momentum together.

On the plots, the part of the pion and proton lines labelled "Fit to reference data" is drawn in black, and the rest, labelled "Extrapolation", in red.