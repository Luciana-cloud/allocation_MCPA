PESTICIDE CARBON ALLOCATION MODEL 

This repository contains a gene-centric model for pesticide degradation in soils affected by different temperature and soil moisture regimes, as well as different initial pesticide concentrations. All related files are written in Matlab.
File Overview

The repository consists of five folders:

1. Data folder. Contains results from a microcosm experiment (three repetitions) of pesticide degradation (2-methyl-4-chlorophenoxyacetic acid - MCPA) under two temperature levels (10°C and 20°C), two soil moisture levels (pF 1.8 and 3.5), and two initial MCPA concentrations (1 and 20 mg/kg). Total experiment duration was 30 days, and mineralization, residual MCPA concentration in soils, related tdfA genes and transcripts were measured. Additionally, carbon use efficiency and the fraction of MCPA-C incorporated into the biomass were measured. 
2. Calibration folder. Contains files for calibration using the gene-centric model and the experimental data from the Data folder. We followed a hierarchical strategy for calibration, and grouped parameters based on assuming i) different bacterial subpopulations under the different MCPA concentrations and ii) physiological and morphological changes under the different soil moisture levels. 
3. Model variants file. Contains the model files organized for the hierarchical calibration. 
4. Calibrated parameters folder. Presents the model's parameter values after calibration, using the simulated annealing tool from Matlab. Only kinetic parameters were calibrated, and sorption parameters were taken from the literature. 
5. Paper_figures file. Shows a visualization of the calibration and simulation results using MATLAB® live scripts.

Contact
In case of questions, please get in touch with Luciana Chavez Rodriguez (lucianagoku@hotmail.com)

License
The modeling code and its documentation itself are licensed under an MIT license (see LICENSE file). The complete experimental dataset can be provided upon request.
