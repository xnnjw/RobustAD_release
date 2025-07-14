Robust Activity Detection for Massive Random Access
Paper Summary
This MATLAB code accompanies the paper "Robust Activity Detection for Massive Random Access" by Xinjue Wang, Esa Ollila, and Sergiy A. Vorobyov, available on arXiv: https://arxiv.org/abs/2505.15555. The paper addresses device activity detection in massive machine-type communications (mMTC) for IoT. Traditional methods assume Gaussian noise, which fails under heavy-tailed or impulsive noise. The authors propose robust algorithms, RCWO (Robust Coordinate-Wise Optimization) and RCL-MP (Robust Covariance Learning-based Matching Pursuit), using robust loss functions to improve detection in non-Gaussian noise environments.
Description
This repository contains MATLAB scripts to simulate activity detection in massive MIMO systems under various noise conditions. The scripts perform Monte Carlo simulations to evaluate performance metrics (e.g., probability of missed detection, PMD) for algorithms including RCWO, RCL-MP, CWO, and CL-MP, across different numbers of antennas (M), pilots (L), and active users (K). Results are visualized in plots.
Files

CL_toolbox/: Contains helper functions for covariance learning and signal processing.
algorithms/: Core algorithm implementations, including RCWO and RCL-MP.
utils/: Utility functions for data processing and simulation setup.
Exp1_synthetic.m: Simulates activity detection with synthetic data, varying M and L.
Exp2_Uplink.m: Simulates a cellular massive MIMO uplink scenario under impulsive noise.

Usage

Ensure the folders CL_toolbox/, algorithms/, and utils/ are in your MATLAB path:addpath 'CL_toolbox/SSR_algorithms/'
addpath 'CL_toolbox/'
addpath 'algorithms/'
addpath 'utils/'


Run the scripts:
For synthetic data experiments: Exp1_synthetic.m
For cellular MIMO uplink simulations: Exp2_Uplink.m


The scripts will simulate and plot:
PMD vs. number of antennas (M) for different K (Exp1).
PMD vs. number of pilots (L) for different K (Exp1).
PMD and runtime vs. SNR for a fixed configuration (Exp2).



General Settings

MC_iters: Number of Monte Carlo iterations (default: 20 for Exp1, 100 for Exp2).
P: Total power scale (default: 1).
N: Number of Machine Type Devices (MTDs) (default: 1000).
fading: Uniform fading with limits from -15 to 0 dB.
pilot: Bernoulli pilot sequence.

Simulation Procedure
1. Experiment 1: Synthetic Data (Exp1_synthetic.m)

Over M (Antenna Count):
Fix number of pilots (L = 30).
Vary number of antennas (M = 30:30:90).
Evaluate PMD for different active users (K = [5, 10, 20, 40]).


Over L (Pilot Count):
Fix number of antennas (M = 30).
Vary number of pilots (L = 30:30:90).
Evaluate PMD for different K.



2. Experiment 2: Cellular MIMO Uplink (Exp2_Uplink.m)

Fix configuration: M = 40, L = 40, K = 20, N = 1000.
Vary SNR from -15 to 0 dB.
Evaluate PMD and runtime for multiple algorithms under impulsive t-noise (ν = 2.5).

Plotting

Exp1_synthetic.m: Generates plots of PMD vs. M and PMD vs. L on logarithmic scales, with legends for algorithms (CWO, CL-MP, RCWO, RCL-MP).
Exp2_Uplink.m: Plots PMD and runtime vs. SNR, comparing nine algorithms (VAMP, SBL, CLMP, CWOpt, ASPG, SNIHT, SOMP, RCWO, RCL-MP).

Outputs

Exp1_synthetic.m: Performance matrices (Pmd, Per) and figures for PMD vs. M and L.
Exp2_Uplink.m: Performance matrices (Pmd, Pfa, tme) and figures for PMD and runtime vs. SNR.

Example Command to Run
>> Exp1_synthetic
>> Exp2_Uplink

Citation
@article{wang2024robust,
  title={Robust Activity Detection for Massive Random Access},
  author={Wang, Xinjue and Ollila, Esa and Vorobyov, Sergiy A.},
  journal={arXiv preprint arXiv:2505.15555},
  year={2024}
}

License
MIT License - see the LICENSE file for details.
