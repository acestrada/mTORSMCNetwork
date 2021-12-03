# mTORSMCNetwork
This github repository contains the MATLAB files for a smooth muscle signaling network model focused on mTOR signaling built by Estrada et al. in
"Roles of mTOR in Thoracic Aortopathy Understood by Complex Intracellular Signaling Interactions", PLOS Computational Biology, 2021.

Included files:

Estrada_etal_SMC_mTOR_network.xlsx: Spreadsheet including the network connections and references. It is compatible with Netflux.

SMC_mTOR_ODE.m: MATLAB code of the network system of ordinary differential equations.

SMC_mTOR_ODE_loadParams.m: MATLAB code that assigns parameter values for the network.

Simulation_run.m: MATLAB code that calls both SMC_mTOR_ODE_loadParams.m and SMC_mTOR_ODE.m to run the Baseline, TSC1/2 KO, and Rapamycin simulations
and qualitative comparison included in the manuscript by Estrada et al.