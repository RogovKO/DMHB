# **Multivariable Harmonic Balance for Delay-Coupled Systems Toolbox (DMHB Toolbox)**  

The **DMHB Toolbox** provides all the necessary scripts and functions to reproduce the results presented in the associated research paper. This toolbox is designed for analyzing delay-coupled systems, focusing on networks exhibiting Hopf bifurcations.  

## **Features**  
- **Pre-built Demos**:  
  Two live scripts are included in the library:  
  1. `DMHB_Supercritical.mlx`  
  2. `DMHB_Subcritical.mlx`  
  Run these scripts to test networks with supercritical and subcritical Hopf bifurcations.  

- **Network Customization**:  
  - Quickly test the included demo networks or experiment with different network sizes and parameter multipliers.  
  - Easily adapt the toolbox to your own network by modifying the **"Individual Dynamics"** section and the **ODE45 model** in the scripts.  

- **Output Metrics**:  
  - Oscillatory profiles of the system, including:
    - **Frequency**
    - **Offset**
    - **Amplitudes**
    - **Phases**  
  - Comparison of theoretical predictions with simulation results.  

## **Getting Started**  

1. **Add the Toolbox to the MATLAB Path**:  
   - Add the folder containing the toolbox to the MATLAB path using the `addpath` command or the MATLAB interface.  

2. **Run the Demo Scripts**:  
   - Open and run `DMHB_Supercritical.mlx` or `DMHB_Subcritical.mlx` to explore the provided demo networks.  

3. **Customize for Your Network**:  
   - Replace the individual dynamics of the node and update the **ODE45 model** section to adapt the toolbox to your custom system.  

The DMHB Toolbox allows you to explore the behavior of delay-coupled systems efficiently and provides a robust framework for studying Hopf bifurcations in networks.
For more details, see [Research paper](https://research.tue.nl/files/182724857/5.0022610.pdf)


