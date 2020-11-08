# DMHB
The Multivariable Harmonic Balance for delay-coupled systems (DEMO)

This library containts all necessary sripts and functions in ordar to reproduce results presented in
related paper.

The two life-scripts are found in the librariry. (DMHB_Supercritical.mlx and DMHB_Subcritical.mlx)
Run these scripts to test one of two networks with different types of a Hopf bifurcation.

In order to use the software add the folder in the MATLAB path.
The demo networks are encoded in scripts, so you can run it to test immediately. Try different zises of the network and parameter multipliers.
In order to use the software with your own network replace the individual dynamics of the node
in the section "individual dynamics" and ODE45 model. 
As output, you get the oscillatory profile including values of frequency, offset, amplitudes, and phases
and comparison of the prediction with results of simulations.  
