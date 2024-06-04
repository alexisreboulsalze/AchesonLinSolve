# AchesonLinSolve
Solve numerically the dispersion relation of Acheson 1978 to study the stability of an axisymmetric toroidal magnetic field in stars. In particular, it solves it for a hypermassive neutron star after a binary neutron merger. 

Necessary Python Pmodules: SymPy, NumPy, Matplotlib

Manual: 
- Run derive_equations.py to create python files with the polynomial coefficients of Acheson dispersion relation under different hypotheses (Original, Spruit 1999 or diffusionless Acheson with a polytropic gas)
- Run Dissipation_impact.py to Reproduce Figure 8 of Reboul-Salze et al. 2024

Caveat: 
With numerical solving, it is important to test the results against theoretical predictions as some roots may be wrong due to catastrophic substractions.
