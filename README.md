# n2pzd1d
N2PZD marine plankton food web model in 1-dimensional, water column setting

README file for N2PZD1D model
June 23, 2025
Katsumi Matsumoto

The N2PZD1D model describes food web interaction between two functional types of phytoplankton and one functional type for zooplankton in a 1-dimensional, upper water column setting. Unique features are flexible phytoplankton C:N:P stoichiometry, homeostatic zooplankton stoichiometry, and grazing based on the stoichiometric imbalance between predator/zooplankton and prey/phytoplankton. N2PZD1D model consists of 10 MATLAB scripts:
•	N2PZD1D.m: main script
•	N2PZD1D_eqs.m: a function called y N2PZD1D.m that contains the ODEs
•	advection.m: calculates particles sinking
•	dailyAverages.m: makes inputs
•	diffprofile.m: calculates the vertical profile of Kz
•	diffusion_explicit.m: calculates the diffusive flux (i.e., vertical mixing)
•	diffusion_i.m: calculates the diffusive flux (i.e., vertical mixing), implicit, preferred
•	ode4.m: ode solver with user-defined dt
•	sincurve.m: makes inputs
•	plot_km.m: makes figures

These scripts have been modified after Sergio Vallina at the Spanish Institute of Oceanography. He developed these scripts as part of his educational/outreach activities.
