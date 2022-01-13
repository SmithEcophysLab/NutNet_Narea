# NutNet_Narea

## Repository Description
Data and code for NutNet Narea analyses

## Summary of main files and folders
The folder [analysis](analysis) contains the analysis script, [narea_analysis](narea_analysis), and the folder [optimal_vcmax_R](optimal_vcmax_R), which contains the R functions necessary for calculating optimal vcmax in R.

The script [calc_optimal_vcmax.R](calc_optimal_vcmax.R) contains code for the optimal vcmax function. The function will calculate optimal vcmax, optimal jmax, and the optimal jmax/vcmax ratio. The inputs required are temperature, PAR, VPD, elevation, and an estimate for the quantum efficiency of photosynthetic electron transport, and the curvature of the light response curve of photosynthetic electron transport.

The folder [functions](functions) contains the functions necessary to run the [calc_optimal_vcmax.R](calc_optimal_vcmax.R) script.

The folder [data](data) contains the processed data files: 

- [leaf_traits](leaf_traits)
- [site_metadata](site_metadata)
- [site_traits](site_traits)

## Contact
Any questions or issues can be submitted via GitHub or directed to Nick Smith (nick.smith@ttu.edu).
