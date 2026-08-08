Make sure the following R packages are installed (and up to date) before running simulations:\
- rms
- dplyr
- survey
- foreach
- readxl
- sandwich
- ggplot2

Brief detail of various R files:
- simulationList.xlsx: An excel file that includes various simulations parameters (and fed into sim.R).
- sim.R: Runs and saves results as .rds files to be fed in getFiguresAndTables.R
- sim_supp_A1.R: Runs and saves results as .rds files for Supplemental Figure A1
- sim_prob.R: Runs and saves results as .rds files for Supplemental Table A1
- getFiguresandTables.R: Reads in .rds files and outputs either plots or numerical values used in tables.