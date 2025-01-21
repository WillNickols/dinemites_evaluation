# DINEMITES evaluation

## Installation

First, clone this GitHub repository, navigate to the cloned folder, and create a conda environment from the `yml` file.
```
conda env create -f dinemites.yml
conda activate dinemites
R
```

Then, in this R, install the following packages:
```
install.packages(c('mvtnorm', 'posterior', 'dplyr', 'tidyr', 'optparse', 'devtools', 'pkgmaker'))
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()
devtools::install_github("WillNickols/dinemites")
```

## Real data analysis

The real data analysis can be run with:
```
Rscript real_data_analysis_scripts/CIS43LS_analysis.R
Rscript real_data_analysis_scripts/mali_analysis.R
Rscript real_data_analysis_scripts/mali_analysis_drop_out.R
Rscript real_data_analysis_scripts/uganda_analysis.R
```

## Synthetic data analysis

The synthetic data analysis can be run with:
```
python simulation_scripts/all_model_testing.py \
  -o dinemites_evaluation/simulation_output/qpcr/  \
  --evaluation qpcr
  
python simulation_scripts/all_model_testing.py \
  -o dinemites_evaluation/simulation_output/other/  \
  --evaluation other
  
python simulation_scripts/all_model_testing.py \
  -o dinemites_evaluation/simulation_output/locus/  \
  --evaluation locus
```

More `--local-jobs` or `--grid-jobs` can be specified if multiple CPUs
or cores are available:
```
python simulation_scripts/all_model_testing.py \
  -o dinemites_evaluation/simulation_output/qpcr/  \
  --evaluation qpcr --grid-jobs 2000 --cores 1 --time 180 --mem 40000 \
  --grid-partition PARTITION_NAME
```