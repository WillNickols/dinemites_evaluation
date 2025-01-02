# novel_infection

Predicting novel malaria plasmodium infections using constrained Bayesian linear models.

## Installation

```
conda env create -f environment.yml
conda env create -f dinemites.yml
```

```
install.packages(c('mvtnorm', 'posterior', 'dplyr', 'tidyr', 'optparse', 'devtools', 'pkgmaker'))
devtools::install_github('WillNickols/dinemites', type = 'source')
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()
```

# Running

```
salloc -p hsph --mem 40000 -t 0-01:00 -c 10 --account=neafsey_lab
conda activate dinemites
python scripts/all_model_testing.py \
  -o /n/holylabs/LABS/neafsey_lab/Lab/wnickols/dinemites_evaluation/qpcr/  \
  --evaluation qpcr --local-jobs 10
  
python scripts/all_model_testing.py \
  -o /n/holylabs/LABS/neafsey_lab/Lab/wnickols/dinemites_evaluation/other/  \
  --evaluation other --local-jobs 10
  
python scripts/all_model_testing.py \
  -o /n/holylabs/LABS/neafsey_lab/Lab/wnickols/dinemites_evaluation/locus/  \
  --evaluation locus --local-jobs 10
  
python scripts/all_model_testing.py \
  -o /n/holylabs/LABS/neafsey_lab/Lab/wnickols/dinemites_evaluation/qpcr/  \
  --evaluation qpcr --grid-jobs 2000 --cores 1 --time 180 --mem 40000 \
  --grid-partition sapphire --grid-options="--account=neafsey_lab"
  
python scripts/all_model_testing.py \
  -o /n/holylabs/LABS/neafsey_lab/Lab/wnickols/dinemites_evaluation/other/  \
  --evaluation other --grid-jobs 2000 --cores 1 --time 180 --mem 40000 \
  --grid-partition sapphire --grid-options="--account=neafsey_lab"
  
python scripts/all_model_testing.py \
  -o /n/holylabs/LABS/neafsey_lab/Lab/wnickols/dinemites_evaluation/locus/  \
  --evaluation locus --grid-jobs 2000 --cores 1 --time 180 --mem 40000 \
  --grid-partition sapphire --grid-options="--account=neafsey_lab"
```