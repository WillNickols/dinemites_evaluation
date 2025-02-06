library(dplyr)
library(doParallel)
library(foreach)
library(dinemites)
library(doRNG)
set.seed(1)
n_cores <- 10

####################
# Data preparation #
####################

# Read in infections
metadata_original <- read.csv('data/uganda/full_meta_6mo.tab', sep='\t')
seq_data_original <- read.csv('data/uganda/jess-ama1_sampInfo.tab', sep='\t')
metadata <- metadata_original %>%
    mutate(merge_id = paste0(date, '-', cohortid)) %>%
    select(cohortid, gender, enrolldate, date, ageyrs, qPCRdich, merge_id, malariacat, malaria) %>%
    mutate(time = as.numeric(as.Date(date) - as.Date(enrolldate)))
seq_data <- seq_data_original %>%
    select(s_Sample, h_popUID)

final_df <- full_join(metadata, seq_data, by = c('merge_id' = 's_Sample'))
final_df <- final_df %>%
    dplyr::rename(subject = cohortid, allele = h_popUID) %>%
    dplyr::select(-merge_id) %>%
    dplyr::filter(!is.na(subject) & !is.na(time) & (!is.na(qPCRdich) | !is.na(allele)))

# Read in treatments
treatments <- metadata %>%
    dplyr::filter(malaria != "no malaria diagnosed today") %>%
    select(cohortid, time) %>%
    rename(subject = cohortid)
dir.create('real_data_outputs', showWarnings = F)
write.table(treatments, 'real_data_outputs/treatments_uganda.tsv', sep = '\t', row.names = F)

# Get precipitation and create season data frame
dir.create('real_data_outputs', showWarnings = F)
if (!file.exists('real_data_outputs/precip_df_uganda.RDS')) {
    url_link <- paste0("https://archive-api.open-meteo.com/v1/era5?latitude=0.7757695437808487&longitude=34.014663992271494&start_date=",
                       as.character(min(final_df$date)),"&end_date=",as.character(max(final_df$date)),"&daily=precipitation_sum")

    precip_list <- jsonlite::fromJSON(url(url_link))
    precip_df <- data.frame(date = as.Date(precip_list$daily$time), precip = precip_list$daily$precipitation_sum)

    rainy_season_dates <- precip_df %>%
        mutate(year = format(date, "%Y")) %>%
        filter(precip > 1) %>%
        group_by(year) %>%
        summarize(
            first_date = min(date),
            last_date = max(date),
            .groups = "drop"
        )

    rainy_season_dates$first_date[rainy_season_dates$year == 2017] <- as.Date("2017-10-04")
    rainy_season_dates$last_date[rainy_season_dates$year == 2017] <- as.Date("2017-12-18")
    rainy_season_dates$first_date[rainy_season_dates$year == 2018] <- as.Date("2018-02-19")
    rainy_season_dates$last_date[rainy_season_dates$year == 2018] <- as.Date("2018-12-31")
    rainy_season_dates$first_date[rainy_season_dates$year == 2019] <- as.Date("2019-03-02")
    rainy_season_dates$last_date[rainy_season_dates$year == 2019] <- as.Date("2019-08-01")

    precip_df <- precip_df %>%
        mutate(year = format(date, "%Y")) %>%
        left_join(rainy_season_dates, by = "year") %>%
        mutate(
            rainy_season = if_else(date >= first_date & date <= last_date, 1, 0, missing = 0)
        ) %>%
        select(-first_date, -last_date, -year)

    saveRDS(precip_df, 'real_data_outputs/precip_df_uganda.RDS')
} else {
    precip_df <- readRDS('real_data_outputs/precip_df_uganda.RDS')
}
final_df <- final_df %>%
    mutate(date = as.Date(date))
final_df <- left_join(final_df, precip_df, by = c("date" = "date"))
dataset <- fill_in_dataset(final_df)

dataset <- dataset %>%
    add_present_infection() %>%
    dplyr::group_by(subject, time) %>%
    dplyr::mutate(present = ifelse(all(present == 0) & qPCRdich == 'Positive', 2, present),
                  present_infection = ifelse(all(present_infection == 0) & qPCRdich == 'Positive', 1, present_infection)) %>%
    dplyr::ungroup()

# Dataset summary stats
print(length(unique(dataset$subject)))
print(length(unique(dataset$allele)))
print(length(unique(interaction(dataset$subject, dataset$time))))
print(estimate_drop_out(dataset))
print(dataset %>% 
          group_by(subject, time) %>% 
          summarize(qpcr = any(present == 2), present = any(present == 1), .groups = 'drop') %>%
          summarize(sum(qpcr), sum(present)))

# Check datasets
if (!all(seq_data_original$s_Sample %in% paste0(dataset$date, "-", dataset$subject) | 
        gsub('.*-', '', seq_data_original$s_Sample) %in% 
        c("3595", "3661", "3604", "3441", "3502", "3860"))) { # These are known to not show up in the metadata
    stop("Preprocessing error: seq_data_original")
}

if (!all(paste0(metadata_original$date, '-', metadata_original$cohortid) %in% 
         paste0(dataset$date, "-", dataset$subject) | is.na(metadata_original$qPCRdich))) { # These are known to not show up in the metadata
    stop("Preprocessing error: metadata_original")
}

# Create imputations
n_imputations <- 50
imputed_dataset <- impute_dataset(dataset, n_imputations, n_cores = floor(n_cores / 2))
dataset <- add_probability_present(dataset, imputation_mat = imputed_dataset)

# Use time points after the first 2 for fitting
dataset <- dataset %>%
    mutate(rownum = seq(nrow(dataset)))

dataset_for_fitting <- dataset %>%
    group_by(subject, allele) %>%
    filter(n() >= 3) %>%
    mutate(second_smallest_time = sort(time)[2]) %>%
    filter(time > second_smallest_time) %>%
    select(-second_smallest_time) %>%
    ungroup()

dataset_not_for_fitting <- anti_join(dataset, dataset_for_fitting)

###################
# Bayesian method #
###################
message("Starting Bayesian section")

# Number of cores to use for parallel processing
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Prepare matrix for output
probabilities_out_mat_bayesian <- matrix(0, nrow = nrow(dataset_for_fitting), ncol = n_imputations)

# Parallel processing with foreach
results <- foreach(i = 1:n_imputations, .packages = c('dplyr', 'dinemites')) %dorng% {
    message("Imputation ", i, " of ", n_imputations)
    dataset_tmp <- dataset
    dataset_tmp$present <- imputed_dataset[,i]

    dataset_tmp <- add_persistent_column(dataset_tmp)
    dataset_tmp <- add_persistent_infection(dataset_tmp)

    dataset_tmp <- add_lag_column(dataset_tmp, lag_time = 30)
    dataset_tmp <- add_lag_infection(dataset_tmp, lag_time = 30)
    dataset_tmp <- add_lag_column(dataset_tmp, lag_time = 60)
    dataset_tmp <- add_lag_infection(dataset_tmp, lag_time = 60)
    dataset_tmp <- add_lag_column(dataset_tmp, lag_time = 90)
    dataset_tmp <- add_lag_infection(dataset_tmp, lag_time = 90)

    dataset_tmp <- add_treatment_column(dataset_tmp, treatments)
    dataset_tmp <- add_treatment_infection(dataset_tmp, treatments)

    dataset_for_fitting_tmp <- dataset_tmp %>%
        group_by(subject, allele) %>%
        filter(n() >= 3) %>%
        mutate(second_smallest_time = sort(time)[2]) %>%
        filter(time > second_smallest_time) %>%
        select(-second_smallest_time) %>%
        ungroup()

    dataset_not_for_fitting_tmp <- anti_join(dataset_tmp, dataset_for_fitting_tmp)

    probabilities_out <- determine_probabilities_bayesian(
        dataset = dataset_for_fitting_tmp,
        infection_persistence_covariates = c('persistent_infection', 'lag_infection_30', 'lag_infection_60',
                                             'lag_infection_90', 'treatment_acute_infection', 'treatment_longitudinal_infection'),
        infection_general_covariates = c("rainy_season"),
        alleles_persistence_covariates = c('persistent', 'lag_30', 'lag_60', 'lag_90', 'treatment_acute', 'treatment_longitudinal'),
        refresh = 1)

    # Return both probability matrices for this imputation
    list(probability_new = probabilities_out$probability_new, 
         fit_bayesian = probabilities_out$fit, 
         dataset_tmp = dataset_tmp)
}

for (i in 1:n_imputations) {
    probabilities_out_mat_bayesian[, i] <- results[[i]]$probability_new
}

saveRDS(results, 'real_data_outputs/models_uganda_bayesian.RDS')

# Stop the cluster
stopCluster(cl)

##################
# Simple version #
##################
message("Starting simple section")

cl <- makeCluster(n_cores)
registerDoParallel(cl)

results <- foreach(i = 1:n_imputations, .packages = c('dplyr', 'dinemites')) %dorng% {
    dataset_tmp <- dataset
    dataset_tmp$present <- imputed_dataset[,i]

    dataset_tmp$probability_new <- determine_probabilities_simple(dataset_tmp)$probability_new
    
    dataset_for_fitting_tmp <- dataset_tmp %>%
        group_by(subject, allele) %>%
        filter(n() >= 3) %>%
        mutate(second_smallest_time = sort(time)[2]) %>%
        filter(time > second_smallest_time) %>%
        select(-second_smallest_time) %>%
        ungroup()

    list(probability_new = dataset_for_fitting_tmp$probability_new, 
         dataset_tmp = dataset_tmp)
}

# Stop the cluster
stopCluster(cl)

probabilities_out_mat_simple <- matrix(0, nrow = nrow(dataset_for_fitting), ncol = n_imputations)
for (i in 1:n_imputations) {
    probabilities_out_mat_simple[, i] <- results[[i]]$probability_new
}

saveRDS(results, 'real_data_outputs/models_uganda_simple.RDS')

##################################
# Wrap-up non-clustering methods #
##################################

estimated_new_infections_bayesian <- 
    estimate_new_infections(dataset_for_fitting, 
                            imputation_mat = imputed_dataset[dataset_for_fitting$rownum,], 
                            probability_mat = probabilities_out_mat_bayesian)
write.table(estimated_new_infections_bayesian, 'real_data_outputs/uganda_infections_bayesian.tsv', sep = '\t', row.names = T)

dataset_for_fitting_bayesian <- add_probability_new(dataset_for_fitting, probabilities_out_mat_bayesian)
final_df <- rbind(cbind(dataset_not_for_fitting, probability_new = NA),
                  dataset_for_fitting_bayesian)
write.table(final_df, 'real_data_outputs/uganda_probabilities_bayesian.tsv', sep = '\t', row.names = F)

estimated_new_infections_simple <- 
    estimate_new_infections(dataset_for_fitting, 
                            imputation_mat = imputed_dataset[dataset_for_fitting$rownum,], 
                            probability_mat = probabilities_out_mat_simple)
write.table(estimated_new_infections_simple, 'real_data_outputs/uganda_infections_simple.tsv', sep = '\t', row.names = T)

dataset_for_fitting_simple <- add_probability_new(dataset_for_fitting, probabilities_out_mat_simple)
final_df <- rbind(cbind(dataset_not_for_fitting, probability_new = NA),
                  dataset_for_fitting_simple)
write.table(final_df, 'real_data_outputs/uganda_probabilities_simple.tsv', sep = '\t', row.names = F)

#####################
# Clustering method #
#####################
message("Starting clustering section")

cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Parallel processing
results <- foreach(i = 1:n_imputations, .packages = c('dplyr', 'dinemites')) %dorng% {
    dataset_tmp <- dataset
    dataset_tmp$present <- imputed_dataset[,i]
    
    dataset_tmp$probability_new <- determine_probabilities_clustering(dataset = dataset_tmp)$probability_new
    
    dataset_for_fitting_tmp <- dataset_tmp %>%
        group_by(subject, allele) %>%
        filter(n() >= 3) %>%
        mutate(second_smallest_time = sort(time)[2]) %>%
        filter(time > second_smallest_time) %>%
        select(-second_smallest_time) %>%
        ungroup()
    
    # Return both probability matrices for this imputation
    list(probability_new = dataset_for_fitting_tmp$probability_new, 
         dataset_tmp = dataset_tmp)
}

stopCluster(cl)

# Prepare matrices for output
probabilities_out_mat_clustering <- matrix(0, nrow = nrow(dataset_for_fitting), ncol = n_imputations)
for (i in 1:n_imputations) {
    probabilities_out_mat_clustering[, i] <- results[[i]]$probability_new
}

saveRDS(results, 'real_data_outputs/models_uganda_clustering.RDS')

estimated_new_infections_clustering <- 
    estimate_new_infections(dataset_for_fitting, 
                            imputation_mat = imputed_dataset[dataset_for_fitting$rownum,], 
                            probability_mat = probabilities_out_mat_clustering)
write.table(estimated_new_infections_clustering, 'real_data_outputs/uganda_infections_clustering.tsv', sep = '\t', row.names = T)

dataset_for_fitting_clustering <- add_probability_new(dataset_for_fitting, probabilities_out_mat_clustering)
final_df <- rbind(cbind(dataset_not_for_fitting, probability_new = NA),
                  dataset_for_fitting_clustering)
write.table(final_df, 'real_data_outputs/uganda_probabilities_clustering.tsv', sep = '\t', row.names = F)
