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

# Read in sequencing
ampseq_in_original <- read.csv('data/CIS43LS/SIX_1372_genotypes_jan29.txt', sep='\t')
ampseq_in <- ampseq_in_original
ampseq_in <- reshape2::melt(ampseq_in, id.vars = 'sample')
ampseq_in <- ampseq_in %>%
    dplyr::mutate(value = strsplit(ampseq_in$value, '_'),
                  sample = gsub("^(.*?-.+?)-.*$", "\\1", sample) %>%
                      gsub(pattern = 'Enroll', replacement = 'enroll')) %>%
    tidyr::unnest(value) %>%
    dplyr::mutate(allele = gsub('\\:.*', '', value)) %>%
    dplyr::mutate(allele = ifelse(allele == '.', paste0(variable, '_reference'), allele)) %>%
    dplyr::filter(!grepl("Positive", sample)) # Remove controls

# Read in and merge metadata
metadata_1_original <- read.csv('data/CIS43LS/CIS43LS PCR Data Listing_DFExplore_Exported25Mar2024_Harvard.txt', sep = '\t', check.names = F)
metadata_1 <- metadata_1_original
metadata_2_original <- read.csv('data/CIS43LS/CIS43LS Smear Data Listing DFExplore JRS dates 27Mar2024.txt', sep = '\t', check.names = F)
metadata_2 <- metadata_2_original
metadata_2$merge_id_piece_1 <- gsub("MAB-|-", "", metadata_2$EnrolIDcap)
metadata_2$merge_id_piece_2 <- metadata_2$Day %>%
    gsub(pattern = "Day 0 pre", replacement = "Day0") %>%
    gsub(pattern = "Day ", replacement = "Day") %>%
    gsub(pattern = "Enroll", replacement = "enroll") %>%
    gsub(pattern = "Illness Visit ", replacement = "UV")
metadata_2$merge_id <- paste0(metadata_2$merge_id_piece_1, "-", metadata_2$merge_id_piece_2)
metadata_2 <- metadata_2 %>%
    dplyr::group_by(EnrolIDcap) %>%
    dplyr::mutate(first_visit = ifelse(length(VISITDATE[Day == "Day 0 pre"]) > 0, 
                                       VISITDATE[Day == "Day 0 pre"], 
                                       NA)) %>%
    dplyr::ungroup()
metadata_2 <- full_join(metadata_2, metadata_1, by =
                            c('EnrolIDcap' = 'Complete Enrollment ID',
                              'Visit Number' = 'Visit Number',
                              'Complete Screening ID' = 'Complete Screening ID')) %>%
    filter(!is.na(merge_id))

# All samples are in metadata_2, some are in metadata_1, some are in ampseq
metadata_2_for_merge <- metadata_2 %>%
    select(`Visit Number`, Day, merge_id, first_visit, VISITDATE, `6. P.falciparum detected?`)

ampseq_in <- full_join(ampseq_in, metadata_2_for_merge, by = c("sample" = "merge_id")) %>%
    dplyr::mutate(first_visit = as.Date(first_visit, format = "%d%b%Y")) %>%
    dplyr::mutate(VISITDATE = as.Date(VISITDATE, format = "%d%b%Y")) %>%
    dplyr::mutate(time = as.numeric(VISITDATE) - as.numeric(first_visit)) %>%
    dplyr::mutate(subject = gsub("\\-.*", "", sample))

# AL treatments at enrollment time
treat_at_enroll_df <- unique(ampseq_in[grepl("Reenroll|Enroll", ampseq_in$Day), 
                                       c("subject", "time", "Day")]) %>%
    dplyr::filter(!is.na(time)) %>%
    dplyr::filter(time > -16)

# Other treatments during the study
treatments_original <- read.csv('data/CIS43LS/VisitsAE130CMplus 26Nov2024 clean with notes.csv')
treatments <- treatments_original
treatments <- treatments %>%
    dplyr::filter(grepl("ART", DRUGNAME)) %>%
    dplyr::mutate(subject = gsub("MAb\\-|\\-", "", EnrolID),
                  VISITDATE = as.Date(VISITDATE, format = "%d%b%Y")) %>%
    left_join(ampseq_in %>% 
                  select(subject, first_visit, VISITDATE) %>% 
                  unique(), 
              by = c("subject", "VISITDATE")) %>%
    dplyr::select(subject, first_visit, VISITDATE) %>%
    dplyr::mutate(time = as.numeric(VISITDATE) - as.numeric(first_visit)) %>%
    dplyr::select(subject, time) %>%
    dplyr::filter(!is.na(time))
treatments <- rbind(treatments, treat_at_enroll_df[,c("subject", "time")])

# Use metadata to determine PCR-only infections
ampseq_in <- ampseq_in %>%
    mutate(detected = case_when(`6. P.falciparum detected?` == 'Low-positive Non-Pf detected' ~ 'PCR_no_infection',
                                `6. P.falciparum detected?` == 'Low-positive Pf detected' ~ 'PCR_no_infection',
                                `6. P.falciparum detected?` == 'Non-Pf detected' ~ 'PCR_no_infection',
                                `6. P.falciparum detected?` == 'Not detected' ~ 'PCR_no_infection',
                                `6. P.falciparum detected?` == 'Pf detected' ~ 'PCR_infection',
                                `6. P.falciparum detected?` == 'Potentially Non-Pf and Pf mixed' ~ 'PCR_infection',
                                is.na(`6. P.falciparum detected?`) ~ 'no_PCR'))

# Get precipitation data and define rainy season
dir.create('real_data_outputs', showWarnings = F)
if (!file.exists('real_data_outputs/precip_df_CIS43LS.RDS')) {
    url_link <- paste0("https://archive-api.open-meteo.com/v1/era5?latitude=12.9522274&longitude=-8.1756305&start_date=",
                       as.character(min(ampseq_in$VISITDATE)),"&end_date=",as.character(max(ampseq_in$VISITDATE)),"&daily=precipitation_sum")

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

    precip_df <- precip_df %>%
        mutate(year = format(date, "%Y")) %>%
        left_join(rainy_season_dates, by = "year") %>%
        mutate(
            rainy_season = if_else(date >= first_date & date <= last_date, 1, 0, missing = 0)
        ) %>%
        select(-first_date, -last_date, -year)

    saveRDS(precip_df, 'real_data_outputs/precip_df_CIS43LS.RDS')
} else {
    precip_df <- readRDS('real_data_outputs/precip_df_CIS43LS.RDS')
}
ampseq_in <- left_join(ampseq_in, precip_df, by = c("VISITDATE" = "date"))

# Exclude individuals enrolled but with no first visit (these are ambiguous new/persistent infections anyway)
ampseq_in <- ampseq_in %>%
    filter(!is.na(time))
ampseq_in <- ampseq_in %>%
    select(sample, variable, allele, time, subject, detected, rainy_season)

ampseq_in <- ampseq_in %>%
    filter(detected != 'no_PCR' | !is.na(allele)) %>%
    rename(locus = variable) %>%
    mutate(locus = as.character(locus))

if (!(max(rowSums(table(ampseq_in$allele[!is.na(ampseq_in$allele)],
                        ampseq_in$locus[!is.na(ampseq_in$allele)]) > 1)) == 1)) {
    stop("Alleles have the same names in >1 genes")
}

dataset <- fill_in_dataset(ampseq_in)
dataset$locus <- plyr::mapvalues(dataset$allele,
                                      ampseq_in$allele[!is.na(ampseq_in$allele)],
                                      ampseq_in$locus[!is.na(ampseq_in$allele)], warn_missing = F)
write.table(treatments, 'real_data_outputs/treatments_CIS43LS.tsv', sep = '\t', row.names = F)

dataset <- dataset %>%
    add_present_infection() %>%
    group_by(subject, time) %>%
    mutate(present = ifelse(all(present == 0) & detected == 'PCR_infection', 2, present),
           present_infection = ifelse(detected == 'PCR_infection', 1, present_infection)) %>%
    ungroup() %>%
    mutate(initial_phase = ifelse(time >= 1 & time <= 9, 1, 0),
           mid_phase = ifelse(time >= 10 & time <= 35, 1, 0))

# Dataset summary stats
print(length(unique(dataset$subject)))
print(length(unique(dataset$allele)))
print(length(unique(interaction(dataset$subject, dataset$time))))
print(estimate_drop_out(dataset))
print(dataset %>% 
          group_by(subject, time) %>% 
          summarize(qpcr = any(present == 2), present = any(present == 1), .groups = 'drop') %>%
          summarize(sum(qpcr), sum(present)))

# Check that all original samples are included

# Unscheduled visit marked for the same day
if (mean(tolower(ampseq_in$sample) %in% 
         c(tolower(dataset$sample), '0106a-day126', '0230x-uv1')) != 1){
    stop('Preprocessing error: ampseq')
}

if (mean(metadata_2_for_merge$merge_id[
    !is.na(metadata_2_for_merge$first_visit) & 
    !is.na(metadata_2_for_merge$`6. P.falciparum detected?`)] %in% 
    c(dataset$sample, '0106A-Day126', '0230X-UV1')) != 1) {
    stop('Preprocessing error: metadata')
}

if (mean(gsub("MAB-|-", "", treatments_original$EnrollCap) %in%
         c(dataset$subject, "")) != 1) {
    stop("Preprocessing error: treatments")
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

if (!file.exists('real_data_outputs/CIS43LS_probabilities_bayesian.tsv')) {
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
            infection_general_covariates = c("rainy_season", "initial_phase", "mid_phase"),
            alleles_persistence_covariates = c('persistent', 'lag_30', 'lag_60', 'lag_90', 'treatment_acute', 'treatment_longitudinal'))
    
        # Return both probability matrices for this imputation
        list(probability_new = probabilities_out$probability_new, 
             fit_bayesian = probabilities_out$fit, 
             dataset_tmp = dataset_tmp)
    }
    
    for (i in 1:n_imputations) {
        probabilities_out_mat_bayesian[, i] <- results[[i]]$probability_new
    }
    
    saveRDS(results, 'real_data_outputs/models_CIS43LS_bayesian.RDS')
    
    # Stop the cluster
    stopCluster(cl)
    
    estimated_new_infections_bayesian <- 
        estimate_new_infections(dataset_for_fitting, 
                                imputation_mat = imputed_dataset[dataset_for_fitting$rownum,], 
                                probability_mat = probabilities_out_mat_bayesian)
    write.table(estimated_new_infections_bayesian, 'real_data_outputs/CIS43LS_infections_bayesian.tsv', sep = '\t', row.names = T)
    
    dataset_for_fitting_bayesian <- add_probability_new(dataset_for_fitting, probabilities_out_mat_bayesian)
    final_df <- rbind(cbind(dataset_not_for_fitting, probability_new = NA),
                      dataset_for_fitting_bayesian)
    write.table(final_df, 'real_data_outputs/CIS43LS_probabilities_bayesian.tsv', sep = '\t', row.names = F)
}

################
# Simple model #
################

if (!file.exists('real_data_outputs/CIS43LS_probabilities_simple.tsv')) {
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
    
    saveRDS(results, 'real_data_outputs/models_CIS43LS_simple.RDS')

    estimated_new_infections_simple <- 
        estimate_new_infections(dataset_for_fitting, 
                                imputation_mat = imputed_dataset[dataset_for_fitting$rownum,], 
                                probability_mat = probabilities_out_mat_simple)
    write.table(estimated_new_infections_simple, 'real_data_outputs/CIS43LS_infections_simple.tsv', sep = '\t', row.names = T)
    
    dataset_for_fitting_simple <- add_probability_new(dataset_for_fitting, probabilities_out_mat_simple)
    final_df <- rbind(cbind(dataset_not_for_fitting, probability_new = NA),
                      dataset_for_fitting_simple)
    write.table(final_df, 'real_data_outputs/CIS43LS_probabilities_simple.tsv', sep = '\t', row.names = F)
}

#####################
# Clustering method #
#####################

if (!file.exists('real_data_outputs/CIS43LS_probabilities_clustering.tsv')) {
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
    
    saveRDS(results, 'real_data_outputs/models_CIS43LS_clustering.RDS')
    
    estimated_new_infections_clustering <- 
        estimate_new_infections(dataset_for_fitting, 
                                imputation_mat = imputed_dataset[dataset_for_fitting$rownum,], 
                                probability_mat = probabilities_out_mat_clustering)
    write.table(estimated_new_infections_clustering, 'real_data_outputs/CIS43LS_infections_clustering.tsv', sep = '\t', row.names = T)
    
    dataset_for_fitting_clustering <- add_probability_new(dataset_for_fitting, probabilities_out_mat_clustering)
    final_df <- rbind(cbind(dataset_not_for_fitting, probability_new = NA),
                      dataset_for_fitting_clustering)
    write.table(final_df, 'real_data_outputs/CIS43LS_probabilities_clustering.tsv', sep = '\t', row.names = F)
}

##########################
# Bayesian single allele #
##########################

if (!file.exists('real_data_outputs/CIS43LS_probabilities_bayesian-single-locus.tsv')) {
    message("Starting Bayesian single locus section")
    
    unique_alleles_each_locus <- rowSums(table(dataset$locus, dataset$allele) > 0)
    best_locus <- names(unique_alleles_each_locus)[which.max(unique_alleles_each_locus)]
    
    dataset <- dataset %>% filter(locus == best_locus)
    dataset$present_infection <- NULL
    
    dataset <- dataset %>%
        add_present_infection() %>%
        mutate(rownum = seq(nrow(dataset)))
    
    # Create imputations
    n_imputations <- 50
    imputed_dataset <- impute_dataset(dataset, n_imputations, n_cores = floor(n_cores / 2))
    dataset$probability_present <- NULL
    dataset <- add_probability_present(dataset, imputation_mat = imputed_dataset)
    
    dataset_for_fitting <- dataset %>%
        group_by(subject, allele) %>%
        filter(n() >= 3) %>%
        mutate(second_smallest_time = sort(time)[2]) %>%
        filter(time > second_smallest_time) %>%
        select(-second_smallest_time) %>%
        ungroup()
    
    dataset_not_for_fitting <- anti_join(dataset, dataset_for_fitting)
    
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
            infection_general_covariates = c("rainy_season", "initial_phase", "mid_phase"),
            alleles_persistence_covariates = c('persistent', 'lag_30', 'lag_60', 'lag_90', 'treatment_acute', 'treatment_longitudinal'))
        
        # Return both probability matrices for this imputation
        list(probability_new = probabilities_out$probability_new, 
             fit_bayesian = probabilities_out$fit, 
             dataset_tmp = dataset_tmp)
    }
    
    for (i in 1:n_imputations) {
        probabilities_out_mat_bayesian[, i] <- results[[i]]$probability_new
    }
    
    saveRDS(results, 'real_data_outputs/models_CIS43LS_bayesian-single-locus.RDS')
    
    # Stop the cluster
    stopCluster(cl)
    
    estimated_new_infections_bayesian <- 
        estimate_new_infections(dataset_for_fitting, 
                                imputation_mat = imputed_dataset[dataset_for_fitting$rownum,], 
                                probability_mat = probabilities_out_mat_bayesian)
    write.table(estimated_new_infections_bayesian, 'real_data_outputs/CIS43LS_infections_bayesian-single-locus.tsv', sep = '\t', row.names = T)
    
    dataset_for_fitting_bayesian <- add_probability_new(dataset_for_fitting, probabilities_out_mat_bayesian)
    final_df <- rbind(cbind(dataset_not_for_fitting, probability_new = NA),
                      dataset_for_fitting_bayesian)
    write.table(final_df, 'real_data_outputs/CIS43LS_probabilities_bayesian-single-locus.tsv', sep = '\t', row.names = F)
}

