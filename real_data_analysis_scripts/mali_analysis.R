library(dplyr)
library(dinemites)

####################
# Data preparation #
####################
haps_2011 <- read.csv('data/mali/haplotypes_2011.csv')
samples_2011_original <- read.csv('data/mali/samples_2011.csv')
samples_2011 <- samples_2011_original %>%
    filter(visit_date_year == 2011) %>%
    arrange(factor(seq_status, levels = c('seq_pos', 'seq_neg'))) %>%
    distinct(subject_id, visit_date_year, visit_date_month, visit_date_day, .keep_all = TRUE)

haps_joined <- haps_2011 %>% select(extraction_id, locus, haplotype)
samps_joined <- unique(samples_2011)

ampseq_in <- full_join(haps_joined, samps_joined, by = c("extraction_id")) %>%
    dplyr::filter(!is.na(subject_id)) %>%
    dplyr::mutate(visit_date = as.Date(paste0(visit_date_day, '-', visit_date_month,'-', visit_date_year), format = "%d-%m-%Y")) %>%
    dplyr::group_by(subject_id) %>%
    dplyr::mutate(first_visit = min(visit_date)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(time = as.numeric(visit_date) - as.numeric(first_visit)) %>%
    dplyr::rename(subject = subject_id, date = visit_date, allele = haplotype, sample = extraction_id) %>%
    dplyr::mutate(allele = ifelse(allele == '.', paste0(locus, '_reference'), allele))

# Get treatments
treatments_original <- readRDS('data/mali/treatment_metadata_2011.rds')
treatments <- treatments_original %>%
    dplyr::filter(subject_id %in% ampseq_in$subject) %>%
    dplyr::mutate(subject = subject_id,
                  date = as.Date(VisitDate, format = "%d%b%Y")) %>%
    left_join(ampseq_in %>% select(subject, first_visit) %>% unique(), by = c("subject")) %>%
    dplyr::select(subject, first_visit, date) %>%
    dplyr::mutate(time = as.numeric(date) - as.numeric(first_visit)) %>%
    dplyr::select(subject, time)

dir.create('real_data_outputs', showWarnings = F)
write.table(treatments, 'real_data_outputs/treatments_mali.tsv', sep = '\t', row.names = F)

# Get precipitation and set rainy season
if (!file.exists('real_data_outputs/precip_df_mali.RDS')) {
    url_link <- paste0("https://archive-api.open-meteo.com/v1/era5?latitude=12.9522274&longitude=-8.1756305&start_date=",
                       as.character(min(ampseq_in$date)),"&end_date=",as.character(max(ampseq_in$date)),"&daily=precipitation_sum")

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

    rainy_season_dates$first_date[rainy_season_dates$year == 2011] <- as.Date("2011-05-12")
    rainy_season_dates$first_date[rainy_season_dates$year == 2016] <- as.Date("2015-05-14")

    precip_df <- precip_df %>%
        mutate(year = format(date, "%Y")) %>%
        left_join(rainy_season_dates, by = "year") %>%
        mutate(
            rainy_season = if_else(date >= first_date & date <= last_date, 1, 0, missing = 0)
        ) %>%
        select(-first_date, -last_date, -year)

    saveRDS(precip_df, 'real_data_outputs/precip_df_mali.RDS')
} else {
    precip_df <- readRDS('real_data_outputs/precip_df_mali.RDS')
}

ampseq_in <- left_join(ampseq_in, precip_df, by = c("date" = "date"))
ampseq_in <- ampseq_in %>%
    select(sample, locus, allele, time, subject, seq_status, rainy_season)

dataset <- fill_in_dataset(ampseq_in)

if (!(max(rowSums(table(ampseq_in$allele[!is.na(ampseq_in$allele)],
                       ampseq_in$locus[!is.na(ampseq_in$allele)]) > 1)) == 1)) {
    stop("Alleles have the same names in >1 genes")
}

dataset$locus <- plyr::mapvalues(dataset$allele,
                                ampseq_in$allele[!is.na(ampseq_in$allele)],
                                ampseq_in$locus[!is.na(ampseq_in$allele)], 
                                warn_missing = F)

dataset <- dataset %>%
    add_present_infection() %>%
    select(-seq_status)

# Dataset summary stats
print(length(unique(dataset$subject)))
print(length(unique(dataset$allele)))
print(length(unique(interaction(dataset$subject, dataset$time))))
print(estimate_drop_out(dataset))

# PST053_F09 is marked sequencing positive, but there are no sequencing results for it.
# All others are accounted for.

if (!all(samples_2011$extraction_id %in% dataset$sample)) {
    stop("Preprocessing error: samples_2011")
}

samples_2012 <- read.csv('data/mali/samples_2012_2016.csv')
if (!all(haps_2011$extraction_id %in% 
         c(samples_2011_original$extraction_id, 
           gsub("-", "_", samples_2012$extraction_id)))) {
    stop("Preprocessing error: haps2011")
}

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

if (!file.exists('real_data_outputs/mali_probabilities_bayesian.tsv')) {
    message('Running full Bayesian model')
    dataset_tmp <- dataset
    
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
        infection_general_covariates = c('rainy_season'),
        alleles_persistence_covariates = c('persistent', 'lag_30', 'lag_60', 'lag_90', 'treatment_acute', 'treatment_longitudinal'),
        drop_out = F)
    
    probabilities_bayesian <- probabilities_out$probability_new
    results <- list(probability_new = probabilities_out$probability_new, 
                    fit_bayesian = probabilities_out$fit, 
                    dataset_tmp = dataset_tmp)
    saveRDS(results, 'real_data_outputs/models_mali_bayesian.RDS')
    
    dataset_for_fitting_bayesian <- dataset_for_fitting
    dataset_for_fitting_bayesian$probability_new <- probabilities_bayesian
    
    estimated_new_infections_bayesian <- estimate_new_infections(dataset_for_fitting_bayesian)
    write.table(estimated_new_infections_bayesian, 'real_data_outputs/mali_infections_bayesian.tsv', sep = '\t', row.names = T)
    
    final_df <- rbind(cbind(dataset_not_for_fitting, probability_new = NA), dataset_for_fitting_bayesian)
    write.table(final_df, 'real_data_outputs/mali_probabilities_bayesian.tsv', sep = '\t', row.names = F)
}

##################
# Simple version #
##################

if (!file.exists('real_data_outputs/mali_probabilities_simple.tsv')) {
    message("Starting simple section")
    dataset_tmp <- dataset
    
    # Run simple division
    dataset_tmp$probability_new <- determine_probabilities_simple(dataset_tmp)$probability_new
    
    dataset_for_fitting_tmp <- dataset_tmp %>%
        group_by(subject, allele) %>%
        filter(n() >= 3) %>%
        mutate(second_smallest_time = sort(time)[2]) %>%
        filter(time > second_smallest_time) %>%
        select(-second_smallest_time) %>%
        ungroup()
    
    probabilities_simple <- dataset_for_fitting_tmp$probability_new
    results <- list(probability_new = dataset_for_fitting_tmp$probability_new, 
                    dataset_tmp = dataset_tmp)
    
    probabilities_simple = dataset_for_fitting_tmp$probability_new
    
    saveRDS(results, 'real_data_outputs/models_mali_simple.RDS')
    
    dataset_for_fitting_simple <- dataset_for_fitting
    dataset_for_fitting_simple$probability_new <- probabilities_simple
    
    estimated_new_infections_simple <- estimate_new_infections(dataset_for_fitting_simple)
    write.table(estimated_new_infections_simple, 'real_data_outputs/mali_infections_simple.tsv', sep = '\t', row.names = T)
    
    final_df <- rbind(cbind(dataset_not_for_fitting, probability_new = NA), dataset_for_fitting_simple)
    write.table(final_df, 'real_data_outputs/mali_probabilities_simple.tsv', sep = '\t', row.names = F)
}

#####################
# Clustering method #
#####################

if (!file.exists('real_data_outputs/mali_probabilities_clustering.tsv')) {
    message("Starting clustering section")
    dataset_tmp <- dataset
    
    probabilities_out <- determine_probabilities_clustering(dataset = dataset_tmp)
    probabilities_clustering <- probabilities_out$probability_new
    
    results <- list(probability_new = probabilities_clustering, 
                    dataset_tmp = dataset_tmp)
    
    saveRDS(results, 'real_data_outputs/models_mali_clustering.RDS')
    
    dataset_for_fitting_clustering <- dataset_for_fitting
    dataset_for_fitting_clustering$probability_new <- probabilities_clustering[dataset_for_fitting$rownum]
    
    estimated_new_infections_clustering <- estimate_new_infections(dataset_for_fitting_clustering)
    write.table(estimated_new_infections_clustering, 'real_data_outputs/mali_infections_clustering.tsv', sep = '\t', row.names = T)
    
    final_df <- rbind(cbind(dataset_not_for_fitting, probability_new = NA), dataset_for_fitting_clustering)
    write.table(final_df, 'real_data_outputs/mali_probabilities_clustering.tsv', sep = '\t', row.names = F)
}

##########################
# Bayesian single allele #
##########################

if (!file.exists('real_data_outputs/mali_probabilities_bayesian-single-locus.tsv')) {
    message('Running Bayesian model with one locus')
    
    unique_alleles_each_locus <- rowSums(table(dataset$locus, dataset$allele) > 0)
    best_locus <- names(unique_alleles_each_locus)[which.max(unique_alleles_each_locus)]
    
    dataset <- dataset %>% filter(locus == best_locus)
    dataset$present_infection <- NULL
    
    dataset <- dataset %>%
        add_present_infection() %>%
        mutate(rownum = seq(nrow(dataset)))
    
    dataset_for_fitting <- dataset %>%
        group_by(subject, allele) %>%
        filter(n() >= 3) %>%
        mutate(second_smallest_time = sort(time)[2]) %>%
        filter(time > second_smallest_time) %>%
        select(-second_smallest_time) %>%
        ungroup()
    
    dataset_not_for_fitting <- anti_join(dataset, dataset_for_fitting)
    
    dataset_tmp <- dataset
    
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
        infection_general_covariates = c('rainy_season'),
        alleles_persistence_covariates = c('persistent', 'lag_30', 'lag_60', 'lag_90', 'treatment_acute', 'treatment_longitudinal'),
        drop_out = F)
    
    probabilities_bayesian <- probabilities_out$probability_new
    results <- list(probability_new = probabilities_out$probability_new, 
                    fit_bayesian = probabilities_out$fit, 
                    dataset_tmp = dataset_tmp)
    saveRDS(results, 'real_data_outputs/models_mali_bayesian-single-locus.RDS')
    
    dataset_for_fitting_bayesian <- dataset_for_fitting
    dataset_for_fitting_bayesian$probability_new <- probabilities_bayesian
    
    estimated_new_infections_bayesian <- estimate_new_infections(dataset_for_fitting_bayesian)
    write.table(estimated_new_infections_bayesian, 'real_data_outputs/mali_infections_bayesian-single-locus.tsv', sep = '\t', row.names = T)
    
    final_df <- rbind(cbind(dataset_not_for_fitting, probability_new = NA), dataset_for_fitting_bayesian)
    write.table(final_df, 'real_data_outputs/mali_probabilities_bayesian-single-locus.tsv', sep = '\t', row.names = F)
}




