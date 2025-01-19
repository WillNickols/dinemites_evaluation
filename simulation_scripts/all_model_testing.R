package_vec = c("optparse", "tidyr", "mvtnorm", "devtools", "dplyr")
invisible(suppressPackageStartupMessages(lapply(package_vec, require, character.only = TRUE)))
library(dinemites)
set.seed(1)

# Command Line Usage
option_list = list(
    make_option(
        c("--output"),
        type = "character"),
    make_option(
        c("--simulation_type"), default = "persistent",
        type = "character"),
    make_option(
        c("--nsubjects"), default = 100,
        type = "integer"),
    make_option(
        c("--nalleles"), default = 100,
        type = "integer"),
    make_option(
        c("--nappointments"), default = 10,
        type = "integer"),
    make_option(
        c("--total_times"), default = 200,
        type = "integer"),
    make_option(
        c("--gene_interaction"), default = FALSE,
        type = "logical"),
    make_option(
        c("--qPCR_only"), default = 0,
        type = "double"),
    make_option(
        c("--drop_out"), default = 0.2,
        type = "double"),
    make_option(
        c("--pois_drop_out"), default = 0.2,
        type = "double"),
    make_option(
        c("--loci"), default = 1,
        type = "integer"),
    make_option(
        c("--multiple_imputation"), default = FALSE,
        type = "logical"),
    make_option(
        c("--synthetic_type"), default = 'bin-present',
        type = "character"),
    make_option(
        c("--model"), default = 'bayesian',
        type = "character"),
    make_option(
        c("--iteration"), default = 1,
        type = "integer"))
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments = TRUE)

simulation_type <- unlist(strsplit(opt$options$simulation_type, "-"))

input_directory <- file.path(opt$options$output, 'input')
dir.create(input_directory, showWarnings = F, recursive = T)

file_name<-paste(opt$options$simulation_type, 
                 opt$options$nsubjects, 
                 opt$options$nalleles, 
                 opt$options$nappointments,
                 opt$options$total_times,
                 opt$options$gene_interaction,
                 opt$options$qPCR_only,
                 opt$options$drop_out,
                 opt$options$loci,
                 opt$options$synthetic_type,
                 opt$options$iteration,
                 sep='_')

synth_data <- read.csv(file = paste0(input_directory, "/", file_name, ".tsv"),
                       sep = '\t')

if ('treatment' %in% simulation_type) {
    treatments <- synth_data[,c("subject", "time", "treatment")]
    treatments <- unique(treatments[treatments$treatment,])
    treatments$treatment <- NULL
}

alleles_persistence_covariates <- c("lag_30", "persistent")
alleles_persistence_covariates <- alleles_persistence_covariates[
    alleles_persistence_covariates %in% simulation_type]
infection_general_covariates <- c("season", "fixed_covariate", "time_varying_covariate", "prevention_covariate")
infection_general_covariates <- infection_general_covariates[
    infection_general_covariates %in% simulation_type]

infection_persistence_covariates <- paste0(alleles_persistence_covariates, '_infection')
infection_persistence_covariates[infection_persistence_covariates == 'lag_30_infection'] <- 'lag_infection_30'

if ('treatment' %in% simulation_type) {
    alleles_persistence_covariates <- c(alleles_persistence_covariates, "treatment_acute", "treatment_longitudinal")
    infection_persistence_covariates <- c(infection_persistence_covariates, "treatment_acute_infection", "treatment_longitudinal_infection")
}

if (any(!setdiff(simulation_type, 'treatment') %in% c(alleles_persistence_covariates, infection_general_covariates))) {
    stop("simulation_type not in persistence_covariates or general_covariates")
}
if (any(!colnames(synth_data) %in% c("allele", "novelty", "time", "subject", "present", "actually_present", "prevention_covariate", "locus", "infection_event","treatment",
                                     alleles_persistence_covariates, infection_persistence_covariates, infection_general_covariates))) {
    stop("extra columns, refusing to run")
}
if (!opt$options$model %in% c('bayesian', 'clustering', 'simple')) {
    stop("model must be bayesian, clustering, or simple")
}
if (length(alleles_persistence_covariates) == 0) {
    alleles_persistence_covariates <- NULL
}
if (length(infection_general_covariates) == 0) {
    infection_general_covariates <- NULL
}
if (length(infection_persistence_covariates) == 0) {
    infection_persistence_covariates <- NULL
}

synth_data <- synth_data %>%
    dplyr::mutate(allele = factor(as.character(allele))) %>%
    dplyr::ungroup()

# Run imputation
if (opt$options$multiple_imputation) {
    n_imputations <- 10
    imputed_dataset <- impute_dataset(synth_data, n_imputations)
    
    probability_new_mat <- matrix(ncol = n_imputations, nrow = nrow(synth_data))
    for (i in 1:n_imputations) {
        set.seed(i)
        message(paste0("Working on imputation: ", i))
        synth_data_tmp <- synth_data
        synth_data_tmp$present <- imputed_dataset[,i]
        
        synth_data_tmp <- add_time_gap(synth_data_tmp, default = 0)
        synth_data_tmp <- add_present_infection(synth_data_tmp)
        
        if ('persistent' %in% simulation_type) {
            synth_data_tmp <- add_persistent_column(synth_data_tmp)
            synth_data_tmp <- add_persistent_infection(synth_data_tmp)
        }
        
        if ('lag_30' %in% simulation_type) {
            synth_data_tmp <- add_lag_column(synth_data_tmp, 30)
            synth_data_tmp <- add_lag_infection(synth_data_tmp, 30)
        }
        
        if ('treatment' %in% simulation_type) {
            synth_data_tmp <- add_treatment_column(synth_data_tmp, treatments)
            synth_data_tmp <- add_treatment_infection(synth_data_tmp, treatments)
        }

        if (opt$options$model == 'bayesian') {
            probability_new <- determine_probabilities_bayesian(synth_data_tmp, 
                                              infection_persistence_covariates,
                                              infection_general_covariates,
                                              alleles_persistence_covariates)
        } else if (opt$options$model == 'clustering') {
            probability_new <- determine_probabilities_clustering(synth_data_tmp)
        } else if (opt$options$model == 'simple') {
            probability_new <- determine_probabilities_simple(synth_data_tmp)
        }
        
        probability_new_mat[,i] <- probability_new$probability_new
    }
    
    synth_data <- add_probability_present(synth_data, imputed_dataset)
    synth_data <- add_probability_new(synth_data, probability_new_mat)
    estimated_new_infections <- estimate_new_infections(
        synth_data, imputation_mat = imputed_dataset, probability_mat = probability_new_mat)
} else {
    synth_data <- add_time_gap(synth_data, default = 0)
    synth_data <- add_present_infection(synth_data)
    synth_data$present[synth_data$present == 2] <- 0
    
    if ('persistent' %in% simulation_type) {
        synth_data <- add_persistent_column(synth_data)
        synth_data <- add_persistent_infection(synth_data)
    }
    
    if ('lag_30' %in% simulation_type) {
        synth_data <- add_lag_column(synth_data, 30)
        synth_data <- add_lag_infection(synth_data, 30)
    }
    
    if ('treatment' %in% simulation_type) {
        synth_data <- add_treatment_column(synth_data, treatments)
        synth_data <- add_treatment_infection(synth_data, treatments)
    }

    if (opt$options$model == 'bayesian') {
        probability_new <- determine_probabilities_bayesian(synth_data, 
                                          infection_persistence_covariates,
                                          infection_general_covariates,
                                          alleles_persistence_covariates)
    } else if (opt$options$model == 'clustering') {
        probability_new <- determine_probabilities_clustering(synth_data)
    } else if (opt$options$model == 'simple') {
        probability_new <- determine_probabilities_simple(synth_data)
    }
    
    synth_data$probability_present <- synth_data$present
    synth_data$probability_new <- probability_new$probability_new
    estimated_new_infections <- estimate_new_infections(synth_data)
}

output_directory <- file.path(opt$options$output, 'output')
dir.create(output_directory, showWarnings = F, recursive = T)

file_name<-paste(opt$options$simulation_type, 
                 opt$options$nsubjects, 
                 opt$options$nalleles, 
                 opt$options$nappointments,
                 opt$options$total_times,
                 opt$options$gene_interaction,
                 opt$options$qPCR_only,
                 opt$options$drop_out,
                 opt$options$loci,
                 opt$options$synthetic_type,
                 opt$options$multiple_imputation,
                 opt$options$model,
                 opt$options$iteration, 
                 sep='_')

write.table(synth_data,
            file = paste0(output_directory, "/", file_name, "_probabilities.tsv"),
            sep = '\t',
            row.names = F)

write.table(estimated_new_infections,
            file = paste0(output_directory, "/", file_name, "_infections.tsv"),
            sep = '\t',
            row.names = T)


