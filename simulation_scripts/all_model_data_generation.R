package_vec = c("optparse", "tidyr", "mvtnorm", "devtools", "dplyr")
invisible(suppressPackageStartupMessages(lapply(package_vec, require, character.only = TRUE)))
library(dinemites)

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
        c("--loci"), default = 1,
        type = "integer"),
    make_option(
        c("--multiple_imputation"), default = FALSE,
        type = "logical"),
    make_option(
        c("--synthetic_type"), default = 'bin-present',
        type = "character"),
    make_option(
        c("--iteration"), default = 1,
        type = "integer"))
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments = TRUE)

simulation_type <- unlist(strsplit(opt$options$simulation_type, "-"))

qPCR_only_fun <- function(dataset, rate = 0.5) {
    unique_person_times <- unique(interaction(dataset$time[dataset$present == 1],
                                       dataset$subject[dataset$present == 1]))
    if (length(unique_person_times) > 0) {
        selected_person_times <- sample(unique_person_times, floor(length(unique_person_times) * rate))
        dataset <- dataset %>%
            mutate(actually_present = if ("actually_present" %in% names(.)) {
                actually_present
            } else {
                present
            }) %>%
            mutate(present = ifelse(interaction(time, subject) %in% selected_person_times,
                                    2, present))
    }
    return(dataset)
}

generate_synthetic_data_bin_present <- function(simulation_type,
                                                nsubjects,
                                                nalleles,
                                                nappointments,
                                                total_times,
                                                gene_interaction = TRUE,
                                                qPCR_only = 0.5,
                                                drop_out = 0.2,
                                                seed = 1) {
    set.seed(seed)

    if (total_times < nappointments) {
        stop("total_times must be greater than nappointments")
    }

    allowed_simulation_types <-
        c("persistent",
          "lag_30",
          "season",
          "fixed_covariate",
          "prevention_covariate",
          "time_varying_covariate",
          "treatment")
    if (any(!simulation_type %in% allowed_simulation_types)) {
        stop(paste0(c("simulation_type must only include:", allowed_simulation_types), sep = ' '))
    }

    allele_set <- 1:nalleles
    abeta <- 1
    bbeta <- 50

    if (gene_interaction) {
        infection_multiplier <- 5
        allele_probs <- rbeta(nalleles, abeta / infection_multiplier, bbeta)
    } else {
        allele_probs <- rbeta(nalleles, abeta, bbeta)
    }

    if ("persistent" %in% simulation_type) {
        persistent_probs <- rbeta(nalleles, 1, 5)
    }

    if ("lag_30" %in% simulation_type) {
        lag_probs <- rbeta(nalleles, 1, 1)
    }

    # Season applies to everyone
    if ("season" %in% simulation_type) {
        Sigma <- matrix(0.4, nrow = nalleles, ncol = nalleles)
        diag(Sigma) <- 1
        season_modifier <- c(t(mvtnorm::rmvnorm(n = 1, mean = rep(0, nalleles), Sigma)))
    }

    # Covariates apply to everyone
    if ("fixed_covariate" %in% simulation_type) {
        Sigma <- matrix(0.25, nrow = nalleles, ncol = nalleles)
        diag(Sigma) <- 0.5
        fixed_covariate_modifier <- c(t(mvtnorm::rmvnorm(n = 1, mean = rep(0, nalleles), Sigma)))
        fixed_covariate <- rbinom(nsubjects, 1, 0.5)
    }

    if ("time_varying_covariate" %in% simulation_type) {
        Sigma <- matrix(0.15, nrow = nalleles, ncol = nalleles)
        diag(Sigma) <- 0.5
        varying_covariate_modifier <- c(t(mvtnorm::rmvnorm(n = 1, mean = rep(0, nalleles), Sigma)))
    }

    # Prevention covariate
    if ("prevention_covariate" %in% simulation_type) {
        prevention_covariate <- rbinom(nsubjects, 1, 0.5)
    }

    # dataset needs to have allele, subject, time, new,
    # and optionally: persistent, lag, season, covariate

    # First generate allele, subject, time, new, season, covariate
    # Others can be calculated afterwards
    dataset <- data.frame()
    for (i in 1:nsubjects) {
        time_points <- floor(seq(1, total_times, length.out = sample(max(floor(nappointments * 3/4),  1):max(floor(nappointments * 5/4),  1), 1)))
        time_varying_covariate <- rbinom(time_points, 1, 0.5)

        infected_list <- list()
        final_infected_list <- list()
        novelty_list <- list()
        time_list <- list()
        final_time_list <- list()
        time_varying_covariate_list <- list()
        treatment_list <- list()
        infection_event_list <- list()
        infected_alleles_next_time <- c()
        for (j in seq(length(time_points))) {
            time_point <- time_points[j]

            adjusted_allele_probs <- allele_probs
            if ("season" %in% simulation_type) {
                if (time_point > total_times / 2) {
                    adjusted_allele_probs <- adjusted_allele_probs * exp(season_modifier)
                }
            }

            if ("fixed_covariate" %in% simulation_type) {
                adjusted_allele_probs <- adjusted_allele_probs *
                    exp(fixed_covariate[i] * fixed_covariate_modifier)
            }

            if ("time_varying_covariate" %in% simulation_type) {
                adjusted_allele_probs <- adjusted_allele_probs *
                    exp(varying_covariate_modifier * time_varying_covariate[j])
            }

            if ("prevention_covariate" %in% simulation_type) {
                adjusted_allele_probs <- adjusted_allele_probs *
                    ifelse(prevention_covariate[i] == 1, 1/2, 1)
            }

            adjusted_allele_probs <- pmax(pmin(adjusted_allele_probs, 1), 0)

            if (j > 1 && any(infection_event_list[[j-1]] == 1)) {
                adjusted_allele_probs <- rep(0, length(adjusted_allele_probs))
            }

            # Convert to time-specfic probability and sample
            infected_alleles <- allele_set[rbinom(nalleles, 1, adjusted_allele_probs) == 1]
            n_new_infections <- length(infected_alleles)

            if (gene_interaction & n_new_infections > 0) {
                if ("prevention_covariate" %in% simulation_type) {
                    additional_infections <- min(rpois(1, (n_new_infections * infection_multiplier * ifelse(prevention_covariate[i] == 1, 1/2, 1))),
                                                 nalleles,
                                                 sum(adjusted_allele_probs[!allele_set %in% infected_alleles] > 0))
                } else {
                    additional_infections <- min(rpois(1, (n_new_infections * infection_multiplier)),
                                                 nalleles,
                                                 sum(adjusted_allele_probs[!allele_set %in% infected_alleles] > 0))
                }

                if (additional_infections > 0) {
                    alleles_still_available <- allele_set[!allele_set %in% infected_alleles]
                    infected_alleles <- c(infected_alleles, sample(alleles_still_available, additional_infections,
                                                                   prob = adjusted_allele_probs[!allele_set %in% infected_alleles],
                                                                   replace = F))
                    n_new_infections <- length(infected_alleles)
                }
            }
            infection_event <- ifelse(n_new_infections > 0, 1, 0)

            infected_alleles_next_time_tmp <- infected_alleles_next_time

            infected_alleles_next_time <- infected_alleles[rbinom(n_new_infections, 1, drop_out) == 1]
            infected_alleles <- setdiff(infected_alleles, infected_alleles_next_time)
            infected_alleles <- c(infected_alleles, infected_alleles_next_time_tmp)
            n_new_infections <- length(infected_alleles)

            if ("persistent" %in% simulation_type) {
                if (j > 1) {
                    prev_infections <- sort(unique(unlist(infected_list[1:(j-1)])))
                } else {
                    prev_infections <- c()
                }

                tmp_probs <- rep(0, nalleles)
                tmp_probs[prev_infections] <- persistent_probs[prev_infections]
                infected_alleles <- c(infected_alleles,
                                      allele_set[rbinom(nalleles, 1, tmp_probs) == 1])
            }

            if ("lag_30" %in% simulation_type) {
                if (j > 1) {
                    prev_infections <- sort(unique(unlist(infected_list[1:(j-1)])[
                        unlist(time_list[1:(j-1)]) >= time_point - 30]))
                } else {
                    prev_infections <- c()
                }

                tmp_probs <- rep(0, nalleles)
                tmp_probs[prev_infections] <- lag_probs[prev_infections]
                infected_alleles <- c(infected_alleles,
                                      allele_set[rbinom(nalleles, 1, tmp_probs) == 1])
            }

            infected_alleles <- unique(infected_alleles)
            infected_list[[j]] <- infected_alleles
            final_infected_list[[j]] <- infected_alleles
            novelty_list[[j]] <- c(rep("yes", n_new_infections),
                                   rep("no", length(infected_alleles) - n_new_infections))
            time_list[[j]] <- rep(time_point, length(infected_alleles))
            final_time_list[[j]] <- rep(time_point, length(infected_alleles))
            time_varying_covariate_list[[j]] <- rep(time_varying_covariate[j], length(infected_alleles))
            infection_event_list[[j]] <- infection_event

            # Treatment
            if ("treatment" %in% simulation_type) {
                treatment_efficacy <- 0.9
                is_treated <- length(unique(infected_alleles)) > 0 & runif(1) < 0.1
                if (is_treated) {
                    all_infections <- unique(unlist(infected_list))
                    breakthrough_infections <- all_infections[rbinom(length(all_infections), 1, 1 - treatment_efficacy) == 1]
                    for (list_iter in 1:length(infected_list)) {
                        time_list[[list_iter]] <- time_list[[list_iter]][infected_list[[list_iter]] %in% breakthrough_infections]
                        infected_list[[list_iter]] <- infected_list[[list_iter]][infected_list[[list_iter]] %in% breakthrough_infections]
                    }
                }
                treatment_list[[j]] <- rep(is_treated, length(infected_alleles))
            }
        }

        df <- data.frame("allele" = unlist(final_infected_list),
                         "novelty" = unlist(novelty_list),
                         "time" = unlist(final_time_list),
                         "subject" = rep(i, length(unlist(final_time_list))),
                         "infection_event" = unlist(lapply(1:length(time_points), FUN = function(x){rep(infection_event_list[[x]], length(final_infected_list[[x]]))})))

        extra_times <- setdiff(time_points, unlist(final_time_list))
        extra_df <- data.frame("allele" = rep(0, length(extra_times)),
                               "novelty" = rep("yes", length(extra_times)),
                               "time" = extra_times,
                               "subject" = rep(i, length(extra_times)),
                               "infection_event" = unlist(lapply(which(time_points %in% extra_times), FUN = function(x){rep(infection_event_list[[x]], 1)})))

        if ("season" %in% simulation_type) {
            df$season <- ifelse(df$time > total_times / 2, 1, 0)
            extra_df$season <- ifelse(extra_df$time > total_times / 2, 1, 0)
        }
        if ("fixed_covariate" %in% simulation_type) {
            df$fixed_covariate <- rep(fixed_covariate[i], length(unlist(final_time_list)))
            extra_df$fixed_covariate <- rep(fixed_covariate[i], length(extra_times))
        }
        if ("time_varying_covariate" %in% simulation_type) {
            df$time_varying_covariate <- unlist(time_varying_covariate_list)
            extra_df$time_varying_covariate <- time_varying_covariate[which(!time_points %in% unlist(final_time_list))]
        }
        if ("prevention_covariate" %in% simulation_type) {
            df$prevention_covariate <- rep(prevention_covariate[i], length(unlist(final_time_list)))
            extra_df$prevention_covariate <- rep(prevention_covariate[i], length(extra_times))
        }
        if ("treatment" %in% simulation_type) {
            df$treatment <- unlist(treatment_list)
            extra_df$treatment <- rep(FALSE, nrow(extra_df))
        }

        dataset <- rbind(dataset, df, extra_df)
    }
    dataset$allele <- factor(dataset$allele)

    # Expand the dataframe so that for each subject it contains all combinations
    dataset <- fill_in_dataset(dataset)

    dataset$novelty <- factor(dataset$novelty, levels = c("yes", "no"))

    dataset <- dataset %>%
        dplyr::filter(allele != "0")

    if (qPCR_only > 0) {
        dataset <- qPCR_only_fun(dataset, qPCR_only)
    }

    dataset <- dataset %>%
        dplyr::arrange(time, subject, allele)

    dataset$locus <- 1

    dataset <- dataset %>%
        dplyr::mutate(allele = as.character(allele))

    if ('locus' %in% colnames(dataset)) {
        dataset <- dataset %>%
            dplyr::mutate(locus = as.character(locus))
    }

    return(dataset)
}

generate_synthetic_data_pois_time <- function(simulation_type,
                                              nsubjects,
                                              nalleles,
                                              nappointments,
                                              total_times,
                                              gene_interaction = TRUE,
                                              qPCR_only = 0.5,
                                              drop_out = 0.2,
                                              loci = 1,
                                              seed = 1) {
    set.seed(seed)

    if (total_times < nappointments) {
        stop("total_times must be greater than nappointments")
    }

    allowed_simulation_types <-
        c("persistent",
          "lag_30",
          "season",
          "fixed_covariate",
          "prevention_covariate",
          "time_varying_covariate",
          "treatment")
    if (any(!simulation_type %in% allowed_simulation_types)) {
        stop(paste0(c("simulation_type must only include:", allowed_simulation_types), sep = ' '))
    }

    allele_set <- 1:nalleles
    abeta <- 1
    bbeta <- 50

    loci_corresponding <- rep(1:loci, each = nalleles / loci)

    if (gene_interaction) {
        infection_multiplier <- ifelse(loci == 1, 5, 2)
        allele_probs <- rbeta(nalleles, abeta / infection_multiplier * loci, bbeta)
    } else {
        allele_probs <- rbeta(nalleles, abeta * loci, bbeta)
    }

    # Season applies to everyone
    if ("season" %in% simulation_type) {
        Sigma <- matrix(0.4, nrow = nalleles, ncol = nalleles)
        diag(Sigma) <- 1
        season_modifier <- c(t(mvtnorm::rmvnorm(n = 1, mean = rep(0, nalleles), Sigma)))
    }

    # Covariates apply to everyone
    if ("fixed_covariate" %in% simulation_type) {
        Sigma <- matrix(0.25, nrow = nalleles, ncol = nalleles)
        diag(Sigma) <- 0.5
        fixed_covariate_modifier <- c(t(mvtnorm::rmvnorm(n = 1, mean = rep(0, nalleles), Sigma)))
        fixed_covariate <- rbinom(nsubjects, 1, 0.5)
    }

    if ("time_varying_covariate" %in% simulation_type) {
        Sigma <- matrix(0.15, nrow = nalleles, ncol = nalleles)
        diag(Sigma) <- 0.5
        varying_covariate_modifier <- c(t(mvtnorm::rmvnorm(n = 1, mean = rep(0, nalleles), Sigma)))
    }

    # Prevention covariate
    if ("prevention_covariate" %in% simulation_type) {
        prevention_covariate <- rbinom(nsubjects, 1, 0.5)
    }

    # dataset needs to have allele, subject, time, new,
    # and optionally: persistent, lag, season, covariate

    # First generate allele, subject, time, new, season, covariate
    # Others can be calculated afterwards
    dataset <- data.frame()
    for (i in 1:nsubjects) {
        # Generate appointments
        time_points <- floor(seq(1, total_times, length.out = sample(max(floor(nappointments * 3/4),  1):max(floor(nappointments * 5/4),  1), 1)))
        time_varying_covariate <- rbinom(time_points, 1, 0.5)

        infected_list <- list()
        final_infected_list <- list()
        novelty_list <- list()
        time_list <- list()
        final_time_list <- list()
        time_varying_covariate_list <- list()
        treatment_list <- list()
        infection_event_list <- list()
        infected_alleles_next_time <- c()

        # To store infections
        infection_events <- list()

        for (j in seq(length(time_points))) {
            time_point <- time_points[j]

            adjusted_allele_probs <- allele_probs
            if ("season" %in% simulation_type) {
                if (time_point > total_times / 2) {
                    adjusted_allele_probs <- adjusted_allele_probs * exp(season_modifier)
                }
            }

            if ("fixed_covariate" %in% simulation_type) {
                adjusted_allele_probs <- adjusted_allele_probs *
                    exp(fixed_covariate[i] * fixed_covariate_modifier)
            }

            if ("time_varying_covariate" %in% simulation_type) {
                adjusted_allele_probs <- adjusted_allele_probs * exp(varying_covariate_modifier)
            }

            if ("prevention_covariate" %in% simulation_type) {
                adjusted_allele_probs <- adjusted_allele_probs *
                    ifelse(prevention_covariate[i] == 1, 1/2, 1)
            }

            adjusted_allele_probs <- pmax(pmin(adjusted_allele_probs, 1), 0)

            if (j > 1 && any(infection_event_list[[j-1]] == 1)) {
                adjusted_allele_probs <- rep(0, length(adjusted_allele_probs))
            }
            
            key_locus <- sample(1:loci, 1)
            infected_alleles <- allele_set[rbinom(nalleles, 1, adjusted_allele_probs) == 1 & loci_corresponding == key_locus]
            n_new_infections <- length(infected_alleles)

            if (gene_interaction & n_new_infections > 0) {
                if ("prevention_covariate" %in% simulation_type) {
                    additional_infections <- min(rpois(1, (n_new_infections * infection_multiplier * ifelse(prevention_covariate[i] == 1, 1/2, 1))),
                                                 nalleles,
                                                 sum(adjusted_allele_probs[!allele_set %in% infected_alleles & loci_corresponding == key_locus] > 0))
                } else {
                    additional_infections <- min(rpois(1, (n_new_infections * infection_multiplier)),
                                                 nalleles,
                                                 sum(adjusted_allele_probs[!allele_set %in% infected_alleles & loci_corresponding == key_locus] > 0))
                }

                if (additional_infections > 0) {
                    alleles_still_available <- allele_set[!allele_set %in% infected_alleles & loci_corresponding == key_locus]
                    infected_alleles <- c(infected_alleles, sample(alleles_still_available, additional_infections,
                                                                   prob = adjusted_allele_probs[!allele_set %in% infected_alleles & loci_corresponding == key_locus],
                                                                   replace = F))
                    n_new_infections <- length(infected_alleles)
                }
            }

            infected_alleles_augmented <- c()
            if (loci > 1) {
                for (loci_val in setdiff(1:loci, key_locus)) {
                    if (n_new_infections > 0) {
                        infected_alleles_augmented <- c(infected_alleles_augmented,
                                                        sample(allele_set[loci_corresponding == loci_val], n_new_infections,
                                                               prob = adjusted_allele_probs[loci_corresponding == loci_val], replace = T))
                    }
                }
            }
            infected_alleles <- c(infected_alleles, unique(infected_alleles_augmented))
            n_new_infections <- length(infected_alleles)

            infection_event <- ifelse(n_new_infections > 0, 1, 0)
            infected_alleles_next_time_tmp <- infected_alleles_next_time
            infected_alleles_next_time <- infected_alleles[rbinom(n_new_infections, 1, drop_out) == 1]
            infected_alleles <- setdiff(infected_alleles, infected_alleles_next_time)
            infected_alleles <- c(infected_alleles, infected_alleles_next_time_tmp)
            n_new_infections <- length(infected_alleles)

            if (n_new_infections > 0) {
                if ("lag" %in% simulation_type & "persistent" %in% simulation_type) {
                    if (runif(1) < 0.2) {
                        additional_time_points <- rpois(1, length(time_points))
                    } else {
                        additional_time_points <- rpois(1, 2)
                    }
                } else if ("persistent" %in% simulation_type) {
                    additional_time_points <- rpois(1, length(time_points))
                } else if ("lag" %in% simulation_type) {
                    additional_time_points <- rpois(1, 2)
                }

                if (additional_time_points > 0) {
                    additional_time_points <- min(additional_time_points, max(length(time_points) - additional_time_points, 1))
                    infection_events <- append(infection_events,
                                               list(list(time_points = time_points[j + 1:additional_time_points],
                                                         infected_alleles = infected_alleles)))
                }
            }

            if ("persistent" %in% simulation_type | "lag" %in% simulation_type) {
                infected_alleles_at_time_point <- unlist(
                    lapply(infection_events, function(entry) {
                        if (time_point %in% entry$time_points) {
                            return(entry$infected_alleles)
                        } else {
                            return(NULL)
                        }
                    })
                )

                if (!is.null(infected_alleles_at_time_point)) {
                    tmp_probs <- rep(0, nalleles)
                    tmp_probs[infected_alleles_at_time_point] <- (1 - drop_out)
                    infected_alleles <- c(infected_alleles,
                                          allele_set[rbinom(nalleles, 1, tmp_probs) == 1])
                }
            }

            infected_alleles <- unique(infected_alleles)
            infected_list[[j]] <- infected_alleles
            final_infected_list[[j]] <- infected_alleles
            novelty_list[[j]] <- c(rep("yes", n_new_infections),
                                   rep("no", length(infected_alleles) - n_new_infections))
            time_list[[j]] <- rep(time_point, length(infected_alleles))
            final_time_list[[j]] <- rep(time_point, length(infected_alleles))
            time_varying_covariate_list[[j]] <- rep(time_varying_covariate[j], length(infected_alleles))
            infection_event_list[[j]] <- infection_event

            if ("treatment" %in% simulation_type) {
                treatment_efficacy <- 0.9
                is_treated <- length(unique(infected_alleles)) > 0 & runif(1) < 0.1
                if (is_treated) {
                    for (entry in seq(length(infection_events))) {
                        if (any(infection_events[[entry]]$time_points > time_point, na.rm = TRUE)) {
                            all_infections <- unique(unlist(infection_events[[entry]]$infected_alleles))
                            infection_events[[entry]]$infected_alleles <-
                                all_infections[rbinom(length(all_infections), 1, 1 - treatment_efficacy) == 1]
                        }
                    }
                }
                treatment_list[[j]] <- rep(is_treated, length(infected_alleles))
            }
        }

        df <- data.frame("allele" = unlist(final_infected_list),
                         "novelty" = unlist(novelty_list),
                         "time" = unlist(final_time_list),
                         "subject" = rep(i, length(unlist(final_time_list))),
                         "infection_event" = unlist(lapply(1:length(time_points), FUN = function(x){rep(infection_event_list[[x]], length(final_infected_list[[x]]))})))

        extra_times <- setdiff(time_points, unlist(final_time_list))
        extra_df <- data.frame("allele" = rep(0, length(extra_times)),
                               "novelty" = rep("yes", length(extra_times)),
                               "time" = extra_times,
                               "subject" = rep(i, length(extra_times)),
                               "infection_event" = unlist(lapply(which(time_points %in% extra_times), FUN = function(x){rep(infection_event_list[[x]], 1)})))

        if ("season" %in% simulation_type) {
            df$season <- ifelse(df$time > total_times / 2, 1, 0)
            extra_df$season <- ifelse(extra_df$time > total_times / 2, 1, 0)
        }
        if ("fixed_covariate" %in% simulation_type) {
            df$fixed_covariate <- rep(fixed_covariate[i], length(unlist(time_list)))
            extra_df$fixed_covariate <- rep(fixed_covariate[i], length(extra_times))
        }
        if ("time_varying_covariate" %in% simulation_type) {
            df$time_varying_covariate <- unlist(time_varying_covariate_list)
            extra_df$time_varying_covariate <- time_varying_covariate[which(!time_points %in% unlist(time_list))]
        }
        if ("prevention_covariate" %in% simulation_type) {
            df$prevention_covariate <- rep(prevention_covariate[i], length(unlist(time_list)))
            extra_df$prevention_covariate <- rep(prevention_covariate[i], length(extra_times))
        }
        if ("treatment" %in% simulation_type) {
            df$treatment <- unlist(treatment_list)
            extra_df$treatment <- rep(FALSE, nrow(extra_df))
        }

        dataset <- rbind(dataset, df, extra_df)
    }

    # Expand the dataframe so that for each subject it contains all combinations
    dataset <- fill_in_dataset(dataset)
    dataset$allele <- factor(dataset$allele)

    dataset$novelty <- factor(dataset$novelty, levels = c("yes", "no"))

    dataset <- dataset %>%
        dplyr::filter(allele != "0")

    if (qPCR_only) {
        dataset <- qPCR_only_fun(dataset, qPCR_only)
    }

    dataset <- dataset %>%
        dplyr::arrange(time, subject, allele)

    dataset$locus <- loci_corresponding[as.numeric(as.character(dataset$allele))]

    dataset <- dataset %>%
        dplyr::mutate(allele = as.character(allele))

    if ('locus' %in% colnames(dataset)) {
        dataset <- dataset %>%
            dplyr::mutate(locus = as.character(locus))
    }

    return(dataset)
}

if (opt$options$synthetic_type == 'bin-present') {
    synth_data <- generate_synthetic_data_bin_present(simulation_type,
                                          opt$options$nsubjects,
                                          opt$options$nalleles,
                                          opt$options$nappointments,
                                          opt$options$total_times,
                                          opt$options$gene_interaction,
                                          opt$options$qPCR_only,
                                          opt$options$drop_out,
                                          opt$options$iteration)
} else if (opt$options$synthetic_type == 'pois-time') {
    synth_data <- generate_synthetic_data_pois_time(simulation_type,
                                                    opt$options$nsubjects,
                                                    opt$options$nalleles,
                                                    opt$options$nappointments,
                                                    opt$options$total_times,
                                                    opt$options$gene_interaction,
                                                    opt$options$qPCR_only,
                                                    opt$options$drop_out,
                                                    opt$options$loci,
                                                    opt$options$iteration)
}

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

write.table(synth_data,
            file = paste0(input_directory, "/", file_name, ".tsv"),
            sep = '\t',
            row.names = F)










