library(ggplot2)
library(dplyr)
library(gridExtra)
library(parallel)
library(doParallel)
library(dinemites)

parse_parameters <- function(lines) {
    param_list <- list()

    for (line in lines) {
        # Split by ':'
        split_line <- strsplit(line, ':', fixed = TRUE)[[1]]
        key <- trimws(split_line[1])
        value <- trimws(split_line[c(2,3)])

        # Further split value by ':' to separate method and options
        combine_method <- value[1]
        options <- strsplit(trimws(value[2]), ' ')[[1]]  # Options after the last colon

        # Store the options in the list
        param_list[[key]] <- list(method = combine_method, options = options)
    }

    return(param_list)
}

generate_combinations <- function(param_list) {
    combinations <- list()

    one_options <- list()
    all_options <- list()

    for (key in names(param_list)) {
        method <- param_list[[key]]$method
        options <- param_list[[key]]$options

        if (method == 'one') {
            # Use the first option by default
            one_options[[key]] <- options[1]  # Only take the first option
        } else if (method == 'all') {
            # Use all options
            all_options[[key]] <- options
        }
    }

    # If there are no 'all' parameters, return just one combination
    if (length(all_options) == 0) {
        combinations <- append(combinations, list(one_options))
        return(combinations)
    }

    # Generate combinations by iterating through 'all' parameters
    all_keys <- names(all_options)
    all_values <- expand.grid(all_options, stringsAsFactors = FALSE)

    # For each combination of 'all' parameters
    for (row_idx in 1:nrow(all_values)) {
        for (key in names(one_options)) {
            # Hold other "one" parameters constant and vary the current one
            current_combination <- one_options
            current_combination[[key]] <- one_options[[key]]  # Hold the first option constant

            # Iterate through all options for each 'one' parameter
            for (option in param_list[[key]]$options) {
                current_combination[[key]] <- option  # Change to the current option

                # Add combinations of 'all' parameters
                for (i in seq_along(all_keys)) {
                    all_key <- all_keys[i]
                    current_combination[[all_key]] <- all_values[row_idx, i]
                }

                combinations <- append(combinations, list(current_combination))
            }
        }
    }

    return(combinations)
}

combinations_to_dataframe <- function(combinations) {
    # Convert list of combinations (each as a named list) into a dataframe
    df <- do.call(rbind, lapply(combinations, as.data.frame))

    # Ensure it's a proper dataframe
    df <- as.data.frame(df, stringsAsFactors = FALSE)

    df <- df[!duplicated(df),]

    return(df)
}

parameters <- 'simulation_scripts/config/locus_testing.txt'
parameter_text <- readLines(parameters)
parsed_params <- parse_parameters(parameter_text)
params_for_files <- generate_combinations(parsed_params)
params_df <- combinations_to_dataframe(params_for_files)
params_df <- params_df[names(parsed_params)]
run_only_args <- c('multiple_imputation', 'model')
params_df <- params_df[c(names(parsed_params)[!names(parsed_params) %in% run_only_args], run_only_args)]

rownames(params_df) <- paste0("sim", 1:nrow(params_df))

input_path <- 'locus/output/'
files_in <- list.files(input_path, full.names = TRUE)

prepare_error <- function(row_num) {
    expanded_df <- tryCatch({
        inputString <- paste0(unlist(c(params_df[row_num,])), collapse = '_')

        if (!(length(files_in[grepl(inputString, files_in)]) > 1)) {
            return(NULL)
        }
        print(cbind(params_df[row_num,], 'count' = length(files_in[grepl(inputString, files_in)])))

        total_input <- data.frame()
        # total_new_infections <- data.frame()
        total_new_COI <- data.frame()
        for (file_in in files_in[grepl(inputString, files_in)]) {
            example_input <- read.csv(file_in, sep = '\t')
            if (!'actually_present' %in% colnames(example_input)) {
                example_input$actually_present <- example_input$present
            }

            if ('probability_present' %in% colnames(example_input)) {
                if (2 %in% example_input$probability_present) {
                    example_input$probability_present <- 1 * (example_input$probability_present == 1)
                }
            }

            example_input$replicate <- as.numeric(gsub(".*_|\\.tsv", "", file_in))

            new_infections <- example_input %>%
                dplyr::group_by(subject, replicate, time) %>%
                dplyr::summarise(new_infections_tmp = 1 * any(infection_event == 1), .groups = 'drop') %>%
                dplyr::group_by(subject, replicate) %>%
                dplyr::summarise(new_infections_true = sum(new_infections_tmp), .groups = 'drop')

            # estimated_new_infections <- estimate_new_infections(example_input)

            # new_infection_comparison <- full_join(new_infections, estimated_new_infections, by = c("subject"))

            total_new_COI_stm <- compute_total_new_COI(example_input, method = 'sum_then_max')
            colnames(total_new_COI_stm) <- c("subject", "stm")
            total_new_COI_mts <- compute_total_new_COI(example_input, method = 'max_then_sum')
            colnames(total_new_COI_mts) <- c("subject", "mts")

            this_new_COI <- full_join(total_new_COI_stm, total_new_COI_mts, by = c('subject'))
            this_new_COI$replicate <- as.numeric(gsub(".*_|\\.tsv", "", file_in))

            total_input <- rbind(total_input, example_input)
            # total_new_infections <- rbind(total_new_infections, new_infection_comparison)
            total_new_COI <- rbind(total_new_COI, this_new_COI)
        }

        binned_props <- total_input %>%
            dplyr::filter(present == 1) %>%
            dplyr::filter(probability_new <= 1) %>%
            dplyr::mutate(probability_bin = cut(probability_new, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)) %>%
            dplyr::group_by(probability_bin, novelty) %>%
            dplyr::summarise(count = n(), .groups = 'drop') %>%
            dplyr::group_by(probability_bin) %>%
            dplyr::mutate(total_count = sum(count)) %>%
            dplyr::mutate(proportion = count / total_count) %>%
            dplyr::filter(novelty == 'yes') %>%
            dplyr::group_by(novelty, probability_bin) %>%
            dplyr::summarise(avg_proportion = mean(proportion, na.rm = TRUE),
                  count = total_count)

        levels_to_include <- levels(cut((0:20)/10, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE))

        bool_indicator <- !levels_to_include %in% binned_props$probability_bin
        binned_props <- rbind(binned_props,
                              data.frame(probability_bin = levels_to_include[bool_indicator],
                                         novelty = rep('yes', sum(bool_indicator)),
                                         'avg_proportion' = rep(0, sum(bool_indicator)),
                                         'count' = rep(0, sum(bool_indicator))))
        binned_props$probability_bin <- factor(binned_props$probability_bin, levels =
                                                   levels(cut((0:20)/10, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)))

        binned_props <- binned_props %>%
            dplyr::mutate(error = 2 * sqrt(avg_proportion * (1 - avg_proportion) / count))
        binned_props <- binned_props %>%
            dplyr::group_by(probability_bin) %>%
            dplyr::mutate(ub = ifelse(is.na(error), NA, binom.test(x = avg_proportion * count, n = count, conf.level = 0.95)$conf.int[2]),
                   lb = ifelse(is.na(error), NA, binom.test(x = avg_proportion * count, n = count, conf.level = 0.95)$conf.int[1])) %>%
            dplyr::mutate(color = ifelse(count > 10 & (ub - lb) < 0.25, 'black', 'black'))

        true_new_COI_tmp <- total_input %>%
            dplyr::group_by(subject, replicate, locus, time) %>%
            dplyr::summarise(true_novelty_imputed = sum(actually_present == 1 & novelty == 'yes'),
                             infection_event = any(infection_event == 1)) %>%
            dplyr::ungroup()

        # Use the fact that infection events and time point immediately afterwards
        # constitute the same infection to be de-duplicated for new_COI
        true_new_COI <- true_new_COI_tmp %>%
            group_by(subject, replicate, locus) %>%
            arrange(time) %>%
            mutate(true_novelty_imputed = if_else(infection_event,
                                               lag(true_novelty_imputed, default = 0) + true_novelty_imputed,
                                               0)) %>%
            ungroup() %>%
            filter(infection_event == 1)

        true_new_COI <- true_new_COI %>%
            dplyr::group_by(subject, replicate, time) %>%
            dplyr::summarize(new_COI_tmp = max(true_novelty_imputed)) %>%
            dplyr::group_by(subject, replicate) %>%
            dplyr::summarize(new_COI = sum(new_COI_tmp))

        total_new_COI <- full_join(true_new_COI, total_new_COI, by = c("subject", "replicate")) %>%
            mutate(new_COI = ifelse(is.na(new_COI), 0, new_COI))

        if ('actually_present' %in% colnames(total_input)) {
            total_input <- total_input %>%
            dplyr::group_by(subject, replicate) %>%
            dplyr::summarise(true_novelty_observed = sum(actually_present[present == 1] == 1 & novelty[present == 1] == 'yes'),
                             predicted_novelty_observed = sum(probability_new[present == 1] * probability_present[present == 1]),
                             true_novelty_imputed = sum(actually_present == 1 & novelty == 'yes'),
                             predicted_novelty_imputed = sum(probability_new * probability_present, na.rm=T),
                             prevention_covariate = ifelse('prevention_covariate' %in% colnames(total_input),
                                                           prevention_covariate[1], NA))
        } else {
            total_input <- total_input %>%
            dplyr::group_by(subject, replicate) %>%
            dplyr::summarise(true_novelty_observed = sum(present[present == 1] == 1 & novelty[present == 1] == 'yes'),
                             predicted_novelty_observed = sum(probability_new[present == 1] * probability_present[present == 1]),
                             true_novelty_imputed = sum(present == 1 & novelty == 'yes'),
                             predicted_novelty_imputed = sum(probability_new * probability_present),
                             prevention_covariate = ifelse('prevention_covariate' %in% colnames(total_input),
                                                           prevention_covariate[1], NA))
        }

        binned_props <- binned_props %>%
            mutate(int_lower = as.numeric(gsub("\\[|\\(|,.*", "", probability_bin)),
                   int_upper = as.numeric(gsub("\\]|\\)|.*,", "", probability_bin))) %>%
            mutate(int_mid = (int_lower + int_upper) / 2)

        total_weighted_probability_error <- sum((binned_props$int_mid - binned_props$avg_proportion) * binned_props$count) /
                                sum(binned_props$count)
        total_weighted_probability_abs_error <- sum(abs(binned_props$int_mid - binned_props$avg_proportion) * binned_props$count) /
                                sum(binned_props$count)

        total_new_infections_error <- mean((total_new_infections$new_infections - total_new_infections$new_infections_true))
        total_new_infections_abs_error <- mean(abs(total_new_infections$new_infections - total_new_infections$new_infections_true))
        total_new_infections_sd <- sd((total_new_infections$new_infections - total_new_infections$new_infections_true))
        total_new_infections_abs_sd <- sd(abs(total_new_infections$new_infections - total_new_infections$new_infections_true))

        total_new_COI_error_stm <- mean(total_new_COI$stm - total_new_COI$new_COI)
        total_new_COI_error_sd_stm <- sd(total_new_COI$stm - total_new_COI$new_COI)
        total_new_COI_abs_error_stm <- mean(abs(total_new_COI$stm - total_new_COI$new_COI))
        total_new_COI_abs_error_sd_stm <- sd(abs(total_new_COI$stm - total_new_COI$new_COI))
        total_new_COI_error_mts <- mean(total_new_COI$mts - total_new_COI$new_COI)
        total_new_COI_error_sd_mts <- sd(total_new_COI$mts - total_new_COI$new_COI)
        total_new_COI_abs_error_mts <- mean(abs(total_new_COI$mts - total_new_COI$new_COI))
        total_new_COI_abs_error_sd_mts <- sd(abs(total_new_COI$mts - total_new_COI$new_COI))

        expanded_df <- data.frame(params_df[row_num,],
        total_weighted_probability_error,
        total_weighted_probability_abs_error,
        total_new_COI_error_stm,
        total_new_COI_error_sd_stm,
        total_new_COI_abs_error_stm,
        total_new_COI_abs_error_sd_stm,
        total_new_COI_error_mts,
        total_new_COI_error_sd_mts,
        total_new_COI_abs_error_mts,
        total_new_COI_abs_error_sd_mts,
        total_new_infections_error,
        total_new_infections_abs_error,
        total_new_infections_sd,
        total_new_infections_abs_sd)

        expanded_df
    }, error = function(e){e})
    return(expanded_df)
}

if (!file.exists('evaluation_intermediates/locus_evaluation.RDS')) {
    # num_cores <- 1
    # cl <- makeCluster(num_cores)
    # registerDoParallel(cl)

    growing_df <- list()
    for (row_num in 1:nrow(params_df)) {
        suppressMessages(growing_df[[row_num]] <- prepare_error(row_num))
    }

    # stopCluster(cl)

    growing_df_save <- dplyr::bind_rows(growing_df)
    saveRDS(growing_df_save, file = 'evaluation_intermediates/locus_evaluation.RDS')
} else {
    growing_df_save <- readRDS('evaluation_intermediates/locus_evaluation.RDS')
}

growing_df <- growing_df_save
growing_df <- growing_df %>%
    filter(model != 'allele-specific') %>%
    mutate(model = case_when(model == 'bayesian' ~ 'Bayesian',
                             model == 'clustering' ~ 'Clustering',
                             model == 'simple' ~ 'Simple'),
           multiple_imputation = ifelse(multiple_imputation, 'Yes', 'No')) %>%
    mutate(synthetic_type = case_when(synthetic_type == 'bin-present' ~ 'Rolling presence probability',
                                      synthetic_type == 'pois-time' ~ 'Poisson time to clearance'))

growing_df <- reshape2::melt(growing_df, measure.vars =
                                 c("total_new_COI_error_stm",
                                 "total_new_COI_error_sd_stm",
                                 "total_new_COI_abs_error_stm",
                                 "total_new_COI_abs_error_sd_stm",
                                 "total_new_COI_error_mts",
                                 "total_new_COI_error_sd_mts",
                                 "total_new_COI_abs_error_mts",
                                 "total_new_COI_abs_error_sd_mts"))

growing_df$aggregation_method <- gsub(".*_", "", growing_df$variable)
growing_df$variable <- gsub("_stm|_mts", "", growing_df$variable)
growing_df <- growing_df %>%
    tidyr::pivot_wider(names_from = variable, values_from = value) %>%
    dplyr::mutate(aggregation_method = ifelse(aggregation_method == 'stm', "Sum then max", "Max then sum"))

########
# Loci #
########

current_df <- growing_df %>%
    dplyr::filter(drop_out == '0.2')

plot1 <- ggplot(current_df,
                aes(x = loci, y = total_weighted_probability_error, fill = model, shape = multiple_imputation, group = interaction(model, multiple_imputation))) +
        geom_hline(yintercept = 0, color = 'black') +
    geom_line(aes(color = model), linewidth = 1) +
    scale_color_manual(
        values = c(
            "Bayesian" = "orange",
            "Clustering" = "maroon",
            "Simple" = "darkblue"
        ), guide = "none"
    ) +
    ggnewscale::new_scale_fill() +
    geom_point(aes(fill = model), size = 3, alpha=1) +
    labs(
        x = "Loci",
        y = "Predicted probability - True probability",
        fill = "Model",
        shape = "Multiple Imputation",
        title = "Probability an infection is new"
    ) +
    theme_bw() +
    scale_fill_manual(
        values = c(
            "Bayesian" = "orange",
            "Clustering" = "maroon",
            "Simple" = "darkblue"
        )
    ) +
    scale_shape_manual(values = c("Yes" = 21, "No" = 22), guide = 'none') +
    scale_y_continuous(labels = scales::percent, breaks = seq(-0.25, 0.1, 0.05), limits = c(-0.25, 0.1)) +
    guides(fill = guide_legend(override.aes = list(shape=21))) +
    facet_wrap( ~ synthetic_type)

ggsave('~/Documents/Harvard University/Rotations/Neafsey/figures/testing/probability_error_loci.png', plot1, width = 6, height = 4)

plot2 <- ggplot(current_df,
                aes(x = loci, y = total_new_COI_error, shape = aggregation_method, group = interaction(aggregation_method, model))) +
        geom_hline(yintercept = 0, color = 'black') +
    geom_line(aes(color = model), linewidth = 1, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = total_new_COI_error - total_new_COI_error_sd,
        ymax = total_new_COI_error + total_new_COI_error_sd,
        color = model),
    position = position_dodge(width = 0.5), width = 0.5) +
    scale_color_manual(
        values = c(
            "Bayesian" = "orange",
            "Clustering" = "maroon",
            "Simple" = "darkblue"
        ), guide = "none"
    ) +
    ggnewscale::new_scale_fill() +
    geom_point(aes(fill = model), size = 3, alpha=1, position = position_dodge(width = 0.5)) +
    labs(
        title = "Cross-locus molFOI",
        x = "Loci",
        y = "Predicted molFOI - True molFOI",
        fill = "Model",
        shape = "Aggregation method"
    ) +
    theme_bw() +
    scale_fill_manual(
        values = c(
            "Bayesian" = "orange",
            "Clustering" = "maroon",
            "Simple" = "darkblue"
        )
    ) +
    scale_shape_manual(values = c("Sum then max" = 22, "Max then sum" = 21)) +
    guides(fill = guide_legend(override.aes = list(shape=21))) +
    facet_wrap( ~ synthetic_type)
ggsave('~/Documents/Harvard University/Rotations/Neafsey/figures/testing/molFOI_obs_error_loci.png', plot2, width = 6, height = 4)

############
# Drop-out #
############

current_df <- growing_df %>%
    dplyr::filter(loci == 5) %>%
    dplyr::mutate(drop_out = as.numeric(drop_out))

plot1 <- ggplot(current_df,
                aes(x = drop_out, y = total_weighted_probability_error, fill = model, shape = multiple_imputation, group = interaction(model, multiple_imputation))) +
    geom_hline(yintercept = 0, color = 'black') +
    geom_line(aes(color = model), linewidth = 1) +
    scale_color_manual(
        values = c(
            "Bayesian" = "orange",
            "Clustering" = "maroon",
            "Simple" = "darkblue"
        ), guide = "none"
    ) +
    ggnewscale::new_scale_fill() +
    geom_point(aes(fill = model), size = 3, alpha=1) +
    labs(
        x = "Drop-out",
        y = "Predicted probability - True probability",
        fill = "Model",
        shape = "Multiple Imputation",
        title = "Probability an infection is new"
    ) +
    theme_bw() +
    scale_fill_manual(
        values = c(
            "Bayesian" = "orange",
            "Clustering" = "maroon",
            "Simple" = "darkblue"
        )
    ) +
    scale_shape_manual(values = c("Yes" = 21, "No" = 22), guide = 'none') +
    scale_y_continuous(labels = scales::percent, breaks = seq(-0.4, 0.1, 0.05), limits = c(-0.4, 0.1)) +
    scale_x_continuous(labels = scales::percent) +
    guides(fill = guide_legend(override.aes = list(shape=21))) +
    facet_wrap( ~ synthetic_type)

ggsave('~/Documents/Harvard University/Rotations/Neafsey/figures/testing/probability_error_drop_out.png', plot1, width = 6, height = 4)

plot2 <- ggplot(current_df,
                aes(x = drop_out, y = total_new_COI_error, shape = aggregation_method, group = interaction(aggregation_method, model))) +
    geom_hline(yintercept = 0, color = 'black') +
    geom_line(aes(color = model), linewidth = 1, position = position_dodge(width = 0.05)) +
    geom_errorbar(aes(ymin = total_new_COI_error - total_new_COI_error_sd,
                      ymax = total_new_COI_error + total_new_COI_error_sd,
                      color = model),
                  position = position_dodge(width = 0.05), width = 0.05) +
    scale_color_manual(
        values = c(
            "Bayesian" = "orange",
            "Clustering" = "maroon",
            "Simple" = "darkblue"
        ), guide = "none"
    ) +
    ggnewscale::new_scale_fill() +
    geom_point(aes(fill = model), size = 3, alpha=1, position = position_dodge(width = 0.05)) +
    labs(
        title = "Cross-locus molFOI",
        x = "Drop-out",
        y = "Predicted molFOI - True molFOI",
        fill = "Model",
        shape = "Aggregation method"
    ) +
    theme_bw() +
    scale_fill_manual(
        values = c(
            "Bayesian" = "orange",
            "Clustering" = "maroon",
            "Simple" = "darkblue"
        )
    ) +
    scale_shape_manual(values = c("Sum then max" = 22, "Max then sum" = 21)) +
    guides(fill = guide_legend(override.aes = list(shape=21))) +
    scale_x_continuous(labels = scales::percent, breaks = c(0, 0.2, 0.4, 0.6, 0.8)) +
    facet_wrap( ~ synthetic_type)
ggsave('~/Documents/Harvard University/Rotations/Neafsey/figures/testing/molFOI_obs_error_drop_out.png', plot2, width = 6, height = 4)





