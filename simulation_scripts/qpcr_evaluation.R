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

parameters <- 'simulation_scripts/config/qpcr_testing.txt'
parameter_text <- readLines(parameters)
parsed_params <- parse_parameters(parameter_text)
params_for_files <- generate_combinations(parsed_params)
params_df <- combinations_to_dataframe(params_for_files)
params_df <- params_df[names(parsed_params)]
run_only_args <- c('multiple_imputation', 'model')
params_df <- params_df[c(names(parsed_params)[!names(parsed_params) %in% run_only_args], run_only_args)]

rownames(params_df) <- paste0("sim", 1:nrow(params_df))

input_path <- 'simulation_output/qpcr/output/'
files_in <- list.files(input_path, full.names = TRUE)

prepare_error <- function(row_num) {
    expanded_df <- tryCatch({
        inputString <- paste0(unlist(c(params_df[row_num,])), collapse = '_')

        if (!(length(files_in[grepl(inputString, files_in)]) > 1)) {
            return(NULL)
        }
        print(cbind(params_df[row_num,], 'count' = length(files_in[grepl(inputString, files_in)])))

        total_input <- data.frame()
        total_new_infections <- data.frame()
        for (file_in in files_in[grepl(inputString, files_in)]) {
            if (grepl('_probabilities.tsv', file_in)) {
                example_input <- read.csv(file_in, sep = '\t')
                if (!'actually_present' %in% colnames(example_input)) {
                    example_input$actually_present <- example_input$present
                }
                
                example_input$replicate <- as.numeric(sub(".*_(.*?)_probabilities\\.tsv$", "\\1", file_in))
                
                new_infections <- example_input %>%
                    dplyr::group_by(subject, replicate, time) %>%
                    dplyr::summarise(new_infections_tmp = 1 * any(infection_event == 1), .groups = 'drop') %>%
                    dplyr::group_by(subject, replicate) %>%
                    dplyr::summarise(new_infections_true = sum(new_infections_tmp), .groups = 'drop') %>%
                    dplyr::mutate(subject = as.character(subject))
                
                estimated_new_infections <- data.frame(new_infections = rowMeans(read.csv(gsub("probabilities", "infections", file_in), sep = '\t')),
                                                       subject = as.character(rownames(read.csv(gsub("probabilities", "infections", file_in), sep = '\t'))))
                
                new_infection_comparison <- full_join(new_infections, estimated_new_infections, by = c("subject"))
                
                total_input <- rbind(total_input, example_input)
                total_new_infections <- rbind(total_new_infections, new_infection_comparison)
            }
        }

        binned_props <- total_input %>%
            dplyr::filter(present == 1) %>%
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

        total_novelty_molFOI_error_observed <- mean((total_input$predicted_novelty_observed - total_input$true_novelty_observed))
        total_novelty_molFOI_abs_error_observed <- mean(abs(total_input$predicted_novelty_observed - total_input$true_novelty_observed))
        total_novelty_molFOI_error_sd_observed <- sd((total_input$predicted_novelty_observed - total_input$true_novelty_observed))
        total_novelty_molFOI_abs_error_sd_observed <- sd(abs(total_input$predicted_novelty_observed - total_input$true_novelty_observed))
        total_novelty_molFOI_error_imputed <- mean((total_input$predicted_novelty_imputed - total_input$true_novelty_imputed))
        total_novelty_molFOI_abs_error_imputed <- mean(abs(total_input$predicted_novelty_imputed - total_input$true_novelty_imputed))
        total_novelty_molFOI_error_sd_imputed <- sd((total_input$predicted_novelty_imputed - total_input$true_novelty_imputed))
        total_novelty_molFOI_abs_error_sd_imputed <- sd(abs(total_input$predicted_novelty_imputed - total_input$true_novelty_imputed))

        total_new_infections_error <- mean((total_new_infections$new_infections - total_new_infections$new_infections_true))
        total_new_infections_abs_error <- mean(abs(total_new_infections$new_infections - total_new_infections$new_infections_true))
        total_new_infections_sd <- sd((total_new_infections$new_infections - total_new_infections$new_infections_true))
        total_new_infections_abs_sd <- sd(abs(total_new_infections$new_infections - total_new_infections$new_infections_true))

        expanded_df <- data.frame(params_df[row_num,],
        total_weighted_probability_error,
        total_weighted_probability_abs_error,
        total_novelty_molFOI_error_observed,
        total_novelty_molFOI_abs_error_observed,
        total_novelty_molFOI_error_sd_observed,
        total_novelty_molFOI_abs_error_sd_observed,
        total_novelty_molFOI_error_imputed,
        total_novelty_molFOI_abs_error_imputed,
        total_novelty_molFOI_error_sd_imputed,
        total_novelty_molFOI_abs_error_sd_imputed,
        total_new_infections_error,
        total_new_infections_abs_error,
        total_new_infections_sd,
        total_new_infections_abs_sd)

        expanded_df
    }, error = function(e){e})
    return(expanded_df)
}

if (!file.exists('evaluation_intermediates/qpcr_evaluation.RDS')) {
    # num_cores <- 1
    # cl <- makeCluster(num_cores)
    # registerDoParallel(cl)

    growing_df <- list()
    for (row_num in 1:nrow(params_df)) {
        suppressMessages(growing_df[[row_num]] <- prepare_error(row_num))
    }

    # stopCluster(cl)

    growing_df_save <- dplyr::bind_rows(growing_df)
    dir.create('evaluation_intermediates', showWarnings = FALSE)
    saveRDS(growing_df_save, file = 'evaluation_intermediates/qpcr_evaluation.RDS')
} else {
    growing_df_save <- readRDS('evaluation_intermediates/qpcr_evaluation.RDS')
}

# Change here
growing_df <- growing_df_save
growing_df <- growing_df %>%
    filter(model != 'allele-specific') %>%
    mutate(model = case_when(model == 'bayesian' ~ 'Bayesian',
                             model == 'clustering' ~ 'Clustering',
                             model == 'simple' ~ 'Simple'),
           multiple_imputation = ifelse(multiple_imputation, 'Yes', 'No')) %>%
    mutate(synthetic_type = case_when(synthetic_type == 'bin-present' ~ 'Rolling presence probability',
                                      synthetic_type == 'pois-time' ~ 'Poisson time to clearance'))

growing_df <- rbind(growing_df %>% dplyr::filter(!(multiple_imputation == 'Yes' & qPCR_only == 0)),
                    growing_df %>% dplyr::filter(multiple_imputation == 'No', qPCR_only == 0) %>%
                        mutate(multiple_imputation = 'Yes'))

# First plot: total_novelty_error vs model
plot1 <- ggplot(growing_df,
                aes(x = factor(qPCR_only), y = total_weighted_probability_error, fill = model, shape = multiple_imputation, group = interaction(model, multiple_imputation))) +
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
        x = "qPCR Only Proportion",
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
    scale_shape_manual(values = c("Yes" = 21, "No" = 22)) +
    scale_y_continuous(labels = scales::percent, breaks = seq(-0.15, 0.25, 0.05)) +
    guides(fill = guide_legend(override.aes = list(shape=21))) +
    facet_wrap( ~ synthetic_type)

dir.create('figures/testing', showWarnings = F, recursive = T)
ggsave('figures/testing/probability_error.png', plot1, width = 6, height = 4)

# Second plot: mean total_novelty_error vs qPCR_only with error bars
plot2 <- ggplot(growing_df,
                aes(x = factor(qPCR_only), y = total_novelty_molFOI_error_observed, shape = multiple_imputation, group = interaction(multiple_imputation, model))) +
        geom_hline(yintercept = 0, color = 'black') +
    geom_line(aes(color = model), linewidth = 1, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = total_novelty_molFOI_error_observed - total_novelty_molFOI_error_sd_observed,
        ymax = total_novelty_molFOI_error_observed + total_novelty_molFOI_error_sd_observed,
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
        title = "Counting new infections only for sequenced points",
        x = "qPCR Only Proportion",
        y = "Predicted molFOI - True molFOI",
        fill = "Model",
        shape = "Multiple Imputation"
    ) +
    theme_bw() +
    scale_fill_manual(
        values = c(
            "Bayesian" = "orange",
            "Clustering" = "maroon",
            "Simple" = "darkblue"
        )
    ) +
    scale_shape_manual(values = c("Yes" = 21, "No" = 22)) +
    guides(fill = guide_legend(override.aes = list(shape=21))) +
    facet_wrap( ~ synthetic_type)
ggsave('figures/testing/molFOI_obs_error.png', plot2, width = 6, height = 4)

plot3 <- ggplot(growing_df,
                aes(x = factor(qPCR_only), y = total_novelty_molFOI_error_imputed, fill = model, shape = multiple_imputation, group = interaction(multiple_imputation, model))) +
        geom_hline(yintercept = 0, color = 'black') +
    geom_line(aes(color = model), linewidth = 1, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = total_novelty_molFOI_error_imputed - total_novelty_molFOI_error_sd_imputed,
            ymax = total_novelty_molFOI_error_imputed + total_novelty_molFOI_error_sd_imputed,
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
        title = "Counting new infections for sequenced and imputed points",
        x = "qPCR Only Proportion",
        y = "Predicted molFOI - True molFOI",
        fill = "Model",
        shape = "Multiple Imputation"
    ) +
    theme_bw() +
    scale_fill_manual(
        values = c(
            "Bayesian" = "orange",
            "Clustering" = "maroon",
            "Simple" = "darkblue"
        )
    ) +
    scale_shape_manual(values = c("Yes" = 21, "No" = 22)) +
    guides(fill = guide_legend(override.aes = list(shape=21))) +
    facet_wrap( ~ synthetic_type)

ggsave('figures/testing/molFOI_inf_error.png', plot3, width = 6, height = 4)

plot4 <- ggplot(growing_df,
                aes(x = factor(qPCR_only), y = total_new_infections_error, fill = model, shape = multiple_imputation, group = interaction(multiple_imputation, model))) +
        geom_hline(yintercept = 0, color = 'black') +
    geom_line(aes(color = model), linewidth = 1, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = total_new_infections_error - total_new_infections_sd,
                ymax = total_new_infections_error + total_new_infections_sd,
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
        title = "Counting new infections for sequenced and imputed points",
        x = "qPCR Only Proportion",
        y = "Predicted infections - True infections",
        fill = "Model",
        shape = "Multiple Imputation"
    ) +
    theme_bw() +
    scale_fill_manual(
        values = c(
            "Bayesian" = "orange",
            "Clustering" = "maroon",
            "Simple" = "darkblue"
        )
    ) +
    scale_shape_manual(values = c("Yes" = 21, "No" = 22)) +
    guides(fill = guide_legend(override.aes = list(shape=21))) +
    facet_wrap( ~ synthetic_type)
ggsave('figures/testing/infections_error.png', plot4, width = 6, height = 4)


plot5 <- ggplot(growing_df,
                aes(x = factor(qPCR_only), y = total_weighted_probability_abs_error, fill = model, shape = multiple_imputation, group = interaction(model, multiple_imputation))) +
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
        x = "qPCR Only Proportion",
        y = "|Predicted probability - True probability|",
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
    scale_shape_manual(values = c("Yes" = 21, "No" = 22)) +
    scale_y_continuous(labels = scales::percent, breaks = seq(-0.15, 0.25, 0.05)) +
    guides(fill = guide_legend(override.aes = list(shape=21))) +
    facet_wrap( ~ synthetic_type)

ggsave('figures/testing/probability_abs_error.png', plot5, width = 6, height = 4)

#########################
# Combined plot for all #
#########################

p1 <- ggplot(growing_df %>% filter(synthetic_type == 'Rolling presence probability'),
                aes(x = factor(qPCR_only), y = total_weighted_probability_error, fill = model, shape = multiple_imputation, group = interaction(model, multiple_imputation))) +
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
        x = "qPCR Only Proportion",
        y = "Predicted probability - True probability",
        fill = "Model",
        shape = "Multiple Imputation",
        title = "Rolling presence probability"
    ) +
    theme_bw() +
    scale_fill_manual(
        values = c(
            "Bayesian" = "orange",
            "Clustering" = "maroon",
            "Simple" = "darkblue"
        )
    ) +
    scale_shape_manual(values = c("Yes" = 21, "No" = 22)) +
    scale_y_continuous(labels = scales::percent, breaks = seq(-0.25, 0.25, 0.05), limits = c(-0.25, 0.25)) +
    guides(fill = guide_legend(override.aes = list(shape=21)))

p2 <- ggplot(growing_df %>% filter(synthetic_type == 'Rolling presence probability'),
                aes(x = factor(qPCR_only), y = total_novelty_molFOI_error_imputed, fill = model, shape = multiple_imputation, group = interaction(multiple_imputation, model))) +
    geom_hline(yintercept = 0, color = 'black') +
    geom_line(aes(color = model), linewidth = 1, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = total_novelty_molFOI_error_imputed - total_novelty_molFOI_error_sd_imputed,
                      ymax = total_novelty_molFOI_error_imputed + total_novelty_molFOI_error_sd_imputed,
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
        title = "Rolling presence probability",
        x = "qPCR Only Proportion",
        y = "Predicted molFOI - True molFOI",
        fill = "Model",
        shape = "Multiple Imputation"
    ) +
    theme_bw() +
    scale_fill_manual(
        values = c(
            "Bayesian" = "orange",
            "Clustering" = "maroon",
            "Simple" = "darkblue"
        )
    ) +
    scale_shape_manual(values = c("Yes" = 21, "No" = 22)) +
    scale_y_continuous(limits = c(-40, 40), breaks = seq(-40, 40, 10)) +
    guides(fill = guide_legend(override.aes = list(shape=21)))

p3 <- ggplot(growing_df %>% filter(synthetic_type == 'Rolling presence probability'),
                aes(x = factor(qPCR_only), y = total_new_infections_error, fill = model, shape = multiple_imputation, group = interaction(multiple_imputation, model))) +
    geom_hline(yintercept = 0, color = 'black') +
    geom_line(aes(color = model), linewidth = 1, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = total_new_infections_error - total_new_infections_sd,
                      ymax = total_new_infections_error + total_new_infections_sd,
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
        title = "Rolling presence probability",
        x = "qPCR Only Proportion",
        y = "Predicted infections - True infections",
        fill = "Model",
        shape = "Multiple Imputation"
    ) +
    theme_bw() +
    scale_fill_manual(
        values = c(
            "Bayesian" = "orange",
            "Clustering" = "maroon",
            "Simple" = "darkblue"
        )
    ) +
    scale_shape_manual(values = c("Yes" = 21, "No" = 22)) +
    scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
    guides(fill = guide_legend(override.aes = list(shape=21)))

p4 <- ggplot(growing_df %>% filter(synthetic_type == 'Poisson time to clearance'),
             aes(x = factor(qPCR_only), y = total_weighted_probability_error, fill = model, shape = multiple_imputation, group = interaction(model, multiple_imputation))) +
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
        x = "qPCR Only Proportion",
        y = "Predicted probability - True probability",
        fill = "Model",
        shape = "Multiple Imputation",
        title = "Poisson time to clearance"
    ) +
    theme_bw() +
    scale_fill_manual(
        values = c(
            "Bayesian" = "orange",
            "Clustering" = "maroon",
            "Simple" = "darkblue"
        )
    ) +
    scale_shape_manual(values = c("Yes" = 21, "No" = 22)) +
    scale_y_continuous(labels = scales::percent, breaks = seq(-0.25, 0.25, 0.05), limits = c(-0.25, 0.25)) +
    guides(fill = guide_legend(override.aes = list(shape=21)))

p5 <- ggplot(growing_df %>% filter(synthetic_type == 'Poisson time to clearance'),
             aes(x = factor(qPCR_only), y = total_novelty_molFOI_error_imputed, fill = model, shape = multiple_imputation, group = interaction(multiple_imputation, model))) +
    geom_hline(yintercept = 0, color = 'black') +
    geom_line(aes(color = model), linewidth = 1, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = total_novelty_molFOI_error_imputed - total_novelty_molFOI_error_sd_imputed,
                      ymax = total_novelty_molFOI_error_imputed + total_novelty_molFOI_error_sd_imputed,
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
        title = "Poisson time to clearance",
        x = "qPCR Only Proportion",
        y = "Predicted molFOI - True molFOI",
        fill = "Model",
        shape = "Multiple Imputation"
    ) +
    theme_bw() +
    scale_fill_manual(
        values = c(
            "Bayesian" = "orange",
            "Clustering" = "maroon",
            "Simple" = "darkblue"
        )
    ) +
    scale_shape_manual(values = c("Yes" = 21, "No" = 22)) +
    scale_y_continuous(limits = c(-40, 40), breaks = seq(-40, 40, 10)) +
    guides(fill = guide_legend(override.aes = list(shape=21)))

p6 <- ggplot(growing_df %>% filter(synthetic_type == 'Poisson time to clearance'),
             aes(x = factor(qPCR_only), y = total_new_infections_error, fill = model, shape = multiple_imputation, group = interaction(multiple_imputation, model))) +
    geom_hline(yintercept = 0, color = 'black') +
    geom_line(aes(color = model), linewidth = 1, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = total_new_infections_error - total_new_infections_sd,
                      ymax = total_new_infections_error + total_new_infections_sd,
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
        title = "Poisson time to clearance",
        x = "qPCR Only Proportion",
        y = "Predicted infections - True infections",
        fill = "Model",
        shape = "Multiple Imputation"
    ) +
    theme_bw() +
    scale_fill_manual(
        values = c(
            "Bayesian" = "orange",
            "Clustering" = "maroon",
            "Simple" = "darkblue"
        )
    ) +
    scale_shape_manual(values = c("Yes" = 21, "No" = 22)) +
    scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
    guides(fill = guide_legend(override.aes = list(shape=21)))

plot_out <- patchwork::wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 3, guides = 'collect', axis_titles = 'collect')
ggsave('figures/testing/all_metrics.png', plot_out, width = 10, height = 6)






###############################
# Combined for absolute scale #
###############################

p1 <- ggplot(growing_df %>% filter(synthetic_type == 'Rolling presence probability'),
             aes(x = factor(qPCR_only), y = total_weighted_probability_abs_error, fill = model, shape = multiple_imputation, group = interaction(model, multiple_imputation))) +
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
        x = "qPCR Only Proportion",
        y = "|Predicted probability - True probability|",
        fill = "Model",
        shape = "Multiple Imputation",
        title = "Rolling presence probability"
    ) +
    theme_bw() +
    scale_fill_manual(
        values = c(
            "Bayesian" = "orange",
            "Clustering" = "maroon",
            "Simple" = "darkblue"
        )
    ) +
    scale_shape_manual(values = c("Yes" = 21, "No" = 22)) +
    scale_y_continuous(labels = scales::percent, breaks = seq(-0.25, 0.25, 0.05), limits = c(0, 0.25)) +
    guides(fill = guide_legend(override.aes = list(shape=21)))

p2 <- ggplot(growing_df %>% filter(synthetic_type == 'Rolling presence probability'),
             aes(x = factor(qPCR_only), y = total_novelty_molFOI_abs_error_imputed, fill = model, shape = multiple_imputation, group = interaction(multiple_imputation, model))) +
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
        title = "Rolling presence probability",
        x = "qPCR Only Proportion",
        y = "|Predicted molFOI - True molFOI|",
        fill = "Model",
        shape = "Multiple Imputation"
    ) +
    theme_bw() +
    scale_fill_manual(
        values = c(
            "Bayesian" = "orange",
            "Clustering" = "maroon",
            "Simple" = "darkblue"
        )
    ) +
    scale_shape_manual(values = c("Yes" = 21, "No" = 22)) +
    scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
    guides(fill = guide_legend(override.aes = list(shape=21)))

p3 <- ggplot(growing_df %>% filter(synthetic_type == 'Rolling presence probability'),
             aes(x = factor(qPCR_only), y = total_new_infections_abs_error, fill = model, shape = multiple_imputation, group = interaction(multiple_imputation, model))) +
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
        title = "Rolling presence probability",
        x = "qPCR Only Proportion",
        y = "|Predicted infections - True infections|",
        fill = "Model",
        shape = "Multiple Imputation"
    ) +
    theme_bw() +
    scale_fill_manual(
        values = c(
            "Bayesian" = "orange",
            "Clustering" = "maroon",
            "Simple" = "darkblue"
        )
    ) +
    scale_shape_manual(values = c("Yes" = 21, "No" = 22)) +
    scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.3)) +
    guides(fill = guide_legend(override.aes = list(shape=21)))

p4 <- ggplot(growing_df %>% filter(synthetic_type == 'Poisson time to clearance'),
             aes(x = factor(qPCR_only), y = total_weighted_probability_abs_error, fill = model, shape = multiple_imputation, group = interaction(model, multiple_imputation))) +
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
        x = "qPCR Only Proportion",
        y = "|Predicted probability - True probability|",
        fill = "Model",
        shape = "Multiple Imputation",
        title = "Poisson time to clearance"
    ) +
    theme_bw() +
    scale_fill_manual(
        values = c(
            "Bayesian" = "orange",
            "Clustering" = "maroon",
            "Simple" = "darkblue"
        )
    ) +
    scale_shape_manual(values = c("Yes" = 21, "No" = 22)) +
    scale_y_continuous(labels = scales::percent, breaks = seq(-0.25, 0.25, 0.05), limits = c(0, 0.25)) +
    guides(fill = guide_legend(override.aes = list(shape=21)))

p5 <- ggplot(growing_df %>% filter(synthetic_type == 'Poisson time to clearance'),
             aes(x = factor(qPCR_only), y = total_novelty_molFOI_abs_error_imputed, fill = model, shape = multiple_imputation, group = interaction(multiple_imputation, model))) +
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
        title = "Poisson time to clearance",
        x = "qPCR Only Proportion",
        y = "|Predicted molFOI - True molFOI|",
        fill = "Model",
        shape = "Multiple Imputation"
    ) +
    theme_bw() +
    scale_fill_manual(
        values = c(
            "Bayesian" = "orange",
            "Clustering" = "maroon",
            "Simple" = "darkblue"
        )
    ) +
    scale_shape_manual(values = c("Yes" = 21, "No" = 22)) +
    scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
    guides(fill = guide_legend(override.aes = list(shape=21)))

p6 <- ggplot(growing_df %>% filter(synthetic_type == 'Poisson time to clearance'),
             aes(x = factor(qPCR_only), y = total_new_infections_abs_error, fill = model, shape = multiple_imputation, group = interaction(multiple_imputation, model))) +
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
        title = "Poisson time to clearance",
        x = "qPCR Only Proportion",
        y = "|Predicted infections - True infections|",
        fill = "Model",
        shape = "Multiple Imputation"
    ) +
    theme_bw() +
    scale_fill_manual(
        values = c(
            "Bayesian" = "orange",
            "Clustering" = "maroon",
            "Simple" = "darkblue"
        )
    ) +
    scale_shape_manual(values = c("Yes" = 21, "No" = 22)) +
    scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.3)) +
    guides(fill = guide_legend(override.aes = list(shape=21)))

plot_out <- patchwork::wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 3, guides = 'collect', axis_titles = 'collect')
ggsave('figures/testing/all_metrics_abs.png', plot_out, width = 10, height = 6)




##################
# In-text values #
##################

growing_df_save %>%
    filter(qPCR_only == 0) %>%
    select(synthetic_type, model, total_weighted_probability_error)

rbind(growing_df_save %>% dplyr::filter(!(multiple_imputation == TRUE & qPCR_only == 0)),
      growing_df_save %>% dplyr::filter(multiple_imputation == FALSE, qPCR_only == 0) %>%
          mutate(multiple_imputation = TRUE)) %>%
    filter(qPCR_only %in% c(0, 0.5)) %>%
    select(synthetic_type, model, qPCR_only, multiple_imputation, total_weighted_probability_error) %>%
    group_by(synthetic_type, model, multiple_imputation) %>%
    arrange(qPCR_only) %>%
    summarise(diff(total_weighted_probability_error)[1])

growing_df_save %>%
    filter(qPCR_only == 0) %>%
    select(synthetic_type, model, total_new_infections_error)

growing_df_save %>%
    filter(qPCR_only == 0) %>%
    select(synthetic_type, model, total_new_infections_abs_error)

growing_df_save %>%
    filter(qPCR_only == 0) %>%
    select(synthetic_type, model, total_novelty_molFOI_error_imputed)















