library(dinemites)
library(dplyr)

treatments <- read.csv('real_data_outputs/treatments_CIS43LS.tsv', sep = '\t')
for (model in c('bayesian', 'clustering', 'simple')) {
    if (model == 'bayesian') {
        figures_folder <- 'figures/CIS43LS/bayesian/'
        final_df <- read.csv('real_data_outputs/CIS43LS_probabilities_bayesian.tsv', sep = '\t')
        infections_df <- read.csv('real_data_outputs/CIS43LS_infections_bayesian.tsv', sep = '\t')
    } else if (model == 'clustering') {
        figures_folder <- 'figures/CIS43LS/clustering/'
        final_df <- read.csv('real_data_outputs/CIS43LS_probabilities_clustering.tsv', sep = '\t')
        infections_df <- read.csv('real_data_outputs/CIS43LS_infections_clustering.tsv', sep = '\t')
    } else if (model == 'simple') {
        figures_folder <- 'figures/CIS43LS/simple/'
        final_df <- read.csv('real_data_outputs/CIS43LS_probabilities_simple.tsv', sep = '\t')
        infections_df <- read.csv('real_data_outputs/CIS43LS_infections_simple.tsv', sep = '\t')
    }
    dir.create(figures_folder, showWarnings = F, recursive = T)
    
    final_df <- final_df %>%
        mutate(locus = case_when(locus == 'PfCSP_1' ~ 'PF3D7_0304600',
                                 locus == 'PfTRAP_1' ~ 'PF3D7_1335900',
                                 TRUE ~ locus))

    plot_dataset(dataset = final_df, 
                 treatments = treatments,
                 estimated_new_infections = infections_df, 
                 output = figures_folder)
    if (model == 'bayesian') {
        freq_df <- final_df %>%
            dplyr::group_by(.data$subject, .data$time) %>%
            dplyr::filter(any(.data$present == 1)) %>%
            dplyr::ungroup() %>%
            dplyr::group_by(.data$allele) %>%
            dplyr::summarise(mean(.data$present == 1),
                             .groups = "drop")
        
        colnames(freq_df) <- c("allele", "prevalence")
        
        freq_table <- freq_df$prevalence
        names(freq_table) <- freq_df$allele
        
        final_df <- final_df %>%
            mutate(prevalence =
                       paste0(round(freq_table[as.character(.data$allele)] * 100, 1),
                              "%")) %>%
            mutate(prevalence =
                       ifelse(.data$prevalence == 'NA%', "", .data$prevalence)) %>%
            ungroup()
        
        
        output_file <- 'figures/CIS43LS/bayesian/subject_0266R.png'
        plot_out <- plot_single_subject('0266R', 
                            final_df %>% dplyr::mutate(locus = gsub('_', ' ', locus)), 
                            estimated_new_infections = infections_df, 
                            output = output_file,
                            treatments = treatments, 
                            height = 7,
                            width = 8)
        plot_out[[1]] <- plot_out[[1]] + ggplot2::ggtitle(paste0('One individual from a CIS43LS trial\nEstimated new infection events: ', round(mean(as.numeric(infections_df['0266R',])), 1)))
        ggplot2::ggsave(filename = output_file, plot_out, height = 6, width = 8)
        
        final_df$probability_new <- NULL
        final_df$prevalence <- NULL
        dir.create('figures/CIS43LS/example', showWarnings = FALSE)
        output_file <- 'figures/CIS43LS/example/subject_0266R.png'
        plot_out <- plot_single_subject('0266R', 
                            final_df %>% dplyr::mutate(locus = gsub('_', ' ', locus)), 
                            output = output_file,
                            treatments = treatments, 
                            height = 4,
                            width = 9, 
                            no_imputation = TRUE)
        plot_out <- plot_out + 
            ggplot2::ggtitle('One individual from a CIS43LS trial')
        ggplot2::ggsave(filename = output_file, plot_out, height = 4, width = 9)
    }
}

