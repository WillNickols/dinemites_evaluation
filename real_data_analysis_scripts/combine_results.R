library(dinemites)
library(dplyr)
library(ggplot2)

dataset_names <- c('CIS43LS', 'mali', 'uganda')

##########################
# Probability comparison #
##########################
plot_list <- list()
for (i in 1:length(dataset_names)) {
    dataset_name <- dataset_names[i]
    model <- 'bayesian'
    final_df <- read.csv(paste0('real_data_outputs/', dataset_name, '_probabilities_', model, '.tsv'), sep = '\t')
    final_df <- final_df %>% dplyr::rename(bayesian = probability_new)
    for (model in c('clustering', 'simple')) {
        final_df_tmp <- read.csv(paste0('real_data_outputs/', dataset_name, '_probabilities_', model, '.tsv'), sep = '\t')
        colnames(final_df_tmp)[colnames(final_df_tmp) == 'probability_new'] <- model
        final_df <- full_join(final_df, final_df_tmp)
    }
    final_df <- final_df %>%
        rename('Bayesian' = 'bayesian', 'Clustering' = 'clustering', 'Simple' = 'simple')

    cormat <- cor(final_df[final_df$present == 1, c("Bayesian", "Clustering", "Simple")],
                  use = 'pairwise.complete.obs')
    
    cormat_inner <- cor(final_df[final_df$present == 1 & final_df$Bayesian > 0.05 & final_df$Bayesian < 0.95 &
                               final_df$Clustering > 0.05 & final_df$Clustering < 0.95, 
                           c("Bayesian", "Clustering", "Simple")],
                  use = 'pairwise.complete.obs')
    
    # In-text value
    print(cormat_inner)
    print(sum(final_df$present == 1 & final_df$Bayesian > 0.05 & final_df$Bayesian < 0.95 &
                  final_df$Clustering > 0.05 & final_df$Clustering < 0.95, na.rm=T))
    print(sum(final_df$present == 1, na.rm=T))
    
    plotting_df <- reshape2::melt(cormat) %>%
        mutate(Var1 = factor(Var1, levels = c("Bayesian", "Clustering", "Simple")),
               Var2 = factor(Var2, levels = c("Bayesian", "Clustering", "Simple"))) %>%
        filter(as.numeric(Var1) >= as.numeric(Var2))

    dataset_name_nice <- case_when(dataset_name == 'CIS43LS' ~ 'CIS43LS cohort',
                                   dataset_name == 'mali' ~ '2011 Malian cohort',
                                   dataset_name == 'uganda' ~ 'Ugandan cohort')

    plot_list[[i]] <- ggplot(plotting_df, aes(x = Var1, y = Var2, fill = value)) +
        geom_tile(color = "black", size = 0.5) +
        geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 4) +
        coord_fixed() +
        labs(fill = 'Correlation', title = dataset_name_nice) +
        theme_minimal() +
        theme(
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
        ) +
        scale_fill_gradient(low = "#FFCCCB", high = "#FF0000", limits = c(0.5, 1)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none')
}

plot_out <- patchwork::wrap_plots(plot_list[[1]], plot_list[[2]], plot_list[[3]], ncol = 3, guides = 'collect')
dir.create('figures/general', showWarnings = F)
ggplot2::ggsave('figures/general/probability_correlation.png', plot = plot_out, height = 3, width = 8, dpi = 1000)

#####################
# molFOI comparison #
#####################

plot_list <- list()
plot_list2 <- list()
for (i in 1:length(dataset_names)) {
    dataset_name <- dataset_names[i]
    model <- 'bayesian'
    final_df <- read.csv(paste0('real_data_outputs/', dataset_name, '_probabilities_', model, '.tsv'), sep = '\t')
    if (!'locus' %in% colnames(final_df)) {
        final_df$locus <- 1
    }
    final_df <- final_df %>%
        group_by(subject, allele) %>%
        filter(n() >= 3) %>%
        mutate(second_smallest_time = sort(time)[2]) %>%
        filter(time > second_smallest_time) %>%
        select(-second_smallest_time) %>%
        ungroup()
    
    final_df <- compute_molFOI(final_df, method = 'sum_then_max')
    final_df <- final_df %>% dplyr::rename(bayesian = molFOI)

    for (model in c('clustering', 'simple')) {
        final_df_tmp <- read.csv(paste0('real_data_outputs/', dataset_name, '_probabilities_', model, '.tsv'), sep = '\t')
        if (!'locus' %in% colnames(final_df_tmp)) {
            final_df_tmp$locus <- 1
        }
        final_df_tmp <- final_df_tmp %>%
            group_by(subject, allele) %>%
            filter(n() >= 3) %>%
            mutate(second_smallest_time = sort(time)[2]) %>%
            filter(time > second_smallest_time) %>%
            select(-second_smallest_time) %>%
            ungroup()
        
        final_df_tmp <- compute_molFOI(final_df_tmp, method = 'sum_then_max')
        colnames(final_df_tmp)[colnames(final_df_tmp) == 'molFOI'] <- model
        final_df <- full_join(final_df, final_df_tmp)
    }
    final_df <- final_df %>%
        rename('Bayesian' = 'bayesian', 'Clustering' = 'clustering', 'Simple' = 'simple')

    cormat <- cor(final_df[, c("Bayesian", "Clustering", "Simple")],
                  use = 'pairwise.complete.obs')
    plotting_df <- reshape2::melt(cormat) %>%
        mutate(Var1 = factor(Var1, levels = c("Bayesian", "Clustering", "Simple")),
               Var2 = factor(Var2, levels = c("Bayesian", "Clustering", "Simple"))) %>%
        filter(as.numeric(Var1) >= as.numeric(Var2))

    dataset_name_nice <- case_when(dataset_name == 'CIS43LS' ~ 'CIS43LS cohort',
                                   dataset_name == 'mali' ~ '2011 Malian cohort',
                                   dataset_name == 'uganda' ~ 'Ugandan cohort')
    
    plot_list[[i]] <- ggplot(plotting_df %>%
                                 mutate(Var1 = case_when(Var1 == 'Bayesian' ~ 'B',
                                                         Var1 == 'Clustering' ~ 'C',
                                                         Var1 == 'Simple' ~ 'S'),
                                        Var2 = case_when(Var2 == 'Bayesian' ~ 'B',
                                                         Var2 == 'Clustering' ~ 'C',
                                                         Var2 == 'Simple' ~ 'S')), 
                             aes(x = Var1, y = Var2, fill = value)) +
        geom_tile(color = "black", size = 0.5) +
        geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 8) +
        coord_fixed() +
        labs(fill = 'Correlation', title = dataset_name_nice) +
        theme_minimal() +
        theme(
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            legend.position = 'none',
            text = element_text(size = 30)
        ) +
        scale_fill_gradient(low = "#FFCCCB", high = "#FF0000", limits = c(0.5, 1))

    plotting_df <- reshape2::melt(final_df, measure.vars = c("Bayesian", "Clustering", "Simple"))
    plotting_df <- plotting_df %>%
        mutate(subject = factor(as.character(subject), levels = (plotting_df %>% group_by(subject) %>%
                                                                     summarize(all_sum = sum(value)) %>%
                                                                     arrange(all_sum))$subject)) %>%
        mutate(subject = as.numeric(subject))

    plotting_df <- plotting_df[plotting_df$subject %in% plotting_df$subject[plotting_df$value > 0],]
    plot_list2[[i]] <- ggplot(plotting_df[sample(1:nrow(plotting_df)),], aes(x = subject, y = value, color = variable)) +
        geom_point(alpha = 0.5) +
        theme_bw() +
        scale_color_manual(
            values = c(
                "Bayesian" = "orange",
                "Clustering" = "maroon",
                "Simple" = "darkblue"
            )
        ) +
        labs(color = 'Model', x = 'Subject number', y = 'molFOI', title = dataset_name_nice) + 
        theme(legend.position = 'none')

    plot_list2[[i]] <- ggExtra::ggMarginal(plot_list2[[i]], type = "histogram", margins = "y", groupFill = T)
}

# plot_out <- patchwork::wrap_plots(plot_list[[1]], plot_list[[2]], plot_list[[3]], ncol = 3, guides = 'collect')
ggplot2::ggsave('figures/general/molFOI_correlation_1.png', plot = plot_list[[1]] + theme(plot.title = element_blank()), height = 4, width = 4, bg = 'white', dpi = 1000)
ggplot2::ggsave('figures/general/molFOI_correlation_2.png', plot = plot_list[[2]] + theme(plot.title = element_blank()), height = 4, width = 4, bg = 'white', dpi = 1000)
ggplot2::ggsave('figures/general/molFOI_correlation_3.png', plot = plot_list[[3]] + theme(plot.title = element_blank()), height = 4, width = 4, bg = 'white', dpi = 1000)

plot_out2 <- patchwork::wrap_plots(plot_list2[[1]], plot_list2[[2]], plot_list2[[3]], ncol = 3, guides = 'collect', axis_titles = 'collect')
ggplot2::ggsave('figures/general/molFOI_count.png', plot = plot_out2, height = 3.5, width = 10, dpi = 1000)

#########################
# Infections comparison #
#########################

plot_list <- list()
plot_list2 <- list()
for (i in 1:length(dataset_names)) {
    dataset_name <- dataset_names[i]
    model <- 'bayesian'
    final_df <- read.csv(paste0('real_data_outputs/', dataset_name, '_infections_', model, '.tsv'), sep = '\t')
    final_df <- data.frame(subject = rownames(final_df), bayesian = rowMeans(final_df))

    for (model in c('clustering', 'simple')) {
        final_df_tmp <- read.csv(paste0('real_data_outputs/', dataset_name, '_infections_', model, '.tsv'), sep = '\t')
        final_df_tmp <- data.frame(subject = rownames(final_df_tmp), new_infections = rowMeans(final_df_tmp))
        colnames(final_df_tmp)[colnames(final_df_tmp) == 'new_infections'] <- model
        final_df <- full_join(final_df, final_df_tmp, by = c('subject'))
    }
    final_df <- final_df %>%
        rename('Bayesian' = 'bayesian', 'Clustering' = 'clustering', 'Simple' = 'simple')
    
    # Mean infections
    print(mean(final_df$Bayesian))

    cormat <- cor(final_df[, c("Bayesian", "Clustering", "Simple")],
                  use = 'pairwise.complete.obs')
    plotting_df <- reshape2::melt(cormat) %>%
        mutate(Var1 = factor(Var1, levels = c("Bayesian", "Clustering", "Simple")),
               Var2 = factor(Var2, levels = c("Bayesian", "Clustering", "Simple"))) %>%
        filter(as.numeric(Var1) >= as.numeric(Var2))

    dataset_name_nice <- case_when(dataset_name == 'CIS43LS' ~ 'CIS43LS cohort',
                                   dataset_name == 'mali' ~ '2011 Malian cohort',
                                   dataset_name == 'uganda' ~ 'Ugandan cohort')
    
    plot_list[[i]] <- ggplot(plotting_df %>%
                                 mutate(Var1 = case_when(Var1 == 'Bayesian' ~ 'B',
                                                         Var1 == 'Clustering' ~ 'C',
                                                         Var1 == 'Simple' ~ 'S'),
                                        Var2 = case_when(Var2 == 'Bayesian' ~ 'B',
                                                         Var2 == 'Clustering' ~ 'C',
                                                         Var2 == 'Simple' ~ 'S')), 
                             aes(x = Var1, y = Var2, fill = value)) +
        geom_tile(color = "black", size = 0.5) +
        geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 8) +
        coord_fixed() +
        labs(fill = 'Correlation', title = dataset_name_nice) +
        theme_minimal() +
        theme(
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            legend.position = 'none',
            text = element_text(size = 30)
        ) +
        scale_fill_gradient(low = "#FFCCCB", high = "#FF0000", limits = c(0.5, 1))
    
    plotting_df <- reshape2::melt(final_df, measure.vars = c("Bayesian", "Clustering", "Simple"))
    plotting_df <- plotting_df %>%
        mutate(subject = factor(as.character(subject), levels = (plotting_df %>% group_by(subject) %>%
                                                         summarize(all_sum = sum(value)) %>%
                                                         arrange(all_sum))$subject)) %>%
        mutate(subject = as.numeric(subject))

    plotting_df <- plotting_df[plotting_df$subject %in% plotting_df$subject[plotting_df$value > 0],]
    plot_list2[[i]] <- ggplot(plotting_df, aes(x = subject, y = value, color = variable)) +
        geom_point(alpha = 0.5) +
        theme_bw() +
        scale_color_manual(
            values = c(
                "Bayesian" = "orange",
                "Clustering" = "maroon",
                "Simple" = "darkblue"
            )
        ) +
        labs(color = 'Model', x = 'Subject number', y = 'New infections', title = dataset_name_nice) +
        scale_y_continuous(breaks = 0:100) + 
        theme(legend.position = 'none')

    plot_list2[[i]] <- ggExtra::ggMarginal(plot_list2[[i]], type = "histogram", margins = "y", groupFill = T)
}

# plot_out <- patchwork::wrap_plots(plot_list[[1]], plot_list[[2]], plot_list[[3]], ncol = 3, guides = 'collect')
ggplot2::ggsave('figures/general/new_infections_correlation_1.png', plot = plot_list[[1]] + theme(plot.title = element_blank()), height = 4, width = 4, bg = 'white', dpi = 1000)
ggplot2::ggsave('figures/general/new_infections_correlation_2.png', plot = plot_list[[2]] + theme(plot.title = element_blank()), height = 4, width = 4, bg = 'white', dpi = 1000)
ggplot2::ggsave('figures/general/new_infections_correlation_3.png', plot = plot_list[[3]] + theme(plot.title = element_blank()), height = 4, width = 4, bg = 'white', dpi = 1000)
plot_out2 <- patchwork::wrap_plots(plot_list2[[1]], plot_list2[[2]], plot_list2[[3]], ncol = 3, guides = 'collect', axis_titles = 'collect')
ggplot2::ggsave('figures/general/new_infections_count.png', plot = plot_out2, height = 3.5, width = 10, dpi = 1000)

##################################
# Simple vs. sequencing analysis #
##################################

plot_list <- list()
for (i in 1:length(dataset_names)) {
    dataset_name <- dataset_names[i]
    model <- 'bayesian'
    final_df <- read.csv(paste0('real_data_outputs/', dataset_name, '_probabilities_', model, '.tsv'), sep = '\t')

    if (!'locus' %in% colnames(final_df)) {
        final_df$locus <- 1
    }

    final_df_sub <- final_df %>%
        dplyr::group_by(subject, time) %>%
        dplyr::summarize(present = 1 * any(present == 1 | present == 2), .groups = 'drop') %>%
        dplyr::mutate(allele = 1) %>%
        dplyr::arrange(time, subject, allele)

    probabilities_new_simple <- final_df_sub %>%
        determine_probabilities_simple()

    final_df_sub$probability_new <- probabilities_new_simple$probability_new

    final_df_sub <- final_df_sub %>%
        group_by(subject, allele) %>%
        filter(n() >= 3) %>%
        mutate(second_smallest_time = sort(time)[2]) %>%
        filter(time > second_smallest_time) %>%
        select(-second_smallest_time) %>%
        ungroup() %>%
        group_by(subject) %>%
        summarize(total_infections = sum(probability_new & present == 1)) %>%
        mutate(subject = as.character(subject))
    colnames(final_df_sub) <- c("subject", "infection_only")

    if (i != 3) {
        final_df_sub_2 <- read.csv(paste0('real_data_outputs/', dataset_name, '_infections_', model, '-single-locus.tsv'), sep = '\t')
        final_df_sub_2 <- data.frame(subject = rownames(final_df_sub_2), bayesian = rowMeans(final_df_sub_2))
        colnames(final_df_sub_2) <- c("subject", "single_locus")
    }

    final_df_sub_3 <- read.csv(paste0('real_data_outputs/', dataset_name, '_infections_', model, '.tsv'), sep = '\t')
    final_df_sub_3 <- data.frame(subject = rownames(final_df_sub_3), bayesian = rowMeans(final_df_sub_3))
    
    if (i != 3) {
        colnames(final_df_sub_3) <- c("subject", "multi_locus")
    } else {
        colnames(final_df_sub_3) <- c("subject", "single_locus")
    }
    
    if (i != 3) {
        final_df <- full_join(final_df_sub, final_df_sub_2, by = c('subject')) %>%
            full_join(final_df_sub_3, by = c('subject'))
    } else {
        final_df <- full_join(final_df_sub, final_df_sub_3, by = c('subject'))
    }

    final_df <- reshape2::melt(final_df, id.vars = c('subject', 'infection_only')) %>%
        mutate(variable = case_when(variable == 'single_locus' ~ "Only the most\nvariable locus",
                                    variable == 'multi_locus' ~ "All loci")) %>%
        mutate(variable = factor(variable, levels = c('Only the most\nvariable locus', 'All loci')))
    
    dataset_name_nice <- case_when(dataset_name == 'CIS43LS' ~ 'CIS43LS cohort',
                                   dataset_name == 'mali' ~ '2011 Malian cohort',
                                   dataset_name == 'uganda' ~ 'Ugandan cohort')
    
    final_df <- final_df %>%
        dplyr::mutate(infection_only = as.factor(as.character(infection_only)))

    p1 <- ggplot(final_df[sample(1:nrow(final_df)),], 
                 aes(x = variable,
                     y = value, 
                     fill = variable, 
                     group = interaction(infection_only, variable))) +
        ggbeeswarm::geom_beeswarm(method = "compactswarm", alpha = 1, cex = 2, shape = 21, size = 2, corral = 'wrap', corral.width = 0.8) +
        theme_bw() +
        labs(x = "Infections with no sequencing", 
             y = "Infection events", 
             fill = 'Sequencing type', 
             title = dataset_name_nice) +
        scale_fill_manual(
            values = c(
                "Only the most\nvariable locus" = "orangered",
                "All loci" = "maroon"
            )
        ) +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(hjust = 0.5),
            panel.spacing = unit(0,'lines')
        ) +
        geom_hline(aes(yintercept = as.numeric(infection_only) - 1), linetype = "11", color = "black", linewidth = 1) +
        scale_x_discrete(
            expand = c(0,0)) + 
        scale_y_continuous(breaks = seq(0, 10000)) + 
        facet_wrap(~infection_only, ncol = 100, strip.position = "bottom")
    
    if (i == 3) {
        p1 <- p1 + theme(legend.position = 'none')
    }

    print("Infections without sequencing:")
    print(final_df %>%
              dplyr::filter(value > 0 | as.numeric(as.character(infection_only)) > 0) %>%
              dplyr::group_by(variable) %>%
              dplyr::summarise(mean(as.numeric(as.character(infection_only)))))

    print("Infections with sequencing")
    print(final_df %>%
              dplyr::filter(value > 0 | as.numeric(as.character(infection_only)) > 0) %>%
              dplyr::group_by(variable) %>%
              dplyr::summarise(mean(value)))

    plot_list[[i]] <- p1
}
plot_out <- patchwork::wrap_plots(plot_list[[1]], plot_list[[2]], plot_list[[3]], ncol = 3, guides = 'collect', axis_titles = 'collect')
ggplot2::ggsave('figures/general/sequencing_comparison.png', plot = plot_out, height = 3.5, width = 10, dpi = 1000)

(1.92-0.920)/0.920

(1.64-0.920)/0.920

(3.45 - 0.751) / 0.751

(3.23 - 0.754) / 0.754

(0.774 - 0.583) / 0.583

#######################
# Drop-out comparison #
#######################

dataset_name <- 'mali'
df_bayesian <- read.csv(paste0('real_data_outputs/', dataset_name, '_probabilities_bayesian.tsv'), sep = '\t')
df_bayesian_drop_out <- read.csv(paste0('real_data_outputs/', dataset_name, '_probabilities_bayesian-drop-out.tsv'), sep = '\t')

df_bayesian <- df_bayesian %>%
    group_by(subject, allele) %>%
    filter(n() >= 3) %>%
    mutate(second_smallest_time = sort(time)[2]) %>%
    filter(time > second_smallest_time) %>%
    select(-second_smallest_time) %>%
    ungroup()

df_bayesian_drop_out <- df_bayesian_drop_out %>%
    group_by(subject, allele) %>%
    filter(n() >= 3) %>%
    mutate(second_smallest_time = sort(time)[2]) %>%
    filter(time > second_smallest_time) %>%
    select(-second_smallest_time) %>%
    ungroup()

joined_df <- full_join(df_bayesian %>% 
              select(locus, allele, time, subject, present, probability_new) %>% 
              rename(no_drop_out = probability_new), 
          df_bayesian_drop_out %>% 
              select(locus, allele, time, subject, present, probability_new) %>% 
              rename(drop_out = probability_new)) %>% 
    filter(present == 1)

p1 <- ggplot(joined_df, aes(x = no_drop_out, y = drop_out)) + 
    geom_point(alpha = 0.5) + 
    labs(x = 'Probability new with drop-out mode off',
         y = 'Probability new with drop-out mode on') + 
    theme_bw() + 
    scale_x_continuous(breaks = seq(0, 1, 0.2)) + 
    scale_y_continuous(breaks = seq(0, 1, 0.2)) + 
    geom_abline(intercept = 0, slope = 1, color = 'red')
p1 <- ggExtra::ggMarginal(p1, type = "histogram", margins = 'both')

joined_df <- full_join(compute_molFOI(df_bayesian, method = 'sum_then_max') %>% 
                           rename(no_drop_out = molFOI), 
                       compute_molFOI(df_bayesian_drop_out, method = 'sum_then_max') %>% 
                           rename(drop_out = molFOI))

mean(joined_df$no_drop_out)
mean(joined_df$drop_out)

p2 <- ggplot(joined_df, aes(x = no_drop_out, y = drop_out)) + 
    geom_point(alpha = 0.5) + 
    labs(x = 'molFOI with drop-out mode off',
         y = 'molFOI with drop-out mode on') + 
    theme_bw() + 
    scale_x_continuous(limits = c(0, 32)) + 
    scale_y_continuous(limits = c(0, 32)) + 
    geom_abline(intercept = 0, slope = 1, color = 'red') + 
    coord_fixed()
p2 <- ggExtra::ggMarginal(p2, type = "histogram", margins = 'both')

df_bayesian <- read.csv(paste0('real_data_outputs/', dataset_name, '_infections_bayesian.tsv'), sep = '\t')
df_bayesian <- data.frame(subject = rownames(df_bayesian), no_drop_out = rowMeans(df_bayesian))

df_bayesian_drop_out <- read.csv(paste0('real_data_outputs/', dataset_name, '_infections_bayesian-drop-out.tsv'), sep = '\t')
df_bayesian_drop_out <- data.frame(subject = rownames(df_bayesian_drop_out), drop_out = rowMeans(df_bayesian_drop_out))

joined_df <- full_join(df_bayesian, df_bayesian_drop_out)

mean(joined_df$no_drop_out)
mean(joined_df$drop_out)

p3 <- ggplot(joined_df, aes(x = no_drop_out, y = drop_out)) + 
    geom_point(alpha = 0.5) + 
    labs(x = 'New infections with drop-out mode off',
         y = 'New infections with drop-out mode on') + 
    theme_bw() + 
    scale_x_continuous(limits = c(0, 9), breaks = 0:9) + 
    scale_y_continuous(limits = c(0, 9), breaks = 0:9) + 
    geom_abline(intercept = 0, slope = 1, color = 'red') + 
    coord_fixed()
p3 <- ggExtra::ggMarginal(p3, type = "histogram", margins = 'both')

plot_out <- patchwork::wrap_plots(p1, p2, p3, ncol = 3, guides = 'collect', axis_titles = 'collect')
ggplot2::ggsave('figures/general/drop_out.png', plot = plot_out, height = 3.5, width = 10, dpi = 1000)

###################
# In-text numbers #
###################

dataset_name <- 'uganda'
model <- 'bayesian'
final_df <- read.csv(paste0('real_data_outputs/', dataset_name, '_probabilities_', model, '.tsv'), sep = '\t')
final_df %>%
    arrange(subject, allele, time) %>%
    group_by(subject, allele) %>%
    summarize(seq_course = grepl("1(0{3,})1", paste0(present, collapse = ''))) %>%
    ungroup() %>%
    summarize(output = sum(seq_course)) %>%
    select(output)

final_df %>%
    arrange(subject, allele, time) %>%
    group_by(subject, allele) %>%
    summarize(seq_course = grepl("1([01]{3,})1", paste0(present, collapse = ''))) %>%
    ungroup() %>%
    summarize(output = sum(seq_course)) %>%
    select(output)

final_df %>%
    arrange(subject, allele, time) %>%
    group_by(subject, allele) %>%
    summarize(seq_course = grepl("1([02]{3,})1", paste0(present, collapse = ''))) %>%
    ungroup() %>%
    summarize(output = sum(seq_course)) %>%
    select(output)

final_df %>% group_by(subject) %>% summarize(time_max = max(time) - min(time)) %>% summarize(mean(time_max))



