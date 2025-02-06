library(dinemites)

treatments <- read.csv('real_data_outputs/treatments_uganda.tsv', sep = '\t')
for (model in c('bayesian', 'clustering', 'simple')) {
    if (model == 'bayesian') {
        figures_folder <- 'figures/uganda/bayesian/'
        final_df <- read.csv('real_data_outputs/uganda_probabilities_bayesian.tsv', sep = '\t')
        infections_df <- read.csv('real_data_outputs/uganda_infections_bayesian.tsv', sep = '\t')
    } else if (model == 'clustering') {
        figures_folder <- 'figures/uganda/clustering/'
        final_df <- read.csv('real_data_outputs/uganda_probabilities_clustering.tsv', sep = '\t')
        infections_df <- read.csv('real_data_outputs/uganda_infections_clustering.tsv', sep = '\t')
    } else if (model == 'simple') {
        figures_folder <- 'figures/uganda/simple/'
        final_df <- read.csv('real_data_outputs/uganda_probabilities_simple.tsv', sep = '\t')
        infections_df <- read.csv('real_data_outputs/uganda_infections_simple.tsv', sep = '\t')
    }
    dir.create(figures_folder, showWarnings = F, recursive = T)
    
    plot_dataset(dataset = final_df, 
                 treatments = treatments,
                 estimated_new_infections = infections_df, 
                 output = figures_folder)
}
