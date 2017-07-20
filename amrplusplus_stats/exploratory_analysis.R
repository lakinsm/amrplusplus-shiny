require(data.table)

source('graphing.R')

alpha_rarefaction_workflow <- function(experiment_data_lists,
                                       sample_column_id,
                                       exploratory_analyses,
                                       graph_output_dir) {
    for(v in 1:length(exploratory_analyses)) {
        # AMR
        if(!is.na(experiment_data_lists[[1]])) {
            meg_alpha_rarefaction(data_list=experiment_data_lists[[2]],
                                  data_names=experiment_data_lists[[5]],
                                  metadata=experiment_data_lists[[11]],
                                  sample_var=sample_column_id,
                                  group_var=exploratory_analyses[[v]]$exploratory_var,
                                  analysis_subset=exploratory_analyses[[v]]$subsets,
                                  outdir=paste(graph_output_dir, 'AMR',
                                               exploratory_analyses[[v]]$name,
                                               sep='/', collapse=''),
                                  data_type='AMR')
        }
        
        # Microbiome
        if(!is.na(experiment_data_lists[[6]])) {
            meg_alpha_rarefaction(data_list=experiment_data_lists[[7]],
                                  data_names=experiment_data_lists[[10]],
                                  metadata=experiment_data_lists[[11]],
                                  sample_var=sample_column_id,
                                  group_var=exploratory_analyses[[v]]$exploratory_var,
                                  analysis_subset=exploratory_analyses[[v]]$subsets,
                                  outdir=paste(graph_output_dir, 'Microbiome',
                                               exploratory_analyses[[v]]$name,
                                               sep='/', collapse=''),
                                  data_type='Microbiome')
        }
    }
}


ordination_workflow <- function(experiment_data_lists,
                                sample_column_id,
                                exploratory_analyses,
                                graph_output_dir) {
    if(!is.na(experiment_data_lists[[1]])) {
        # AMR NMDS
        meg_ordination(data_list=experiment_data_lists[[1]],
                       data_names=experiment_data_lists[[5]],
                       metadata=experiment_data_lists[[11]],
                       sample_var=sample_column_id,
                       hull_var=exploratory_analyses[[v]]$exploratory_var,
                       analysis_subset=exploratory_analyses[[v]]$subsets,
                       outdir=paste(graph_output_dir, 'AMR',
                                    exploratory_analyses[[v]]$name,
                                    sep='/', collapse=''),
                       data_type='AMR',
                       method='NMDS')
        
        # AMR PCA
        meg_ordination(data_list=experiment_data_lists[[1]],
                       data_names=experiment_data_lists[[5]],
                       metadata=experiment_data_lists[[11]],
                       sample_var=sample_column_id,
                       hull_var=exploratory_analyses[[v]]$exploratory_var,
                       analysis_subset=exploratory_analyses[[v]]$subsets,
                       outdir=paste(graph_output_dir, 'AMR',
                                    exploratory_analyses[[v]]$name,
                                    sep='/', collapse=''),
                       data_type='AMR',
                       method='PCA')
    }
    
    if(!is.na(experiment_data_lists[[6]])) {
        # Microbiome NMDS
        meg_ordination(data_list=experiment_data_lists[[6]],
                       data_names=experiment_data_lists[[10]],
                       metadata=experiment_data_lists[[11]],
                       sample_var=sample_column_id,
                       hull_var=exploratory_analyses[[v]]$exploratory_var,
                       analysis_subset=exploratory_analyses[[v]]$subsets,
                       outdir=paste(graph_output_dir, 'Microbiome',
                                    exploratory_analyses[[v]]$name,
                                    sep='/', collapse=''),
                       data_type='Microbiome',
                       method='NMDS')
        
        # Microbiome PCA
        meg_ordination(data_list=experiment_data_lists[[6]],
                       data_names=experiment_data_lists[[10]],
                       metadata=experiment_data_lists[[11]],
                       sample_var=sample_column_id,
                       hull_var=exploratory_analyses[[v]]$exploratory_var,
                       analysis_subset=exploratory_analyses[[v]]$subsets,
                       outdir=paste(graph_output_dir, 'AMR',
                                    exploratory_analyses[[v]]$name,
                                    sep='/', collapse=''),
                       data_type='Microbiome',
                       method='PCA')
    }
}


heatmap_workflow <- function(experiment_data_lists,
                             sample_column_id,
                             exploratory_analyses,
                             graph_output_dir) {
    if(!is.na(experiment_data_lists[[1]])) {
        # AMR Heatmaps for each level
        for( v in 1:length(exploratory_analyses) ) {
            for( l in 1:length(experiment_data_lists[[5]]) ) {
                meg_heatmap(melted_data=experiment_data_lists[[3]],
                            metadata=metadata,
                            sample_var=sample_column_id,
                            group_var=exploratory_analyses[[v]]$exploratory_var,
                            level_var=experiment_data_lists[[5]][l],
                            analysis_subset=exploratory_analyses[[v]]$subsets,
                            outdir=paste(graph_output_dir, 'AMR', exploratory_analyses[[v]]$name,
                                         sep='/', collapse=''),
                            data_type='AMR')
            }
        }
    }
    
    if(!is.na(experiment_data_lists[[6]])) {
        # Microbiome Heatmaps for each level
        for( v in 1:length(exploratory_analyses) ) {
            for( l in 1:length(experiment_data_lists[[10]]) ) {
                meg_heatmap(melted_data=experiment_data_lists[[8]],
                            metadata=metadata,
                            sample_var=sample_column_id,
                            group_var=exploratory_analyses[[v]]$exploratory_var,
                            level_var=experiment_data_lists[[10]][l],
                            analysis_subset=exploratory_analyses[[v]]$subsets,
                            outdir=paste(graph_output_dir, 'Microbiome', exploratory_analyses[[v]]$name,
                                         sep='/', collapse=''),
                            data_type='Microbiome')
            }
        }
    }
}

barplot_workflow <- function(experiment_data_lists,
                             sample_column_id,
                             exploratory_analyses,
                             graph_output_dir) {
    if(!is.na(experiment_data_lists[[1]])) {
        # AMR barplots for each level
        for( v in 1:length(exploratory_analyses) ) {
            for( l in 1:length(experiment_data_lists[[5]]) ) {
                suppressWarnings(
                    meg_barplot(melted_data=experiment_data_lists[[3]],
                                metadata=metadata,
                                sample_var=sample_column_id,
                                group_var=exploratory_analyses[[v]]$exploratory_var,
                                level_var=experiment_data_lists[[5]][l],
                                analysis_subset=exploratory_analyses[[v]]$subsets,
                                outdir=paste(graph_output_dir, 'AMR', exploratory_analyses[[v]]$name,
                                             sep='/', collapse=''),
                                data_type='AMR')
                )
            }
        }
    }
    
    if(!is.na(experiment_data_lists[[6]])) {
        # Microbiome barplots for each level
        for( v in 1:length(exploratory_analyses) ) {
            for( l in 1:length(experiment_data_lists[[10]]) ) {
                suppressWarnings(
                    meg_barplot(melted_data=experiment_data_lists[[8]],
                                metadata=metadata,
                                sample_var=sample_column_id,
                                group_var=exploratory_analyses[[v]]$exploratory_var,
                                level_var=experiment_data_lists[[10]][l],
                                analysis_subset=exploratory_analyses[[v]]$subsets,
                                outdir=paste(graph_output_dir, 'AMR', exploratory_analyses[[v]]$name,
                                             sep='/', collapse=''),
                                data_type='AMR')
                )
            }
        }
    }
}


























