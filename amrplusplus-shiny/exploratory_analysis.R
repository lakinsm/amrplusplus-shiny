require(data.table)

source('graphing.R')

alpha_rarefaction_workflow <- function(experiment_data_lists,
                                       sample_column_id,
                                       feature_name,
                                       graph_output_dir,
                                       data_type) {
    if(data_type == 'Resistome') {
        if(!is.na(experiment_data_lists[[1]])) {
            meg_alpha_rarefaction(data_list=experiment_data_lists[[2]],
                                  data_names=experiment_data_lists[[5]],
                                  metadata=experiment_data_lists[[11]],
                                  sample_var=sample_column_id,
                                  group_var=feature_name,
                                  analysis_subset=list(),
                                  outdir=paste(graph_output_dir, data_type,
                                               feature_name,
                                               sep='/', collapse=''),
                                  data_type=data_type)
        }
    }
    else if(data_type == 'Microbiome') {
        if(!is.na(experiment_data_lists[[6]])) {
            meg_alpha_rarefaction(data_list=experiment_data_lists[[7]],
                                  data_names=experiment_data_lists[[10]],
                                  metadata=experiment_data_lists[[11]],
                                  sample_var=sample_column_id,
                                  group_var=feature_name,
                                  analysis_subset=list(),
                                  outdir=paste(graph_output_dir, data_type,
                                               feature_name,
                                               sep='/', collapse=''),
                                  data_type=data_type)
        }
    }
}


ordination_workflow <- function(experiment_data_lists,
                                sample_column_id,
                                feature_name,
                                graph_output_dir,
                                data_type) {
    if(data_type == 'Resistome') {
        if(!is.na(experiment_data_lists[[1]])) {
            # AMR NMDS
            meg_ordination(data_list=experiment_data_lists[[1]],
                           data_names=experiment_data_lists[[5]],
                           metadata=experiment_data_lists[[11]],
                           sample_var=sample_column_id,
                           hull_var=feature_name,
                           analysis_subset=list(),
                           outdir=paste(graph_output_dir, data_type,
                                        feature_name,
                                        sep='/', collapse=''),
                           data_type=data_type,
                           method='NMDS')
            
            # AMR PCA
            meg_ordination(data_list=experiment_data_lists[[1]],
                           data_names=experiment_data_lists[[5]],
                           metadata=experiment_data_lists[[11]],
                           sample_var=sample_column_id,
                           hull_var=feature_name,
                           analysis_subset=list(),
                           outdir=paste(graph_output_dir, data_type,
                                        feature_name,
                                        sep='/', collapse=''),
                           data_type=feature_name,
                           method='PCA')
        }
    }
    else if(data_type == 'Microbiome') {
        if(!is.na(experiment_data_lists[[6]])) {
            # Microbiome NMDS
            meg_ordination(data_list=experiment_data_lists[[6]],
                           data_names=experiment_data_lists[[10]],
                           metadata=experiment_data_lists[[11]],
                           sample_var=sample_column_id,
                           hull_var=feature_name,
                           analysis_subset=list(),
                           outdir=paste(graph_output_dir, data_type,
                                        feature_name,
                                        sep='/', collapse=''),
                           data_type=data_type,
                           method='NMDS')
            
            # Microbiome PCA
            meg_ordination(data_list=experiment_data_lists[[6]],
                           data_names=experiment_data_lists[[10]],
                           metadata=experiment_data_lists[[11]],
                           sample_var=sample_column_id,
                           hull_var=feature_name,
                           analysis_subset=list(),
                           outdir=paste(graph_output_dir, data_type,
                                        feature_name,
                                        sep='/', collapse=''),
                           data_type=data_type,
                           method='PCA')
        }
    }
}


heatmap_workflow <- function(experiment_data_lists,
                             sample_column_id,
                             feature_name,
                             graph_output_dir,
                             data_type) {
    if(data_type == 'Resistome') {
        if(!is.na(experiment_data_lists[[1]])) {
            # AMR Heatmaps for each level
            for( l in 1:length(experiment_data_lists[[5]]) ) {
                meg_heatmap(melted_data=experiment_data_lists[[3]],
                            metadata=experiment_data_lists[[11]],
                            sample_var=sample_column_id,
                            group_var=feature_name,
                            level_var=experiment_data_lists[[5]][l],
                            analysis_subset=list(),
                            outdir=paste(graph_output_dir, data_type,
                                         feature_name,
                                         sep='/', collapse=''),
                            data_type=data_type)
            }
        }
    }
    else if(data_type == 'Microbiome') {
        if(!is.na(experiment_data_lists[[6]])) {
            # Microbiome Heatmaps for each level
            for( l in 1:length(experiment_data_lists[[10]]) ) {
                meg_heatmap(melted_data=experiment_data_lists[[8]],
                            metadata=experiment_data_lists[[11]],
                            sample_var=sample_column_id,
                            group_var=feature_name,
                            level_var=experiment_data_lists[[10]][l],
                            analysis_subset=list(),
                            outdir=paste(graph_output_dir, data_type,
                                         feature_name,
                                         sep='/', collapse=''),
                            data_type=data_type)
            }
        }
    }
}

barplot_workflow <- function(experiment_data_lists,
                             sample_column_id,
                             feature_name,
                             graph_output_dir,
                             data_type) {
    if(data_type == 'Resistome') {
        if(!is.na(experiment_data_lists[[1]])) {
            # AMR barplots for each level
            for( l in 1:length(experiment_data_lists[[5]]) ) {
                suppressWarnings(
                    meg_barplot(melted_data=experiment_data_lists[[3]],
                                metadata=experiment_data_lists[[11]],
                                sample_var=sample_column_id,
                                group_var=feature_name,
                                level_var=experiment_data_lists[[5]][l],
                                analysis_subset=list(),
                                outdir=paste(graph_output_dir, data_type,
                                             feature_name,
                                             sep='/', collapse=''),
                                data_type=data_type)
                )
            }
        }
    }
    else if(data_type == 'Microbiome') {
        if(!is.na(experiment_data_lists[[6]])) {
            # Microbiome barplots for each level
            for( l in 1:length(experiment_data_lists[[10]]) ) {
                suppressWarnings(
                    meg_barplot(melted_data=experiment_data_lists[[8]],
                                metadata=experiment_data_lists[[11]],
                                sample_var=sample_column_id,
                                group_var=feature_name,
                                level_var=experiment_data_lists[[10]][l],
                                analysis_subset=list(),
                                outdir=paste(graph_output_dir, data_type,
                                             feature_name,
                                             sep='/', collapse=''),
                                data_type=data_type)
                )
            }
        }
    }
}
