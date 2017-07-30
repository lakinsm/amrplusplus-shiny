# Data preview code for exploratory preview section


microbiome_lookup <- c(
    'Domain' = 1,
    'Phylum' = 2,
    'Class' = 3,
    'Order' = 4,
    'Family' = 5,
    'Genus' = 6,
    'Species' = 7
)

amr_lookup <- c(
    'Class' = 1,
    'Mechanism' = 2,
    'Group' = 3,
    'Gene' = 4
)

generate_exploratory_preview <- function(data,
                                         data_type,
                                         metadata,
                                         annotation_level,
                                         analysis_type,
                                         metadata_feature,
                                         low_pass_filter_threshold,
                                         norm_method,
                                         sample_depth,
                                         subset_strings) {
    out_plot <- character(0)
    
    if(data_type == 'Resistome') {
        
        # Normalization
        if(norm_method == 'Rarefaction') {
            amr_data_list <- rarefy_normalize_and_extract(sampling_depth=sample_depth,
                                                          filtering_quantile=low_pass_filter_threshold,
                                                          amr_MRexp=data[[1]])
            amr_data_list[[7]] <- data[[2]]
        }
        else if(norm_method == 'Cumulative Sum Scaling') {
            amr_data_list <- CSS_normalize_and_extract(filtering_quantile=low_pass_filter_threshold,
                                                       amr_MRexp=data[[1]])
            amr_data_list[[7]] <- data[[2]]
        }
        
        # Aggregation and data set creation
        analytic_list <- aggregate_and_filter(amr_data_list, data.table(metadata))
        
        # Plot generation
        if(analysis_type == 'NMDS' || analysis_type == 'PCA') {
            out_plot <- exploratory_ordination(analytic_MRexp=analytic_list[[1]][[amr_lookup[annotation_level]]],
                                               metadata=analytic_list[[11]],
                                               annotation_level=annotation_level,
                                               feature_var=colnames(analytic_list[[11]])[1],
                                               hull_var=metadata_feature,
                                               analysis_subset=subset_strings,
                                               data_type='Resistome',
                                               method=analysis_type)
        }
        else if(analysis_type == 'Bar Graph') {
            out_plot <- exploratory_barplot(melted_data=analytic_list[[3]],
                                            metadata=analytic_list[[11]],
                                            sample_var=colnames(analytic_list[[11]])[1],
                                            group_var=metadata_feature,
                                            level_var=annotation_level,
                                            analysis_subset=subset_strings,
                                            data_type='Resistome')
        }
        else if(analysis_type == 'Heatmap') {
            out_plot <- exploratory_heatmap(melted_data=analytic_list[[3]],
                                            metadata=analytic_list[[11]],
                                            sample_var=colnames(analytic_list[[11]])[1],
                                            group_var=metadata_feature,
                                            level_var=annotation_level,
                                            analysis_subset=subset_strings,
                                            data_type='Resistome')
        }
    }
    else if(data_type == 'Microbiome') {
        
        # Normalization
        if(norm_method == 'Rarefaction') {
            microbiome_data_list <- rarefy_normalize_and_extract(sampling_depth=sample_depth,
                                                                 filtering_quantile=low_pass_filter_threshold,
                                                                 kraken_MRexp=data[[1]])
            microbiome_data_list[[7]] <- NA
        }
        else if(norm_method == 'Cumulative Sum Scaling') {
            microbiome_data_list <- CSS_normalize_and_extract(filtering_quantile=low_pass_filter_threshold,
                                                              kraken_MRexp=data[[1]])
            microbiome_data_list[[7]] <- NA
        }
        
        # Aggregation and data set creation
        analytic_list <- aggregate_and_filter(microbiome_data_list, data.table(metadata))
        
        # Plot generation
        if(analysis_type == 'NMDS' || analysis_type == 'PCA') {
            out_plot <- exploratory_ordination(analytic_MRexp=analytic_list[[6]][[microbiome_lookup[annotation_level]]],
                                               metadata=analytic_list[[11]],
                                               annotation_level=annotation_level,
                                               feature_var=colnames(analytic_list[[11]])[1],
                                               hull_var=metadata_feature,
                                               analysis_subset=subset_strings,
                                               data_type='Microbiome',
                                               method=analysis_type)
        }
        else if(analysis_type == 'Bar Graph') {
            out_plot <- exploratory_barplot(melted_data=analytic_list[[8]],
                                            metadata=analytic_list[[11]],
                                            sample_var=colnames(analytic_list[[11]])[1],
                                            group_var=metadata_feature,
                                            level_var=annotation_level,
                                            analysis_subset=subset_strings,
                                            data_type='Microbiome')
        }
        else if(analysis_type == 'Heatmap') {
            out_plot <- exploratory_heatmap(melted_data=analytic_list[[8]],
                                            metadata=analytic_list[[11]],
                                            sample_var=colnames(analytic_list[[11]])[1],
                                            group_var=metadata_feature,
                                            level_var=annotation_level,
                                            analysis_subset=subset_strings,
                                            data_type='Microbiome')
        }
    }
    
    return(out_plot)
}


generate_statistical_preview <- function(data,
                                         data_type,
                                         metadata,
                                         annotation_level,
                                         analysis_type,
                                         metadata_feature,
                                         low_pass_filter_threshold,
                                         norm_method,
                                         sample_depth,
                                         subset_strings,
                                         model_string,
                                         random_effect,
                                         pval_threshold,
                                         num_tophits,
                                         sort_by) {
    if(data_type == 'Resistome') {
        
        # Normalization
        if(norm_method == 'Rarefaction') {
            amr_data_list <- rarefy_normalize_and_extract(sampling_depth=sample_depth,
                                                          filtering_quantile=low_pass_filter_threshold,
                                                          amr_MRexp=data[[1]])
            amr_data_list[[7]] <- data[[2]]
        }
        else if(norm_method == 'Cumulative Sum Scaling') {
            amr_data_list <- CSS_normalize_and_extract(filtering_quantile=low_pass_filter_threshold,
                                                       amr_MRexp=data[[1]])
            amr_data_list[[7]] <- data[[2]]
        }
        
        # Aggregation and data set creation
        analytic_list <- aggregate_and_filter(amr_data_list, data.table(metadata))
        
        # Regression
        if(analysis_type == 'ZIG Regression') {
            out_dt <- exploratory_zero_inflated_gaussian_regression(analytic_MRexp=analytic_list[[1]][[amr_lookup[annotation_level]]],
                                                                    metadata=analytic_list[[11]],
                                                                    annotation_level=annotation_level,
                                                                    metadata_feature=metadata_feature,
                                                                    zero_mod=model.matrix(~1 + log(libSize(data[[1]]))),
                                                                    data_mod=model_string,
                                                                    random_effect_var=random_effect,
                                                                    data_type=data_type,
                                                                    pval=pval_threshold,
                                                                    top_hits=num_tophits,
                                                                    sort_by=sort_by)
        }
        else if(analysis_type == 'Elastic Net Regression') {
            
        }
    }
    else if(data_type == 'Microbiome') {
        
        # Normalization
        if(norm_method == 'Rarefaction') {
            microbiome_data_list <- rarefy_normalize_and_extract(sampling_depth=sample_depth,
                                                                 filtering_quantile=low_pass_filter_threshold,
                                                                 kraken_MRexp=data[[1]])
            microbiome_data_list[[7]] <- NA
        }
        else if(norm_method == 'Cumulative Sum Scaling') {
            microbiome_data_list <- CSS_normalize_and_extract(filtering_quantile=low_pass_filter_threshold,
                                                              kraken_MRexp=data[[1]])
            microbiome_data_list[[7]] <- NA
        }
        
        # Aggregation and data set creation
        analytic_list <- aggregate_and_filter(microbiome_data_list, data.table(metadata))
    }
    
    # Regression
    if(analysis_type == 'ZIG Regression') {
        out_dt <- exploratory_zero_inflated_gaussian_regression(analytic_MRexp=analytic_list[[6]][[microbiome_lookup[annotation_level]]],
                                                                metadata=analytic_list[[11]],
                                                                annotation_level=annotation_level,
                                                                metadata_feature=metadata_feature,
                                                                zero_mod=model.matrix(~1 + log(libSize(data[[1]]))),
                                                                data_mod=model_string,
                                                                random_effect_var=random_effect,
                                                                data_type=data_type,
                                                                pval=pval_threshold,
                                                                top_hits=num_tophits,
                                                                sort_by=sort_by)
    }
    else if(analysis_type == 'Elastic Net Regression') {
        
    }
}
















