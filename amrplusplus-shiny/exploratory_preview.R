# Data preview code for exploratory preview section

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
        print(length(analytic_list))
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
        print(length(analytic_list))
    }
}