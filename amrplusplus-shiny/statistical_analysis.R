require(data.table)
require(metagenomeSeq)

generate_statistical_output <- function(data,
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
}


generate_contrasts <- function(analytic_MRexp,
                               feature) {
    values <- data.table(pData(analytic_MRexp))
    values <- as.character(unique(values[[feature]]))
    combs <- utils::combn(values, 2)
    strs <- apply(combs, 2, function(x) {
        paste(paste(feature, x[1], sep='', collapse=''),
              ' - ',
              paste(feature, x[2], sep='', collapse=''),
              sep='', collapse='')
    })
    return(strs)
}


exploratory_zero_inflated_gaussian_regression <- function(analytic_MRexp,
                                                          metadata,
                                                          annotation_level,
                                                          metadata_feature,
                                                          zero_mod,
                                                          data_mod,
                                                          random_effect_var,
                                                          data_type,
                                                          pval,
                                                          top_hits,
                                                          sort_by) {
    settings <- zigControl(maxit=50, verbose=F)
    
    local_MRexp <- analytic_MRexp
    
    col_selection <- as.integer(which(colSums(MRcounts(local_MRexp) > 0) > 1))
    local_MRexp <- local_MRexp[, col_selection]
    pData
    
    mod_select <- model.matrix(eval(parse(text=data_mod)), data=pData(local_MRexp))
    zero_mod_select <- zero_mod[col_selection, ]
    
    cumNorm(local_MRexp)  # This is a placeholder for metagenomeSeq; we don't actually use these values
    
    res <- NA
    
    tryCatch(
        {
            if( is.na(random_effect_var) ) {
                res <- fitZig(obj=local_MRexp,
                                   mod=mod_select,
                                   zeroMod=zero_mod_select,
                                   control=settings,
                                   useCSSoffset=F)
            }
            else {
                res <- fitZig(obj=local_MRexp,
                                   mod=mod_select,
                                   zeroMod=zero_mod_select,
                                   control=settings,
                                   useCSSoffset=F,
                                   useMixedModel=T,
                                   block=pData(local_MRexp)[, random_effect_var])
            }
        },
        error=function(e) {
            print(paste('Model failed to converge for ', data_type, ' ', annotation_level,
                        sep='', collapse=''))
        }
    )
    
    contrast_list <- generate_contrasts(local_MRexp,
                                        metadata_feature)
    
    local_contrasts <- as.list(contrast_list)
    local_contrasts[[length(local_contrasts)+1]] <- res$fit$design
    names(local_contrasts)[length(local_contrasts)] <- 'levels'
    
    contrast_matrix <- do.call(makeContrasts, local_contrasts)
    colnames(contrast_matrix) <- make.names(contrast_list)
    
    contrast_fit <- contrasts.fit(res$fit, contrast_matrix)        
    contrast_fit <- eBayes(contrast_fit)
    
    stats_results <- data.table(
        Node.Name=character(),
        Contrast=character(),
        logFC=numeric(),
        CI.L=numeric(),
        CI.R=numeric(),
        AveExpr=numeric(),
        t=numeric(),
        P.Value=numeric(),
        adj.P.Val=numeric(),
        B=numeric()
    )
    
    local_sort <- NA
    # 'P-value', 'Effect Size', 'Abundance', 'F-statistic'
    if(sort_by == 'P-value') {
        local_sort <- 'P'
    }
    else if(sort_by == 'Effect Size') {
        local_sort <- 'M'
    }
    else if(sort_by == 'Abundance') {
        local_sort <- 'A'
    }
    else if(sort_by == 'T-statistic') {
        local_sort <- 't'
    }
    print(pval)
    print(top_hits)
    print(local_sort)
    print(c)
    for( c in 1:ncol(contrast_fit$contrasts) ) {
        tophits <- topTable(contrast_fit, p.value=pval, confint=T,
                            number=top_hits, sort.by=local_sort, coef=c)
        
        if( nrow(tophits) > 0) {
            temp_res <- data.table(
                Node.Name=rownames(tophits),
                Contrast=rep(colnames(contrast_fit$contrasts)[c], nrow(tophits))
            )
            temp_res <- cbind(temp_res, tophits)
            stats_results <- rbind(stats_results, temp_res)
        }
        else {
            print(paste('No significant results for', data_type,
                        annotation_level, colnames(contrast_fit$contrasts)[c],
                        sep=' ', collapse=''))
        }
    }
    
    return(stats_results)
}









