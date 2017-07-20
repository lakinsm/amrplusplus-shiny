require(data.table)
require(metagenomeSeq)

meg_fitZig <- function(data_list,
                       data_names,
                       metadata,
                       zero_mod,
                       data_mod,
                       filter_min_threshold,
                       contrast_list,
                       random_effect_var,
                       outdir,
                       analysis_name,
                       analysis_subset,
                       data_type,
                       pval,
                       top_hits) {
    settings <- zigControl(maxit=50, verbose=F)
    
    local_obj <- data_list
    res <- list()
    for( l in 1:length(local_obj) ) {
        
        filter_threshold <- quantile(rowSums(MRcounts(local_obj[[l]])), 0.15)
        if( filter_threshold > filter_min_threshold ) filter_threshold <- filter_min_threshold
        local_obj[[l]] <- local_obj[[l]][which(rowSums(MRcounts(local_obj[[l]])) >= filter_threshold ), ]
        
        if(length(analysis_subset) > 0) {
            local_obj[[l]] <- data_subset(local_obj[[l]], analysis_subset)
        }
        
        
        col_selection <- as.integer(which(colSums(MRcounts(local_obj[[l]]) > 0) > 1))
        local_obj[[l]] <- local_obj[[l]][, col_selection]
        
        mod_select <- model.matrix(eval(parse(text=data_mod)), data=pData(local_obj[[l]]))
        zero_mod_select <- zero_mod[col_selection, ]
        
        cumNorm(local_obj[[l]])  # This is a placeholder for metagenomeSeq; we don't actually use these values
        
        tryCatch(
            {
                if( is.na(random_effect_var) ) {
                    res[[l]] <- fitZig(obj=local_obj[[l]],
                                       mod=mod_select,
                                       zeroMod=zero_mod_select,
                                       control=settings,
                                       useCSSoffset=F)
                }
                else {
                    res[[l]] <- fitZig(obj=local_obj[[l]],
                                       mod=mod_select,
                                       zeroMod=zero_mod_select,
                                       control=settings,
                                       useCSSoffset=F,
                                       useMixedModel=T,
                                       block=pData(local_obj[[l]])[, random_effect_var])
                }
            },
            error=function(e) {
                print(paste('Model failed to converge for ', data_type, ' ', data_names[l], ' ', analysis_name,
                            sep='', collapse=''))
            },
            finally={
                if( length(res) != l ) {
                    next
                }
            }
        )
        
        local_contrasts <- contrast_list
        local_contrasts[[length(local_contrasts)+1]] <- res[[l]]$fit$design
        names(local_contrasts)[length(local_contrasts)] <- 'levels'
        
        contrast_matrix <- do.call(makeContrasts, local_contrasts)
        colnames(contrast_matrix) <- make.names(contrast_list)
        
        contrast_fit <- contrasts.fit(res[[l]]$fit, contrast_matrix)        
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
        
        for( c in 1:ncol(contrast_fit$contrasts) ) {
            tophits <- topTable(contrast_fit, p.value=pval, confint=T,
                                number=top_hits, sort.by='AveExpr', coef=c)
            
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
                            data_names[l], analysis_name,
                            colnames(contrast_fit$contrasts)[c],
                            sep=' ', collapse=''))
            }
        }
        
        if( nrow(stats_results) > 0 ) {
            write.csv(stats_results,
                      file=paste(outdir, '/', analysis_name, '_', data_type, '_',
                                 data_names[l], '_',
                                 contrast_list[1], '_Model_Contrasts.csv',
                                 sep='', collapse=''),
                      quote=F, row.names=F)
        }
    }
}


zero_inflated_gaussian_regression_workflow <- function(experiment_data_lists,
                                                       original_MRexps,
                                                       sample_column_id,
                                                       statistical_analyses,
                                                       stats_output_dir,
                                                       sig_level,
                                                       top_hits,
                                                       lowpass_filter_threshold) {
    if(!is.na(experiment_data_lists[[1]])) {
        # AMR ZIG Regression
        for( a in 1:length(statistical_analyses) ) {
            meg_fitZig(data_list=experiment_data_lists[[1]],
                       data_names=experiment_data_lists[[5]],
                       metadata=experiment_data_lists[[11]],
                       zero_mod=model.matrix(~1 + log(libSize(original_MRexps[[1]]))),
                       data_mod=statistical_analyses[[a]]$model_matrix,
                       filter_min_threshold=lowpass_filter_threshold,
                       contrast_list=statistical_analyses[[a]]$contrasts,
                       random_effect_var=statistical_analyses[[a]]$random_effect,
                       outdir=paste(stats_output_dir, 'AMR', statistical_analyses[[a]]$name,
                                    sep='/', collapse=''),
                       analysis_name=statistical_analyses[[a]]$name,
                       analysis_subset=statistical_analyses[[a]]$subsets,
                       data_type='AMR',
                       pval=sig_level,
                       top_hits=top_hits)
        }
    }
    
    if(!is.na(experiment_data_lists[[6]])) {
        # Microbiome ZIG Regression
        for( a in 1:length(statistical_analyses) ) {
            meg_fitZig(data_list=experiment_data_lists[[6]],
                       data_names=experiment_data_lists[[10]],
                       metadata=experiment_data_lists[[11]],
                       zero_mod=model.matrix(~1 + log(libSize(original_MRexps[[2]]))),
                       data_mod=statistical_analyses[[a]]$model_matrix,
                       filter_min_threshold=lowpass_filter_threshold,
                       contrast_list=statistical_analyses[[a]]$contrasts,
                       random_effect_var=statistical_analyses[[a]]$random_effect,
                       outdir=paste(stats_output_dir, 'Microbiome', statistical_analyses[[a]]$name,
                                    sep='/', collapse=''),
                       analysis_name=statistical_analyses[[a]]$name,
                       analysis_subset=statistical_analyses[[a]]$subsets,
                       data_type='Microbiome',
                       pval=sig_level,
                       top_hits=top_hits)
        }
    }
}


elastic_net_regression_workflow <- function(experiment_data_lists,
                                            sample_column_id,
                                            statistical_analyses,
                                            stats_output_dir) {
    
}























