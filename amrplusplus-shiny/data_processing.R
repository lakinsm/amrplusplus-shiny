require(data.table)
require(metagenomeSeq)

snp_regex = c('ACRR',
              'CATB',
              'CLS',
              'DFRC',
              'DHFR',
              'DHFRIII',
              'DHFRIX',
              'EMBA',
              'embB',
              'EMBB',
              'EMBC',
              'EMBR',
              'ETHA',
              'FOLP',
              'GIDB',
              'GYRA',
              'gyrB',
              'GYRB',
              'INHA',
              'INIA',
              'INIC',
              'KASA',
              'LIAFSR',
              'LMRA',
              'MARR',
              'MEXR',
              'MEXZ',
              'mprF',
              'MPRF',
              'NDH',
              'omp36',
              'OMP36',
              'OMPF',
              'OPRD',
              'PARC',
              'parE',
              'PARE',
              'PGSA',
              'phoP',
              'PHOP',
              'PNCA',
              'POR',
              'PORB',
              'RAMR',
              'rpoB',
              'RPOB',
              'RPOC',
              'RPSL',
              'SOXS',
              'tetR',
              'TETR',
              'TLYA',
              'TUFAB')


validate_data_inputs <- function(amr_counts,
                                 amr_annotations,
                                 microbiome_data,
                                 metadata) {
    
    if(is.null(amr_counts) && is.null(amr_annotations) &&
       is.null(microbiome_data)) {
        cat('Error: Please select AMR or Microbiome files for upload')
        return()
    }
    
    # AMR validation logic
    if(!is.null(amr_counts) && !is.null(amr_annotations)) {
        
        if(!all(rownames(MRcounts(amr_counts)) %in% amr_annotations[['header']])) {
            cat('Error: AMR count matrix gene names do not match the gene',
                ' names in the annotation file. Make sure your annotation file',
                ' matches the version of the MEGARes database used in the',
                ' AMR++ nextflow pipeline.\n',
                sep='')
        }
        
        if(length(colnames(amr_annotations)) != 4) {
            cat('Error: AMR annotation file has an improper number of columns\n',
                'Check file integrity and proper CSV format.\n',
                sep='')
        }
        
        if(is.null(metadata)) {
            cat('Error: Metadata file not present.\n')
            return()
        }
        else {
            
            sample_column_id <- colnames(metadata)[1]
            
            if(!all(colnames(amr_counts) %in% metadata[[sample_column_id]])) {
                cat('Error: AMR count matrix sample names do not match the',
                    ' sample names in the metadata file (first column).',
                    ' Check file integrity and proper CSV format.\n',
                    sep='')
            }
            cat('AMR files validated\n')
        }
    }
    else if(!is.null(amr_counts) && is.null(amr_annotations)) {
        cat('Error: AMR count matrix detected but AMR annotations missing.\n')
    }
    else if(is.null(amr_counts) && !is.null(amr_annotations)) {
        cat('Error: AMR annotations detected but AMR count matrix missing.\n')
    }
    else {
        cat('No AMR data files detected: ',
            'resistance analysis will not be performed.\n',
            sep='')
    }
    
    
    # Microbiome validation logic
    if(!is.null(microbiome_data)) {
        
        if(is.null(metadata)) {
            cat('Error: Metadata file not present.\n')
            return()
        }
        else {
            sample_column_id <- colnames(metadata)[1]
            
            if(!all(colnames(microbiome_data) %in% metadata[[sample_column_id]])) {
                cat('Error: Microbiome count matrix sample names do not match the',
                    ' sample names in the metadata file (first column).',
                    ' Check file integrity and proper CSV format.\n',
                    sep='')
            }
            cat('Microbiome file validated\n')
        }
    }
    else {
        cat('No microbiome data file detected: microbiome analysis',
            ' will not be performed.\n',
            sep='')
    }
}


create_working_directories <- function(graph_output_dir,
                                       stats_output_dir,
                                       statistical_analyses,
                                       exploratory_analyses,
                                       data_types) {
    # If subdirs for stats and exploratory variables don't exist, create them
    ifelse(!dir.exists(file.path(graph_output_dir)),
           dir.create(file.path(graph_output_dir), mode='777'), FALSE)
    ifelse(!dir.exists(file.path(stats_output_dir)),
           dir.create(file.path(stats_output_dir), mode='777'), FALSE)
    
    for( dtype in data_types ) {
        ifelse(!dir.exists(file.path(graph_output_dir, dtype)),
               dir.create(file.path(graph_output_dir, dtype),
                          mode='777'), FALSE)
        
        for( v in 1:length(exploratory_analyses) ) {
            ifelse(!dir.exists(file.path(graph_output_dir, dtype,
                                         exploratory_analyses[[v]]$name)),
                   dir.create(file.path(graph_output_dir, dtype,
                                        exploratory_analyses[[v]]$name),
                              mode='777'), FALSE)
        }
        
        ifelse(!dir.exists(file.path(stats_output_dir, dtype)),
               dir.create(file.path(stats_output_dir, dtype),
                          mode='777'), FALSE)
        
        for( a in 1:length(statistical_analyses) ) {
            ifelse(!dir.exists(file.path(stats_output_dir, dtype,
                                         statistical_analyses[[a]]$name)),
                   dir.create(file.path(stats_output_dir, dtype,
                                        statistical_analyses[[a]]$name),
                              mode='777'), FALSE)
        }
    }
    
    if('AMR' %in% data_types) {
        ifelse(!dir.exists(file.path('amr_matrices')),
               dir.create(file.path('amr_matrices'), mode='777'), FALSE)
        ifelse(!dir.exists(file.path('amr_matrices/sparse_normalized')),
               dir.create(file.path('amr_matrices/sparse_normalized'),
                          mode='777'), FALSE)
        ifelse(!dir.exists(file.path('amr_matrices/normalized')),
               dir.create(file.path('amr_matrices/normalized'),
                          mode='777'), FALSE)
        ifelse(!dir.exists(file.path('amr_matrices/raw')),
               dir.create(file.path('amr_matrices/raw'), mode='777'), FALSE)
    }
    
    if('Microbiome' %in% data_types) {
        ifelse(!dir.exists(file.path('kraken_matrices')),
               dir.create(file.path('kraken_matrices'), mode='777'), FALSE)
        ifelse(!dir.exists(file.path('kraken_matrices/sparse_normalized')),
               dir.create(file.path('kraken_matrices/sparse_normalized'),
                          mode='777'), FALSE)
        ifelse(!dir.exists(file.path('kraken_matrices/normalized')),
               dir.create(file.path('kraken_matrices/normalized'),
                          mode='777'), FALSE)
        ifelse(!dir.exists(file.path('kraken_matrices/raw')),
               dir.create(file.path('kraken_matrices/raw'), mode='777'), FALSE)
    }
}


CSS_normalize_and_extract <- function(filtering_quantile,
                                      kraken_MRexp=NA,
                                      amr_MRexp=NA) {
    amr_norm <- NA
    amr_raw <- NA
    kraken_norm <- NA
    kraken_raw <- NA
    if(!is.na(kraken_MRexp)) {
        valid_samples <- apply(MRcounts(kraken_MRexp), 2, function(x) {sum(x > 0)}) > 1
        kraken_MRexp <- kraken_MRexp[, valid_samples]
        cumNorm(kraken_MRexp)
        kraken_norm <- data.table(MRcounts(kraken_MRexp, norm=T))
        kraken_norm_filtering_threshold <- quantile(rowSums(kraken_norm), filtering_quantile)
        filter_indices <- rowSums(kraken_norm) >= kraken_norm_filtering_threshold
        kraken_norm <- kraken_norm[filter_indices]
        kraken_norm[, id := ( rownames(kraken_MRexp)[filter_indices] )]
        
        kraken_raw <- data.table(MRcounts(kraken_MRexp, norm=F))
        kraken_raw_filtering_threshold <- quantile(rowSums(kraken_raw), filtering_quantile)
        filter_indices <- rowSums(kraken_raw) >= kraken_raw_filtering_threshold
        kraken_raw <- kraken_raw[filter_indices]
        kraken_raw[, id := ( rownames(kraken_MRexp)[filter_indices] )]
    }
    if(!is.na(amr_MRexp)) {
        valid_samples <- apply(MRcounts(amr_MRexp), 2, function(x) {sum(x > 0)}) > 1
        amr_MRexp <- amr_MRexp[, valid_samples]
        cumNorm(amr_MRexp)
        amr_norm <- data.table(MRcounts(amr_MRexp, norm=T))
        amr_norm_filtering_threshold <- quantile(rowSums(amr_norm), filtering_quantile)
        filter_indices <- which(rowSums(amr_norm) >= amr_norm_filtering_threshold)
        amr_norm <- amr_norm[filter_indices]
        amr_norm[, header :=( rownames(amr_MRexp)[filter_indices] )]
        
        amr_raw <- data.table(MRcounts(amr_MRexp, norm=F))
        amr_raw_filtering_threshold <- quantile(rowSums(amr_raw), filtering_quantile)
        filter_indices <- rowSums(amr_raw) >= amr_raw_filtering_threshold
        amr_raw <- amr_raw[filter_indices]
        amr_raw[, header :=( rownames(amr_MRexp)[filter_indices] )]
    }
    return(list(kraken_MRexp, kraken_norm, kraken_raw, amr_MRexp, amr_norm, amr_raw))
}


rarefy_normalize_and_extract <- function(sampling_depth,
                                         filtering_quantile,
                                         kraken_MRexp=NA,
                                         amr_MRexp=NA) {
    amr_norm <- NA
    amr_raw <- NA
    kraken_norm <- NA
    kraken_raw <- NA
    if(!is.na(kraken_MRexp)) {
        kraken_sample_ranking <- colSums(MRcounts(kraken_MRexp))[order(colSums(MRcounts(kraken_MRexp)))]
        kraken_rare_threshold <- kraken_sample_ranking[sampling_depth]
        
        kraken_norm <- data.table(rrarefy(MRcounts(kraken_MRexp), sample=kraken_rare_threshold))
        kraken_norm_filtering_threshold <- quantile(rowSums(kraken_norm), filtering_quantile)
        filter_indices <- rowSums(kraken_norm) >= kraken_norm_filtering_threshold
        kraken_norm <- kraken_norm[filter_indices]
        kraken_norm[, id :=( rownames(kraken_MRexp)[filter_indices] )]
        
        kraken_raw <- data.table(MRcounts(kraken_MRexp, norm=F))
        kraken_raw_filtering_threshold <- quantile(rowSums(kraken_raw), filtering_quantile)
        filter_indices <- rowSums(kraken_raw) >= kraken_raw_filtering_threshold
        kraken_raw <- kraken_raw[filter_indices]
        kraken_raw[, id :=( rownames(kraken_MRexp)[filter_indices] )]
    }
    if(!is.na(amr_MRexp)) {
        amr_sample_ranking <- colSums(MRcounts(amr_MRexp))[order(colSums(MRcounts(amr_MRexp)))]
        amr_rare_threshold <- amr_sample_ranking[sampling_depth]
        
        amr_norm <- data.table(rrarefy(MRcounts(amr_MRexp), sample=amr_rare_threshold))
        amr_norm_filtering_threshold <- quantile(rowSums(amr_norm), filtering_quantile)
        filter_indices <- rowSums(amr_norm) >= amr_norm_filtering_threshold
        amr_norm <- amr_norm[filter_indices]
        amr_norm[, header :=( rownames(amr_MRexp)[filter_indices] )]
        
        amr_raw <- data.table(MRcounts(amr_MRexp, norm=F))
        amr_raw_filtering_threshold <- quantile(rowSums(amr_raw), filtering_quantile)
        filter_indices <- rowSums(amr_raw) >= amr_raw_filtering_threshold
        amr_raw <- amr_raw[filter_indices]
        amr_raw[, header :=( rownames(amr_MRexp)[filter_indices] )]
    }
    return(list(kraken_MRexp, kraken_norm, kraken_raw, amr_MRexp, amr_norm, amr_raw))
}


aggregate_and_filter <- function(dt_list,
                                 metadata,
                                 filter=TRUE) {
    kraken <- dt_list[[1]]
    kraken_norm <- dt_list[[2]]
    kraken_raw <- dt_list[[3]]
    amr <- dt_list[[4]]
    amr_norm <- dt_list[[5]]
    amr_raw <- dt_list[[6]]
    annotations <- dt_list[[7]]
    
    sample_column_id <- colnames(metadata)[1]
    
    # Outputs
    AMR_analytic_data <- NA
    AMR_raw_analytic_data <- NA
    amr_melted_analytic <- NA
    amr_melted_raw_analytic <- NA
    kraken_analytic_data <- NA
    kraken_raw_analytic_data <- NA
    kraken_melted_analytic <- NA
    kraken_melted_raw_analytic <- NA
    AMR_analytic_names <- NA
    kraken_analytic_names <- NA
    
    if(!is.na(amr_norm)) {
        setkey(amr_norm, header)
        amr_norm <- annotations[amr_norm]
        
        setkey(amr_raw, header)
        amr_raw <- annotations[amr_raw]
        
        if(filter == TRUE) {
            # Remove groups that correspond to potentially wild-type genes
            amr_raw <- amr_raw[!(group %in% snp_regex), ]
            amr_norm<- amr_norm[!(group %in% snp_regex), ]
        }
        
        # Group the AMR data by level for analysis
        amr_class <- amr_norm[, lapply(.SD, sum), by='class', .SDcols=!c('header', 'mechanism', 'group')]
        amr_class_analytic <- newMRexperiment(counts=amr_class[, .SD, .SDcols=!'class'])
        rownames(amr_class_analytic) <- amr_class$class
        
        amr_class_raw <- amr_raw[, lapply(.SD, sum), by='class', .SDcols=!c('header', 'mechanism', 'group')]
        amr_class_raw_analytic <- newMRexperiment(counts=amr_class_raw[, .SD, .SDcols=!'class'])
        rownames(amr_class_raw_analytic) <- amr_class_raw$class
        
        amr_mech <- amr_norm[, lapply(.SD, sum), by='mechanism', .SDcols=!c('header', 'class', 'group')]
        amr_mech_analytic <- newMRexperiment(counts=amr_mech[, .SD, .SDcols=!'mechanism'])
        rownames(amr_mech_analytic) <- amr_mech$mechanism
        
        amr_mech_raw <- amr_raw[, lapply(.SD, sum), by='mechanism', .SDcols=!c('header', 'class', 'group')]
        amr_mech_raw_analytic <- newMRexperiment(counts=amr_mech_raw[, .SD, .SDcols=!'mechanism'])
        rownames(amr_mech_raw_analytic) <- amr_mech_raw$mechanism
        
        amr_group <- amr_norm[, lapply(.SD, sum), by='group', .SDcols=!c('header', 'mechanism', 'class')]
        amr_group_analytic <- newMRexperiment(counts=amr_group[, .SD, .SDcols=!'group'])
        rownames(amr_group_analytic) <- amr_group$group
        
        amr_group_raw <- amr_raw[, lapply(.SD, sum), by='group', .SDcols=!c('header', 'mechanism', 'class')]
        amr_group_raw_analytic <- newMRexperiment(counts=amr_group_raw[, .SD, .SDcols=!'group'])
        rownames(amr_group_raw_analytic) <- amr_group_raw$group
        
        amr_gene_analytic <- newMRexperiment(
            counts=amr_norm[!(group %in% snp_regex),
                            .SD, .SDcols=!c('header', 'class', 'mechanism', 'group')])
        amr_gene_raw_analytic <- newMRexperiment(
            counts=amr_raw[!(group %in% snp_regex),
                           .SD, .SDcols=!c('header', 'class', 'mechanism', 'group')])
        
        rownames(amr_gene_analytic) <- amr_norm$header
        rownames(amr_gene_raw_analytic) <- amr_raw$header
        
        
        # Make long data frame for plotting with ggplot2
        amr_melted_analytic <- rbind(melt_dt(MRcounts(amr_class_analytic), 'Class'),
                                     melt_dt(MRcounts(amr_mech_analytic), 'Mechanism'),
                                     melt_dt(MRcounts(amr_group_analytic), 'Group'),
                                     melt_dt(MRcounts(amr_gene_analytic), 'Gene'))
        amr_melted_raw_analytic <- rbind(melt_dt(MRcounts(amr_class_raw_analytic), 'Class'),
                                         melt_dt(MRcounts(amr_mech_raw_analytic), 'Mechanism'),
                                         melt_dt(MRcounts(amr_group_raw_analytic), 'Group'),
                                         melt_dt(MRcounts(amr_gene_raw_analytic), 'Gene'))
        
        # Vector of objects for iteration and their names
        AMR_analytic_data <- c(amr_class_analytic,
                               amr_mech_analytic,
                               amr_group_analytic,
                               amr_gene_analytic)
        AMR_analytic_names <- c('Class', 'Mechanism', 'Group', 'Gene')
        AMR_raw_analytic_data <- c(amr_class_raw_analytic,
                                   amr_mech_raw_analytic,
                                   amr_group_raw_analytic,
                                   amr_gene_raw_analytic)
        
        for( l in 1:length(AMR_analytic_data) ) {
            sample_idx <- match(colnames(MRcounts(AMR_analytic_data[[l]])), metadata[[sample_column_id]])
            pData(AMR_analytic_data[[l]]) <- data.frame(
                metadata[sample_idx, .SD, .SDcols=!sample_column_id])
            rownames(pData(AMR_analytic_data[[l]])) <- metadata[sample_idx, .SD, .SDcols=sample_column_id][[sample_column_id]]
            fData(AMR_analytic_data[[l]]) <- data.frame(Feature=rownames(MRcounts(AMR_analytic_data[[l]])))
            rownames(fData(AMR_analytic_data[[l]])) <- rownames(MRcounts(AMR_analytic_data[[l]]))
        }
        
        for( l in 1:length(AMR_raw_analytic_data) ) {
            sample_idx <- match(colnames(MRcounts(AMR_raw_analytic_data[[l]])), metadata[[sample_column_id]])
            pData(AMR_raw_analytic_data[[l]]) <- data.frame(
                metadata[sample_idx, .SD, .SDcols=!sample_column_id])
            rownames(pData(AMR_raw_analytic_data[[l]])) <- metadata[sample_idx, .SD, .SDcols=sample_column_id][[sample_column_id]]
            fData(AMR_raw_analytic_data[[l]]) <- data.frame(Feature=rownames(MRcounts(AMR_raw_analytic_data[[l]])))
            rownames(fData(AMR_raw_analytic_data[[l]])) <- rownames(MRcounts(AMR_raw_analytic_data[[l]]))
        }
        
        metadata <- data.table(metadata[match(colnames(MRcounts(amr_class_analytic)), metadata[[sample_column_id]])])
        setkeyv(metadata, sample_column_id)
        
        # write.csv(make_sparse(amr_class, 'class', c('class')), 'amr_matrices/sparse_normalized/AMR_Class_Sparse_Normalized.csv',
        #           row.names=T)
        # write.table(amr_class, 'amr_matrices/normalized/AMR_Class_Normalized.csv', sep=',', row.names = F, col.names = T)
        # write.table(amr_class_raw, 'amr_matrices/raw/AMR_Class_Raw.csv', sep=',', row.names = F, col.names = T)
        # 
        # 
        # write.csv(make_sparse(amr_mech, 'mechanism', c('mechanism')), 'amr_matrices/sparse_normalized/AMR_Mechanism_Sparse_Normalized.csv',
        #           row.names=T)
        # write.table(amr_mech, 'amr_matrices/normalized/AMR_Mechanism_Normalized.csv', sep=',', row.names = F, col.names = T)
        # write.table(amr_mech_raw, 'amr_matrices/raw/AMR_Mechanism_Raw.csv', sep=',', row.names = F, col.names = T)
        # 
        # write.csv(make_sparse(amr_group, 'group', c('group')), 'amr_matrices/sparse_normalized/AMR_Group_Sparse_Normalized.csv',
        #           row.names=T)
        # write.table(amr_group, 'amr_matrices/normalized/AMR_Group_Normalized.csv', sep=',', row.names = F, col.names = T)
        # write.table(amr_mech_raw, 'amr_matrices/raw/AMR_Group_Raw.csv', sep=',', row.names = F, col.names = T)
        # 
        # write.csv(make_sparse(amr_norm, 'header', c('header', 'class', 'mechanism', 'group')),
        #           'amr_matrices/sparse_normalized/AMR_Gene_Sparse_Normalized.csv',
        #           row.names=T)
        # write.table(amr_norm, 'amr_matrices/normalized/AMR_Gene_Normalized.csv', sep=',', row.names = F, col.names = T)
        # write.table(amr_raw, 'amr_matrices/raw/AMR_Gene_Raw.csv', sep=',', row.names = F, col.names = T)
    }
    
    if(!is.na(kraken_norm)) {
        kraken_taxonomy <- data.table(id=rownames(kraken))
        num_sep = max(sapply(kraken_taxonomy$id, function(x) {sum(unlist(strsplit(x, '')) == '|')}))
        
        if(num_sep <= 6) {
            setDT(kraken_taxonomy)[, c('Domain',
                                       'Phylum',
                                       'Class',
                                       'Order',
                                       'Family',
                                       'Genus',
                                       'Species') := tstrsplit(id, '|', type.convert = TRUE, fixed = TRUE)]
        } else if(num_sep == 7) {
            setDT(kraken_taxonomy)[, c('Domain',
                                       'Kingdom',
                                       'Phylum',
                                       'Class',
                                       'Order',
                                       'Family',
                                       'Genus',
                                       'Species') := tstrsplit(id, '|', type.convert = TRUE, fixed = TRUE)]
        }
        
        kraken_taxonomy[, Kingdom:=NULL]
        
        setkey(kraken_taxonomy, id)
        setkey(kraken_norm, id)
        kraken_norm <- kraken_taxonomy[kraken_norm]  # left outer join
        
        setkey(kraken_raw, id)
        kraken_raw <- kraken_taxonomy[kraken_raw]  # left outer join
        
        
        # Group the kraken data by level for analysis, removing NA entries
        kraken_domain <- kraken_norm[!is.na(Domain) & Domain != 'NA', lapply(.SD, sum), by='Domain', .SDcols=!1:8]
        kraken_domain_analytic <- newMRexperiment(counts=kraken_domain[, .SD, .SDcols=!'Domain'])
        rownames(kraken_domain_analytic) <- kraken_domain$Domain
        
        kraken_domain_raw <- kraken_raw[!is.na(Domain) & Domain != 'NA', lapply(.SD, sum), by='Domain', .SDcols=!1:8]
        kraken_domain_raw_analytic <- newMRexperiment(counts=kraken_domain_raw[, .SD, .SDcols=!'Domain'])
        rownames(kraken_domain_raw_analytic) <- kraken_domain_raw$Domain
        
        kraken_phylum <- kraken_norm[!is.na(Phylum) & Phylum != 'NA', lapply(.SD, sum), by='Phylum', .SDcols=!1:8]
        kraken_phylum_analytic <- newMRexperiment(counts=kraken_phylum[, .SD, .SDcols=!'Phylum'])
        rownames(kraken_phylum_analytic) <- kraken_phylum$Phylum
        
        kraken_phylum_raw <- kraken_raw[!is.na(Phylum) & Phylum != 'NA', lapply(.SD, sum), by='Phylum', .SDcols=!1:8]
        kraken_phylum_raw_analytic <- newMRexperiment(counts=kraken_phylum_raw[, .SD, .SDcols=!'Phylum'])
        rownames(kraken_phylum_raw_analytic) <- kraken_phylum_raw$Phylum
        
        kraken_class <- kraken_norm[!is.na(Class) & Class != 'NA', lapply(.SD, sum), by='Class', .SDcols=!1:8]
        kraken_class_analytic <- newMRexperiment(counts=kraken_class[, .SD, .SDcols=!'Class'])
        rownames(kraken_class_analytic) <- kraken_class$Class
        
        kraken_class_raw <- kraken_raw[!is.na(Class) & Class != 'NA', lapply(.SD, sum), by='Class', .SDcols=!1:8]
        kraken_class_raw_analytic <- newMRexperiment(counts=kraken_class_raw[, .SD, .SDcols=!'Class'])
        rownames(kraken_class_raw_analytic) <- kraken_class_raw$Class
        
        kraken_order <- kraken_norm[!is.na(Order) & Order != 'NA', lapply(.SD, sum), by='Order', .SDcols=!1:8]
        kraken_order_analytic <- newMRexperiment(counts=kraken_order[, .SD, .SDcols=!'Order'])
        rownames(kraken_order_analytic) <- kraken_order$Order
        
        kraken_order_raw <- kraken_raw[!is.na(Order) & Order != 'NA', lapply(.SD, sum), by='Order', .SDcols=!1:8]
        kraken_order_raw_analytic <- newMRexperiment(counts=kraken_order_raw[, .SD, .SDcols=!'Order'])
        rownames(kraken_order_raw_analytic) <- kraken_order_raw$Order
        
        kraken_family <- kraken_norm[!is.na(Family) & Family != 'NA', lapply(.SD, sum), by='Family', .SDcols=!1:8]
        kraken_family_analytic <- newMRexperiment(counts=kraken_family[, .SD, .SDcols=!'Family'])
        rownames(kraken_family_analytic) <- kraken_family$Family
        
        kraken_family_raw <- kraken_raw[!is.na(Family) & Family != 'NA', lapply(.SD, sum), by='Family', .SDcols=!1:8]
        kraken_family_raw_analytic <- newMRexperiment(counts=kraken_family_raw[, .SD, .SDcols=!'Family'])
        rownames(kraken_family_raw_analytic) <- kraken_family_raw$Family
        
        kraken_genus <- kraken_norm[!is.na(Genus) & Genus != 'NA', lapply(.SD, sum), by='Genus', .SDcols=!1:8]
        kraken_genus_analytic <- newMRexperiment(counts=kraken_genus[, .SD, .SDcols=!'Genus'])
        rownames(kraken_genus_analytic) <- kraken_genus$Genus
        
        kraken_genus_raw <- kraken_raw[!is.na(Genus) & Genus != 'NA', lapply(.SD, sum), by='Genus', .SDcols=!1:8]
        kraken_genus_raw_analytic <- newMRexperiment(counts=kraken_genus_raw[, .SD, .SDcols=!'Genus'])
        rownames(kraken_genus_raw_analytic) <- kraken_genus_raw$Genus
        
        kraken_species <- kraken_norm[!is.na(Species) & Species != 'NA', lapply(.SD, sum), by='Species', .SDcols=!1:8]
        kraken_species_analytic <- newMRexperiment(counts=kraken_species[, .SD, .SDcols=!'Species'])
        rownames(kraken_species_analytic) <- kraken_species$Species
        
        kraken_species_raw <- kraken_raw[!is.na(Species) & Species != 'NA', lapply(.SD, sum), by='Species', .SDcols=!1:8]
        kraken_species_raw_analytic <- newMRexperiment(counts=kraken_species_raw[, .SD, .SDcols=!'Species'])
        rownames(kraken_species_raw_analytic) <- kraken_species_raw$Species
        
        
        # Make long data frame for plotting with ggplot2
        kraken_melted_analytic <- rbind(melt_dt(MRcounts(kraken_domain_analytic), 'Domain'),
                                        melt_dt(MRcounts(kraken_phylum_analytic), 'Phylum'),
                                        melt_dt(MRcounts(kraken_class_analytic), 'Class'),
                                        melt_dt(MRcounts(kraken_order_analytic), 'Order'),
                                        melt_dt(MRcounts(kraken_family_analytic), 'Family'),
                                        melt_dt(MRcounts(kraken_genus_analytic), 'Genus'),
                                        melt_dt(MRcounts(kraken_species_analytic), 'Species'))
        kraken_melted_raw_analytic <- rbind(melt_dt(MRcounts(kraken_domain_raw_analytic), 'Domain'),
                                            melt_dt(MRcounts(kraken_phylum_raw_analytic), 'Phylum'),
                                            melt_dt(MRcounts(kraken_class_raw_analytic), 'Class'),
                                            melt_dt(MRcounts(kraken_order_raw_analytic), 'Order'),
                                            melt_dt(MRcounts(kraken_family_raw_analytic), 'Family'),
                                            melt_dt(MRcounts(kraken_genus_raw_analytic), 'Genus'),
                                            melt_dt(MRcounts(kraken_species_raw_analytic), 'Species'))
        
        # Vector of objects for iteration and their names
        kraken_analytic_data <- c(kraken_domain_analytic,
                                  kraken_phylum_analytic,
                                  kraken_class_analytic,
                                  kraken_order_analytic,
                                  kraken_family_analytic,
                                  kraken_genus_analytic,
                                  kraken_species_analytic)
        kraken_analytic_names <- c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
        kraken_raw_analytic_data <- c(kraken_domain_raw_analytic,
                                      kraken_phylum_raw_analytic,
                                      kraken_class_raw_analytic,
                                      kraken_order_raw_analytic,
                                      kraken_family_raw_analytic,
                                      kraken_genus_raw_analytic,
                                      kraken_species_raw_analytic)

        for( l in 1:length(kraken_analytic_data) ) {
            sample_idx <- match(colnames(MRcounts(kraken_analytic_data[[l]])), metadata[[sample_column_id]])
            pData(kraken_analytic_data[[l]]) <- data.frame(
                metadata[sample_idx, .SD, .SDcols=!sample_column_id])
            rownames(pData(kraken_analytic_data[[l]])) <- metadata[[sample_column_id]][sample_idx]
            fData(kraken_analytic_data[[l]]) <- data.frame(Feature=rownames(MRcounts(kraken_analytic_data[[l]])))
            rownames(fData(kraken_analytic_data[[l]])) <- rownames(MRcounts(kraken_analytic_data[[l]]))
        }
        
        for( l in 1:length(kraken_raw_analytic_data) ) {
            sample_idx <- match(colnames(MRcounts(kraken_raw_analytic_data[[l]])), metadata[[sample_column_id]])
            pData(kraken_raw_analytic_data[[l]]) <- data.frame(
                metadata[sample_idx, .SD, .SDcols=!sample_column_id])
            rownames(pData(kraken_raw_analytic_data[[l]])) <- metadata[sample_idx, .SD, .SDcols=sample_column_id][[sample_column_id]]
            fData(kraken_raw_analytic_data[[l]]) <- data.frame(Feature=rownames(MRcounts(kraken_raw_analytic_data[[l]])))
            rownames(fData(kraken_raw_analytic_data[[l]])) <- rownames(MRcounts(kraken_raw_analytic_data[[l]]))
        }
        
        # Ensure that the metadata entries match the factor order of the MRexperiments
        metadata <- data.table(metadata[match(colnames(MRcounts(kraken_domain_analytic)), metadata[[sample_column_id]])])
        setkeyv(metadata, sample_column_id)
        
        # write.csv(make_sparse(kraken_domain, 'Domain', c('Domain')),
        #           'kraken_matrices/sparse_normalized/kraken_Domain_Sparse_Normalized.csv',
        #           row.names=T)
        # write.table(kraken_domain, 'kraken_matrices/normalized/kraken_Domain_Normalized.csv', sep=',', row.names=F, col.names=T)
        # write.table(kraken_domain_raw, 'kraken_matrices/raw/kraken_Domain_Raw.csv', sep=',', row.names=F, col.names=T)
        # 
        # write.csv(make_sparse(kraken_phylum, 'Phylum', c('Phylum')),
        #           'kraken_matrices/sparse_normalized/kraken_Phylum_Sparse_Normalized.csv',
        #           row.names=T)
        # write.table(kraken_phylum, 'kraken_matrices/normalized/kraken_Phylum_Normalized.csv', sep=',', row.names=F, col.names=T)
        # write.table(kraken_phylum_raw, 'kraken_matrices/raw/kraken_Phylum_Raw.csv', sep=',', row.names=F, col.names=T)
        # 
        # write.csv(make_sparse(kraken_class, 'Class', c('Class')),
        #           'kraken_matrices/sparse_normalized/kraken_Class_Sparse_Normalized.csv',
        #           row.names=T)
        # write.table(kraken_class, 'kraken_matrices/normalized/kraken_Class_Normalized.csv', sep=',', row.names=F, col.names=T)
        # write.table(kraken_class_raw, 'kraken_matrices/raw/kraken_Class_Raw.csv', sep=',', row.names=F, col.names=T)
        # 
        # write.csv(make_sparse(kraken_order, 'Order', c('Order')),
        #           'kraken_matrices/sparse_normalized/kraken_Order_Sparse_Normalized.csv',
        #           row.names=T)
        # write.table(kraken_order, 'kraken_matrices/normalized/kraken_Order_Normalized.csv', sep=',', row.names=F, col.names=T)
        # write.table(kraken_order_raw, 'kraken_matrices/raw/kraken_Order_Raw.csv', sep=',', row.names=F, col.names=T)
        # 
        # write.csv(make_sparse(kraken_family, 'Family', c('Family')),
        #           'kraken_matrices/sparse_normalized/kraken_Family_Sparse_Normalized.csv',
        #           row.names=T)
        # write.table(kraken_family, 'kraken_matrices/normalized/kraken_Family_Normalized.csv', sep=',', row.names=F, col.names=T)
        # write.table(kraken_family_raw, 'kraken_matrices/raw/kraken_Family_Raw.csv', sep=',', row.names=F, col.names=T)
        # 
        # write.csv(make_sparse(kraken_genus, 'Genus', c('Genus')),
        #           'kraken_matrices/sparse_normalized/kraken_Genus_Sparse_Normalized.csv',
        #           row.names=T)
        # write.table(kraken_genus, 'kraken_matrices/normalized/kraken_Genus_Normalized.csv', sep=',', row.names=F, col.names=T)
        # write.table(kraken_genus_raw, 'kraken_matrices/raw/kraken_Genus_Raw.csv', sep=',', row.names=F, col.names=T)
        # 
        # write.csv(make_sparse(kraken_species, 'Species', c('Species')),
        #           'kraken_matrices/sparse_normalized/kraken_Species_Sparse_Normalized.csv',
        #           row.names=T)
        # write.table(kraken_species, 'kraken_matrices/normalized/kraken_Species_Normalized.csv', sep=',', row.names=F, col.names=T)
        # write.table(kraken_species_raw, 'kraken_matrices/raw/kraken_Species_Raw.csv', sep=',', row.names=F, col.names=T)
    }
    
    
    
    return(list(AMR_analytic_data,
                AMR_raw_analytic_data,
                amr_melted_analytic,
                amr_melted_raw_analytic,
                AMR_analytic_names,
                kraken_analytic_data,
                kraken_raw_analytic_data,
                kraken_melted_analytic,
                kraken_melted_raw_analytic,
                kraken_analytic_names,
                metadata))
}
