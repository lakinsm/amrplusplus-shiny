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


CSS_normalize_and_extract <- function(kraken_MRexp=NA,
                                      amr_MRexp=NA) {
    amr_norm <- NA
    amr_raw <- NA
    kraken_norm <- NA
    kraken_raw <- NA
    if(!is.na(kraken_MRexp)) {
        cumNorm(kraken_MRexp)
        kraken_norm <- data.table(MRcounts(kraken, norm=T))
        kraken_raw <- data.table(MRcounts(kraken, norm=F))
    }
    if(!is.na(amr_MRexp)) {
        cumNorm(amr_MRexp)
        amr_norm <- data.table(MRcounts(amr, norm=T))
        amr_raw <- data.table(MRcounts(amr, norm=F))
    }
    return(list(kraken_norm, kraken_raw, amr_norm, amr_raw))
}


aggregate_and_filter <- function(dt_list,
                                 metadata,
                                 filter=TRUE) {
    kraken_norm <- dt_list[[1]]
    kraken_raw <- dt_list[[2]]
    amr_norm <- dt_list[[3]]
    amr_raw <- dt_list[[4]]
    
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
        amr_norm[, header :=( rownames(amr) ), ]
        setkey(amr_norm, header)
        amr_norm <- annotations[amr_norm]
        
        amr_raw[, header :=( rownames(amr) ), ]
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
    }
    
    if(!is.na(kraken_norm)) {
        kraken_taxonomy <- data.table(id=rownames(kraken))
        setDT(kraken_taxonomy)[, c('Domain',
                                   'Phylum',
                                   'Class',
                                   'Order',
                                   'Family',
                                   'Genus',
                                   'Species') := tstrsplit(id, '|', type.convert = TRUE, fixed = TRUE)]
        setkey(kraken_taxonomy, id)
        kraken_norm[, id :=(rownames(kraken)), ]
        setkey(kraken_norm, id)
        kraken_norm <- kraken_taxonomy[kraken_norm]  # left outer join
        
        kraken_raw[, id :=(rownames(kraken)), ]
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
            rownames(pData(kraken_analytic_data[[l]])) <- metadata[sample_idx, .SD, .SDcols=sample_column_id][[sample_column_id]]
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
    }
    
    # Ensure that the metadata entries match the factor order of the MRexperiments
    metadata <- data.table(metadata[match(colnames(MRcounts(amr_class_analytic)), metadata[, sample_column_id]), ])
    setkeyv(metadata, sample_column_id)
    
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










