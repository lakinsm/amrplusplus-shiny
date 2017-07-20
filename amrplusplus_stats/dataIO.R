require(data.table)
require(metagenomeSeq)

load_amr_data <- function(amr_filepath,
                          annotations_filepath,
                          metadata_filepath) {
    # Load the data, MEGARes annotations, and metadata
    amr <- newMRexperiment(read.table('AMR_analytic_matrix.csv', header=T, row.names=1, sep=','))
    annotations <- data.table(read.csv(megares_annotation_filename, header=T))
    setkey(annotations, header)
    return(list(amr, annotations))
}


load_kraken_data <- function(kraken_filepath) {
    temp_kraken <- read.table('kraken_analytic_matrix.csv', header=T, row.names=1, sep=',')
    kraken <- newMRexperiment(temp_kraken[rowSums(temp_kraken) > 0, ])
    return(kraken)
}


load_metadata <- function(metadata_filepath,
                          sample_column_id) {
    metadata <- read.csv(metadata_filepath, header=T)
    metadata[, sample_column_id] <- make.names(metadata[, sample_column_id])
    return(metadata)
}