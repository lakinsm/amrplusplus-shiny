require(data.table)
require(metagenomeSeq)

load_amr_counts <- function(amr_filepath) {
    # Load the data, MEGARes annotations, and metadata
    out <- tryCatch(
        {
            newMRexperiment(read.table(amr_filepath, header=T, row.names=1, sep=','))
        },
        error=function(e) {
            print('Error in loading of the AMR count file\nCheck file integrity')
            print(e)
            return(NULL)
        },
        warning=function(w) {
            print('Warning encountered in loading of the AMR count file')
            print(w)
            return(NULL)
        }
    )
    return(out)
}


load_amr_annotations <- function(annotations_filepath) {
    out <- tryCatch(
        {
            annotations <- data.table(read.csv(annotations_filepath, header=T))
            setkey(annotations, header)
            return(annotations)
        },
        error=function(e) {
            print('Error in loading of the AMR annotation file\nCheck file integrity')
            print(e)
            return(NULL)
        },
        warning=function(w) {
            print('Warning encountered in loading of the AMR annotation file')
            print(w)
            return(NULL)
        }
    )
    return(out)
}


load_kraken_data <- function(kraken_filepath) {
    out <- tryCatch(
        {
            temp_kraken <- read.table(kraken_filepath, header=T, row.names=1, sep=',')
            kraken <- newMRexperiment(temp_kraken[rowSums(temp_kraken) > 0, ])
        },
        error=function(e) {
            print('Error in loading of the Kraken count file\nCheck file integrity')
            print(e)
            return(NULL)
        },
        warning=function(w) {
            print('Warning encountered in loading of the Kraken count file')
            print(w)
            return(NULL)
        }
    )
    return(out)
}


load_metadata <- function(metadata_filepath) {
    out <- tryCatch(
        {
            metadata <- read.csv(metadata_filepath, header=T)
            sample_column_id <- colnames(metadata)[1]
            metadata[, sample_column_id] <- make.names(metadata[, sample_column_id])
            return(metadata)
        },
        error=function(e) {
            print('Error in loading of the metadata file\nCheck file integrity')
            print(e)
            return(NULL)
        },
        warning=function(w) {
            print('Warning encountered in loading of the metadata file')
            print(w)
            return(NULL)
        }
    )
    return(out)
}
