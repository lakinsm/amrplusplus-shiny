server <- function(input, output) {
    set.seed(154)
    
    # Load AMR count data
    amr_counts <- reactive({
        counts <- input$amr_count_matrix
        if(is.null(counts)) {
            return(NULL)
        }
        load_amr_counts(counts$datapath)
    })
    
    # Load AMR annotations
    amr_annotations <- reactive({
        annotations <- input$amr_annotation_file
        if(is.null(annotations)) {
            return(NULL)
        }
        load_amr_annotations(annotations$datapath)
    })
    
    # Load microbiome count data
    microbiome_data <- reactive({
        counts <- input$kraken_count_matrix
        if(is.null(counts)) {
            return(NULL)
        }
        load_kraken_data(counts$datapath)
    })
    
    # Load experiment metadata
    metadata <- reactive({
        infile <- input$metadata_file
        if(is.null(infile)) {
            return(NULL)
        }
        load_metadata(infile$datapath)
    })
    
    observeEvent(input$validate_inputs, {
        output$validation_output <- renderPrint({
            isolate({
                invisible(
                    validate_data_inputs(amr_counts(),
                                     amr_annotations(),
                                     microbiome_data(),
                                     metadata())
                )
            })
        })
    })
}
