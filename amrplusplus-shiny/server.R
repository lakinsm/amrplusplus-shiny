server <- function(input, output, session) {
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
    
    observe({
        if(!is.null(metadata())) {
            choices <- c(colnames(metadata())[2:ncol(metadata())])
            updateCheckboxGroupInput(session = session,
                              inputId = 'metadata_fields',
                              choices = choices)
        }
    })
    
    active_metadata_fields <- reactiveValues(data=c())
    invisible(
        observe({
            vals <- input$metadata_fields
            names(vals) <- input$metadata_fields
            if(length(vals) > 0) {
                active_metadata_fields$data <- vals
                output$metadata_types <- renderUI({
                    lapply(names(active_metadata_fields$data), function(x) {
                        selectInput(inputId = as.character(x),
                                    label = paste(x, ' - Select Data Type', collapse=''),
                                    choices = c('Categorical', 'Numeric', 'Ordinal'))
                    })
                })
            }
        })
    )
    
    metadata_updated <- reactiveValues(data=data.frame())
    observeEvent(input$apply_metadata_types, {
        isolate({
            metadata_updated$data <- metadata()
            local_active_fields <- active_metadata_fields$data
            invisible({
                for(field in names(local_active_fields)) {
                    dtype <- input[[field]]
                    if(dtype == 'Numeric' || dtype == 'Ordinal') {
                        metadata_updated$data[, field] <- as.numeric(metadata_updated$data[, field])
                    }
                    else if(dtype == 'Categorical') {
                        metadata_updated$data[, field] <- as.factor(metadata_updated$data[, field])
                    }
                }
            })
        })
        invisible({
            output$metadata_console <- renderPrint({
                cat('Metadata Summary:\n')
                str(metadata_updated$data)
            })
        })
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
