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
    
    # Validate data inputs
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
    
    # Load experiment metadata
    metadata <- reactive({
        infile <- input$metadata_file
        if(is.null(infile)) {
            return(NULL)
        }
        load_metadata(infile$datapath)
    })
    
    # Create checkboxes for metadata column features
    observe({
        if(!is.null(metadata())) {
            choices <- c(colnames(metadata())[2:ncol(metadata())])
            updateCheckboxGroupInput(session = session,
                                     inputId = 'metadata_fields',
                                     choices = choices)
        }
    })
    
    # Populate selectInputs for checked features
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
            else {
                active_metadata_fields$data <- c()
            }
        })
    )
    
    # Apply the selected data types to the metadata table for analysis
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
    
    # Populate or remove an additional sliderInput if rarefaction is active/inactive
    observe({
        if(input$normalization_select == 'Cumulative Sum Scaling') {
            output$exploratory_rarefaction_sample_threshold <- renderUI({
                character(0)
            })
        }
        else if(input$normalization_select == 'Rarefaction') {
            min_samples <- c()
            isolate({
                if(!is.null(metadata()))  {
                    min_samples <- c(min_samples, nrow(metadata()))
                }
                if(!is.null(amr_counts())) {
                    min_samples <- c(min_samples, ncol(amr_counts()))
                }
                if(!is.null(microbiome_data())) {
                    min_samples <- c(min_samples, ncol(microbiome_data()))
                }
            })
            if(length(min_samples) > 0) {
                output$exploratory_rarefaction_sample_threshold <- renderUI({
                    sliderInput(inputId = 'rarefaction_slider',
                                label = 'Rarefy to Ranked Sample # (lowest to highest):',
                                min = 1,
                                max = min(min_samples),
                                step = 1,
                                value = 1)
                })
            }
        }
    })
    
    # Populate the options for Data Type exploratory preview selectInput
    observe({
        choices <- c()
        if(!is.null(amr_counts()) && !is.null(amr_annotations())) {
            choices <- c(choices, 'Resistome')
        }
        if(!is.null(microbiome_data())) {
            choices <- c(choices, 'Microbiome')
        }
        if(length(choices) > 0) {
            updateSelectInput(session = session,
                              inputId = 'exploratory_data_type_preview_select',
                              choices = choices)
        }
    })
    
    # Populate the options for exploratory preview annotation level dependent on data type
    observe({
        if(input$exploratory_data_type_preview_select == 'Resistome') {
            updateSelectInput(session = session,
                              inputId = 'exploratory_annotation_level_preview_select',
                              choices = c('Class',
                                          'Mechanism',
                                          'Group',
                                          'Gene'))
        }
        else if(input$exploratory_data_type_preview_select == 'Microbiome') {
            updateSelectInput(session = session,
                              inputId = 'exploratory_annotation_level_preview_select',
                              choices = c('Domain',
                                          'Phylum',
                                          'Class',
                                          'Order',
                                          'Family',
                                          'Genus',
                                          'Species'))
        }
    })
    
    # Populate the metadata features for exploratory preview dependent on metadata loading
    observe({
        if(length(active_metadata_fields$data) > 0) {
            updateSelectInput(session = session,
                              inputId = 'exploratory_feature_select',
                              choices = active_metadata_fields$data)
        }
        else {
            updateSelectInput(session = session,
                              inputId = 'exploratory_feature_select',
                              choices = character(0))
        }
    })
    
    # Reactive variables for exploratory subsetting
    number_exploratory_subsets <- reactiveVal(value = 0)
    exploratory_subset_obj_names <- reactiveValues(data=c())
    
    # Add an exploratory subset rule for a given feature type
    observeEvent(input$exploratory_add_subset, {
        if(length(input$exploratory_feature_select) > 0) {
            newval <- number_exploratory_subsets() + 1
            number_exploratory_subsets(newval)
            old_names <- exploratory_subset_obj_names$data
            exploratory_subset_obj_names$data <- c(old_names,
                                                   paste('value_', number_exploratory_subsets(), sep='', collapse=''),
                                                   paste('logic_', number_exploratory_subsets(), sep='', collapse=''))
            dtype <- input[[input$exploratory_feature_select]]
            if(dtype == 'Categorical') {
                output$exploratory_subset_rules <- renderUI({
                    lapply(1:number_exploratory_subsets(), function(x) {
                        splitLayout(cellWidths = c('25%', '75%'),
                                    selectInput(inputId = paste('logic_', x, sep='', collapse=''),
                                                label = NULL,
                                                choices = c('==', '!=')),
                                    selectInput(inputId = paste('value_', x, sep='', collapse=''),
                                                label = NULL,
                                                choices = as.character(unique(metadata_updated$data[, input$exploratory_feature_select])))
                        )
                    })
                })
            }
            else if(dtype == 'Numeric' || dtype == 'Ordinal') {
                output$exploratory_subset_rules <- renderUI({
                    lapply(1:number_exploratory_subsets(), function(x) {
                        splitLayout(cellWidths = c('25%', '75%'),
                                    selectInput(inputId = paste('logic_', x, sep='', collapse=''),
                                                label = NULL,
                                                choices = c('==', '!=', '>', '<', '>=', '<=')),
                                    textInput(inputId = paste('value_', x, sep='', collapse=''),
                                              label = NULL)
                        )
                    })
                })
            }
        }
    })
    
    # Remove an exploratory subset rule
    observeEvent(input$exploratory_remove_subset, {
        if(length(input$exploratory_feature_select) > 0) {
            if(number_exploratory_subsets() > 1) {
                newval <- number_exploratory_subsets() - 1
                number_exploratory_subsets(newval)
                old_names <- exploratory_subset_obj_names$data
                exploratory_subset_obj_names$data <- old_names[c(-1, -2)]
                dtype <- input[[input$exploratory_feature_select]]
                if(dtype == 'Categorical') {
                    output$exploratory_subset_rules <- renderUI({
                        lapply(1:number_exploratory_subsets(), function(x) {
                            splitLayout(cellWidths = c('25%', '75%'),
                                        selectInput(inputId = paste('logic_', x, sep='', collapse=''),
                                                    label = NULL,
                                                    choices = c('==', '!=')),
                                        selectInput(inputId = paste('value_', x, sep='', collapse=''),
                                                    label = NULL,
                                                    choices = as.character(unique(metadata_updated$data[, input$exploratory_feature_select])))
                            )
                        })
                    })
                }
                else if(dtype == 'Numeric' || dtype == 'Ordinal') {
                    output$exploratory_subset_rules <- renderUI({
                        lapply(1:number_exploratory_subsets(), function(x) {
                            splitLayout(cellWidths = c('25%', '75%'),
                                        selectInput(inputId = paste('logic_', x, sep='', collapse=''),
                                                    label = NULL,
                                                    choices = c('==', '!=', '>', '<', '>=', '<=')),
                                        textInput(inputId = paste('value_', x, sep='', collapse=''),
                                                  label = NULL)
                            )
                        })
                    })
                }
            }
            else if(number_exploratory_subsets() == 1){
                newval <- number_exploratory_subsets() - 1
                number_exploratory_subsets(newval)
                exploratory_subset_obj_names$data <- c()
                output$exploratory_subset_rules <- renderUI({
                    character(0)
                })
            }
        }
        
    })
    
    # Make sure to clear all exploratory subset rules on feature change
    observe({
        x <- input$exploratory_feature_select
        output$exploratory_subset_rules <- renderUI({
            character(0)
        })
        number_exploratory_subsets(0)
        exploratory_subset_obj_names$data <- c()
    })
    
    # Run preview code for exploratory analysis on button click
    observeEvent(input$exploratory_update_preview, {
        isolate({
            subset_string <- c()
            if(length(exploratory_subset_obj_names$data) > 0) {
                for(i in 1:(length(exploratory_subset_obj_names$data) / 2)) {
                    logic_name <- exploratory_subset_obj_names$data[(2*(i-1))+2]
                    value_name <- exploratory_subset_obj_names$data[(2*(i-1))+1]
                    logic_str <- input[[logic_name]]
                    value_str <- input[[value_name]]
                    combined_str <- paste(input$exploratory_feature_select, logic_str, value_str, sep=' ', collapse='')
                    subset_string <- c(subset_string, combined_str)
                }
            }
            
            
            sample_depth <- 0
            if(input$normalization_select == 'Rarefaction') {
                sample_depth <- as.numeric(input$rarefaction_slider)
            }
            
            filtering_value <- as.numeric(input$filtering_slider) / 100
            
            if(input$exploratory_data_type_preview_select == 'Resistome') {
                g <- generate_exploratory_preview(data=list(amr_counts(), amr_annotations()),
                                             data_type='Resistome',
                                             metadata=metadata_updated$data,
                                             annotation_level=input$exploratory_annotation_level_preview_select,
                                             analysis_type=input$exploratory_analysis_type_select,
                                             metadata_feature=input$exploratory_feature_select,
                                             low_pass_filter_threshold=filtering_value,
                                             norm_method=input$normalization_select,
                                             sample_depth=sample_depth,
                                             subset_string=subset_string)
            }
            else if(input$exploratory_data_type_preview_select == 'Microbiome') {
                g <- generate_exploratory_preview(data=list(microbiome_data()),
                                             data_type='Microbiome',
                                             metadata=metadata_updated$data,
                                             annotation_level=input$exploratory_annotation_level_preview_select,
                                             analysis_type=input$exploratory_analysis_type_select,
                                             metadata_feature=input$exploratory_feature_select,
                                             low_pass_filter_threshold=filtering_value,
                                             norm_method=input$normalization_select,
                                             sample_depth=sample_depth,
                                             subset_string=subset_string)
            }
            
            if('gg' %in% class(g)) {
                output$exploratory_preview <- renderPlot({
                    print(g)
                }, width=1024, height=720)
            }
        })
    })
}

