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
            output$statistical_rarefaction_sample_threshold <- renderUI({
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
            # ||
            #     input$statistical_data_type_preview_select == 'Microbiome') {
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
                                                   paste('exploratory_value_', number_exploratory_subsets(), sep='', collapse=''),
                                                   paste('exploratory_logic_', number_exploratory_subsets(), sep='', collapse=''))
            dtype <- input[[input$exploratory_feature_select]]
            if(dtype == 'Categorical') {
                output$exploratory_subset_rules <- renderUI({
                    lapply(1:number_exploratory_subsets(), function(x) {
                        splitLayout(cellWidths = c('25%', '75%'),
                                    selectInput(inputId = paste('exploratory_logic_', x, sep='', collapse=''),
                                                label = NULL,
                                                choices = c('==', '!=')),
                                    selectInput(inputId = paste('exploratory_value_', x, sep='', collapse=''),
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
                                    selectInput(inputId = paste('exploratory_logic_', x, sep='', collapse=''),
                                                label = NULL,
                                                choices = c('==', '!=', '>', '<', '>=', '<=')),
                                    textInput(inputId = paste('exploratory_value_', x, sep='', collapse=''),
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
                                        selectInput(inputId = paste('exploratory_logic_', x, sep='', collapse=''),
                                                    label = NULL,
                                                    choices = c('==', '!=')),
                                        selectInput(inputId = paste('exploratory_value_', x, sep='', collapse=''),
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
                                        selectInput(inputId = paste('exploratory_logic_', x, sep='', collapse=''),
                                                    label = NULL,
                                                    choices = c('==', '!=', '>', '<', '>=', '<=')),
                                        textInput(inputId = paste('exploratory_value_', x, sep='', collapse=''),
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
    
    # Make sure to clear all subset rules on exploratory feature change
    observe({
        x <- input$exploratory_feature_select
        output$exploratory_subset_rules <- renderUI({
            character(0)
        })
        number_exploratory_subsets(0)
        exploratory_subset_obj_names$data <- c()
    })
    
    # Add statistical parameter menus if a statistical method is selected
    observe({
        x <- input$exploratory_analysis_type_select
        isolate({
            if(x %in% c('ZIG Regression', 'Elastic Net Regression')) {
                if(input$exploratory_analysis_type_select == 'ZIG Regression') {
                    output$exploratory_statistical_parameters <- renderUI({
                        tagList(
                            sliderInput(inputId = 'exploratory_statistical_pvalue_slider',
                                        label = 'Minimum P-value Threshold',
                                        min = 0.01,
                                        max = 0.2,
                                        step=0.01,
                                        value=0.1),
                            sliderInput(inputId = 'exploratory_statistical_number_slider',
                                        label = '# Significant Results to Display',
                                        min = 1,
                                        max = 100,
                                        step=1,
                                        value=20),
                            selectInput(inputId = 'exploratory_statistical_sortby_select',
                                        label = 'Sort Results by Field',
                                        choices = c('P-value', 'Effect Size', 'Abundance', 'T-statistic')),
                            checkboxGroupInput(inputId = 'exploratory_statistical_feature_checkboxes',
                                               label = 'Features to Include in Regression',
                                               choices = active_metadata_fields$data,
                                               selected = active_metadata_fields$data),
                            selectInput(inputId = 'exploratory_statistical_random_effect_select',
                                        label = 'Random Effect to Include in Regression',
                                        choices = c('None', active_metadata_fields$data[!(active_metadata_fields$data %in% 
                                                                                              input$exploratory_feature_select)]),
                                        selected='None')
                        )
                    })
                }
                else if(input$exploratory_analysis_type_select == 'Elastic Net Regression') {
                    output$exploratory_statistical_parameters <- renderUI({
                        tagList(
                            sliderInput(inputId = 'exploratory_statistical_pvalue_slider',
                                        label = 'Minimum P-value Threshold',
                                        min = 0.01,
                                        max = 0.2,
                                        step=0.01,
                                        value=0.1),
                            sliderInput(inputId = 'exploratory_statistical_number_slider',
                                        label = '# Significant Results to Display',
                                        min = 1,
                                        max = 100,
                                        step=1,
                                        value=20),
                            selectInput(inputId = 'exploratory_statistical_sortby_select',
                                        label = 'Sort Results by Field',
                                        choices = c('P-value', 'Effect Size', 'Abundance', 'T-statistic')),
                            checkboxGroupInput(inputId = 'exploratory_statistical_feature_checkboxes',
                                               label = 'Fixed Effects to Include in Regression',
                                               choices = active_metadata_fields$data,
                                               selected = active_metadata_fields$data),
                            selectInput(inputId = 'exploratory_statistical_random_effect_select',
                                        label = 'Random Effect to Include in Regression',
                                        choices = c('None', active_metadata_fields$data[!(active_metadata_fields$data %in% 
                                                                                              input$exploratory_feature_select) &&
                                                                                            !(active_metadata_fields$data %in% 
                                                                                                  input$exploratory_statistical_feature_checkboxes)]),
                                        selected='None')
                        )
                    })
                }
            }
            else {
                output$exploratory_statistical_parameters <- renderUI({
                    character(0)
                })
            }
        })
    })
    
    # Update random effect choices if feature select changes
    observe({
        if(input$exploratory_analysis_type_select %in% c('ZIG Regression', 'Elastic Net Regression')) {
            x <- input$exploratory_feature_select
            updateSelectInput(session = session,
                              inputId = 'exploratory_statistical_random_effect_select',
                              choices = c('None', active_metadata_fields$data[!(active_metadata_fields$data %in% 
                                                                                    input$exploratory_feature_select)]),
                              selected='None')
        }
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
            print('check1')
            
            
            sample_depth <- 0
            if(input$normalization_select == 'Rarefaction') {
                sample_depth <- as.numeric(input$rarefaction_slider)
            }
            
            print('check2')
            
            filtering_value <- as.numeric(input$filtering_slider) / 100
            
            if(input$exploratory_analysis_type_select %in% c('NMDS', 'PCA', 'Bar Graph', 'Heatmap')) {
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
                
                output$exploratory_preview <- renderUI({
                    plotOutput(outputId = 'exploratory_output',
                               width = '100%',
                               height='720px')
                })
                if('gg' %in% class(g)) {
                    output$exploratory_output <- renderPlot({
                        print(g)
                    }, width=1024, height=720)
                }
            }
            else if(input$exploratory_analysis_type_select %in% c('ZIG Regression', 'Elastic Net Regression')) {
                if(length(input$exploratory_statistical_feature_checkboxes) > 0 &&
                   input$exploratory_feature_select %in% input$exploratory_statistical_feature_checkboxes) {
                    covars <- input$exploratory_statistical_feature_checkboxes[!(input$exploratory_statistical_feature_checkboxes %in% 
                                                                                   input$exploratory_statistical_random_effect_select)]
                    primary_effect_idx <- which(covars == input$exploratory_feature_select)
                    print(primary_effect_idx)
                    other_idx <- which(1:length(covars) != primary_effect_idx)
                    covars <- covars[c(primary_effect_idx, other_idx)]
                    covars <- c('~ 0', covars)
                    covar_str <- paste(covars, sep='', collapse=' + ')

                    rand_effect <- NA
                    if(input$exploratory_statistical_random_effect_select != 'None') {
                        rand_effect <- input$exploratory_statistical_random_effect_select
                    }
                    
                    if(input$exploratory_data_type_preview_select == 'Resistome') {
                        dt <- generate_statistical_preview(data=list(amr_counts(), amr_annotations()),
                                                           data_type='Resistome',
                                                           metadata=metadata_updated$data,
                                                           annotation_level=input$exploratory_annotation_level_preview_select,
                                                           analysis_type=input$exploratory_analysis_type_select,
                                                           metadata_feature=input$exploratory_feature_select,
                                                           low_pass_filter_threshold=filtering_value,
                                                           norm_method=input$normalization_select,
                                                           sample_depth=sample_depth,
                                                           subset_strings=subset_string,
                                                           model_string=covar_str,
                                                           random_effect=rand_effect,
                                                           pval_threshold=as.numeric(input$exploratory_statistical_pvalue_slider),
                                                           num_tophits=as.numeric(input$exploratory_statistical_number_slider),
                                                           sort_by=input$exploratory_statistical_sortby_select)
                    }
                    else if(input$exploratory_data_type_preview_select == 'Microbiome') {
                        dt <- generate_statistical_preview(data=list(microbiome_data()),
                                                           data_type='Microbiome',
                                                           metadata=metadata_updated$data,
                                                           annotation_level=input$exploratory_annotation_level_preview_select,
                                                           analysis_type=input$exploratory_analysis_type_select,
                                                           metadata_feature=input$exploratory_feature_select,
                                                           low_pass_filter_threshold=filtering_value,
                                                           norm_method=input$normalization_select,
                                                           sample_depth=sample_depth,
                                                           subset_strings=subset_string,
                                                           model_string=covar_str,
                                                           random_effect=rand_effect,
                                                           pval_threshold=as.numeric(input$exploratory_statistical_pvalue_slider),
                                                           num_tophits=as.numeric(input$exploratory_statistical_number_slider),
                                                           sort_by=input$exploratory_statistical_sortby_select)
                    }
                    
                    output$exploratory_preview <- renderUI({
                        dataTableOutput(outputId = 'exploratory_output')
                    })
                    output$exploratory_output <- renderDataTable({
                        dt
                    })
                }
            }
        })
    })

    # Reactive data structure for storing experimental designs
    experimental_designs <- reactiveValues(names=list(),
                                           experiments=list(),
                                           activity=list(),
                                           taglist=list())
    
    
    # Function to make observers dynamically for the experimental designs
    # We will need a primary observer for remove event, and observers for 
    # the other values as well to update the data structure for persistent
    # values
    makeObservers <- eventReactive(input$experimental_design_add, {
        local_activity <- experimental_designs$activity
        
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
        
        res <- lapply(1:length(local_activity), function(i) {
            if(local_activity[[i]]) {
                # Create an observer for updating random effect field on change of primary feature
                observe({
                    x <- input[[paste(experimental_designs$names[[i]], '_feature_select', sep='', collapse='')]]
                    updateSelectInput(session = session,
                                      inputId = paste(experimental_designs$names[[i]], '_random_effect_select', sep='', collapse=''),
                                      choices = c('None', input[[paste(experimental_designs$names[[i]], '_selected_features_boxes', sep='', collapse='')]][
                                          input[[paste(experimental_designs$names[[i]], '_selected_features_boxes', sep='', collapse='')]] != x
                                      ]))
                })
                
                # Create an observer for removal
                observeEvent(input[[paste(experimental_designs$names[[i]], '_delete', sep='', collapse='')]], {
                    experimental_designs$activity[[i]] <- FALSE
                    experimental_designs$taglist[[i]] <- character(0)
                    
                    # Update current values into data structure
                    exp_list <- experimental_designs$taglist
                    if(length(experimental_designs$activity) > 0) {
                        for(i in 1:(length(experimental_designs$activity))) {
                            if(experimental_designs$activity[[i]]) {
                                exp_list[[i]] <- tagList(
                                    box(
                                        splitLayout(cellWidths = c('50%', '50%'),
                                                    textInput(inputId = paste(experimental_designs$names[[i]], '_name', sep='', collapse=''),
                                                              label = 'Experiment Name',
                                                              value = input[[paste(experimental_designs$names[[i]], '_name', sep='', collapse='')]]),
                                                    actionButton(inputId = paste(experimental_designs$names[[i]], '_delete', sep='', collapse=''),
                                                                 label = 'Remove Experiment',
                                                                 width = '70%')
                                        ),
                                        selectInput(inputId = paste(experimental_designs$names[[i]], '_feature_select', sep='', collapse=''),
                                                    label = 'Select a Primary Feature',
                                                    choices = experimental_designs$experiments[[i]][['active_features']],
                                                    selected = input[[paste(experimental_designs$names[[i]], '_feature_select', sep='', collapse='')]]),
                                        selectInput(inputId = paste(experimental_designs$names[[i]], '_analysis_select', sep='', collapse=''),
                                                    label = 'Select Regression Type',
                                                    choices = c('ZIG Regression'),
                                                    selected = input[[paste(experimental_designs$names[[i]], '_analysis_select', sep='', collapse='')]]),
                                        selectInput(inputId = paste(experimental_designs$names[[i]], '_normalization_select', sep='', collapse=''),
                                                    label = 'Select Normalization Method',
                                                    choices = c('Cumulative Sum Scaling', 'Rarefaction'),
                                                    selected = input[[paste(experimental_designs$names[[i]], '_normalization_select', sep='', collapse='')]]),
                                        sliderInput(inputId = paste(experimental_designs$names[[i]], '_sample_depth_slider', sep='', collapse = ''),
                                                    label = 'Rarefaction Depth (Ranked Low to High)',
                                                    min = 1,
                                                    max = min(min_samples),
                                                    step = 1,
                                                    value = input[[paste(experimental_designs$names[[i]], '_sample_depth_slider', sep='', collapse = '')]]),
                                        sliderInput(inputId = paste(experimental_designs$names[[i]], '_filter_slider', sep='', collapse=''),
                                                    label = 'Filtering Threshold (Quantile)',
                                                    min = 0,
                                                    max = 100,
                                                    step = 1,
                                                    value = input[[paste(experimental_designs$names[[i]], '_filter_slider', sep='', collapse='')]]),
                                        sliderInput(inputId = paste(experimental_designs$names[[i]], '_pval', sep='', collapse=''),
                                                    label = 'P-Value for Significance',
                                                    min = 0.01,
                                                    max = 0.2,
                                                    step = 0.01,
                                                    value = input[[paste(experimental_designs$names[[i]], '_pval', sep='', collapse='')]]),
                                        checkboxGroupInput(inputId = paste(experimental_designs$names[[i]], '_selected_features_boxes', sep='', collapse=''),
                                                           label = 'Features to Include in This Experiment',
                                                           choices = experimental_designs$experiments[[i]][['active_features']],
                                                           selected = input[[paste(experimental_designs$names[[i]], '_selected_features_boxes', sep='', collapse='')]]),
                                        selectInput(inputId = paste(experimental_designs$names[[i]], '_random_effect_select', sep='', collapse=''),
                                                    label = 'Random Effect for This Experiment',
                                                    choices = c('None', input[[paste(experimental_designs$names[[i]], '_selected_features_boxes', sep='', collapse='')]][
                                                        input[[paste(experimental_designs$names[[i]], '_feature_select', sep='', collapse='')]] !=
                                                            input[[paste(experimental_designs$names[[i]], '_selected_features_boxes', sep='', collapse='')]]
                                                        ]),
                                                    selected = input[[paste(experimental_designs$names[[i]], '_random_effect_select', sep='', collapse='')]])
                                    )
                                )
                            }
                        }
                    }
                    experimental_designs$taglist <- exp_list
                    output$experimental_design_boxes <- renderUI({
                        exp_list
                    })
                })
            }
        })
        res
    })
    
    
    # Add an experiment button
    observeEvent(input$experimental_design_add, {
        # Default experiment object
        def_experiment <- list(name='',
                               analysis_type='ZIG Regression',
                               filter_threshold=15,
                               normalization='Cumulative Sum Scaling',
                               sample_depth=1,
                               primary_feature=active_metadata_fields$data[1],
                               active_features=active_metadata_fields$data,
                               selected_features=active_metadata_fields$data,
                               random_effect='None',
                               pval=0.1,
                               tophits=100)
        
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
        
        experimental_designs$names[[length(experimental_designs$names) + 1]] <- paste('exp_',
                                                                                      length(experimental_designs$names) + 1,
                                                                                      sep='', collapse='')
        experimental_designs$experiments[[length(experimental_designs$experiments) + 1]] <- def_experiment
        experimental_designs$activity[[length(experimental_designs$activity) + 1]] <- TRUE
        
        exp_list <- experimental_designs$taglist
        i <- length(experimental_designs$taglist) + 1
        exp_list[[i]] <- tagList(
            box(
                splitLayout(cellWidths = c('50%', '50%'),
                            textInput(inputId = paste(experimental_designs$names[[i]], '_name', sep='', collapse=''),
                                      label = 'Experiment Name',
                                      value = paste('Experiment ', i, sep='', collapse='')),
                            actionButton(inputId = paste(experimental_designs$names[[i]], '_delete', sep='', collapse=''),
                                         label = 'Remove Experiment',
                                         width = '70%')
                ),
                selectInput(inputId = paste(experimental_designs$names[[i]], '_feature_select', sep='', collapse=''),
                            label = 'Select a Primary Feature',
                            choices = experimental_designs$experiments[[i]][['active_features']],
                            selected = experimental_designs$experiments[[i]][['primary_feature']]),
                selectInput(inputId = paste(experimental_designs$names[[i]], '_analysis_select', sep='', collapse=''),
                            label = 'Select Regression Type',
                            choices = c('ZIG Regression'),
                            selected = experimental_designs$experiments[[i]][['analysis_type']]),
                selectInput(inputId = paste(experimental_designs$names[[i]], '_normalization_select', sep='', collapse=''),
                            label = 'Select Normalization Method',
                            choices = c('Cumulative Sum Scaling', 'Rarefaction'),
                            selected = experimental_designs$experiments[[i]][['normalization']]),
                sliderInput(inputId = paste(experimental_designs$names[[i]], '_sample_depth_slider', sep='', collapse = ''),
                            label = 'Rarefaction Depth (Ranked Low to High)',
                            min = 1,
                            max = min(min_samples),
                            step = 1,
                            value = experimental_designs$experiments[[i]][['sample_depth']]),
                sliderInput(inputId = paste(experimental_designs$names[[i]], '_filter_slider', sep='', collapse=''),
                            label = 'Filtering Threshold (Quantile)',
                            min = 0,
                            max = 100,
                            step = 1,
                            value = experimental_designs$experiments[[i]][['filter_threshold']]),
                sliderInput(inputId = paste(experimental_designs$names[[i]], '_pval', sep='', collapse=''),
                            label = 'P-Value for Significance',
                            min = 0.01,
                            max = 0.2,
                            step = 0.01,
                            value = experimental_designs$experiments[[i]][['pval']]),
                checkboxGroupInput(inputId = paste(experimental_designs$names[[i]], '_selected_features_boxes', sep='', collapse=''),
                                   label = 'Features to Include in This Experiment',
                                   choices = experimental_designs$experiments[[i]][['active_features']],
                                   selected = experimental_designs$experiments[[i]][['selected_features']]),
                selectInput(inputId = paste(experimental_designs$names[[i]], '_random_effect_select', sep='', collapse=''),
                            label = 'Random Effect for This Experiment',
                            choices = c('None', experimental_designs$experiments[[i]][['active_features']][
                                experimental_designs$experiments[[i]][['primary_feature']] !=
                                    experimental_designs$experiments[[i]][['active_features']]
                                ]),
                            selected = 'None')
            )
        )
        
        # Update current values into data structure
        if(length(experimental_designs$activity) > 1) {
            for(i in 1:(length(experimental_designs$activity) - 1)) {
                if(experimental_designs$activity[[i]]) {
                    exp_list[[i]] <- tagList(
                        box(
                            splitLayout(cellWidths = c('50%', '50%'),
                                        textInput(inputId = paste(experimental_designs$names[[i]], '_name', sep='', collapse=''),
                                                  label = 'Experiment Name',
                                                  value = input[[paste(experimental_designs$names[[i]], '_name', sep='', collapse='')]]),
                                        actionButton(inputId = paste(experimental_designs$names[[i]], '_delete', sep='', collapse=''),
                                                     label = 'Remove Experiment',
                                                     width = '70%')
                            ),
                            selectInput(inputId = paste(experimental_designs$names[[i]], '_feature_select', sep='', collapse=''),
                                        label = 'Select a Primary Feature',
                                        choices = experimental_designs$experiments[[i]][['active_features']],
                                        selected = input[[paste(experimental_designs$names[[i]], '_feature_select', sep='', collapse='')]]),
                            selectInput(inputId = paste(experimental_designs$names[[i]], '_analysis_select', sep='', collapse=''),
                                        label = 'Select Regression Type',
                                        choices = c('ZIG Regression'),
                                        selected = input[[paste(experimental_designs$names[[i]], '_analysis_select', sep='', collapse='')]]),
                            selectInput(inputId = paste(experimental_designs$names[[i]], '_normalization_select', sep='', collapse=''),
                                        label = 'Select Normalization Method',
                                        choices = c('Cumulative Sum Scaling', 'Rarefaction'),
                                        selected = input[[paste(experimental_designs$names[[i]], '_normalization_select', sep='', collapse='')]]),
                            sliderInput(inputId = paste(experimental_designs$names[[i]], '_sample_depth_slider', sep='', collapse = ''),
                                        label = 'Rarefaction Depth (Ranked Low to High)',
                                        min = 1,
                                        max = min(min_samples),
                                        step = 1,
                                        value = input[[paste(experimental_designs$names[[i]], '_sample_depth_slider', sep='', collapse = '')]]),
                            sliderInput(inputId = paste(experimental_designs$names[[i]], '_filter_slider', sep='', collapse=''),
                                        label = 'Filtering Threshold (Quantile)',
                                        min = 0,
                                        max = 100,
                                        step = 1,
                                        value = input[[paste(experimental_designs$names[[i]], '_filter_slider', sep='', collapse='')]]),
                            sliderInput(inputId = paste(experimental_designs$names[[i]], '_pval', sep='', collapse=''),
                                        label = 'P-Value for Significance',
                                        min = 0.01,
                                        max = 0.2,
                                        step = 0.01,
                                        value = input[[paste(experimental_designs$names[[i]], '_pval', sep='', collapse='')]]),
                            checkboxGroupInput(inputId = paste(experimental_designs$names[[i]], '_selected_features_boxes', sep='', collapse=''),
                                               label = 'Features to Include in This Experiment',
                                               choices = experimental_designs$experiments[[i]][['active_features']],
                                               selected = input[[paste(experimental_designs$names[[i]], '_selected_features_boxes', sep='', collapse='')]]),
                            selectInput(inputId = paste(experimental_designs$names[[i]], '_random_effect_select', sep='', collapse=''),
                                        label = 'Random Effect for This Experiment',
                                        choices = c('None', input[[paste(experimental_designs$names[[i]], '_selected_features_boxes', sep='', collapse='')]][
                                            input[[paste(experimental_designs$names[[i]], '_feature_select', sep='', collapse='')]] !=
                                                input[[paste(experimental_designs$names[[i]], '_selected_features_boxes', sep='', collapse='')]]
                                            ]),
                                        selected = input[[paste(experimental_designs$names[[i]], '_random_effect_select', sep='', collapse='')]])
                        )
                    )
                }
            }
        }
        
        
        experimental_designs$taglist <- exp_list
        output$experimental_design_boxes <- renderUI({
            exp_list
        })
        makeObservers()
    })
    
    # Observer for rarefaction as a choice in experimental design exploratory parameters
    observe({
        x <- input$experimental_design_exploratory_norm_select
        if(x == 'Rarefaction') {
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
                output$experimental_design_exploratory_sample_threshold <- renderUI({
                    sliderInput(inputId = 'experimental_design_exploratory_rarefaction_slider',
                                label = 'Rarefy to Ranked Sample # (lowest to highest):',
                                min = 1,
                                max = min(min_samples),
                                step = 1,
                                value = 1)
                })
            }
        }
        else if(x == 'Cumulative Sum Scaling') {
            output$experimental_design_exploratory_sample_threshold <- renderUI({
                character(0)
            })
        }
    })
    
    # Persistent data structure for tempdir handling
    temp_reactive <- reactiveVal(value = character(0))
    
    # Run analyses button
    observeEvent(input$experimental_design_run, {
        temp_dir <- tempdir()
        temp_reactive(temp_dir)
        
        data_present <- c()
        if(!is.null(microbiome_data())) {
            data_present <- c(data_present, 'Microbiome')
        }
        if(!is.null(amr_counts()) && !is.null(amr_annotations())) {
            data_present <- c(data_present, 'Resistome')
        }
        
        isolate({
            # Generate output directories
            output_prefix <- create_output_directories(temp_dir=temp_dir,
                                                       project_name=input$project_name,
                                                       data_types=data_present,
                                                       active_features=active_metadata_fields$data,
                                                       experiment_names=experimental_designs$names,
                                                       experiment_activity=experimental_designs$activity)
            
            filter_quantile <- input$experimental_design_exploratory_filter_slider / 100
            if('Resistome' %in% data_present) {
                withProgress(message = 'Running Resistome Analyses: ', value = 0, {
                    n <- (length(active_metadata_fields$data) * 4) + sum(unlist(experimental_designs$activity)) + 1
                    incProgress(amount = 1/n, detail = 'Normalization and Aggregation')
                    if(input$experimental_design_exploratory_norm_select == 'Rarefaction') {
                        amr_data_list <- rarefy_normalize_and_extract(sampling_depth=input$experimental_design_exploratory_rarefaction_slider,
                                                                      filtering_quantile=filter_quantile,
                                                                      amr_MRexp=amr_counts())
                        amr_data_list[[7]] <- amr_annotations()
                    }
                    else if(input$experimental_design_exploratory_norm_select == 'Cumulative Sum Scaling') {
                        amr_data_list <- CSS_normalize_and_extract(filtering_quantile=filter_quantile,
                                                                   amr_MRexp=amr_counts())
                        amr_data_list[[7]] <- amr_annotations()
                    }
                    
                    # Aggregation and data set creation
                    analytic_list <- aggregate_and_filter(amr_data_list, data.table(metadata_updated$data))
                    
                    # Exploratory analyses
                    for(i in 1:length(active_metadata_fields$data)) {
                        # Ordination
                        incProgress(amount = 1/n,
                                    detail = paste('Ordination', active_metadata_fields$data[[i]], sep=' '))
                        ordination_workflow(experiment_data_lists=analytic_list,
                                            sample_column_id=colnames(metadata_updated$data)[1],
                                            feature_name=active_metadata_fields$data[[i]],
                                            graph_output_dir=paste(output_prefix, 'Graphs', sep='/'),
                                            data_type='Resistome')
                        # Barplot
                        incProgress(amount = 1/n,
                                    detail = paste('Barplots', active_metadata_fields$data[[i]], sep=' '))
                        barplot_workflow(experiment_data_lists=analytic_list,
                                         sample_column_id=colnames(metadata_updated$data)[1],
                                         feature_name=active_metadata_fields$data[[i]],
                                         graph_output_dir=paste(output_prefix, 'Graphs', sep='/'),
                                         data_type='Resistome')
                        # Heatmap
                        incProgress(amount = 1/n,
                                    detail = paste('Heatmaps', active_metadata_fields$data[[i]], sep=' '))
                        heatmap_workflow(experiment_data_lists=analytic_list,
                                         sample_column_id=colnames(metadata_updated$data)[1],
                                         feature_name=active_metadata_fields$data[[i]],
                                         graph_output_dir=paste(output_prefix, 'Graphs', sep='/'),
                                         data_type='Resistome')
                        # Alpha rarefaction
                        incProgress(amount = 1/n,
                                    detail = paste('Alpha Rarefaction', active_metadata_fields$data[[i]], sep=' '))
                        alpha_rarefaction_workflow(experiment_data_lists=analytic_list,
                                                   sample_column_id=colnames(metadata_updated$data)[1],
                                                   feature_name=active_metadata_fields$data[[i]],
                                                   graph_output_dir=paste(output_prefix, 'Graphs', sep='/'),
                                                   data_type='Resistome')
                    }
                    
                    # Statistical analyses
                    if(length(experimental_designs$activity) > 0) {
                        for(i in 1:length(experimental_designs$activity)) {
                            if(experimental_designs$activity[[i]]) {
                                # Regression
                                exp_name <- input[[paste(experimental_designs$names[[i]], '_name', sep='', collapse='')]]
                                analysis_type <- input[[paste(experimental_designs$names[[i]], '_analysis_select', sep='', collapse='')]]
                                feature <- input[[paste(experimental_designs$names[[i]], '_feature_select', sep='', collapse='')]]
                                filtering_value <- input[[paste(experimental_designs$names[[i]], '_filter_slider', sep='', collapse='')]] / 100
                                norm_method <- input[[paste(experimental_designs$names[[i]], '_normalization_select', sep='', collapse='')]]
                                sample_depth <- input[[paste(experimental_designs$names[[i]], '_sample_depth_slider', sep='', collapse = '')]]
                                subset_string <- c()
                                checkboxes <- input[[paste(experimental_designs$names[[i]], '_selected_features_boxes', sep='', collapse='')]]
                                rand_select <- input[[paste(experimental_designs$names[[i]], '_random_effect_select', sep='', collapse='')]]
                                
                                covars <- checkboxes[!(checkboxes %in% rand_select)]
                                primary_effect_idx <- which(covars == feature)
                                other_idx <- which(1:length(covars) != primary_effect_idx)
                                covars <- covars[c(primary_effect_idx, other_idx)]
                                covars <- c('~ 0', covars)
                                covar_str <- paste(covars, sep='', collapse=' + ')
                                
                                rand_effect <- NA
                                if(rand_select != 'None') {
                                    rand_effect <- rand_select
                                }
                                
                                pval_threshold <- as.numeric(input[[paste(experimental_designs$names[[i]], '_pval', sep='', collapse='')]])
                                
                                incProgress(amount = 1/n,
                                            detail = paste('Regression', experimental_designs$names[[i]], sep=' '))
                                for(l in names(amr_lookup)) {
                                    dt <- generate_statistical_preview(data=list(amr_counts(), amr_annotations()),
                                                                      data_type='Resistome',
                                                                      metadata=metadata_updated$data,
                                                                      annotation_level=l,
                                                                      analysis_type=analysis_type,
                                                                      metadata_feature=feature,
                                                                      low_pass_filter_threshold=filtering_value,
                                                                      norm_method=norm_method,
                                                                      sample_depth=sample_depth,
                                                                      subset_strings=subset_string,
                                                                      model_string=covar_str,
                                                                      random_effect=rand_effect,
                                                                      pval_threshold=pval_threshold,
                                                                      num_tophits=9000,
                                                                      sort_by='P-value')
                                    if(nrow(dt) > 0) {
                                        fname <- paste(gsub(' ', '_', exp_name), l, gsub(' ', '_', analysis_type), sep='_')
                                        write.csv(dt, file=paste(paste(output_prefix, 'Statistics', 'Resistome',
                                                                       experimental_designs$names[[i]],
                                                                       fname, sep='/'), '.csv', sep='', collapse=''))
                                    }
                                }
                                
                            }
                        }
                    }
                })
            }
            if('Microbiome' %in% data_present) {
                withProgress(message = 'Running Microbiome Analyses: ', value = 0, {
                    n <- (length(active_metadata_fields$data) * 4) + sum(unlist(experimental_designs$activity)) + 1
                    incProgress(amount = 1/n, detail = 'Normalization and Aggregation')
                    # Normalization
                    if(input$experimental_design_exploratory_norm_select == 'Rarefaction') {
                        microbiome_data_list <- rarefy_normalize_and_extract(sampling_depth=input$experimental_design_exploratory_rarefaction_slider,
                                                                             filtering_quantile=filter_quantile,
                                                                             kraken_MRexp=microbiome_data())
                        microbiome_data_list[[7]] <- NA
                    }
                    else if(input$experimental_design_exploratory_norm_select == 'Cumulative Sum Scaling') {
                        microbiome_data_list <- CSS_normalize_and_extract(filtering_quantile=filter_quantile,
                                                                          kraken_MRexp=microbiome_data())
                        microbiome_data_list[[7]] <- NA
                    }
                    
                    # Aggregation and data set creation
                    analytic_list <- aggregate_and_filter(microbiome_data_list, data.table(metadata_updated$data))
                    
                    # Exploratory analyses
                    for(i in 1:length(active_metadata_fields$data)) {
                        # Ordination
                        incProgress(amount = 1/n,
                                    detail = paste('Ordination', active_metadata_fields$data[[i]], sep=' '))
                        ordination_workflow(experiment_data_lists=analytic_list,
                                            sample_column_id=colnames(metadata_updated$data)[1],
                                            feature_name=active_metadata_fields$data[[i]],
                                            graph_output_dir=paste(output_prefix, 'Graphs', sep='/'),
                                            data_type='Microbiome')
                        # Barplot
                        incProgress(amount = 1/n,
                                    detail = paste('Barplots', active_metadata_fields$data[[i]], sep=' '))
                        barplot_workflow(experiment_data_lists=analytic_list,
                                         sample_column_id=colnames(metadata_updated$data)[1],
                                         feature_name=active_metadata_fields$data[[i]],
                                         graph_output_dir=paste(output_prefix, 'Graphs', sep='/'),
                                         data_type='Microbiome')
                        # Heatmap
                        incProgress(amount = 1/n,
                                    detail = paste('Heatmaps', active_metadata_fields$data[[i]], sep=' '))
                        heatmap_workflow(experiment_data_lists=analytic_list,
                                         sample_column_id=colnames(metadata_updated$data)[1],
                                         feature_name=active_metadata_fields$data[[i]],
                                         graph_output_dir=paste(output_prefix, 'Graphs', sep='/'),
                                         data_type='Microbiome')
                        # Alpha rarefaction
                        incProgress(amount = 1/n,
                                    detail = paste('Alpha Rarefaction', active_metadata_fields$data[[i]], sep=' '))
                        alpha_rarefaction_workflow(experiment_data_lists=analytic_list,
                                                   sample_column_id=colnames(metadata_updated$data)[1],
                                                   feature_name=active_metadata_fields$data[[i]],
                                                   graph_output_dir=paste(output_prefix, 'Graphs', sep='/'),
                                                   data_type='Microbiome')
                    }
                    
                    # Statistical analyses
                    if(length(experimental_designs$activity) > 0) {
                        for(i in 1:length(experimental_designs$activity)) {
                            if(experimental_designs$activity[[i]]) {
                                # Regression
                                exp_name <- input[[paste(experimental_designs$names[[i]], '_name', sep='', collapse='')]]
                                analysis_type <- input[[paste(experimental_designs$names[[i]], '_analysis_select', sep='', collapse='')]]
                                feature <- input[[paste(experimental_designs$names[[i]], '_feature_select', sep='', collapse='')]]
                                filtering_value <- input[[paste(experimental_designs$names[[i]], '_filter_slider', sep='', collapse='')]] / 100
                                norm_method <- input[[paste(experimental_designs$names[[i]], '_normalization_select', sep='', collapse='')]]
                                sample_depth <- input[[paste(experimental_designs$names[[i]], '_sample_depth_slider', sep='', collapse = '')]]
                                subset_string <- c()
                                checkboxes <- input[[paste(experimental_designs$names[[i]], '_selected_features_boxes', sep='', collapse='')]]
                                rand_select <- input[[paste(experimental_designs$names[[i]], '_random_effect_select', sep='', collapse='')]]
                                
                                covars <- checkboxes[!(checkboxes %in% rand_select)]
                                primary_effect_idx <- which(covars == feature)
                                other_idx <- which(1:length(covars) != primary_effect_idx)
                                covars <- covars[c(primary_effect_idx, other_idx)]
                                covars <- c('~ 0', covars)
                                covar_str <- paste(covars, sep='', collapse=' + ')

                                rand_effect <- NA
                                if(rand_select != 'None') {
                                    rand_effect <- rand_select
                                }
                                
                                pval_threshold <- as.numeric(input[[paste(experimental_designs$names[[i]], '_pval', sep='', collapse='')]])
                                
                                incProgress(amount = 1/n,
                                            detail = paste('Regression', experimental_designs$names[[i]], sep=' '))
                                for(l in names(microbiome_lookup)) {
                                    dt <- generate_statistical_preview(data=list(microbiome_data()),
                                                                      data_type='Microbiome',
                                                                      metadata=metadata_updated$data,
                                                                      annotation_level=l,
                                                                      analysis_type=analysis_type,
                                                                      metadata_feature=feature,
                                                                      low_pass_filter_threshold=filtering_value,
                                                                      norm_method=norm_method,
                                                                      sample_depth=sample_depth,
                                                                      subset_strings=subset_string,
                                                                      model_string=covar_str,
                                                                      random_effect=rand_effect,
                                                                      pval_threshold=pval_threshold,
                                                                      num_tophits=9000,
                                                                      sort_by='P-value')
                                    if(nrow(dt) > 0) {
                                        fname <- paste(gsub(' ', '_', exp_name), l, gsub(' ', '_', analysis_type), sep='_')
                                        write.csv(dt, file=paste(paste(output_prefix, 'Statistics', 'Microbiome',
                                                                       experimental_designs$names[[i]],
                                                                       fname, sep='/'), '.csv', sep='', collapse=''))
                                    }
                                }
                            }
                        }
                    }
                })
            }
        })
    })
    
    # Download resulting data
    output$experimental_design_download <- downloadHandler(
        filename = function() {
            paste(gsub(' ', '_', input$project_name), '.zip', sep='')
        },
        content = function(con) {
            utils::zip(zipfile=con,
                       files=list.files(paste(temp_reactive(), input$project_name, sep='/'),
                                        recursive = TRUE,
                                        full.names = TRUE))
        },
        contentType = 'application/zip'
    )
    
    # # Adonis PERMANOVA parameter code
    # # Data type
    # observe({
    #     data_present <- c()
    #     if(!is.null(amr_counts()) && !is.null(amr_annotations())) {
    #         data_present <- c(data_present, 'Resistome')
    #     }
    #     if(!is.null(microbiome_data())) {
    #         data_present <- c(data_present, 'Microbiome')
    #     }
    #     
    #     updateSelectInput(session = session,
    #                       inputId = 'adonis_data_type',
    #                       choices = data_present)
    # })
    # 
    # # Annotation Level
    # observe({
    #     x <- input$adonis_data_type
    #     if(x == 'Resistome') {
    #         updateSelectInput(session = session,
    #                           inputId = 'adonis_annotation_level',
    #                           choices = names(amr_lookup))
    #     }
    #     else if(x == 'Microbiome') {
    #         updateSelectInput(session = session,
    #                           inputId = 'adonis_annotation_level',
    #                           choices = names(microbiome_lookup))
    #     }
    # })
    # 
    # # Rarefaction sampling depth
    # observe({
    #     x <- input$adonis_normalization
    #     if(x == 'Cumulative Sum Scaling') {
    #         output$adonis_rarefaction_slider <- renderUI({
    #             character(0)
    #         })
    #     }
    #     else if(x == 'Rarefaction') {
    #         min_samples <- c()
    #         isolate({
    #             if(!is.null(metadata()))  {
    #                 min_samples <- c(min_samples, nrow(metadata()))
    #             }
    #             if(!is.null(amr_counts())) {
    #                 min_samples <- c(min_samples, ncol(amr_counts()))
    #             }
    #             if(!is.null(microbiome_data())) {
    #                 min_samples <- c(min_samples, ncol(microbiome_data()))
    #             }
    #         })
    #         if(length(min_samples) > 0) {
    #             output$adonis_rarefaction_slider <- renderUI({
    #                 sliderInput(inputId = 'adonis_rarefaction_sampling_depth',
    #                             label = 'Rarefy to Ranked Sample # (lowest to highest):',
    #                             min = 1,
    #                             max = min(min_samples),
    #                             step = 1,
    #                             value = 1)
    #             })
    #         }
    #     }
    # })
    # 
    # # Metadata features
    # observe({
    #     x <- active_metadata_fields$data
    #     updateSelectInput(session = session,
    #                       inputId = 'adonis_feature_select',
    #                       choices = x)
    #     updateCheckboxGroupInput(session = session,
    #                              inputId = 'adonis_checkbox_features',
    #                              choices = x,
    #                              selected = x)
    #     updateSelectInput(session = session,
    #                       inputId = 'adonis_strata_feature',
    #                       choices = c('None', x),
    #                       selected = 'None')
    # })
    # 
    # # Feature types
    # observe({
    #     x <- input$adonis_checkbox_features
    #     y <- input$adonis_feature_select
    #     z <- input$adonis_strata_feature
    #     
    #     choices <- x[!(x %in% z)]
    #     
    #     if(length(choices) > 0 && (z != 'None')) {
    #         output$adonis_nested_choices <- renderUI({
    #             checkboxGroupInput(inputId = 'adonis_nested_choices',
    #                                label = 'Factors Nested Under Strata Choice',
    #                                choices = choices)
    #         })
    #     }
    #     else {
    #         output$adonis_nested_choices <- renderUI({
    #             character(0)
    #         })
    #     }
    # })
    # 
    # # PERMANOVA run analysis
    # observeEvent(input$adonis_run_model, {
    #     isolate({
    #         filtering_threshold <- as.numeric(input$adonis_filter_threshold) / 100
    #         
    #         if(input$adonis_normalization == 'Rarefaction') {
    #             depth <- input$adonis_rarefaction_sampling_depth
    #         }
    #         
    #         strata <- input$adonis_strata_feature
    #         if(strata != 'None') {
    #             nested_features <- input$adonis_nested_choices
    #             fixed_features <- input$adonis_checkbox_features[!(input$adonis_checkbox_features %in% strata) &&
    #                                                                  !(input$adonis_checkbox_features %in% nested_features)]
    #         }
    #         
    #         if(input$adonis_data_type == 'Resistome') {
    #             dt <- generate_permanova_table(data=list(microbiome_data()),
    #                                            data_type='Resistome',
    #                                            metadata=metadata_updated$data,
    #                                            annotation_level=input$adonis_annotation_level,
    #                                            metadata_feature=input$adonis_feature_select,
    #                                            low_pass_filter_threshold=filtering_threshold,
    #                                            norm_method=input$adonis_normalization,
    #                                            sample_depth=depth,
    #                                            strata=strata,
    #                                            nested_features=nested_features,
    #                                            fixed_features=fixed_features)
    #         }
    #     })
    #     
    #     
    # })
}

