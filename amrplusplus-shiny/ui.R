ui <- dashboardPage(
    dashboardHeader(title = 'AMR++ Shiny'),
    dashboardSidebar(
        sidebarMenu(
            menuItem('1. Load Data', tabName = 'load_data', icon = icon('upload')),
            menuItem('2. Format Metadata', tabName = 'format_metadata', icon = icon('th')),
            menuItem('3. Data Exploration', tabName = 'data_exploration', icon = icon('cog')),
            menuItem('4. Experimental Designs', tabName = 'experimental_designs', icon = icon('dashboard')),
            menuItem('5. PERMANOVA', tabName = 'permanova', icon = icon('bar-chart')),
            menuItem('6. Procrustes', tabName = 'procrustes', icon = icon('bar-chart'))
        )
    ),
    dashboardBody(
        # Boxes need to be put in a row (or column)
        tabItems(
            tabItem(tabName = 'load_data',
                fluidRow(
                    box(
                        textInput(
                            inputId = 'project_name',
                            label = 'Choose a Project Name',
                            value = 'My Project'
                        ),
                        fileInput(
                            inputId = 'amr_count_matrix',
                            label = 'Select AMR Count Matrix',
                            accept = c(
                                'text/csv',
                                'text/comma-separated-values,text/plain',
                                '.csv'
                            )
                        ),
                        fileInput(
                            inputId = 'amr_annotation_file',
                            label = 'Select AMR Annotation File',
                            accept = c(
                                'text/csv',
                                'text/comma-separated-values,text/plain',
                                '.csv'
                            )
                        ),
                        fileInput(
                            inputId = 'kraken_count_matrix',
                            label = 'Select Kraken Count Matrix',
                            accept = c(
                                'text/csv',
                                'text/comma-separated-values,text/plain',
                                '.csv'
                            )
                        ),
                        fileInput(
                            inputId = 'metadata_file',
                            label = 'Select Metadata File',
                            accept = c(
                                'text/csv',
                                'text/comma-separated-values,text/plain',
                                '.csv'
                            )
                        )
                    ),
                    box(
                        tags$h2('Loading Data into AMR++ Shiny'),
                        'AMR++ Shiny requires a metadata file describing your
                        samples as well as the following outputs from the AMR++
                        nextflow pipeline:',
                        tags$ol(
                            tags$li('an antimicrobial resistance count matrix'),
                            tags$li('a corresponding MEGARes annotation file'),
                            tags$li('a microbiome count matrix')
                        ),
                        tags$p(),
                        'In this menu, you should first type a name for your
                        project, then select the required files from your
                        computer to be uploaded into the AMR++ Shiny App.
                        The sample IDs in your metadata file (the first column)
                        must match the column names of your count matrices.  
                        Once you have completed data and metadata upload and 
                        validation with the Validate button, you are ready to 
                        proceed to the next step: Format Metadata.',
                        collapsible=TRUE
                    )
                ),
                fluidRow(
                    box(
                        actionButton(inputId = 'validate_inputs',
                                     label = 'Validate Inputs',
                                     width = '150px'),
                        tags$p(),
                        tags$h3('Validation Status:'),
                        verbatimTextOutput(outputId = 'validation_output'),
                        width=12,
                        background='light-blue'
                    )
                )
            ),
            tabItem(tabName = 'format_metadata',
                fluidRow(
                    box(
                        checkboxGroupInput(inputId = 'metadata_fields',
                                           label = 'Features to Use for Analysis',
                                           choices = character(0)),
                        uiOutput(outputId = 'metadata_types')
                    ),
                    box(
                        tags$h2('Choose Metadata Types'),
                        'The AMR++ Shiny App needs to understand how the 
                        features in your metadata are encoded so that 
                        analysis can be done properly.  On the left, select 
                        each feature in your metadata that you would like 
                        to be included in your analyses, then select from 
                        the dropdown menus the data type that best 
                        represents that feature.  Date/time data should be 
                        converted to equally spaced days prior to using 
                        the AMR++ Shiny App for analysis, then ordinal 
                        should be selected as the data type at this step.',
                        tags$p(),
                        'When you have finished your selections, click the 
                        Apply Metadata Types button at the bottom of the 
                        window.  A summary of your metadata table will be 
                        given below in R str() summary format.  You are then 
                        ready to proceed to the next step: Data Processing Setup.',
                        collapsible=TRUE
                    )
                ),
                fluidRow(
                    box(
                        actionButton(inputId = 'apply_metadata_types',
                                     label = 'Apply Metadata Types',
                                     width = '150px'),
                        tags$p(),
                        tags$h3('Metadata Status:'),
                        verbatimTextOutput(outputId = 'metadata_console'),
                        width=12,
                        background='light-blue'
                    )
                )
            ),
            tabItem(tabName = 'data_exploration',
                    fluidRow(
                        box(tags$h2('Explore Your Data in Real Time'),
                            'In this section, you can determine how combinations of data 
                            preprocessing parameters affect the exploratory and statistical
                            analysis of your samples.  You must have completed 
                            the previous steps in order for this section to be 
                            populated appropriately.',
                            tags$p(),
                            'Explore the following parameters from either 
                            resistome or microbiome data at a given annotation 
                            level, stratified by one of your metadata features:',
                            tags$ul(
                                tags$li(tags$b('Sample Normalization Method:'), 
                                        'to control for differences in sequencing 
                                        depth across samples'),
                                tags$li(tags$b('Low Pass Filter Threshold:'),
                                        'to remove low abundance taxa or genes 
                                        that are likely to be false positive 
                                        classifications')
                            ),
                            'Use the selections to choose the optimal parameter 
                            set for your experiment, and update the preview with 
                            the Update Exploratory Preview button.  When the 
                            optimal parameters have been chosen, proceed to the 
                            next section: Statistical Setup.',
                            width = 12,
                            collapsible=TRUE)
                    ),
                    fluidRow(
                        box(
                            tags$h3('Data Preprocessing Parameters'),
                            selectInput(inputId = 'normalization_select',
                                       label = 'Sample Normalization Method',
                                       choices = c('Cumulative Sum Scaling',
                                                   'Rarefaction'),
                                       selected='Cumulative Sum Scaling'),
                            uiOutput(outputId = 'exploratory_rarefaction_sample_threshold'),
                            sliderInput(inputId = 'filtering_slider',
                                        label = 'Low Pass Filter Threshold (Quantile)',
                                        min = 0,
                                        max = 100,
                                        value = 15,
                                        step = 1),
                            uiOutput(outputId = 'exploratory_statistical_parameters'),
                            actionButton(inputId = 'exploratory_update_preview',
                                         label = 'Update Exploratory Preview',
                                         width = '200px')
                        ),
                        box(
                            tags$h3('Exploratory Preview Settings'),
                            selectInput(inputId = 'exploratory_data_type_preview_select',
                                        label = 'Data Type',
                                        choices = character(0)),
                            selectInput(inputId = 'exploratory_annotation_level_preview_select',
                                        label = 'Annotation Level',
                                        choices = character(0)),
                            selectInput(inputId = 'exploratory_analysis_type_select',
                                        label = 'Analysis Type',
                                        choices = c('NMDS',
                                                    'PCA',
                                                    'Bar Graph',
                                                    'Heatmap',
                                                    'ZIG Regression'),
                                        selected='PCA'),
                            selectInput(inputId = 'exploratory_feature_select',
                                        label = 'Metadata Feature',
                                        choices = character(0)),
                            uiOutput(outputId = 'exploratory_subset_rules'),
                            actionButton(inputId = 'exploratory_add_subset',
                                         label = 'Add a Subset Rule',
                                         width = '200px'),
                            actionButton(inputId = 'exploratory_remove_subset',
                                         label = 'Remove a Subset Rule',
                                         width = '200px')
                        )
                    ),
                    fluidRow(
                        column(align='center',
                               width=12,
                               uiOutput(outputId = 'exploratory_preview')
                        )
                    )
                ),
            tabItem(
                tabName = 'experimental_designs',
                fluidRow(
                    box(width = 4,
                        tags$h3('Exploratory Analysis Parameters'),
                        selectInput(inputId = 'experimental_design_exploratory_norm_select',
                                    label = 'Exploratory Graphing Normalization Method',
                                    choices = c('Rarefaction', 'Cumulative Sum Scaling'),
                                    selected = 'Cumulative Sum Scaling'),
                        uiOutput(outputId = 'experimental_design_exploratory_sample_threshold'),
                        sliderInput(inputId = 'experimental_design_exploratory_filter_slider',
                                    label = 'Low Pass Filter Threshold (Quantile)',
                                    min = 0,
                                    max = 100,
                                    step = 1,
                                    value = 15)
                        
                    ),
                    box(width = 8,
                        collapsible = TRUE,
                        tags$h2('Setup and Run Your Experiments'),
                        'Now that you\'ve explored your data, you can set
                        up experiments in this tab to perform ordination, plotting,  
                        and regression-based analyses using different 
                        combinations of parameters.  Click on Add an Experiment to 
                        get started.  You can remove an experiment by clicking on 
                        Remove Experiment #, and it will not affect the status or
                         values of your other experiments.  When experimental setup 
                        is complete, run the experiments by clicking on Run All Experiments.
                          Graphs and tables will be zipped and available for download 
                        following completion of the analyses.'
                    )
                ),
                fluidRow(
                    box(width = 8,
                        tags$h3('Statistical Analysis Parameters'),
                        uiOutput(outputId = 'experimental_design_boxes'),
                        actionButton(inputId = 'experimental_design_add',
                                     label = 'Add an Experiment',
                                     width = '200px')
                    ),
                    box(width = 4,
                        background='light-blue',
                        actionButton(inputId = 'experimental_design_run',
                                     label = 'Run All Experiments',
                                     width = '200px'),
                        uiOutput(outputId = 'experimental_design_run_output')
                    )
                )
            )
        )
    )
)

ui
