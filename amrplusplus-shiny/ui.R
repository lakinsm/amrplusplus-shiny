ui <- dashboardPage(
    dashboardHeader(title = 'AMR++ Shiny'),
    dashboardSidebar(
        sidebarMenu(
            menuItem('1. Load Data', tabName = 'load_data', icon = icon('upload')),
            menuItem('2. Format Metadata', tabName = 'format_metadata', icon = icon('th')),
            menuItem('3. Data Preprocessing Setup', tabName = 'data_processing_setup', icon = icon('cog')),
            menuItem('4. View Results', tabName = 'view_results', icon = icon('bar-chart'))
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
                        proceed to the next step: Format Metadata.'
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
                        ready to proceed to the next step: Data Processing Setup.'
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
            tabItem(tabName = 'data_processing_setup',
                    fluidRow(
                        box(tags$h2('Explore Your Data in Real Time'),
                            'In this section, you can determine how combinations of data 
                            preprocessing parameters affect the exploratory 
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
                            width = 12)
                    ),
                    fluidRow(
                        box(
                            tags$h3('Data Preprocessing Parameters'),
                            selectInput(inputId = 'normalization_select',
                                       label = 'Sample Normalization Method',
                                       choices = c('Cumulative Sum Scaling',
                                                   'Rarefaction')),
                            uiOutput(outputId = 'exploratory_rarefaction_sample_threshold'),
                            sliderInput(inputId = 'filtering_slider',
                                        label = 'Low Pass Filter Threshold (Quantile)',
                                        min = 0,
                                        max = 100,
                                        value = 15,
                                        step = 1),
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
                                                    'Heatmap')),
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
                               plotOutput(outputId = 'exploratory_preview',
                                          width='100%',
                                          height='800px')
                        )
                    )
                )
        )
    )
)

ui
