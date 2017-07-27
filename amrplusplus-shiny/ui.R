ui <- dashboardPage(
    dashboardHeader(title = 'AMR++ Shiny'),
    
    dashboardSidebar(
        sidebarMenu(
            menuItem('Load Data', tabName = 'load_data', icon = icon('upload')),
            menuItem('Format Metadata', tabName = 'format_metadata', icon = icon('th')),
            menuItem('Configure Experiments', tabName = 'configure_experiments', icon = icon('cog')),
            menuItem('View Results', tabName = 'view_results', icon = icon('bar-chart'))
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
                        must match the column names of your count matrices'
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
                        )
                    )
            )
        )
    )
)

ui
