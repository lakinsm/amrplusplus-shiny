## app.R ##
library(shinydashboard)

source('dataIO.R')
source('data_processing.R')
source('exploratory_analysis.R')
source('statistical_analysis.R')

ui <- dashboardPage(
    dashboardHeader(title = "Basic dashboard"),
    
    dashboardSidebar(
        sidebarMenu(
            menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
            menuItem("Load Data", tabName = "load_data", icon = icon("upload")),
            menuItem("Configure Experiments", tabName = "configure_experiments", icon = icon("cog")),
            menuItem("View Results", tabName = "view_results", icon = icon("bar-chart"))
        )
    ),
    
    dashboardBody(
        # Boxes need to be put in a row (or column)
        fluidRow(
            box(plotOutput("plot1", height = 250)),
            
            box(
                title = "Controls",
                sliderInput("slider", "Number of observations:", 1, 100, 50)
            )
        )
    )
)

server <- function(input, output) {
    set.seed(122)
    histdata <- rnorm(500)
    
    output$plot1 <- renderPlot({
        data <- histdata[seq_len(input$slider)]
        hist(data)
    })
}

shinyApp(ui, server)

