#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

library(shinydashboard)

ui <- dashboardPage(
  dashboardHeader(title = "Basic dashboard"),
  dashboardSidebar(sidebarMenu(
    menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
    menuItem("Sample Information", tabName = "Sample Information", icon = icon("th")),
    menuItem("Counts Matrix Exploration", tabName = "Counts Matrix Exploration", icon = icon("th")),
    menuItem("Differential Expression", tabName = "Differential Expression", icon = icon("th")),
    menuItem("Gene Set Enrichment Analysis", tabName = "Gene Set Enrichment Analysis", icon = icon("th"))
    )),
  dashboardBody(
    tabItems(
      # First tab content
      tabItem(tabName = "dashboard",
              fluidRow(
                box(plotOutput("plot1", height = 250)),
                
                box(
                  title = "Controls",
                  sliderInput("slider", "Number of observations:", 1, 100, 50)
                )
              )
      ),
      
      # Second tab content
      tabItem(tabName = "Sample Information",
              h2("Add summary table")
      ),
      tabItem(tabName = "Counts Matrix Exploration",
              h2("count matrix stuff")
      ),
      tabItem(tabName = "Differential Expression",
              h2("Add Differential Expression")
      ),
      tabItem(tabName = "Gene Set Enrichment Analysis",
              h2("Add Gene Set Enrichment Analysis")
      )
    )
),

server <- function(input, output) {
  set.seed(122)
  histdata <- rnorm(500)
  
  output$plot1 <- renderPlot({
    data <- histdata[seq_len(input$slider)]
    hist(data)
  })
}

shinyApp(ui, server)
