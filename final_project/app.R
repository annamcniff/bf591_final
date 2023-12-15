library(shiny)
library(bslib)
library(ggplot2)
library(colourpicker) # you might need to install this
library(shinydashboard)
library(DT)
library(readr)
library(dplyr)

ui <- dashboardPage(
  dashboardHeader(title = "Basic Bioinformatics App"),
  dashboardSidebar(sidebarMenu(
    menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
    menuItem("Sample Information", tabName = "Sample_Information", icon = icon("th")),
    menuItem("Counts Matrix Exploration", tabName = "counts", icon = icon("filter")),
    menuItem("Differential Expression", tabName = "dex", icon = icon("sliders")),
    menuItem("Gene Set Enrichment Analysis", tabName = "gsea", icon = icon("dna")),
    menuItem("Individual Gene Expression", tabName = "ige", icon = icon("dna"))
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
      tabItem(tabName = "Sample_Information",
              tabBox(title = "Sample Information",
                # The id lets us use input$tabset1 on the server to find the current tab
                id = "SIhome", height = "250px", width = "80%",
                tabPanel("input", fileInput("metadata_file", "Choose CSV File",
                                            accept = c(".csv")),
                         # Action button to trigger plot and table generation
                         actionButton("generateBtn", "Generate Plot and Table")),
                tabPanel("table", DTOutput("table",width = "80%")),
                tabPanel("summary",  DTOutput("summary",width = "100%")),
                tabPanel("plots",
                         selectInput("hist_variable", "Select Variable for Histogram",
                                     choices = c("age_of_death", "mrna_seq_reads", "pmi", "RIN",
                                                 "age_of_onset", "cag", "duration", "h_v_cortical_score",
                                                 "h_v_striatal_score", "vonsattel_grade")),
                         plotOutput("hist_plot"))
              )
      ),
      tabItem(tabName = "counts",
              tabBox(title = "Counts Matrix",
                     # The id lets us use input$tabset1 on the server to find the current tab
                     id = "Count_home", height = "250px",
                     tabPanel("input", fileInput("counts_file", "Choose CSV File",
                                                 accept = c(".csv")),
                              
                              sliderInput("var_slider", "Minimum percent variance",
                                          min = 0, max = 100, value = 1),
                              sliderInput("non0slider", "Minimum non-zero samples",
                                          min = 0, max = 20, value = 1), #replace 20 with # samples
                              
                              # Action button to trigger plot and table generation
                              actionButton("countsbutton", "Generate Plots and Tables")),
                     tabPanel("table", tableOutput("counts_table")),
                     tabPanel("summary", dataTableOutput("count_summary")),
                     tabPanel("plots", plotOutput("count_plot"))
      )
      ),
      tabItem(tabName = "dex",
              tabBox(title = "Differential Expression",
                     # The id lets us use input$tabset1 on the server to find the current tab
                     id = "diff_home", height = "250px",
                     tabPanel("input", fileInput("dex_file", "Choose CSV File",
                                                 accept = c(".csv")),
                              # Dynamic UI for radio buttons
                              uiOutput("sort_by"),
                              
                              # Action button to trigger plot and table generation
                              actionButton("dexbutton", "Create Table and Plot")),
                     tabPanel("table", tableOutput("dex_table")),
                     tabPanel("volcano plot",  # Dynamic UI for radio buttons
                              uiOutput("variable1_selector"),
                              uiOutput("variable2_selector"),
                              
                              colourInput("color1", "Select Color for Significant", value = "darkgreen"),
                              colourInput("color2", "Select Color for Non-Significant", value = "lightgrey"),
                              
                              sliderInput("slider", "Magnitude of adjusted p-value",
                                          min = -300, max = 1, value = -200),
                              actionButton("countsbutton", "Generate Plot"),
                              plotOutput("volcano"))
                     )
      ),
      tabItem(tabName = "gsea",
              tabBox(title = "Gene Set Enrichment Analysis",
                     
                     id = "gsea_home", height = "250px",
                     tabPanel("input", fileInput("gsea_file", "Choose CSV File",
                                                 accept = c(".csv")),
                              
                              
                              # Action button to trigger plot and table generation
                              actionButton("gseabutton", "Create Table and Plots")),
                     tabPanel("barplot", plotOutput("boxplot")), #add sliders to each tab 
                     tabPanel("table", tableOutput("gseatable")),
                     tabPanel("scatter plot", tableOutput("scatter"))
              )
      ),
      tabItem(tabName = "ige",
              tabBox(title = "Individual Gene Expression",
                     
                     id = "ige_home", height = "250px",
                     tabPanel("input", fileInput("ige_file1", "Choose CSV File for counts matrix",
                                                 accept = c(".csv")),
                     tabPanel("input", fileInput("ige_file2", "Choose CSV File for metadata",
                                                 accept = c(".csv")),
                              textInput("gene_search", "Search for a Gene", placeholder = "gene symbol"), 
                              radioButtons("plot_type", "Plot type:",
                                           c("Bar Plot" = "bar",
                                             "Box Plot" = "box",
                                             "Violin Plot" = "vio",
                                             "Beeswarm Plot" = "bees")),

                              # Action button to trigger plot and table generation
                              actionButton("ige_button", "Create Plot")),
                              plotOutput("distPlot")
                  
              ))
    )
  )
))

server <- function(input, output) {
  load_data <- reactive({
    #for metadata
    req(input$metadata_file)
    data <- read_csv_file(input$metadata_file)
    return(data)
  })
  load_counts <- reactive({
    req(input$counts_file)
    data<-read.csv(input$counts_file$datapath)
    return(data)
  })
  read_csv_file <- function(file) {
    if (is.null(file)) return(NULL)
    read.csv(file$datapath, header = TRUE, stringsAsFactors = FALSE)
  }
  
  # Dynamic UI for radio buttons
  output$variable1_selector <- renderUI({
    choices <- names(load_data())
    radioButtons("variable1", "Choose the column for the x-axis", choices, selected = choices[1])
  })
  output$variable2_selector <- renderUI({
    choices <- names(load_data())
    radioButtons("variable2", "Choose the column for the y-axis", choices, selected = choices[1])
  })

  ige_data1 <- reactive({
    read_csv_file(input$ige_file1)
  })
  
  ige_data2 <- reactive({
    read_csv_file(input$ige_file2)
  })
  
  #front page filler
  set.seed(122)
  histdata <- rnorm(500)
  #front page filler
  output$plot1 <- renderPlot({
    data <- histdata[seq_len(input$slider)]
    hist(data)
  })
  
  volcano_plot <-function(dataf, x_name, y_name, slider, color1, color2) {
    req(x_name, y_name)
    threshold <- 10^slider
    ggplot(dataf, aes(x = !!sym(x_name), y = -log10(!!sym(y_name)), color = padj< threshold)) +
      geom_point(size = 3) +
      scale_color_manual(values = c(color1, color2), breaks = c(TRUE, FALSE), labels = c("significant", "not significant")) +
      labs(title = "Volcano Plot",
           x = x_name,
           y = y_name,
           color = y_name)
  }
  draw_table <- function(dataf) {
    data <- req(input$metadata_file)
    df <- as.data.frame(dataf)
    return(df)
  }
  
  counts_table <- function(dataf, var_slider) {
    req(dataf)
    
    # Calculate the variance for each gene, starting from the second column
    gene_variances <- apply(dataf[, -1], 1, var)
    
    # Find genes with non-zero variance
    nonzero_variance_genes <- gene_variances != 0
    
    # Filter the input tibble to keep only genes with non-zero variance
    filtered_counts <- dataf[nonzero_variance_genes, ]
    return(filtered_counts)
    
  }
  
  get_summary <- function(dataf) {
    data <- dataf
    summary_data <- data %>%
      group_by(diagnosis) %>%
      summarise(
        AvgAgeOfDeath = mean(age_of_death, na.rm = TRUE),
        AvgAgeOfOnset = mean(age_of_onset, na.rm = TRUE),
        UniqueDiagnoses = paste(unique(diagnosis), collapse = ", "),
        DiagnosesCount = as.character(table(diagnosis))
      ) %>%
      as.data.frame()
    return(summary_data)
  }
    
  filtered_data <- reactive({
    gene_search <- input$gene_search
    if (is.null(gene_search) || gene_search == "") {
      return(ige_data1())
    } else {
      ige_data1() %>% filter(Gene == gene_search)
    }
  })
  
  # Reactive expression to create the plot based on user input
  output$distPlot <- renderPlot({
    req(input$ige_button)
    
    # Check if any CSV file is selected
    if (is.null(input$ige_file1) || is.null(input$ige_file2)) {
      return(NULL)
    }
    
    # Check if the selected gene exists in the data
    if (nrow(filtered_data()) == 0) {
      return(NULL)
    }
    
    # Create the plot based on the selected plot type
    plot_type <- input$plot_type
    if (plot_type == "bar") {
      # Code for bar plot
      # Example: barplot(...)
    } else if (plot_type == "box") {
      # Code for box plot
      # Example: boxplot(...)
    } else if (plot_type == "vio") {
      # Code for violin plot
      # Example: violinplot(...)
    } else if (plot_type == "bees") {
      # Code for beeswarm plot
      # Example: beeswarmplot(...)
    }
    
    # Additional code to customize the plot as needed
    
  })
  

#METADATA       
    output$volcano <- renderPlot(volcano_plot(load_data(), input$variable1, input$variable2, input$slider, input$color1, input$color2)) 
    output$table <- renderDT(draw_table(load_data())) 
    
    summary_data <- reactive({
      get_summary(load_data())
    })
    
    output$summary <- renderDT({
      req(summary_data())
      datatable((summary_data()))
    })
    # Render the histogram plot based on user input
    output$hist_plot <- renderPlot({
      req(input$metadata_file, input$hist_variable)
      
      # Read the selected CSV file
      data <- read_csv_file(input$metadata_file)
      
      # Create a histogram based on the selected variable and group by diagnosis
      ggplot(data, aes(x = !!sym(input$hist_variable), fill = diagnosis)) +
        geom_histogram(binwidth = 1, position = "identity", alpha = 0.7) +
        labs(title = paste("Histogram of", input$hist_variable),
             x = input$hist_variable,
             y = "Frequency",
             fill = "Diagnosis")
    })

#COUNTS
    output$counts_table <- renderTable(counts_table(load_counts(), input$var_slider))

}

shinyApp(ui, server)
