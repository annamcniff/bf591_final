library(shiny)
library(bslib)
library(ggplot2)
library(colourpicker) # you might need to install this
library(shinydashboard)
library(DT)
library(readr)
library(dplyr)
library(pheatmap)
library(plotly)

read_csv_file <- function(file) {
  if (is.null(file)) return(NULL)
  read.csv(file$datapath, header = TRUE, stringsAsFactors = FALSE)
}

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
                     tabPanel("table", DTOutput("table", width = "80%")),
                     tabPanel("summary",  DTOutput("summary", width = "100%")),
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
                                          min = 0, max = 70, value = 1), #replace 20 with # samples
                              # Action button to trigger plot and table generation
                              actionButton("countsbutton", "Generate Plots and Tables")),
                     tabPanel("summary", DTOutput("counts_summary_table", width = "100%")),
                     tabPanel("plots", plotOutput("count_plot_variance"), plotOutput("count_plot_zeros")),
                     tabPanel("heatmap", plotOutput("clustered_heatmap")),
                     tabPanel("PCA plot", sliderInput("pca_components", "Top N Principal Components to Plot", min = 1, max = 10, value = 2),
                              actionButton("pca_button", "Generate PCA Plot"),
                              plotOutput("pca_plot"))
              )
      ),
      
      tabItem(tabName = "dex",
              tabBox(title = "Differential Expression",
                     id = "diff_home", height = "250px",
                     tabPanel("input", fileInput("dex_file", "Choose CSV File",
                                                 accept = c(".csv")),
                              # Dynamic UI for radio buttons
                              uiOutput("sort_by"),
                              # Gene symbol search input
                              textInput("gene_symbol_search", "Search by Gene Symbol", ""),
                              # Action button to trigger table and volcano plot generation
                              actionButton("dexbutton", "Create Table and Plot")),
                     tabPanel("table", DTOutput("dex_table")),
                     tabPanel("volcano plot",
                              uiOutput("variable2_selector"),
                              colourInput("color1", "Select Color for Significant", value = "hotpink"),
                              colourInput("color2", "Select Color for Non-Significant", value = "lightgrey"),
                              sliderInput("slider", "Magnitude of adjusted p-value",
                                          min = -300, max = 1, value = -200),
                              actionButton("dexplotbutton", "Generate Plot"),
                              plotOutput("dex_volcano"))
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
                                       # Input control to choose a categorical field
                                       selectInput("categorical_field", "Choose Categorical Field", ""),
                                       textInput("gene_search", "Search for a Gene", placeholder = "gene symbol"), 
                                       radioButtons("plot_type", "Plot type:",
                                                    c("Bar Plot" = "bar",
                                                      "Box Plot" = "box",
                                                      "Violin Plot" = "vio",
                                                      "Beeswarm Plot" = "bees")),
                                       # Action button to trigger plot and table generation
                                       actionButton("ige_button", "Create Plot")),
                              plotOutput("distPlot")
                     )
              )
      )
    ))
)

server <- function(input, output) {
  options(shiny.maxRequestSize = 30*1024^2)
  load_data <- reactive({
    #for metadata
    req(input$metadata_file)
    data <- read_csv_file(input$metadata_file)
    return(data)
  })
  
  #COUNTS
  load_counts <- reactive({
    req(input$counts_file)
    df <- readr::read_csv(input$counts_file$datapath)
    
    # Make 'gene' column character type
    df$gene <- as.character(df$gene)
    
    # Convert all columns except 'gene' to numeric
    numeric_cols <- sapply(df[-1], as.numeric)
    df <- cbind(df["gene"], numeric_cols)
    df <- as.data.frame(df)
    
    #print(dim(df))  # Print dimensions of the resulting data frame
    return(df)
  })
  
  
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
#COUNTS  
  counts_table <- function(dataf, var_slider,non0slider) {
    req(dataf)
    #print(dim(dataf))
    # Calculate the variance for each gene, starting from the second column
    gene_variances <- apply(dataf[, -1], 1, var)
    
    # Find genes with non-zero variance and at least non0slider non-zero samples
    nonzero_variance_genes <- gene_variances != 0
    non_zero_samples <- apply(dataf[, -1] > 0, 1, sum)
    non_zero_samples_pass_filter <- non_zero_samples >= non0slider
    
    # Filter the input tibble to keep only genes with non-zero variance and passing non-zero samples filter
    filtered_counts <- dataf[nonzero_variance_genes & non_zero_samples_pass_filter, ]
    
    # Apply additional filter based on minimum percent variance
    filtered_counts <- filtered_counts[, c(TRUE, apply(filtered_counts[, -1], 2, function(x) var(x) >= (var_slider/100)))]
    
    return(filtered_counts)
  }
  # Function to calculate the median count and number of zeros
  calculate_median_and_zeros <- function(dataf) {
    median_count <- apply(dataf[, -1], 1, median, na.rm = TRUE)
    num_zeros <- apply(dataf[, -1] == 0, 1, sum)
    return(data.frame(Gene = dataf$gene, MedianCount = median_count, NumZeros = num_zeros))
  }
  
  # Add this reactive expression
  heatmap_data <- reactive({
    req(counts_data())
    
    # Log-transform counts for visualization
    log_counts <- log1p(counts_data()[, -1])
    
    # Prepare data for heatmap, you may need to customize this based on your data structure
    heatmap_data <- t(log_counts)
    
    return(heatmap_data)
  })
  
  perform_pca <- function(counts_data) {
    # Extract gene names and counts
    genes <- counts_data$gene
    counts <- as.matrix(counts_data[, -1])  # Exclude the 'gene' column
    
    # Perform PCA
    pca_result <- prcomp(counts, scale. = TRUE, center = TRUE)
    
    # Return a list with PCA results
    return(list(
      x = pca_result$x,      # PCA scores
      sdev = pca_result$sdev  # Standard deviations of principal components
    ))
  }
  
  
#METADATA
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
    #output$counts_table <- renderTable(counts_table(load_counts(), input$var_slider))
    counts_data <- reactive({
      counts_table(load_counts(), input$var_slider,input$non0slider)
    })
    # Perform PCA (replace with your actual PCA function)
    pca_results <- reactive({
      perform_pca(load_counts())
    })
    
    output$counts_summary_table <- renderDT({
      req(counts_data())
      # Calculate summary statistics
      total_genes <- nrow(counts_data())
      passed_filter <- nrow(counts_data())
      not_passed_filter <- total_genes - passed_filter
      total_samples <- ncol(counts_data()) - 1  # excluding the 'gene' column
      
      # Create a summary data frame
      summary_data <- data.frame(
        Metric = c("Total Number of Genes", "Genes Passing Filter", "Genes Not Passing Filter", "Total Number of Samples"),
        Value = c(total_genes, passed_filter, not_passed_filter, total_samples),
        Percentage = c(NA, passed_filter / total_genes * 100, not_passed_filter / total_genes * 100, NA)
      )
      
      # Render the summary table
      datatable(summary_data, options = list(pageLength = 10))
    })
    
    # Plot for median count vs variance
    output$count_plot_variance <- renderPlot({
      req(counts_data())
      
      # Calculate median and zeros
      median_data <- calculate_median_and_zeros(counts_data())
      
      # Create a plot for Median Count vs Variance
      ggplot(median_data, aes(x = MedianCount, y = apply(counts_data()[, -1], 1, var), color = NumZeros > 0)) +
        geom_point() +
        scale_color_manual(values = c("FALSE" = "lightblue", "TRUE" = "darkblue")) +
        labs(title = "Median Count vs Variance",
             x = "Median Count",
             y = "Variance",
             color = "Passes Filters")
    })
    # Plot for median count vs number of zeros
    output$count_plot_zeros <- renderPlot({
      req(counts_data())
      
      # Calculate median and zeros
      median_data <- calculate_median_and_zeros(counts_data())
      
      # Create a plot for Median Count vs Number of Zeros
      ggplot(median_data, aes(x = log(MedianCount), y = NumZeros, color = NumZeros > 0)) +
        geom_point() +
        scale_color_manual(values = c("FALSE" = "lightblue", "TRUE" = "darkblue")) +
        labs(title = "Median Count vs Number of Zeros",
             x = "Log of Median Count",
             y = "Number of Zeros",
             color = "Passes Filters")
    })
    output$clustered_heatmap <- renderPlot({
      req(heatmap_data())
      
      # Create a clustered heatmap
      pheatmap(
        heatmap_data(),
        color = colorRampPalette(c("yellow", "blue"))(100),
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        main = "Clustered Heatmap",
        show_colnames = FALSE,
        fontsize_row = 6,
        fontsize_col = 8,
        legend = TRUE
      )
    })
    pca_results <- reactive({
      # Generate a subset of the counts data
      test_counts_data <- subset(load_counts(), select = c("gene", "C_0002", "C_0003", "C_0004"))
      
      # Perform PCA on the subset
      perform_pca(test_counts_data)
    })
    
    # Render the PCA scatter plot
    output$pca_plot <- renderPlotly({
      req(pca_results())
      
      # Extract the selected number of principal components
      top_components <- pca_results()$x[, 1:input$pca_components]
      
      # Use plot_ly for an interactive plot
      plot_ly(x = ~top_components[, 1], y = ~top_components[, 2], type = "scatter",
              mode = "markers", marker = list(color = "darkblue", size = 10),
              text = ~paste("PC1: ", round(top_components[, 1], 2), "<br>PC2: ", round(top_components[, 2], 2))) %>%
        layout(title = "PCA Scatter Plot",
               xaxis = list(title = paste("PC1 (", round(pca_results()$sdev[1]^2 / sum(pca_results()$sdev^2) * 100, 2), "%)")),
               yaxis = list(title = paste("PC2 (", round(pca_results()$sdev[2]^2 / sum(pca_results()$sdev^2) * 100, 2), "%)")))
    })
    
#DEX 
    dex_data <- reactive({
      req(input$dex_file)
      read_csv_file(input$dex_file)
    })
    
    # Dynamic UI for sorting options
    output$sort_by <- renderUI({
      req(dex_data())
      choices <- names(dex_data())
      radioButtons("sort_by_column", "Sort by:", choices, selected = choices[1])
    })
    
    # Filtered data based on gene symbol search
    filtered_dex_data <- reactive({
      gene_search <- input$gene_symbol_search
      if (is.null(gene_search) || gene_search == "") {
        return(dex_data())
      } else {
        dex_data() %>% filter(symbol == gene_search)
      }
    })
    
    # Render the table based on filtered data
    output$dex_table <- renderDT({
      req(dex_data(), input$sort_by_column)
      datatable(filtered_dex_data() %>% arrange(!!sym(input$sort_by_column)))
    })
    output$dex_volcano <- renderPlot({
      req(input$dexbutton, input$dex_file)
      
      # Read the selected CSV file
      dataf <- read_csv_file(input$dex_file)
      
      # Create the Volcano Plot based on "padj" and "log2FoldChange"
      ggplot(dataf, aes(x = log2FoldChange, y = -log10(padj), color = padj < 10^input$slider)) +
        geom_point(size = 3) +
        scale_color_manual(values = c(input$color1, input$color2), breaks = c(TRUE, FALSE),
                           labels = c("significant", "not significant")) +
        labs(title = "Volcano Plot",
             x = "log2FoldChange",
             y = "padj",
             color = "padj")
    })
#IGE 
    observe({
      req(load_data())  
      # Extract column names from the loaded sample information matrix
      col_names <- colnames(load_data())
      
      # Update the choices for the selectInput
      updateSelectInput(session, "categorical_field", choices = col_names, selected = col_names[1])
    })
    
    output$distPlot <- renderPlot({
      req(input$ige_button, input$ige_file1, input$ige_file2, input$gene_search, input$plot_type)
      
      # Read the counts matrix file
      counts_data <- read_csv_file(input$ige_file1)
      
      # Read the sample information file
      metadata <- read_csv_file(input$ige_file2)
      
      # Filter counts matrix for the selected gene
      selected_gene <- input$gene_search
      filtered_counts <- counts_data %>%
        filter(gene == selected_gene)
      
      if (nrow(filtered_counts) == 0) {
        # Gene not found, return NULL or show an error message
        return(NULL)
      }
      # Extract relevant columns
      sample_col <- input$sample_variable  # replace with your actual variable selection
      plot_type <- input$plot_type
      
      # Perform the requested plot type
      if (plot_type == "bar") {
        # Bar Plot
        ggplot(filtered_counts, aes(x = !!sym(sample_col), y = counts, fill = !!sym(sample_col))) +
          geom_bar(stat = "identity") +
          labs(title = paste("Bar Plot for Gene: ", selected_gene),
               x = sample_col,
               y = "Count",
               fill = sample_col)
      } else if (plot_type == "box") {
        # Box Plot
        ggplot(filtered_counts, aes(x = !!sym(sample_col), y = counts, fill = !!sym(sample_col))) +
          geom_boxplot() +
          labs(title = paste("Box Plot for Gene: ", selected_gene),
               x = sample_col,
               y = "Count",
               fill = sample_col)
      } else if (plot_type == "vio") {
        # Violin Plot
        ggplot(filtered_counts, aes(x = !!sym(sample_col), y = counts, fill = !!sym(sample_col))) +
          geom_violin() +
          labs(title = paste("Violin Plot for Gene: ", selected_gene),
               x = sample_col,
               y = "Count",
               fill = sample_col)
      } else if (plot_type == "bees") {
        # Beeswarm Plot
        ggplot(filtered_counts, aes(x = !!sym(sample_col), y = counts, color = !!sym(sample_col))) +
          geom_beeswarm() +
          labs(title = paste("Beeswarm Plot for Gene: ", selected_gene),
               x = sample_col,
               y = "Count",
               color = sample_col)
      }
    })
    # Additional controls for the Individual Gene Expression tab
    output$sample_variable_selector <- renderUI({
      req(input$ige_file2)
      choices <- names(read_csv_file(input$ige_file2))
      radioButtons("sample_variable", "Choose sample variable", choices, selected = choices[1])
    })
    
    # Updated input section in the UI
    tabPanel(
      "input",
      fileInput("ige_file1", "Choose CSV File for counts matrix", accept = c(".csv")),
      fileInput("ige_file2", "Choose CSV File for metadata", accept = c(".csv")),
      textInput("gene_search", "Search for a Gene", placeholder = "gene symbol"),
      radioButtons("sample_variable_selector", "Choose sample variable", choices = ""),
      radioButtons("plot_type", "Plot type:",
                   c("Bar Plot" = "bar",
                     "Box Plot" = "box",
                     "Violin Plot" = "vio",
                     "Beeswarm Plot" = "bees")),
      actionButton("ige_button", "Create Plot")
    )
    
    # Additional controls for the Individual Gene Expression tab
    output$sample_variable_selector <- renderUI({
      req(input$ige_file2)
      choices <- names(read_csv_file(input$ige_file2))
      radioButtons("sample_variable", "Choose sample variable", choices, selected = choices[1])
    })
}


shinyApp(ui, server)

