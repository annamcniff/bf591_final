# bf591_final
Final Project for BF591 course 

## Project Components: 
### Sample Information Exploration 
Tab with a summary of the table that includes a summary of the type and values in each column, e.g.: Number of rows: X Number of columns: Y
Tab with a data table displaying the sample information, with sortable columns
Tab with histograms of different variables sourced from input csv, with data grouped by diagnosis. 

### Counts Matrix Exploration
Takes as input a normalized counts matrix in CSV format
Input controls that filter out genes based on their statistical properties 
Returns: 
a table summarizing the effect of the filtering, including: number of samples, total number of genes, number and % of genes passing current filter and number and % of genes not passing current filter. 
Tab with diagnostic scatter plots, where genes passing filters are marked in a darker color, and genes filtered out are lighter:
Tab with a clustered heatmap of counts remaining after filtering
Tab with a scatter plot of principal component analysis projections

### Differential Expression
Differential expression identifies which genes, if any, are implicated in a specific biological comparison. This component allows the user to load and explore a differential expression dataset. 
This functionality takes as input the results of a differential expression analysis in CSV format and returns a Tab with sortable table displaying differential expression results with a gene name search functionality to filter rows of the table. 
To visualize this data there is a tab which generates a Volcano plot with user selection of variable to be plotted as well as significance threshold. 

### Individual Gene Expression 
Visualizing individual gene counts is useful for examining or verifying patterns identified by differential expression analysis. There are many different ways of visualizing counts for a single gene. This app allows counts from an arbitrary gene to be selected and visualized broken out by a desired sample information variable.
This functionality takes as it's input a normalized counts matrix in CSV format as well as a sample information matrix in CSV format. 
Input control that allows the user to choose one of the categorical fields found in the sample information matrix file as well as a gene found in the counts matrix. This component also allows the user to choose the type of visualization generated from a bar plot, boxplot, violin plot, or beeswarm plot. 
