#---- Libraries ----

library(plotly)
library(dplyr)
library(tidyquant)


#---- Import data ----

# GC content (<--- TO CHANGE!)
data <- read.delim("/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/tellseq_analysis/gc_content/NTContent/output/AAV_genome_gc_content.txt", header=FALSE, comment.char="#")

# Rename columns
colnames(data) <- c("positions", "gc_content")

# Subset region of interest (example: 5'ITR)
cinq_ITR <- subset(data, positions <= 130)

# Subset region of interest (example: 5'ITR)
trois_ITR <- subset(data, positions >= 4729)


#---- Plotly ----

# Write a function to plot each region (plotly)
## One parameter: region of interest
plotly_gc_content <- function(subset_region) {
  fig <- plot_ly(subset_region, type = 'scatter', mode = 'lines', fill = 'tozeroy') %>%
    add_trace(x = subset_region$positions, y = subset_region$gc_content, name = 'GC content') %>%
    layout(showlegend = F)
  options(warn = -1)
  
  return(fig)
}

# Plot all regions
fig1 <- plotly_gc_content(cinq_ITR)
fig2 <- plotly_gc_content(trois_ITR)


#---- Subplot ----

## Parameters: shareX (link the x axes of subplots) & shareY (link the y axes of subplots)
fig <- subplot(fig1, fig2, nrows=2) %>% 
  layout(title = list(text = "GC content per position"),
         xaxis = list( 
           title = "Positions (AAV vector reference)",
           zerolinecolor = '#ffff', 
           zerolinewidth = 2, 
           gridcolor = 'ffff'), 
         yaxis = list(
           zerolinecolor = '#ffff', 
           zerolinewidth = 2, 
           gridcolor = 'ffff'),
         legend = list(
           title = list(text="<b>Regions of interest</b>")
         ))
fig

