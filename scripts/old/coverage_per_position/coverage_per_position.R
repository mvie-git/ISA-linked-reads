#---- Libraries ----

library(plotly)
library(dplyr)


#---- Import data ----

# Coverage for position
data <- read.csv("/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/tellseq_analysis/coverage_per_position/coverage_per_position_all_sample.csv", header=TRUE)
head(data)


# Test: keep only one sample (149)
#data <- data[data$sample == "149", ]

# List of samples
sample <- list(149, 150, 151, 152, 153, 154)



#---- [OPTIONAL]Subset regions of interest ----

# Use subset() function to extract only coverage in 3'ITR region
#data <- subset(data, position >= 4729)


#---- Plot data ----

# X-axis: all the position of the AAV vector sequence
# Y-axis: coverage value

# Write a function to plot each sample (plotly)
## One parameter: sample
plotly_coverage <- function(sample) {
  # Plot the data
  fig <- plot_ly(x = data[data$sample==sample, "position"], y = data[data$sample==sample, "coverage"], 
                 type = 'scatter', mode = 'lines+markers', name = sample, fill = 'tozeroy')
  
  # Add title, legend
  fig <-  fig %>%
    layout(title = 'Coverage per position (AAV vector)', xaxis = list(title = 'Positions'), 
           yaxis = list(title = 'Number of reads'), legend = list(title=list(text='<b> Sample </b>')))
  
  return(fig)
}


# Plot all samples
fig1 <- plotly_coverage(149)
fig2 <- plotly_coverage(150)
fig3 <- plotly_coverage(151)
fig4 <- plotly_coverage(152)
fig5 <- plotly_coverage(153)
fig6 <- plotly_coverage(154)



#---- (Don't use it) Add a vertical line: zero coverage positions ----

# Write a function to add a vertical line corresponds to zero coverage position
## Two parameters: sample and plot
zero_coverage_line <- function(sample, fig) {
  
  # Subset the data to extract for one sample
  data_subset <- data[data$sample==sample, ]
  
  # Extract positions where coverage is equal to 0
  position_zero_coverage <- data_subset[data_subset$coverage == 0, "position"]
  
  # Loop through the set of position to add a line on the plot at each position where there is zero coverage
  for (i in 1:length(position_zero_coverage)) {
    fig <- fig %>%
      add_trace(x = position_zero_coverage[i], type = 'scatter', mode = 'lines',
                line = list(color = "red", opacity=0.01), text = "Zero coverage region", showlegend = F)
  }
  return(fig)
}


# Add lines to zero coverage regions
#fig1 <- zero_coverage_line(149, fig1)
#fig2 <- zero_coverage_line(150, fig2)
#fig3 <- zero_coverage_line(151, fig3)
#fig4 <- zero_coverage_line(152, fig4)
#fig5 <- zero_coverage_line(153, fig5)
#fig6 <- zero_coverage_line(154, fig6)


#---- Find the range where values equal to zero ----

extract_zero_cov_regions <- function(data, sample) {
  
  # Subset the data to extract only one sample
  data_subset <- data[data$sample==sample, ]
  
  # Extract all zero coverage positions
  position_zero_coverage <- data_subset[data_subset$coverage == 0, "position"]
  
  # Find the range where values equal to zero
  positions <- do.call(rbind,by(position_zero_coverage,cumsum(c(0,diff(position_zero_coverage)!=1)),function(g) data.frame(x0=min(g),x1=max(g))))
  
  # Transpose
  pos <- t(positions)
  
  return(pos)
}

# Extract zero coverage regions (to optimize!)
zero_cov_positions_149 <- extract_zero_cov_regions(data, 149)
zero_cov_positions_150 <- extract_zero_cov_regions(data, 150)
zero_cov_positions_151 <- extract_zero_cov_regions(data, 151)
zero_cov_positions_152 <- extract_zero_cov_regions(data, 152)
zero_cov_positions_153 <- extract_zero_cov_regions(data, 153)
zero_cov_positions_154 <- extract_zero_cov_regions(data, 154)

# Make a list of all the zero coverage regions
zero_cov_positions <- list(zero_cov_positions_149, zero_cov_positions_150, zero_cov_positions_151, zero_cov_positions_152, zero_cov_positions_153, zero_cov_positions_154)

# Rename the elements of the list
names(zero_cov_positions) <- c(149, 150, 151, 152, 153, 154)



#---- Add rectangle corresponds to zero coverage regions (to optimize!) ----
#---- Sample 149 ----

# Function to add a rectangle at positions with zero coverage regions (need optimization!!!)
add_rect_zero_cov_regions <- function (fig, data, sample) {
  
  # Initialize plot
  fig <- plotly_coverage(sample)
  
  # Add rectangles
  fig <- layout(fig, title = 'Coverage per position',
                shapes = 
                  list(
                    
                    # First region
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 1], x1 = zero_cov_positions[[sample]]["x1", 1], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),
                    
                    # Second region
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 2], x1 = zero_cov_positions[[sample]]["x1", 2], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y")
                  ))
  return(fig)
}

# Apply function to samples
ncol(zero_cov_positions[["149"]])
fig1 <- add_rect_zero_cov_regions(fig1, data, "149")
fig1


#---- Sample 150 ----

# Function to add a rectangle at positions with zero coverage regions (need optimization!!!)
add_rect_zero_cov_regions <- function (fig, data, sample) {
  
  # Initialize plot
  fig <- plotly_coverage(sample)
  
  # Add rectangles
  fig <- layout(fig, title = 'Coverage per position',
                shapes = 
                  list(
                    
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 1], x1 = zero_cov_positions[[sample]]["x1", 1], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),
                    
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 2], x1 = zero_cov_positions[[sample]]["x1", 2], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),
                    
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 3], x1 = zero_cov_positions[[sample]]["x1", 3], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),
                    
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 4], x1 = zero_cov_positions[[sample]]["x1", 4], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),
                    
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 5], x1 = zero_cov_positions[[sample]]["x1", 5], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),
                    
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 6], x1 = zero_cov_positions[[sample]]["x1", 6], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),
                    
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 7], x1 = zero_cov_positions[[sample]]["x1", 7], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),
                    
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 8], x1 = zero_cov_positions[[sample]]["x1", 8], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),  
                    
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 9], x1 = zero_cov_positions[[sample]]["x1", 9], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),
                    
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 10], x1 = zero_cov_positions[[sample]]["x1", 10], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),
                    
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 11], x1 = zero_cov_positions[[sample]]["x1", 11], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),
                    
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 12], x1 = zero_cov_positions[[sample]]["x1", 12], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),
                    
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 13], x1 = zero_cov_positions[[sample]]["x1", 13], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),
                    
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 14], x1 = zero_cov_positions[[sample]]["x1", 14], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),
                    
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 15], x1 = zero_cov_positions[[sample]]["x1", 15], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),
                    
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 16], x1 = zero_cov_positions[[sample]]["x1", 16], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"), 
                    
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 17], x1 = zero_cov_positions[[sample]]["x1", 17], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"), 
                    
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 18], x1 = zero_cov_positions[[sample]]["x1", 18], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"), 
                    
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 19], x1 = zero_cov_positions[[sample]]["x1", 19], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y")
                    
                  ))
  return(fig)
}

# Apply function to samples
ncol(zero_cov_positions[["150"]])
fig2 <- add_rect_zero_cov_regions(fig2, data, "150")
fig2


#---- Sample 151 ----

# Function to add a rectangle at positions with zero coverage regions (need optimization!!!)
add_rect_zero_cov_regions <- function (fig, data, sample) {
  
  # Initialize plot
  fig <- plotly_coverage(sample)
  
  # Add rectangles
  fig <- layout(fig, title = 'Coverage per position',
                shapes = 
                  list(
                    
                    # First region
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 1], x1 = zero_cov_positions[[sample]]["x1", 1], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),
                    
                    # Second region
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 2], x1 = zero_cov_positions[[sample]]["x1", 2], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),
                    
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 3], x1 = zero_cov_positions[[sample]]["x1", 3], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),                    
                    
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 4], x1 = zero_cov_positions[[sample]]["x1", 4], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),     
                    
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 5], x1 = zero_cov_positions[[sample]]["x1", 5], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y")                      
                    
                  ))
  return(fig)
}

# Apply function to samples
ncol(zero_cov_positions[["151"]])
fig3 <- add_rect_zero_cov_regions(fig1, data, "151")
fig3



#---- Sample 152 ----

# Function to add a rectangle at positions with zero coverage regions (need optimization!!!)
add_rect_zero_cov_regions <- function (fig, data, sample) {
  
  # Initialize plot
  fig <- plotly_coverage(sample)
  
  # Add rectangles
  fig <- layout(fig, title = 'Coverage per position',
                shapes = 
                  list(
                    
                    # First region
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 1], x1 = zero_cov_positions[[sample]]["x1", 1], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),
                    
                    # Second region
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 2], x1 = zero_cov_positions[[sample]]["x1", 2], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y")
                    
                  ))
  return(fig)
}

# Apply function to samples
ncol(zero_cov_positions[["152"]])
fig4 <- add_rect_zero_cov_regions(fig4, data, "152")
fig4



#---- Sample 153 ----

# Function to add a rectangle at positions with zero coverage regions (need optimization!!!)
add_rect_zero_cov_regions <- function (fig, data, sample) {
  
  # Initialize plot
  fig <- plotly_coverage(sample)
  
  # Add rectangles
  fig <- layout(fig, title = 'Coverage per position',
                shapes = 
                  list(
                    
                    # First region
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 1], x1 = zero_cov_positions[[sample]]["x1", 1], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),
                    
                    # Second region
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 2], x1 = zero_cov_positions[[sample]]["x1", 2], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),
                    
                    # First region
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 3], x1 = zero_cov_positions[[sample]]["x1", 3], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),
                    
                    # Second region
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 4], x1 = zero_cov_positions[[sample]]["x1", 4], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y")                    
                    
                  ))
  return(fig)
}

# Apply function to samples
ncol(zero_cov_positions[["153"]])
fig5 <- add_rect_zero_cov_regions(fig5, data, "153")
fig5



#---- Sample 154 ----

# Function to add a rectangle at positions with zero coverage regions (need optimization!!!)
add_rect_zero_cov_regions <- function (fig, data, sample) {
  
  # Initialize plot
  fig <- plotly_coverage(sample)
  
  # Add rectangles
  fig <- layout(fig, title = 'Coverage per position',
                shapes = 
                  list(
                    
                    # First region
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 1], x1 = zero_cov_positions[[sample]]["x1", 1], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),
                    
                    # Second region
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 2], x1 = zero_cov_positions[[sample]]["x1", 2], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),
                    
                    # First region
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 3], x1 = zero_cov_positions[[sample]]["x1", 3], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),
                    
                    # Second region
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 4], x1 = zero_cov_positions[[sample]]["x1", 4], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),  
                    
                    # First region
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 5], x1 = zero_cov_positions[[sample]]["x1", 5], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),
                    
                    # Second region
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 6], x1 = zero_cov_positions[[sample]]["x1", 6], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y"),
                    
                    list(
                      type = "rect",
                      fillcolor = "grey", line = list(color = "grey"), opacity = 0.3,
                      x0 = zero_cov_positions[[sample]]["x0", 7], x1 = zero_cov_positions[[sample]]["x1", 7], xref = "x",
                      y0 = 0, y1 = max(data[data$sample==sample, "coverage"]), yref = "y")
                    
                  ))
  return(fig)
}

# Apply function to samples
ncol(zero_cov_positions[["154"]])
fig6 <- add_rect_zero_cov_regions(fig6, data, "154")
fig6



#---- Subplot ----

# Change order to display non-tumor samples before tumor samples
## Parameters: shareX (link the x axes of subplots) & shareY (link the y axes of subplots)
fig <- subplot(fig1, fig4, fig2, fig3, fig5, fig6, nrows = length(sample), shareX = TRUE, shareY = TRUE) %>% 
  layout(title = list(text = "Coverage per position"),
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
           title = list(text="<b>Samples</b>")
         )) 
fig


