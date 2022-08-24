# Libraries ----
library(dplyr)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(forcats)
library(readr)


# Example ----

# Create data
names <- c(rep("A", 80) , rep("B", 50) , rep("C", 70))
value <- c( rnorm(80 , mean=10 , sd=9) , rnorm(50 , mean=2 , sd=15) , rnorm(70 , mean=30 , sd=10) )
data <- data.frame(names,value)

# Basic boxplot
boxplot(data$value ~ data$names , col=terrain.colors(4) )

# Add data points
mylevels <- levels(data$names)
levelProportions <- summary(data$names)/nrow(data)
for(i in 1:length(mylevels)) {
  thislevel <- mylevels[i]
  thisvalues <- data[data$names==thislevel, "value"]
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
  points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.9)) 
  
}
 
# Test ----

# Import data
data <- read_delim("/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/ggplot/r_ggplot_boxplot/input/all.csv", 
                   delim = "\t", escape_double = FALSE, 
                   trim_ws = TRUE)

colnames(data) <- c("name", "value")
data$name <- as.character(data$name)

# Boxplot basic
data %>%
  ggplot( aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Length of chimeric molecules for all the samples") +
  xlab("Samples") +
  ylab("Length") +
  theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
  theme(axis.title.x = element_text(size = rel(1.8), angle = 0)) +
  stat_summary(fun=mean, colour="black", geom="text", size=6, show.legend = FALSE, 
               vjust=-0.9, aes( label=round(..y.., digits=0)))

