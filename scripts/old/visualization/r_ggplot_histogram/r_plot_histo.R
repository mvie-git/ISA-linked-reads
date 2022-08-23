# Libraries ----
library(dplyr)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(forcats)
library(readr)


# Example ----

# Load dataset from github
data <- read.table("https://raw.githubusercontent.com/zonination/perceptions/master/probly.csv", header=TRUE, sep=",")
data <- data %>%
  gather(key="text", value="value") %>%
  mutate(text = gsub("\\.", " ",text)) %>%
  mutate(value = round(as.numeric(value),0))

# plot
p <- data %>%
  mutate(text = fct_reorder(text, value)) %>%
  ggplot( aes(x=value, color=text, fill=text)) +
  geom_histogram(alpha=0.6, binwidth = 5) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  ) +
  xlab("") +
  ylab("Assigned Probability (%)") +
  facet_wrap(~text)

p

# Test ----

# Import data
data <- read_delim("phd/project/GSDIa_louisa/ISA/visualization/r_ggplot_histogram/input/all.csv", 
                  delim = "\t", escape_double = FALSE, 
                  trim_ws = TRUE)

colnames(data) <- c("text", "value")
data$text <- as.character(data$text)

# plot
p <- data %>%
  ggplot( aes(x=value, color=text, fill=text)) +
  geom_histogram(alpha=0.6, binwidth = 5) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 14)
  ) +
  ggtitle("Number of reads per molecule for all the samples") +
  xlab("Number of reads") +
  ylab("Number of molecules") +
  facet_wrap(~text) +
  theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
  theme(axis.title.x = element_text(size = rel(1.8), angle = 0))

p
