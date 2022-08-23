#---- Import data ----

# Import library
library(dplyr)

# Data
data <- read.delim("~/Bureau/phd/tell_seq/PROJECT_GSDIa/tellseq_analysis/similarity_vector_mouse_genome/BLAST/blast_output_vector_mouse.txt", header=FALSE)

# Rename the columns
names(data) <- c("title","chromosome","identity","?1","?2","?3","v_start","v_end","m_start","m_end","expect","score")


#---- Set up the elements of the AAV vector sequence ----

cinq_prim_ITR_inf <- 1
cinq_prim_ITR_sup <- 130
hGPE_inf <- 184
hGPE_sup <- 2883
HBB2_inf <- 2958
HBB2_sup <- 3398
hG6PCwt_inf <- 3412
hG6PCwt_sup <- 4485
bGH_inf <- 4513
bGH_sup <- 4675
trois_prim_ITR_inf <- 4729
trois_prim_ITR_sup <- 4856


#---- Set up conditions to add a new column ----

data <- data %>% 
  mutate(v_elements = case_when(
    v_start >= cinq_prim_ITR_inf & v_end <= cinq_prim_ITR_sup ~ "5'ITR",
    v_start >= hGPE_inf & v_end <= hGPE_sup ~ "promoter (hGPE)",
    v_start >= HBB2_inf & v_end <= HBB2_sup ~ "intron (HBB2)",
    v_start >= hG6PCwt_inf & v_end <= hG6PCwt_sup ~ "gene (hG6PCwt)",
    v_start >= bGH_inf & v_end <= bGH_sup ~ "polyA (bGH)",
    v_start >= trois_prim_ITR_inf & v_end <= trois_prim_ITR_sup ~ "3'ITR",
    TRUE ~ "other")
  )


#---- Export data to a table ----

write.table(data, file = "blast_output_vector_mouse_vector_elements.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)


