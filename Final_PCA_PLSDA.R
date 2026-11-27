################### PCA analysis #####################

library(tidyverse)
library(mixOmics)
library(factoextra)
library(FactoMineR)
library(readxl)

rm(list = ls())

# YOU HAVE TO USE TRANSPOSED DATA (SAMPLES = ROWS, FEATURES = COLUMNS)

Brain_metabolites_12m <- read_excel("C:/Users/veera/Documents/Thesis project/Muscle/Muscle Metabolomics/Muscle_12m_skyline_transposed.xlsx")


##let's change the sample names into characters to make life easier in the next steps

Brain_metabolites_12m$Identification <- as.character(Brain_metabolites_12m$Identification)


### median normalization ###

# calculate median for every sample and add it as a new column

Brain_metabolites_12m$Median <- apply(Brain_metabolites_12m[, 2:186], 1, median, na.rm = TRUE)

# calculate the maximum median
maxmedian <- max(Brain_metabolites_12m$Median)

# divide each data point by sample's median and multiply by maximum median value

divided <- Brain_metabolites_12m %>% 
  rowwise() %>%
  mutate(across(2:186, ~.x / Median)) %>%
  ungroup()


median_normalised <- divided %>% mutate(across(2:186, ~ . * maxmedian))

#remove the Median column and rename for later
filtered <- median_normalised[, 1:186]


### Adding genotype and sex columns

#for genotype
genotype <- as.character(c("WT",	"WT",	"ALS",	"WT",	"WT",	"WT",	"WT",	"ALS",	"ALS",	"ALS",	"WT",	"ALS",	"ALS",	"WT",	"ALS",	"ALS",	"ALS",	"ALS",	"WT",	"WT"))

sample_names <- Brain_metabolites_12m$Identification

named_genotypes <- setNames(genotype, sample_names)

filtered$Genotype <- named_genotypes[filtered$Identification]

#for sex
sex <- as.character(c("F",	"F",	"M",	"M",	"M",	"F",	"M",	"F",	"M",	"F",	"M",	"M",	"F",	"M",	"M",	"M",	"F",	"F",	"F",	"F"))

named_sex <- setNames(sex, sample_names)

filtered$Sex <- named_sex[filtered$Identification]


##### PCA analysis #####

## if you use transformed data, then you don't do log2 and scale = FALSE in pca
pca_full <- pca(X = filtered[, 2:186], # the input will exclude the samples column and the label column
                ncomp = 5, # Number of components is set to 5
                scale = FALSE)


## if you use non-transformed data, then you do log2 and scale = TRUE in pca

# Perform PCA on the full dataset
pca_full <- pca(X = log2(filtered[, 2:186]), # the input will exclude the samples column and the label columns
                ncomp = 5, # Number of components is set to 5
                scale = F)


# Plotting PCA

pca_plo <- plotIndiv(pca_full, # First argument is the PCA list object
                     group = filtered$Genotype, # The grouping variable in the pdac_full dataframe
                     ind.names = F, # Labels for each datapoint will not be included
                     rep.space = "XY-variate",
                     comp = c(1, 2), # Plotting components 2 and 3 on the X and Y axis
                     style = "ggplot2", # Produces a ggplot-like plot
                     title = "PCA on Muscle metabolites 12m",
                     col.per.group = c("royalblue", "orange", "red2"), # Colour assignment for the grouping variable
                     ellipse = TRUE, # Draws 95% confidence ellipse around datapoints according to group
                     centroid = FALSE, # Includes an additional datapoint at the centroid of each ellipse
                     legend = TRUE, # Figure legend
                     legend.title = "Group",
                     abline = T, 
                     alpha = 0)

# Extract PCA loadings data from PCA object

pca_loadings <- pca_full$loadings$X %>% # extracts the loadings information
  as.data.frame() %>% # Converts to dataframe
  mutate(features = rownames(.))%>% # Creates a feature column referring to the lipid species
  dplyr::select(PC1, features) %>% 
  {bind_rows(
    slice_max(., order_by = PC1, n = 10),
    slice_min(., order_by = PC1, n = 10)
  )} %>%
  arrange(desc(PC1))


### this is the original code that only extracts the largest positive values
#arrange(-PC1) %>% # Arranging the dataframe by loading values (descending order)
#  slice_head(n = 15) # Selecting the top 20 features - this reduces over plotting but can be changed 

# Plot loadings using ggplot2
ggplot(data = pca_loadings,
       mapping = aes(x = PC1,
                     y = features)) +
  geom_col(alpha = 0.5,
           fill = "firebrick") +
  labs(x = "",
       y = "Metabolites",
       title = "Loadings on Component 1 for Muscle 12m") +
  xlim(-0.25, 0.25) +
  theme_bw() +
  theme(text = element_text(face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)))


##### Sex separated analysis #####

# separating by sex

female_data <- subset(filtered, Sex == "F")

pca_female <- pca(X = log2(female_data[, 2:186]), # the input will exclude the samples column and the label column
                  ncomp = 5, # Number of components is set to 5
                  scale = FALSE)

plot_female <- plotIndiv(pca_female, # First argument is the PCA list object
                         group = female_data$Genotype, # The grouping variable in the pdac_full dataframe
                         ind.names = F, # Labels for each datapoint will not be included
                         rep.space = "XY-variate",
                         comp = c(1, 2), # Plotting components 2 and 3 on the X and Y axis
                         style = "ggplot2", # Produces a ggplot-like plot
                         title = "PCA on female Muscle 12m",
                         col.per.group = c("royalblue", "orange", "red2"), # Colour assignment for the grouping variable
                         ellipse = TRUE, # Draws 95% confidence ellipse around datapoints according to group
                         centroid = FALSE, # Includes an additional datapoint at the centroid of each ellipse
                         legend = TRUE, # Figure legend
                         legend.title = "Group",
                         abline = T, 
                         alpha = 0)


# Extract PCA loadings data from PCA object

pca_loadings_female <- pca_female$loadings$X %>% # extracts the loadings information
  as.data.frame() %>% # Converts to dataframe
  mutate(features = rownames(.))%>% # Creates a feature column referring to the lipid species
  dplyr::select(PC1, features) %>% 
  {bind_rows(
    slice_max(., order_by = PC1, n = 10),
    slice_min(., order_by = PC1, n = 10)
  )} %>%
  arrange(desc(PC1))


### this is the original code that only extracts the largest positive values
#arrange(-PC1) %>% # Arranging the dataframe by loading values (descending order)
#  slice_head(n = 15) # Selecting the top 20 features - this reduces over plotting but can be changed 

# Plot loadings using ggplot2
ggplot(data = pca_loadings_female,
       mapping = aes(x = PC1,
                     y = features)) +
  geom_col(alpha = 0.5,
           fill = "firebrick") +
  labs(x = "",
       y = "Metabolites",
       title = "Loadings on Component 1 for Muscle 12m") +
  xlim(-0.25, 0.25) +
  theme_bw() +
  theme(text = element_text(face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)))


# and then males

male_data <- subset(filtered, Sex == "M")


pca_male <- pca(X = log2(male_data[, 2:186]), # the input will exclude the samples column and the label column
                ncomp = 5, # Number of components is set to 5
                scale = FALSE)

plot_male <- plotIndiv(pca_male, # First argument is the PCA list object
                       group = male_data$Genotype, # The grouping variable in the pdac_full dataframe
                       ind.names = F, # Labels for each datapoint will not be included
                       rep.space = "XY-variate",
                       comp = c(1, 2), # Plotting components 2 and 3 on the X and Y axis
                       style = "ggplot2", # Produces a ggplot-like plot
                       title = "PCA on male Muscle 12m",
                       col.per.group = c("royalblue", "orange", "red2"), # Colour assignment for the grouping variable
                       ellipse = TRUE, # Draws 95% confidence ellipse around datapoints according to group
                       centroid = FALSE, # Includes an additional datapoint at the centroid of each ellipse
                       legend = TRUE, # Figure legend
                       legend.title = "Group",
                       abline = T, 
                       alpha = 0)

# Extract PCA loadings data from PCA object

pca_loadings_male <- pca_male$loadings$X %>% # extracts the loadings information
  as.data.frame() %>% # Converts to dataframe
  mutate(features = rownames(.))%>% # Creates a feature column referring to the lipid species
  dplyr::select(PC1, features) %>% 
  {bind_rows(
    slice_max(., order_by = PC1, n = 10),
    slice_min(., order_by = PC1, n = 10)
  )} %>%
  arrange(desc(PC1))


### this is the original code that only extracts the largest positive values
#arrange(-PC1) %>% # Arranging the dataframe by loading values (descending order)
#  slice_head(n = 15) # Selecting the top 20 features - this reduces over plotting but can be changed 

# Plot loadings using ggplot2
ggplot(data = pca_loadings_male,
       mapping = aes(x = PC1,
                     y = features)) +
  geom_col(alpha = 0.5,
           fill = "firebrick") +
  labs(x = "",
       y = "Metabolites",
       title = "Loadings on Component 1 for Muscle 12m") +
  xlim(-0.25, 0.25) +
  theme_bw() +
  theme(text = element_text(face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)))



################## PLS (DA) analysis ########################

library(tidymodels)
library(caret)
library(tidyverse)
detach("package:mixOmics", unload = TRUE)


# adding log2 transformation

mat <- as.matrix(filtered[, c(-1, -187, -188)])

transformed <- ifelse(mat == 0 | is.na(mat), 0, log2(mat))
# this one is not the right way to make this work - we would have to find a way to remove zeros without making them NAs
# as the PLS-DA does not work if you have NAs

# add the feature names back
normalised <- data.frame(
  Identification = filtered$Identification,
  transformed
)

# now we can create a separate data frame for the PLS analysis
X <- 
  normalised %>% 
  dplyr::select(-Identification)


#adding back genotype and sex columns

#for genotype
genotype <- as.character(c("WT",	"WT",	"ALS",	"WT",	"WT",	"WT",	"WT",	"ALS",	"ALS",	"ALS",	"WT",	"ALS",	"ALS",	"WT",	"ALS",	"ALS",	"ALS",	"ALS",	"WT",	"WT"))

sample_names <- normalised$Identification

named_genotypes <- setNames(genotype, sample_names)

normalised$Genotype <- named_genotypes[normalised$Identification]

#for sex
sex <- as.character(c("F",	"F",	"M",	"M",	"M",	"F",	"M",	"F",	"M",	"F",	"M",	"M",	"F",	"M",	"M",	"M",	"F",	"F",	"F",	"F"))

named_sex <- setNames(sex, sample_names)

normalised$Sex <- named_sex[normalised$Identification]


#### if you're using transformed data ###

#X <-
#  Brain_metabolites_12m %>% 
#  dplyr::select(-Identification, -Median)

Y <- filtered$Genotype # could also just use the genotype vector


### Analysis ###

#you have to find a way to remove NAs for this to work!

# Fitting our model over different tuning parameters:
PLSDA_brain <- train(x = X,            # Our predictors - lipids. Can be matrix or data frame. Samples must be in rows, variables in columns with column names!
                     y = Y, 
                     tuneLength = 10, # Outcome variable: N or T. Must be numeric or factor vector.
                     method = "pls")   # We want PLS algorithm.


## checking the data
PLSDA_brain$results

plot(PLSDA_brain, main = 'PLS-DA model tuning', xlab = 'Number of components')

# Best tune = number of components in the final PLS-DA equation:
PLSDA_brain$bestTune



# Extracting PLS scores:
PLSDA_brain$finalModel$scores

# We turn them into a data frame so the scores plot can be prepared:
scores <- as.matrix.data.frame(PLSDA_brain$finalModel$scores)
scores <- as.data.frame(scores)

# We will add the column with biological groups - so we can color N and T:
scores$Genotype <- filtered$Genotype

# We will also compute the variance explained by the first two components.
# To access the variance explained by first two components run:
PLSDA_brain$finalModel$Xvar

# Extracting total variance explained:
PLSDA_brain$finalModel$Xtotvar


# Variance explained by Component 1 = 4376, Component 2 = 3587, Total variance = 18034
# Relative variance explained by Component 1 
(1721.8095/7531.801)*100 # in [%]
# Relative variance explained by Component 2
(1546.3284/7531.801)*100 # in [%]
# Relative variance explained by Component 3
#(1316.955/8569)*100 # in [%]


#### visualising it ###
# Scores plot from ggplot2:
ggplot(scores, aes(x = V1, y= V2, color = Genotype)) + 
  geom_point(size = 3)+
  scale_color_manual(values = c('red2', 'royalblue'))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 12, 
                                   color = 'black'),
        axis.text.y = element_text(size = 12, 
                                   color = 'black'),
        legend.text = element_text(size = 12,
                                   color = 'black'),
        legend.title = element_text(size = 12, 
                                    color = 'black'))+
  labs(x = "Component 1 (var explained = 22.9%)",
       y = "Component 2 (var explained = 20.5%)"
  )+
  ggtitle('Muscle 12m')
