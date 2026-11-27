##################### Volcano plots and median normalisation for 12m ##############################

rm(list = ls())  # clear environment

library(EnhancedVolcano)
library(tidyverse)
library(rstatix)
library(readxl)
library(mixOmics)
library(factoextra)
library(FactoMineR)
library(writexl)


### YOU HAVE TO USE TRANSPOSED WIDE FORM DATA FOR NORMALISATION

Brain_metabolites_12m <- read_excel("C:/Users/veera/Documents/Thesis project/Liver/Liver Metabolomics/Liver_12m_skyline_transposed.xlsx")
#glimpse(Brain_metabolites_12m)


##let's change the sample names into characters to make life easier in the next steps

Brain_metabolites_12m$Identification <- as.character(Brain_metabolites_12m$Identification)


### median normalization ###

# calculate median for every sample and add it as a new column

Brain_metabolites_12m$Median <- apply(Brain_metabolites_12m[, 2:155], 1, median, na.rm = TRUE)

# calculate the maximum median
maxmedian <- max(Brain_metabolites_12m$Median)

# divide each data point by sample's median and multiply by maximum median value

divided <- Brain_metabolites_12m %>% 
  rowwise() %>%
  mutate(across(2:155, ~.x / Median)) %>%
  ungroup()


median_normalised <- divided %>% mutate(across(2:155, ~ . * maxmedian))

#remove the Median column and rename for later
filtered <- median_normalised[, 1:155]

# for lipids use this
filtered <- Brain_metabolites_12m

# save the normalised values

#write_xlsx(filtered, path = "C:/Users/veera/Documents/Thesis project/Brain/Lipidomics/Lipids_brain_12m_normalised.xlsx")



### t.tests and volcanos ###

# Step 1 - Create 'data.long' tibble:
data.long <- 
  filtered %>%
  pivot_longer(cols = c(-Identification),
               names_to = "Metabolite",
               values_to = "Concentrations")



## Creating metadata columns

#for genotype
genotype <- as.character(c("WT",	"WT",	"ALS",	"WT",	"WT",	"WT",	"WT",	"ALS",	"ALS",	"ALS",	"WT",	"ALS",	"ALS",	"WT",	"ALS",	"ALS",	"ALS",	"ALS",	"WT",	"WT"))

sample_names <- Brain_metabolites_12m$Identification

named_genotypes <- setNames(genotype, sample_names)

data.long$Genotype <- named_genotypes[data.long$Identification]

#for sex
sex <- as.character(c("F",	"F",	"M",	"M",	"M",	"F",	"M",	"F",	"M",	"F",	"M",	"M",	"F",	"M",	"M",	"M",	"F",	"F",	"F",	"F"))

named_sex <- setNames(sex, sample_names)

data.long$Sex <- named_sex[data.long$Identification]


# Step 2 - Compute Welch's t-test and adjust p-values (rstatix):
t.test <- 
  data.long %>%
  group_by(Metabolite) %>%
  rstatix::t_test(Concentrations ~ Genotype,
                  var.equal = T) %>%
  adjust_pvalue(method = "BH")


# Step 3 - Compute mean concentrations for all lipids within healthy controls (N):
metabolites.WT <- 
  data.long %>%
  group_by(Metabolite, Genotype) %>%
  summarise(mean = mean(Concentrations, na.rm = T)) %>%
  arrange(Genotype) %>%
  filter(Genotype == "WT")

# Step 4 - Compute mean concentrations for all lipids within patients with PDAC (T):
metabolites.ALS <- 
  data.long %>%
  group_by(Metabolite, Genotype) %>%
  summarise(mean = mean(Concentrations, na.rm = T)) %>%
  arrange(Genotype) %>%
  filter(Genotype == "ALS")

# Step 5 - Compute fold-change and store it in the "Welch.test" as `Fold change`:
t.test$FC <- metabolites.ALS$mean/metabolites.WT$mean

# Step 6 - Using mutate(), compute log2 of fold change:
t.test <-
  t.test %>%
  mutate(log.2.FC = log2(FC))

# saving the results
#write_xlsx(t.test[, c("Metabolite", "p", "log.2.FC")], path = "C:/Users/veera/Documents/Thesis project/Liver/Lipidomics/Lipids_liver_12m_results.xlsx")

# Step 8 - Create a default volcano plot from the EnhancedVolcano:
#EnhancedVolcano(t.test,
#                lab = t.test$Metabolite,
#                x = "log.2.FC", 
#                y = "p",
#                xlab = "Fold change", 
#                ylab = "p-value(-log10)",
#                title = "Muscle lipids 12m",
#                pCutoff = 0.05,
#                FCcutoff = 0,
#                labSize = 3,
#                drawConnectors = TRUE,
#                arrowheads = F,
#                max.overlaps = 10,
#                col = c("grey","grey", "royalblue", "red4"),
#                ylim = c(0, 4),
#                xlim = c(-2, 2))



# Separating up- and downregulated groups

t.test$new.condition <- ifelse(t.test$p < 0.05 & t.test$log.2.FC > 0, 'Significant upregulation', 
                                   ifelse(t.test$p < 0.05 & t.test$log.2.FC < 0, 'Significant downregulation', 'Insignificant'))


# adding names to colours
colours <- ifelse(t.test$new.condition == "Significant upregulation", "red4", 
                  ifelse(t.test$new.condition == "Significant downregulation", "royalblue", "grey"))


names(colours)[colours == 'grey'] <- 'Insignificant'
names(colours)[colours == 'red4'] <- 'Significant upregulation'
names(colours)[colours == 'royalblue'] <- 'Significant downregulation'

EnhancedVolcano(t.test,
                lab = t.test$Metabolite,
                x = "log.2.FC", 
                y = "p",
                xlab = "log2(Fold change)", 
                ylab = "-log10(p-value)",
                title = "Liver metabolites 12m",
                pCutoff = 0.05,
                FCcutoff = 0,
                labSize = 4,
                drawConnectors = TRUE,
                arrowheads = F,
                max.overlaps = 20,
                colCustom = colours,
                ylim = c(0, 6),
                xlim = c(-3, 3)) 




# male only

male_data <- subset(data.long, Sex == "M")

# Compute Welch's t-test and adjust p-values (rstatix):
t.test.male <- 
  male_data %>%
  group_by(Metabolite) %>%
  rstatix::t_test(Concentrations ~ Genotype, var.equal = T) %>%
  adjust_pvalue(method = "BH")


# Compute mean concentrations for all lipids within healthy controls (N):
male.WT <- 
  male_data %>%
  group_by(Metabolite, Genotype) %>%
  summarise(mean = mean(Concentrations, na.rm = T)) %>%
  arrange(Genotype) %>%
  filter(Genotype == "WT")

# Compute mean concentrations for all lipids within patients with PDAC (T):
male.ALS <- 
  male_data %>%
  group_by(Metabolite, Genotype) %>%
  summarise(mean = mean(Concentrations, na.rm = TRUE)) %>%
  arrange(Genotype) %>%
  filter(Genotype == "ALS")

# Compute fold-change and store it in the "Welch.test" as `Fold change`:
t.test.male$FC <- male.ALS$mean/male.WT$mean

# Using mutate(), compute log2 of fold change:
t.test.male <-
  t.test.male %>%
  mutate(log.2.FC = log2(FC))



# Separating up- and downregulated groups

t.test.male$new.condition <- ifelse(t.test.male$p < 0.05 & t.test.male$log.2.FC > 0, 'Significant upregulation', 
                                    ifelse(t.test.male$p < 0.05 & t.test.male$log.2.FC < 0, 'Significant downregulation', 'Insignificant'))

# adding names to colours
colours.male <- ifelse(t.test.male$new.condition == "Significant upregulation", "red4", 
                       ifelse(t.test.male$new.condition == "Significant downregulation", "royalblue", "grey"))


names(colours.male)[colours.male == 'grey'] <- 'Insignificant'
names(colours.male)[colours.male == 'red4'] <- 'Significant upregulation'
names(colours.male)[colours.male == 'royalblue'] <- 'Significant downregulation'

EnhancedVolcano(t.test.male,
                lab = t.test.male$Metabolite,
                x = "log.2.FC", 
                y = "p",
                xlab = "log2(Fold change)", 
                ylab = "-log10(p-value)",
                title = "Liver metabolites male 12m",
                pCutoff = 0.05,
                FCcutoff = 0,
                labSize = 4,
                drawConnectors = TRUE,
                arrowheads = F,
                max.overlaps = 15,
                colCustom = colours.male,
                ylim = c(0, 4),
                xlim = c(-3, 3)) 



# female only

female_data <- subset(data.long, Sex == "F")

# Compute Welch's t-test and adjust p-values (rstatix):
t.test.female <- 
  female_data %>%
  group_by(Metabolite) %>%
  rstatix::t_test(Concentrations ~ Genotype, var.equal = T) %>%
  adjust_pvalue(method = "BH")

# Compute mean concentrations for all lipids within healthy controls (N):
female.WT <- 
  female_data %>%
  group_by(Metabolite, Genotype) %>%
  summarise(mean = mean(Concentrations, na.rm = T)) %>%
  arrange(Genotype) %>%
  filter(Genotype == "WT")

# Compute mean concentrations for all lipids within patients with PDAC (T):
female.ALS <- 
  female_data %>%
  group_by(Metabolite, Genotype) %>%
  summarise(mean = mean(Concentrations, na.rm = T)) %>%
  arrange(Genotype) %>%
  filter(Genotype == "ALS")

# Compute fold-change and store it in the "Welch.test" as `Fold change`:
t.test.female$FC <- female.ALS$mean/female.WT$mean

# Using mutate(), compute log2 of fold change:
t.test.female <-
  t.test.female %>%
  mutate(log.2.FC = log2(FC))


# Separating up- and downregulated groups

t.test.female$new.condition <- ifelse(t.test.female$p < 0.05 & t.test.female$log.2.FC > 0, 'Significant upregulation', 
                                      ifelse(t.test.female$p < 0.05 & t.test.female$log.2.FC < 0, 'Significant downregulation', 'Insignificant'))

# adding names to colours
colours.female <- ifelse(t.test.female$new.condition == "Significant upregulation", "red4", 
                         ifelse(t.test.female$new.condition == "Significant downregulation", "royalblue", "grey"))


names(colours.female)[colours.female == 'grey'] <- 'Insignificant'
names(colours.female)[colours.female == 'red4'] <- 'Significant upregulation'
names(colours.female)[colours.female == 'royalblue'] <- 'Significant downregulation'

EnhancedVolcano(t.test.female,
                lab = t.test.female$Metabolite,
                x = "log.2.FC", 
                y = "p",
                xlab = "log2(Fold change)", 
                ylab = "-log10(p-value)",
                title = "Liver metabolites female 12m",
                pCutoff = 0.05,
                FCcutoff = 0,
                labSize = 4,
                drawConnectors = TRUE,
                arrowheads = F,
                max.overlaps = 12,
                colCustom = colours.female,
                ylim = c(0, 4),
                xlim = c(-3, 3))



## saving the results

write_xlsx(t.test.male[, c("Metabolite", "p", "log.2.FC")], path = "C:/Users/veera/Documents/Thesis project/Liver/Lipidomics/Lipids_liver_12m_results_male.xlsx")

write_xlsx(t.test.female[, c("Metabolite", "p", "log.2.FC")], path = "C:/Users/veera/Documents/Thesis project/Liver/Lipidomics/Lipids_liver_12m_results_female.xlsx")
