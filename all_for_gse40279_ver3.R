
### Part 1
## Download gse40279 and run analysis by using "MEAL" package (from original code 4.R)
library(MEAL)
library(minfi)
library(limma)
library(ggplot2)

# Install mediation from CRAN
install.packages("mediation")

# Install remotes from CRAN
install.packages("remotes")
# Download modified GEOquery package from my github by using function(install_github()) from 'remotes' package. 
library(remotes)
install_github("curryhank08/GEOquery_with_modifiable_timeout_seconds", force = TRUE)
# Load modified GEOquery
library(GEOquery)
# Setting the max timeout_seconds
options(timeout=100000)
# Check the input timeout_seconds
getOption("timeout")

# Download GSE40279 by a fuction getGEO() from modified GEOquery package.
gse40279 <- getGEO("GSE40279", GSEMatrix = TRUE, AnnotGPL = TRUE)
gse40279_matrix <- gse40279[[1]]

data <- exprs(gse40279_matrix)
#####先試跑到以上

# Create age categories
age <- pData(gse40279_matrix)$characteristics_ch1

# Remove "age (y):" and convert to numeric
age <- sub("^\\s*age \\(y\\): ", "", age)
age <- as.numeric(age)

# The ^ character denotes the start of the string,
# \\s* matches any number of leading whitespace characters,
# and "age \\(y\\): " matches the exact string "age (y): ". 

# Assign age values to a new column in pData of gse40279_matrix
pData(gse40279_matrix)$age <- age

# Define age categories based on specific age ranges
age_categories <- cut(age,
                      breaks = c(0, 30, 65, Inf),
                      labels = c("Young", "Middle", "Old"),
                      include.lowest = TRUE)

# Assign age categories to the pData of gse40279_matrix
pData(gse40279_matrix)$age_category <- age_categories

# Assign seqnames to the pData of gse40279_matrix from CHR
fData(gse40279_matrix)$seqnames <- as.numeric(fData(gse40279_matrix)$CHR)

# Run MEAL pipeline on the categorized data
res <- runPipeline(set = gse40279_matrix,
                   variable_names = "age_category",
                   betas = TRUE,
                   analyses = c("DiffMean", "DiffVar"))

# Extract the result of the DiffMean analysis
result_Meal <- getProbeResults(res, rid = 1, 
                               fNames = c("UCSC_RefGene_Name", "RANGE_START", "CHR", "ID"))

# Remove rows with missing values
result_Meal_clean <- na.omit(result_Meal)

## Subset the data matrix for the filtered samples
# filtered_data <- data[, colnames(data) %in% rownames(young_old_samples)]

# Create a new age category with 2 levels (Old and Young)
new_age_categories <- factor(age_categories, levels = c("Old", "Young"))
new_age_categories <- na.omit(new_age_categories)

# Extract samples belonging to "Young" and "Old" age categories
subset_forMeal_YO <- gse40279_matrix[, age_categories %in% c("Young", "Old")]
subset_forMeal_YO_data <- exprs(subset_forMeal_YO)
# Assign the new age category to pData of the subset
pData(subset_forMeal_YO)$new_age_category <- new_age_categories
# subeset samples info
subset_forMeal_YO_samplesinfo <- pData(subset_forMeal_YO)


# Run MEAL pipeline on the subset
res_sub_YO <- runPipeline(set = subset_forMeal_YO,
                   variable_names = "new_age_category",
                   betas = TRUE,
                   analyses = c("DiffMean", "DiffVar"))



# Extract the result of the DiffMean analysis
result_Meal_sub_YO <- getProbeResults(res_sub_YO, rid = 1, 
                                   fNames = c("UCSC_RefGene_Name", "RANGE_START", "CHR", "ID"))

# Remove rows with missing values
result_Meal_sub_YO_clean <- na.omit(result_Meal_sub_YO)

### Part 1-2
## manhattan plot 
install.packages("qqman") # if you have not installed
library(qqman)

# Extract data for manhattan plot
res_M_manhattan <- result_Meal_sub_YO_clean

# Convert CHR column to numeric
res_M_manhattan$CHR <- as.numeric(res_M_manhattan$CHR)

# Remove rows with missing values
res_M_manhattan_clean <- na.omit(res_M_manhattan)


# function from qqman to plot manhattan 
manhattan(res_M_manhattan_clean, 
          main = "Manhattan Plot for gse40279 (Analysis of DiffMean on MEAL)",
          cex = 0.6,
          ylim = c(0, 25),
          chr="CHR", 
          bp="RANGE_START", 
          snp= "ID", 
          p="P.Value" )

## volcano plot
# Extract data for volcano plot (same as res_manhattan)
res_M_volcano <- result_Meal 
# Remove rows with missing values
res_M_volcano_clean <- na.omit(res_M_volcano)

# Add log-transformed p-value column to res_M_volcano_clean
res_M_volcano_clean$neg_logP <- -log10(res_M_volcano_clean$P.Value)

# Create volcano plot
ggplot(res_M_volcano_clean, aes(x = logFC, y = neg_logP)) +
  geom_point(size = 0.5, color = "BLUE")+
  ggtitle("Volcano plot for gse40279 (Analysis of DiffMean on MEAL)")+
  labs(x="log2(Flod change)", y="-log10(p-value)")

# Way 2 to plot volcano 
plot(res, rid = "DiffMean", type = "volcano", tPV = 14, tFC = 0.1, 
     show.labels = FALSE) + ggtitle("Volcano")

# QQ plot
plot(res, rid = 1, type = "qq") + ggtitle("QQplot")

# Plot the beta values distribution of a CpG
plotFeature(set = gse40279_matrix, feat = "cg16762684", variables = "age_category") + 
  ggtitle("cg16762684 (Lowest P.value on the analysis of DIffMean)") +
  ylab("Methylation (Beta value)")

# Regional plotting (Error) not yet solved
targetRange <- GRanges("18:60000000-100000000")

# Check the current value of fData(res)$Strand
currentStrand <- fData(gse40279_matrix)$Strand

# Transform '+' if the value is 'F', and '-' if the value is 'R'
transformedStrand <- ifelse(currentStrand == 'F', '+', ifelse(currentStrand == 'R', '-', currentStrand))

# Assign the transformed values back to fData(res)$Strand
fData(gse40279_matrix)$Strand <- transformedStrand

# Run MEAL pipeline on the categorized data
res <- runPipeline(set = gse40279_matrix,
                   variable_names = "age_category",
                   betas = TRUE,
                   analyses = c("DiffMean", "DiffVar"))

plotRegion(res, targetRange)
plotRegion(res, targetRange, results = c("DiffMean"), tPV = 24, fNames = c("chromosome", "start", "end"))




### Part 2
## Load gse40279 and run analysis by using "limma" packages (from original code limma.R)
library(GEOquery)
library(limma)
library(umap)

gset <- gse40279
if (length(gset) > 1) idx <- grep("GPL13534", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex <- exprs(gset)

# Assign age values to a new column in pData of gset
pData(gset)$age <- age
# Assign age categories to the pData of gset
pData(gset)$age_category <- age_categories

# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex)}

conditions <- gset$age_category
f <- factor(conditions, levels = c("Young", "Middle", "Old"))
design <- model.matrix(~0+f)
colnames(design) <- c("Young", "Middle", "Old")
fit <- lmFit(gset, design)
contrast.matrix <- makeContrasts(Young-Middle, Young-Old, Middle-Old, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
result_limma_YM <- topTable(fit2, coef=1, number = Inf, adjust.method = "BH")
result_limma_YO <- topTable(fit2, coef=2, number = Inf, adjust.method = "BH")
result_limma_MO <- topTable(fit2, coef=3, number = Inf, adjust.method = "BH")

# Outcome of each hypothesis test
results <- decideTests(fit2)

# Showing numbers of genes significant in each comparison
vennDiagram(results)

# Remove rows with missing values
result_limma_coef1_clean <- na.omit(result_limma_2)


### Part 3
## Merge the results of analysis from MEAL and limma, and then create a scatter plot (from original code 6.R)
library(ggplot2)

# Merge results
res_limma_YO_c_p <- data.frame(name = row.names(result_limma_YO), limma_p = result_limma_YO$P.Value)
res_MEAL_YO_c_p <- data.frame(name = row.names(result_Meal_sub_YO_clean), MEAL_p = result_Meal_sub_YO_clean$P.Value)

limma_Meal_YO_p <- merge(res_limma_YO_c_p,
                      res_MEAL_YO_c_p,
                      by.x = "name",
                      by.y = "name")

# Create a scatter plot with x-axis: p-value from limma and y-axis: p-value from MEAL.
plot(limma_Meal_YO_p$limma_p, limma_Meal_YO_p$MEAL_p, 
     xlab = "limma_P.value", ylab = "MEAL_P.value", 
     main = "P.value from Analysis of DiffMean (Young - Old) on MEAL and limma", 
     pch = 20, col = "#8bc34a", cex = 1)
