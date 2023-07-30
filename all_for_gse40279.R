#### author: 姚博瀚

### Part 1
## Download gse40279 and run analysis by using "MEAL" package (from original code 4.R)
library(MEAL)
library(minfi)
library(limma)
library(ggplot2)

# Install remotes from CRAN (Optional)
install.packages("remotes")
# Download modified GEOquery package from my github 
# by using function(install_github()) from 'remotes' package. 
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
# result_Meal_clean <- na.omit(result_Meal)

### Part 1-2
## manhattan plot for all
library(qqman)
# function from qqman to plot manhattan 
result_Meal$CHR <- as.numeric(result_Meal$CHR)
manhattan(result_Meal, 
          main = "Manhattan Plot for gse40279 (Analysis of DiffMean on MEAL)",
          cex = 0.6,
          ylim = c(0, 25),
          chr="CHR", 
          bp="RANGE_START", 
          snp= "ID", 
          p="P.Value" )

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

# Run MEAL pipeline on the subset YO
res_sub_YO <- runPipeline(set = subset_forMeal_YO,
                          variable_names = "new_age_category",
                          betas = TRUE,
                          analyses = c("DiffMean", "DiffVar"))

# Extract the result of the DiffMean analysis
result_Meal_sub_YO <- getProbeResults(res_sub_YO, rid = 1, 
                                      fNames = c("UCSC_RefGene_Name", "RANGE_START", "CHR", "ID"))

# Remove rows with missing values
result_Meal_sub_YO_clean <- na.omit(result_Meal_sub_YO)



# Create a new age category with 2 levels (Middle and Old)
MO_age_categories <- factor(age_categories, levels = c("Middle", "Old"))
MO_age_categories <- na.omit(MO_age_categories)

# Extract samples belonging to "Middle" and "Old" age categories
subset_forMeal_MO <- gse40279_matrix[, age_categories %in% c("Middle", "Old")]
subset_forMeal_MO_data <- exprs(subset_forMeal_MO)
# Assign the new age category to pData of the subset
pData(subset_forMeal_MO)$MO_age_category <- MO_age_categories

# Run MEAL pipeline on the subset YO
res_sub_MO <- runPipeline(set = subset_forMeal_MO,
                          variable_names = "MO_age_category",
                          betas = TRUE,
                          analyses = c("DiffMean", "DiffVar"))

# Extract the result of the DiffMean analysis
result_Meal_sub_MO <- getProbeResults(res_sub_MO, rid = 1, 
                                      fNames = c("UCSC_RefGene_Name", "RANGE_START", "CHR", "ID"))

# Remove rows with missing values
result_Meal_sub_MO_clean <- na.omit(result_Meal_sub_MO)



## manhattan plot for subset YO 
install.packages("qqman") # if you have not installed
library(qqman)

# Extract data for manhattan plot
res_M_YO_manhattan <- result_Meal_sub_YO_clean

# Convert CHR column to numeric
res_M_YO_manhattan$CHR <- as.numeric(res_M_YO_manhattan$CHR)

# Remove rows with missing values
res_M_YO_manhattan_clean <- na.omit(res_M_YO_manhattan)

# function from qqman to plot manhattan 
manhattan(res_M_YO_manhattan_clean, 
          main = "Manhattan Plot for gse40279 (subset_YO) (Analysis of DiffMean on MEAL) ",
          cex = 0.3,
          ylim = c(0, 50),
          chr="CHR", 
          bp="RANGE_START", 
          snp= "ID", 
          p="P.Value",
          genomewideline = FALSE,
          suggestiveline = -log10(1e-05))

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
     show.labels = FALSE) + ggtitle("Volcano plot for gse40279 (Analysis of DiffMean on MEAL)")

# Way 2 to plot volcano for subset YO
plot(res_sub_YO, rid = "DiffMean", type = "volcano", tPV = 14, tFC = 0.1, 
     show.labels = FALSE) + ggtitle("Volcano plot for gse40279 (subset YO)(Analysis of DiffMean on MEAL)")

# QQ plot for all
plot(res, rid = 1, type = "qq") + ggtitle("QQplot for gse40279 (Analysis of DiffMean on MEAL)")

# QQ plot for subset YO
plot(res_sub_YO, rid = 1, type = "qq") 
+ ggtitle("QQplot for gse40279 (subset YO)(Analysis of DiffMean on MEAL)")

# QQ plot for subset MO
plot(res_sub_MO, rid = 1, type = "qq") + ggtitle("QQplot for gse40279 (subset MO)(Analysis of DiffMean on MEAL)")

# Plot the beta values distribution of a CpG
plotFeature(set = gse40279_matrix, feat = "cg16762684", variables = "age_category") + 
  ggtitle("cg16762684 (Lowest P.value on the analysis of DIffMean)") +
  ylab("Methylation (Beta value)")




## (Error) not yet solved
# Regional plotting 
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
## (Error) not yet solved



### Part 2
## Load gse40279 and run analysis by using "limma" packages (from original code limma.R)
library(limma)
# Load modified GEOquery
library(GEOquery)
# Setting the max timeout_seconds
options(timeout=100000)
# Check the input timeout_seconds
getOption("timeout")
# Download GSE40279 by a fuction getGEO() from modified GEOquery package.
gse40279 <- getGEO("GSE40279", GSEMatrix = TRUE, AnnotGPL = TRUE) # if have not downloaded

gset <- gse40279
if (length(gset) > 1) idx <- grep("GPL13534", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# Create age categories
age <- pData(gset)$characteristics_ch1
# Remove "age (y):" and convert to numeric
age <- sub("^\\s*age \\(y\\): ", "", age)
age <- as.numeric(age)
# Assign age values to a new column in pData of gset
pData(gset)$age <- age

# Define age categories based on specific age ranges
age_categories <- cut(age,
                      breaks = c(0, 30, 65, Inf),
                      labels = c("Young", "Middle", "Old"),
                      include.lowest = TRUE)

# Assign age categories to the pData of gset
pData(gset)$age_category <- age_categories

ex <- exprs(gset)
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

result_limma_YM <- topTable(fit2, coef=1, number = Inf, adjust.method = "BH", sort.by = "logFC")
result_limma_YO <- topTable(fit2, coef=2, number = Inf, adjust.method = "BH", sort.by = "logFC")
result_limma_MO <- topTable(fit2, coef=3, number = Inf, adjust.method = "BH", sort.by = "logFC")

# Outcome of each hypothesis test
results <- decideTests(fit2, method = "separate", adjust.method = "BH", p.value = 1e-5)

# Showing numbers of genes significant in each comparison
vennDiagram(results, show.include = FALSE) + title('Numbers of genes significant (p.value<1e-5) in each comparison')

library(ggplot2)
# Add log-transformed p-value column to result_limma_YM
result_limma_YM$neg_logP <- -log10(result_limma_YM$P.Value)
# Specify the type of the probes in result_limma_YM(from gse40279)
probe_type <- ifelse(result_limma_YM$neg_logP > -log10(0.05), "sig", "ns")
# Add the processed probe types column to result_limma_YM
result_limma_YM$probe_type <- probe_type
# Define the colors for the values(sig or ns) in probe_type
cols <- c("sig" = "red", "ns" = "grey")
# Create a volcano plot with the configurable colors for the points   
ggplot(result_limma_YM, aes(x = logFC, y = neg_logP, color = probe_type)) +
  geom_point(size = 1) +
  scale_color_manual(values = cols) +
  labs(x = "log(Fold Change)", y = "-log10(P.Value)") +
  ggtitle("Volcano Plot: Y vs. M") +
  theme_minimal() +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")

# Remove rows with missing values
# result_limma_coef1_clean <- na.omit(result_limma_2)


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

# Mean-Difference Plot from limma
as.numeric(limma_Meal_YO_p$limma_p)
as.numeric(limma_Meal_YO_p$MEAL_p)
x1 <- log(limma_Meal_YO_p$limma_p)
x2 <- log(limma_Meal_YO_p$MEAL_p)
oldpar <- par(mfrow=c(1,2))
plot(x1,x2,xlab="loge(limma_p)",ylab="loge(MEAL_p)",main="Scatter Plot")
mdplot(cbind(x1,x2),bg.pch=1,bg.cex=1,main = "Mean-Difference Plot")
par(oldpar)
