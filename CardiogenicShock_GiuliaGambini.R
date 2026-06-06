# --- Genomic profiling of Cardiogenic Shock Patients under Extracorporeal Membrane Oxygenation ---

# Is it possible to identify a gene expression signature in cardiogenic shock patients that can
# predict the early clinical outcome (Success vs. Failure) of ECMO support?"
# --- Giulia Gambini --- 


# WEEK 1 & 2: Data Loading & Exploration
# Libraries
install.packages("BiocManager")
BiocManager::install("GEOquery")
library("GEOquery")
BiocManager::install("useful")
library("useful")

# 1. Data set upload
gse <- getGEO("GSE182600", getGPL = TRUE)
gse <- gse[[1]]

ex <- exprs(gse)    # Expression matrix
pheno <- pData(gse)

# 2. Outcome and time informations
outcome <- sub("outcome: ", "", pheno$characteristics_ch1.3)
outcome[outcome == "failure"] <- "Failure" # To uniform the etiquettes
time <- sub("time: ", "", pheno$characteristics_ch1.5)

# 3. Filter the dataset: mantain only time 0 data

idx_0hr <- which(time == "0hr")
ex_0hr <- ex[, idx_0hr]
outcome_0hr <- outcome[idx_0hr]

# 4. Quality control 
boxplot(log(ex_0hr), main = "Boxplot at t=0hr", las=2, col="lightblue")

# NA removal
ex_clean <- na.omit(as.matrix(log(ex_0hr)))

dim(ex)
dim(ex_clean)
colnames(ex_clean)

# 5. PCA (Principal Component Analysis)
pca <- prcomp(t(ex_clean))

# Screeplot
screeplot(pca, main = "pca", col="lightblue")

# 6. PLOT PCA
grpcol <- ifelse(outcome_0hr == "Success", "skyblue", "indianred2")

# set up a 1x3 plotting area to compare all combinations
par(mfrow = c(1,3))

# Plot 1: PC1 vs PC2
plot(pca$x[,1], pca$x[,2],
     xlab = "PC1", ylab = "PC2",
     main = "PCA: PC1 vs PC2",
     pch = 16, col = grpcol)
legend("topright", legend = c("Success", "Failure"), 
       col = c("skyblue", "indianred2"), pch = 16, cex = 0.8)

# Plot 2: PC1 vs PC3
plot(pca$x[,1], pca$x[,3],
     xlab = "PC1", ylab = "PC3",
     main = "PCA: PC1 vs PC3",
     pch = 16, col = grpcol)
legend("topright", legend = c("Success", "Failure"), 
       col = c("skyblue", "indianred2"), pch = 16, cex = 0.8)

# Plot 3: PC2 vs PC3
plot(pca$x[,2], pca$x[,3],
     xlab = "PC2", ylab = "PC3",
     main = "PCA: PC2 vs PC3",
     pch = 16, col = grpcol)

# Legend
legend("topright", legend = c("Success", "Failure"), 
       col = c("skyblue", "indianred2"), pch = 16, cex = 0.8)

# Reset plotting parameters
par(mfrow = c(1, 1))


# WEEK 3: Data Clustering (Unsupervised Learning)
# Do 0hr patients naturally group into Successful/Failure?

library("useful")
install.packages('dendextend')
# 1. K-Means Clustering
# We set k = 2 to test our hypothesis of 2 clinical outcomes (Success vs Failure)
set.seed(123) # Ensures you get the exact same result every time 
k_kmeans <- 2
kmeans_result <- kmeans(t(ex_clean), centers = k_kmeans)

# Compare the k-means clusters to our actual clinical outcomes
cat("--- K-Means (k=2) vs Clinical Outcome ---\n")
print(table(Cluster = kmeans_result$cluster, Actual_Outcome = outcome_0hr))

# Visualize the 2-means clusters
plot(kmeans_result, data = t(ex_clean), main = "K-Means Clustering (k=2) at 0hr")

table(kmeans_result$cluster, time)

# 2. Hierarchical Clustering
# Calculate Euclidean distance between all patients
dist_matrix <- dist(t(ex_clean))

# Perform hierarchical clustering using the "average" linkage method
hc_result <- hclust(dist_matrix, method = "ave")

# Plot the dendrogram of the patients
plot(hc_result, hang = -1, 
     main = "Hierarchical Clustering Dendrogram (0hr)", 
     xlab = "Patients", ylab = "Distance", sub = "")


# Draw red boxes around the 2 main branches of the tree
rect.hclust(hc_result, k = 2, border = "indianred2")

# Extract the group assignments and compare with actual outcomes
hc_groups <- cutree(hc_result, k = 2)
cat("\n--- Hierarchical Clustering (k=2) vs Clinical Outcome ---\n")
print(table(Cluster = hc_groups, Actual_Outcome = outcome_0hr))



# WEEK 5: Random Forest (Supervised Learning)
library(randomForest)

# 1. Prepare and Clean the Data
ex_rf <- na.omit(as.matrix(ex_0hr))
ex_rf <- ex_rf[apply(ex_rf, 1, function(x) all(is.finite(x))), ]

# Ensure the outcome is a factor for classification
group_0hr <- as.factor(outcome_0hr)

# Pre-processing: feature selection
# Evaluate the variance for each gene
varianza_geni <- apply(ex_rf, 1, var)

# Order the genes from the most variable to the least one
geni_ordinati <- order(varianza_geni, decreasing = TRUE)

# Select only the top 5000
top5000_idx <- head(geni_ordinati, 5000)
ex_rf_filtrata <- ex_rf[top5000_idx, ]

cat("Filter applied, genes:", nrow(ex_rf_filtrata), "\n")


# 2. Build the Random Forest Model
set.seed(1234) # Ensures reproducibility
rf <- randomForest(x = t(ex_rf_filtrata), y = group_0hr, ntree = 1000) # 1000 trees

# top importance genes graph
varImpPlot(rf, n.var = 30, main = "Top 30 Important Genes (Random Forest - 0hr)")

ngenes <- 400 

# Extract the importance and order it
importanza_ordinata <- sort(rf$importance, decreasing = TRUE)

# Plot it
plot(importanza_ordinata[1:ngenes], 
     xlab = "Index", 
     ylab = "sort(rf$importance, decreasing = TRUE)[1:ngenes]")

cat("Avvio calcolo statistico sui geni...\n")

# 1. Calculate a p-value for each gene, is it different from Success and Failure patients?
# tryCatch block to ignore the genes that have identical expression
p_values <- apply(ex_rf, 1, function(x) {
  obj <- try(t.test(x ~ group_0hr), silent = TRUE)
  if (is(obj, "try-error")) return(1) else return(obj$p.value)
})

# 2. Order the genes for significance
geni_significativi <- order(p_values, decreasing = FALSE)

# 3. Take top 150 genes
top150_idx <- head(geni_significativi, 150)
ex_rf_mirata <- ex_rf[top150_idx, ]

cat("Filter applied, genes:", nrow(ex_rf_mirata), "\n")

# RANDOM FOREST on restricted data
set.seed(1234)
# Using the new matrix with only 150 top genes
rf_mirato <- randomForest(x = t(ex_rf_mirata), y = group_0hr, ntree = 1000)

# Errors plot
plot(rf_mirato, main = "Random Forest Error")
legend("topright", colnames(rf_mirato$err.rate), col = 1:3, lty = 1:3)
print(rf_mirato) 
# OOB estimate of error rate:18.18%




# Heatmap (Top 25 Genes)
library(RColorBrewer)

# 1. Extract values only for first most significant genes
top25_idx <- head(geni_significativi, 25)
sig_eset <- ex_rf[top25_idx, ] 

# 2. Set up colors for the heatmap (Purple/Orange color palette)
hmcol <- colorRampPalette(brewer.pal(11, "PuOr"))(256)

# 3. Create a color bar for the patients (Skyblue = Success, Indianred2 = Failure)
csc <- rep("indianred2", length(group_0hr)) 
csc[group_0hr == "Success"] <- "skyblue"    

# 4. Draw the Heatmap
heatmap(sig_eset, 
        scale = "row",        
        col = hmcol, 
        ColSideColors = csc,  
        main = "Heatmap: Top 25 Predictive Genes (0hr)",
        margins = c(5, 10))



# WEEK 6 & 7: Linear Discriminant Analysis (LDA)
# Goal: Build a linear classifier to separate Success/Failure at 0h
BiocManager::install("genefilter")
BiocManager::install("pROC")
library("genefilter")
library("MASS")
library("pROC")

# 1. Convert in factor
f_0hr <- factor(outcome_0hr) 

# 2. Pre-filtering for dimensionality reduction (p < 0.1)
tt_0hr <- rowttests(ex_clean, f_0hr)
keepers_0hr <- which(tt_0hr$p.value < 0.1)

# New reduced matrix 
ex_filtered <- ex_clean[keepers_0hr, ]
tex_filtered <- t(ex_filtered)

# Create the dataframe for the LDA
dat_0hr <- data.frame(tex_filtered, COLUMN = f_0hr)

# 3. Test and training sets
idx_success <- which(dat_0hr$COLUMN == "Success")
idx_failure <- which(dat_0hr$COLUMN == "Failure")

# Select randomly the test group
set.seed(123) 
test_success <- sample(idx_success, 5)
test_failure <- sample(idx_failure, 5)
test_idx <- c(test_success, test_failure)

# The remaining ones become the train set
train_idx <- setdiff(1:nrow(dat_0hr), test_idx)

# 4. Training the LDA model
mod_lda <- lda(COLUMN ~ ., data = dat_0hr, prior = c(0.5, 0.5), subset = train_idx)

# 5. Evaluation on the train set (how the model sees the data on which it trained)
preds_train <- predict(mod_lda, dat_0hr[train_idx, ])

# Assign colors based on the outcome
point_colors <- ifelse(dat_0hr[train_idx, "COLUMN"] == "Success", "skyblue", "indianred2")

# Plot of the discriminant scores
plot(preds_train$x[,1], 
     ylab = "LDA Score", 
     xlab = "Patient Index",
     main = "LDA Discriminant Scores (Training Data - 0h)",
     col = point_colors, 
     pch = 16)
legend("topright", legend = c("Success", "Failure"), col = c("skyblue", "indianred2"), pch = 16)

# 6. Validation on the test set
preds_test <- predict(mod_lda, dat_0hr[test_idx, ])

# Confusion matrix
cat("--- LDA Confusion Matrix (Test Set) ---\n")
print(table(Predicted = preds_test$class, Actual = dat_0hr[test_idx, "COLUMN"]))


# 7. ROC curve and AUC
# To measure the diagnostic capacity of the model on the unknown data
roc_obj <- roc(dat_0hr[test_idx, "COLUMN"], preds_test$posterior[, "Success"])

plot.roc(roc_obj, 
         main = "ROC Curve: LDA Performance (0h)", 
         col = "darkblue", 
         lwd = 4, 
         print.auc = TRUE,
         # auc.polygon = TRUE,    
         auc.polygon.col = "lightblue",
         grid = TRUE)



# WEEK 7: Repeated Cross-Validation with Caret
# Goal: Test the LDA model's average performance across multiple splits
library("caret")

# 1. Data preparation
colnames(dat_0hr) <- make.names(colnames(dat_0hr))
colnames(dat_0hr)[ncol(dat_0hr)] <- "COLUMN"

x_dati <- dat_0hr[, -ncol(dat_0hr)]
y_target <- as.factor(dat_0hr$COLUMN)

# 1. 5 folds, 3 repetitions (15 runs total)
control_cv_fast <- trainControl(method = "repeatedcv", number = 5, repeats = 3)

# 2. Block the mtry parameter
parametro_fisso <- expand.grid(mtry = floor(sqrt(ncol(x_dati))))

# Train the LDA
set.seed(123)
cat("Addestramento LDA in corso...\n")
fit_lda_cv <- train(x = x_dati, y = y_target, method = "lda", 
                    metric = "Accuracy", trControl = control_cv_fast)

# Train the Random Forest
set.seed(123)
cat("Addestramento Random Forest in corso (attendi circa 10-30 secondi)...\n")
fit_rf_cv <- train(x = t(ex_rf_mirata), y = y_target, method = "rf", 
                   metric = "Accuracy", trControl = control_cv_fast, 
                   tuneGrid = parametro_fisso, 
                   ntree = 150) 

# Comparison of the results
results_cv <- resamples(list(LDA = fit_lda_cv, RF = fit_rf_cv))

cat("\n--- Average results of the Cross-Validation---\n")
summary(results_cv)

# Plot
ggplot(results_cv) + 
  labs(y = "Accuracy", title = "Model comparison: LDA vs Random Forest") +
  theme_minimal()


# A better plot
bwplot(results_cv, 
       main = "Model Comparison: LDA vs Random Forest",
       ylab = "Accuracy / Kappa",
       scales = list(tck = c(1,0), x = list(cex = 1.2), y = list(cex = 1.2)))



# WEEK 8-1: LASSO Regression (Regularized Learning)
# Goal: Reduce overfitting and select the core gene signature
library("glmnet")
library("caret")

# 1. Prepare data for Lasso
# x must be a matrix, y must be numeric (1 for Success, 0 for Failure)
x_lasso <- t(ex_clean) 
y_lasso <- ifelse(outcome_0hr == "Success", 1, 0)

# 2. Run Cross-Validated LASSO
# This finds the optimal penalty (lambda)
set.seed(123)
cv_fit <- cv.glmnet(x_lasso, y_lasso, family = "binomial", alpha = 1)

# Plot 1: The Error Curve
# Look for the two vertical lines: lambda.min and lambda.1se
plot(cv_fit, main = "LASSO: Binomial Deviance vs. Lambda")

# Plot 2: The Coefficient Paths
# Shows how gene coefficients shrink to zero as penalty increases
plot(cv_fit$glmnet.fit, xvar = "lambda", label = TRUE, main = "LASSO Coefficient Shrinkage")

# 3. Extract the Winning Gene Signature
# We use lambda.1se for a more robust, "parsimonious" model
lasso_coefs <- coef(cv_fit, s = "lambda.1se")
active_genes <- rownames(lasso_coefs)[which(lasso_coefs != 0)]
active_genes <- active_genes[active_genes != "(Intercept)"]
cat("LASSO selected", length(active_genes), "key genes at 0hr.\n")

# Probe ID best gene
cat("Il Probe ID del gene selezionato è:", active_genes, "\n")
# ILMN_1726666



# WEEK 8-2: SCUDO 
library("rScudo")
library("igraph")
library("caret")

# 1. Data preparation
scudo_matrix <- ex_rf_mirata
scudo_labels <- group_0hr
set.seed(123)
inTrain_scudo <- createDataPartition(scudo_labels, p = 0.70, list = FALSE) 
trainData <- scudo_matrix[, inTrain_scudo]
testData  <- scudo_matrix[, -inTrain_scudo]
trainGroups <- scudo_labels[inTrain_scudo]
testGroups  <- scudo_labels[-inTrain_scudo]

# 2. Parameter optimization
model_scudo <- scudoModel(nTop = (2:6)*5, nBottom = (2:6)*5, N = 0.2)
control_scudo <- trainControl(method = "cv", number = 5, summaryFunction = multiClassSummary)
fit_scudo <- train(x = t(trainData), y = trainGroups, method = model_scudo, trControl = control_scudo)

best_nTop <- fit_scudo$bestTune$nTop
best_nBottom <- fit_scudo$bestTune$nBottom

# 3. Final training 
train_res <- scudoTrain(trainData, groups = trainGroups, 
                        nTop = best_nTop, nBottom = best_nBottom, alpha = 0.05)
train_net <- scudoNetwork(train_res, N = 0.2)

# Training visualization (colors corresponding to the outcome)
scudoPlot(train_net, vertex.label = NA, main = paste("SCUDO Training Network"))

# 4. Test and validate
test_res <- scudoTest(train_res, testData, testGroups, 
                      nTop = best_nTop, nBottom = best_nBottom)
test_net <- scudoNetwork(test_res, N = 0.2)

# Test visualization (colors corresponding to the outcome)
scudoPlot(test_net, vertex.label = NA, main = "SCUDO Test Network")

# 5. Clustering (Walktrap) and Final visualization
test_clust <- igraph::cluster_walktrap(test_net)
plot(test_clust, test_net, vertex.label = NA, 
     main = "SCUDO Test Set")

# 6. Classification and Confusion Matrix
class_res <- scudoClassify(trainData, testData, N = 0.2,
                           nTop = best_nTop, nBottom = best_nBottom,
                           trainGroups = trainGroups, alpha = 0.05)
previsioni <- class_res$predicted 
print(table(Predicted = previsioni, Actual = testGroups))
confusionMatrix(previsioni, testGroups)



# G-profiler analysis
library("caret")

# 1. Importance extraction from the LDA model
lda_importance <- varImp(fit_lda_cv, scale = FALSE)
lda_matrix <- as.matrix(lda_importance$importance)

# 2. Selection of the first 500 probe IDs
lda_ordered_probes <- rownames(lda_matrix)[order(lda_matrix[, 1], decreasing = TRUE)]
lda_top500_probes <- head(lda_ordered_probes, 500)

# 3. File writing
write.table(lda_top500_probes, 
            file = "LDA_Top500_Probes_Grezze.txt", 
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE, 
            sep = "\n")

cat("File 'LDA_Top500_Probes.txt' created correctly.\n")


# WEEK 10: pathfindR
install.packages("BiocManager")
BiocManager::install("KEGGREST")
BiocManager::install("KEGGgraph")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
library("KEGGREST")
library("KEGGgraph")
library("AnnotationDbi")
library("org.Hs.eg.db")
install.packages("pathfindR")
library("pathfindR")

# 1. Importance extraction
imp <- varImp(fit_lda_cv, scale = FALSE)
imp_matrix <- as.matrix(imp$importance)

# 2. Creation of the input
my_input_pathfindR <- data.frame(
  Gene.symbol = rownames(imp_matrix),
  Change = as.numeric(imp_matrix[,1]),
  p.value = 0.01                      
)
head(my_input_pathfindR)


# 1. Libraries
library(illuminaHumanv4.db)
library(AnnotationDbi)

# 2. Conversion in genic symbols
sonde_list <- rownames(imp_matrix)

traduzione <- AnnotationDbi::select(illuminaHumanv4.db, 
                                    keys = as.character(sonde_list), 
                                    columns = c("SYMBOL"), 
                                    keytype = "PROBEID")

# Verify
head(traduzione)

# 3. Data preparation
input_pathfindR_finale <- data.frame(
  Gene.symbol = traduzione$SYMBOL,
  Change = as.numeric(imp_matrix[match(traduzione$PROBEID, rownames(imp_matrix)), 1]),
  p.value = 0.01
)

input_pathfindR_finale <- input_pathfindR_finale[!is.na(input_pathfindR_finale$Gene.symbol), ]
input_pathfindR_finale <- aggregate(Change ~ Gene.symbol, data = input_pathfindR_finale, FUN = max)
input_pathfindR_finale$p.value <- 0.01

# Verify
head(input_pathfindR_finale)

# 4. Executing
RA_demo <- run_pathfindR(input_pathfindR_finale, iterations = 1)

term_gene_graph(RA_demo)

# 4. Executing
RA_demo <- run_pathfindR(my_input_pathfindR, iterations = 1)

term_gene_graph(RA_demo)

# 5. Visualization
RA_demo_clu <- cluster_enriched_terms(RA_demo)

RA_demo2 <- run_pathfindR(input_pathfindR_finale,
                          gene_sets = "Reactome",
                          pin_name_path = "STRING",
                          n_processes = 1)

head(RA_demo2)
term_gene_graph(RA_demo2)


# WEEK 11: CCA-SHRINKAGE
BiocManager::install("mixOmics")
library(mixOmics)

# 1. Complexity reduction
varianze <- apply(t(ex_0hr), 2, var)
top_geni_idx <- order(varianze, decreasing = TRUE)[1:500]
X_ridotta <- t(ex_0hr)[, top_geni_idx] # Matrix 33x500
Y <- t(ex_rf_mirata) 

# 2. CCA Analysis
result.rcca <- rcc(X_ridotta, Y, method = 'shrinkage')

# 3. Projection plot
plotIndiv(result.rcca, 
          group = outcome_0hr, 
          legend = TRUE, 
          title = 'Cardiogenic Shock: CCA-shrinkage')

# 4. Plot variable correlation
plotVar(result.rcca, cutoff = 0.5, cex = c(1, 1.5))


# MIXOMICS 
# Reduce X to the top 500 more variable genes for evaluation
varianze <- apply(t(ex_0hr), 2, var)
top_geni_idx <- order(varianze, decreasing = TRUE)[1:500]

X <- t(ex_0hr)[, top_geni_idx] # Matrix 33x500
Y <- t(ex_rf_mirata)          # Matrix 33x150

# 2. CCA analysis
result.rcca.shock <- rcc(Y, X, method = 'shrinkage')

# 3. Plot projection in the canonic variables subspace
colori_tesi <- c("indianred", "steelblue") 

# Plot Indiv (Patients)
plotIndiv(result.rcca.shock, 
          group = outcome_0hr, 
          col = c("indianred", "steelblue"), 
          legend = TRUE, 
          ind.names = TRUE,
          #pch = 16,
          title = 'Cardiogenic Shock: CCA')

# Plot Variables
plotVar(result.rcca.shock, 
        cutoff = 0.2, 
        col = c("forestgreen", "darkorange"), 
        cex = c(1, 1.5))



# WEEK 11: PLS (Partial Least Squares) Analysis

# 1. PLS Analysis
result.pls.shock <- pls(X, Y, ncomp = 2)

# 2. Plot projection in the canonic variables subspace
plotIndiv(result.pls.shock, 
          group = outcome_0hr, 
          col = c("indianred", "steelblue"), 
          legend = TRUE, 
          ind.names = TRUE,
          title = 'Cardiogenic Shock: PLS')

# 3. Plot Variables
plotVar(result.pls.shock, 
        cutoff = 0.5, 
        col = c("forestgreen", "darkorange"), 
        cex = c(1, 1.5),
        title = 'Cardiogenic Shock: PLS Variable Correlation')

