# Using 1000 Genomes population PCs to train a random forest classifier and 
# using the model to classify samples using the ancestry labels in the 1000 genomes data

library(caret)
library(randomForest)
library(pROC)
library(plyr)

setwd("/scratch/gen1/rn180/umap/")

# cohorts+1000G PCs 
dt<- read.table("allcohorts_plus_1000G.pca.eigenvec")
for (i in 1:20){
  names(dt)[i+2] <- paste0("PC", i)
}

# kgp metadata
kgp <-read.csv("/scratch/gen1/rn180/refs/KGP/20130606_sample_info.csv")

# keep unique sample ID and Population
kgp <- kgp[,c(1,3)]
colnames(kgp) <- c("V2", "population")
# join pcs to population by the individual ID
kgp_pcs <- merge(dt, kgp, by="V2")

# make super groups to reduce colour scheme
kgp_pcs$superpop <- "superpop"
kgp_pcs$superpop <- ifelse(kgp_pcs$population=="ACB", "AFR", kgp_pcs$superpop)
kgp_pcs$superpop <- ifelse(kgp_pcs$population=="ASW", "AFR", kgp_pcs$superpop)
kgp_pcs$superpop <- ifelse(kgp_pcs$population=="GWD", "AFR", kgp_pcs$superpop)
kgp_pcs$superpop <- ifelse(kgp_pcs$population=="ESN", "AFR", kgp_pcs$superpop)
kgp_pcs$superpop <- ifelse(kgp_pcs$population=="LWK", "AFR", kgp_pcs$superpop)
kgp_pcs$superpop <- ifelse(kgp_pcs$population=="MSL", "AFR", kgp_pcs$superpop)
kgp_pcs$superpop <- ifelse(kgp_pcs$population=="YRI", "AFR", kgp_pcs$superpop)

kgp_pcs$superpop <- ifelse(kgp_pcs$population=="BEB", "SAS", kgp_pcs$superpop)
kgp_pcs$superpop <- ifelse(kgp_pcs$population=="GIH", "SAS", kgp_pcs$superpop)
kgp_pcs$superpop <- ifelse(kgp_pcs$population=="ITU", "SAS", kgp_pcs$superpop)
kgp_pcs$superpop <- ifelse(kgp_pcs$population=="PJL", "SAS", kgp_pcs$superpop)
kgp_pcs$superpop <- ifelse(kgp_pcs$population=="STU", "SAS", kgp_pcs$superpop)

kgp_pcs$superpop <- ifelse(kgp_pcs$population=="GBR", "EUR", kgp_pcs$superpop)
kgp_pcs$superpop <- ifelse(kgp_pcs$population=="CEU", "EUR", kgp_pcs$superpop)
kgp_pcs$superpop <- ifelse(kgp_pcs$population=="FIN", "EUR", kgp_pcs$superpop)
kgp_pcs$superpop <- ifelse(kgp_pcs$population=="TSI", "EUR", kgp_pcs$superpop)
kgp_pcs$superpop <- ifelse(kgp_pcs$population=="IBS", "EUR", kgp_pcs$superpop)

kgp_pcs$superpop <- ifelse(kgp_pcs$population=="CDX", "EAS", kgp_pcs$superpop)
kgp_pcs$superpop <- ifelse(kgp_pcs$population=="CHB", "EAS", kgp_pcs$superpop)
kgp_pcs$superpop <- ifelse(kgp_pcs$population=="CHS", "EAS", kgp_pcs$superpop)
kgp_pcs$superpop <- ifelse(kgp_pcs$population=="JPT", "EAS", kgp_pcs$superpop)
kgp_pcs$superpop <- ifelse(kgp_pcs$population=="KHV", "EAS", kgp_pcs$superpop)

kgp_pcs$superpop <- ifelse(kgp_pcs$population=="CLM", "AMR", kgp_pcs$superpop)
kgp_pcs$superpop <- ifelse(kgp_pcs$population=="MXL", "AMR", kgp_pcs$superpop)
kgp_pcs$superpop <- ifelse(kgp_pcs$population=="PEL", "AMR", kgp_pcs$superpop)
kgp_pcs$superpop <- ifelse(kgp_pcs$population=="PUR", "AMR", kgp_pcs$superpop)



# Split the data into training and testing sets

set.seed(123)  # Set seed for reproducibility
splitIndex <- createDataPartition(kgp_pcs$superpop, p = 0.8, list = FALSE)
train_data <- kgp_pcs[splitIndex, ]
test_data <- kgp_pcs[-splitIndex, ]
train_labels <- kgp_pcs$superpop[splitIndex]
test_labels <- kgp_pcs$superpop[-splitIndex]


# Define the training control for cross validation
train_control <- trainControl(method = "repeatedcv",
                              number = 5,
                              repeats = 5,
                              classProbs = TRUE)

# Parameter tuning - how many trees to try
tuning_parameters <- expand.grid(mtry = c(1:10))

# Build/train RF model
set.seed(123)
rf_model <- train(superpop ~ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,
                      data = train_data,
                      method = "rf",
                      trControl = train_control,
                      verbose = TRUE,
                      tuneGrid=tuning_parameters,
                      metric="Accuracy")

rf_model
#png(file="rf_model.png",width=600, height=350)
plot(rf_model, metric = "Accuracy")
#dev.off()

# Make predictions on the test set
predictions <- predict(rf_model, test_data, type="prob")

predictions2 <- predict(rf_model, test_data)

# Evaluate the model
conf_matrix <- confusionMatrix(predictions2, as.factor(test_labels))

# Plot Confusion Matrix
conf_matrix$table
conf_matrix
#plot(conf_matrix$table, col = c("darkgreen", "darkred"), main = "Confusion Matrix", cex = 1.2)

plt <- as.data.frame(conf_matrix$table)
plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))

ggplot(plt, aes(Reference,Prediction, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="#009194") +
  labs(x = "Reference",y = "Prediction") +
  scale_x_discrete(labels= colnames(conf_matrix$table)) +
  scale_y_discrete(labels=rev(colnames(conf_matrix$table)))+
  theme_bw()+
  ggtitle("confusion matrix")
ggsave("1000G_confusion1.png")

#================================
# get ipf only
ipf <- dt[!dt$V2 %in% kgp$V2, ]

pred_ipf <- predict(rf_model, ipf)
ipf$pred <- pred_ipf
table(ipf$pred)

pred_ipf2 <- predict(rf_model, ipf, type="prob")

ipf_grouped <- cbind(ipf,pred_ipf2)
ipf_grouped$pred0.9 <- "NA"
ipf_grouped$pred0.9 <- ifelse(ipf_grouped$AFR > 0.9, "AFR", ipf_grouped$pred0.9)
ipf_grouped$pred0.9 <- ifelse(ipf_grouped$EUR > 0.9, "EUR", ipf_grouped$pred0.9)
ipf_grouped$pred0.9 <- ifelse(ipf_grouped$AMR > 0.9, "AMR", ipf_grouped$pred0.9)
ipf_grouped$pred0.9 <- ifelse(ipf_grouped$EAS > 0.9, "EAS", ipf_grouped$pred0.9)
ipf_grouped$pred0.9 <- ifelse(ipf_grouped$SAS > 0.9, "SAS", ipf_grouped$pred0.9)
table(ipf_grouped$pred0.9)

ipf_grouped$pred0.8 <- "NA"
ipf_grouped$pred0.8 <- ifelse(ipf_grouped$AFR > 0.8, "AFR", ipf_grouped$pred0.8)
ipf_grouped$pred0.8 <- ifelse(ipf_grouped$EUR > 0.8, "EUR", ipf_grouped$pred0.8)
ipf_grouped$pred0.8 <- ifelse(ipf_grouped$AMR > 0.8, "AMR", ipf_grouped$pred0.8)
ipf_grouped$pred0.8 <- ifelse(ipf_grouped$EAS > 0.8, "EAS", ipf_grouped$pred0.8)
ipf_grouped$pred0.8 <- ifelse(ipf_grouped$SAS > 0.8, "SAS", ipf_grouped$pred0.8)
table(ipf_grouped$pred0.8)

ipf_grouped$pred0.7 <- "NA"
ipf_grouped$pred0.7 <- ifelse(ipf_grouped$AFR > 0.7, "AFR", ipf_grouped$pred0.7)
ipf_grouped$pred0.7 <- ifelse(ipf_grouped$EUR > 0.7, "EUR", ipf_grouped$pred0.7)
ipf_grouped$pred0.7 <- ifelse(ipf_grouped$AMR > 0.7, "AMR", ipf_grouped$pred0.7)
ipf_grouped$pred0.7 <- ifelse(ipf_grouped$EAS > 0.7, "EAS", ipf_grouped$pred0.7)
ipf_grouped$pred0.7 <- ifelse(ipf_grouped$SAS > 0.7, "SAS", ipf_grouped$pred0.7)
table(ipf_grouped$pred0.7)

ipf_grouped$pred0.5 <- "NA"
ipf_grouped$pred0.5 <- ifelse(ipf_grouped$AFR > 0.5, "AFR", ipf_grouped$pred0.5)
ipf_grouped$pred0.5 <- ifelse(ipf_grouped$EUR > 0.5, "EUR", ipf_grouped$pred0.5)
ipf_grouped$pred0.5 <- ifelse(ipf_grouped$AMR > 0.5, "AMR", ipf_grouped$pred0.5)
ipf_grouped$pred0.5 <- ifelse(ipf_grouped$EAS > 0.5, "EAS", ipf_grouped$pred0.5)
ipf_grouped$pred0.5 <- ifelse(ipf_grouped$SAS > 0.5, "SAS", ipf_grouped$pred0.5)
table(ipf_grouped$pred0.5)
ncol(ipf_grouped)
table_stats <- lapply(ipf_grouped[,c(23,29,30,31)], table)
summary_df <- as.data.frame(do.call(rbind, table_stats))
summary_df$total <- rowSums(summary_df)

#--------------------not important! but i like to keep it here
summary_df$prob <- row.names(summary_df)
library(tidyr)
library(tibble)
summary_df %>% gather(key, values)

df_long <- gather(data = summary_df, key = "Ancestry", value = "samples", -prob)

ggplot(df_long, aes(x = Ancestry, y = samples, color = prob, group = prob)) +
  stat_summary(fun = "mean", geom = "bar", position = "dodge") +
  labs(x = "Ancestry", y = "Number of samples") +
  theme_minimal() +
  scale_color_discrete(name = "Probability")
#-------------------------------

# plot the predictions
pal <- c(
  "AFR" = "#E41A1C",
  "EUR" = "#377EB8", 
  "AMR" = "#4DAF4A", 
  "EAS" = "#FF7F00",
  "SAS" = "#984EA3",
  "NA" = "black"
)

plot_pca <- function(pca, first_PC, second_PC, title, pred) {
  p_pca <- ggplot(pca, aes_string(x=first_PC, y=second_PC, color=pred)) +
    geom_point() +
    theme_classic() +
    scale_color_manual(values = pal,limits = names(pal))+
    guides(color=guide_legend(title="Population")) +
    scale_shape_manual(name = "Samples", values = c("IPF" = 1, "1000G" = 2))+
    theme(text=element_text(size=9))+
    ggtitle(title)
  return(p_pca)
}

plot_pca(ipf_grouped, "PC1", "PC2", "IPF by Highest Prob", "pred")
ggsave("IPF_pred.png")

plot_pca(kgp_pcs, "PC1", "PC2", "1000Genomes", "superpop")
ggsave("1000G_PCA.png")

plot_pca(ipf_grouped, "PC1", "PC2", "IPF by Prob>0.9", "pred0.9")
ggsave("IPF_pred90.png")

plot_pca(ipf_grouped, "PC1", "PC2", "IPF by Prob>0.8", "pred0.8")
ggsave("IPF_pred80.png")

plot_pca(ipf_grouped, "PC1", "PC2", "IPF by Prob>0.7", "pred0.7")
ggsave("IPF_pred70.png")

plot_pca(ipf_grouped, "PC1", "PC2", "IPF by Prob>0.5", "pred0.5")
ggsave("IPF_pred50.png")
