set.seed(12345)
library (bamdit)
library (ggplot2)
library (readr)

library(meta)

# Read your data
serum_data <- read.csv("serum.csv")

# Check the structure of your data
str(serum_data)
head(serum_data)

# Create study labels by combining author and year
study_labels <- paste0(serum_data$author, " (", serum_data$year, ")")

# Display the labels to verify they look correct
print(study_labels)

# SENSITIVITY FOREST PLOT
# =======================

# Calculate events and totals for sensitivity
# Sensitivity = TP / (TP + FN)
sens_events <- serum_data$TP
sens_total <- serum_data$TP + serum_data$FN

# Create meta-analysis object for sensitivity
sens_meta <- metaprop(event = sens_events,
                      n = sens_total,
                      studlab = study_labels,
                      sm = "PRAW",  # Use raw proportions
                      fixed = FALSE,
                      random = TRUE,
                      method.ci = "CP",  # Clopper-Pearson CI
                      title = "Sensitivity Forest Plot")

# Print summary
summary(sens_meta)

# Create sensitivity forest plot
pdf("sensitivity_forest.pdf", width = 12, height = 14)  # Optional: save to PDF
forest(sens_meta,
       sortvar = -TE,  # Sort by effect size (descending)
       xlab = "Sensitivity",
       smlab = "",
       leftcols = "studlab",
       leftlabs = "Study",
       rightcols = c("effect", "ci"),
       rightlabs = c("Sensitivity", "95% CI"),
       xlim = c(0, 1),
       col.diamond = "blue",
       col.diamond.lines = "blue",
       print.I2 = FALSE,  # Don't print I-squared in this context
       print.tau2 = FALSE,  # Don't print tau-squared
       cex = 0.7,  # Overall text size
       cex.studlab = 0.6,  # Study label size
       fontsize = 10)
dev.off()  # Close PDF

# Also display in R
forest(sens_meta,
       sortvar = -TE,
       xlab = "Sensitivity",
       smlab = "",
       leftcols = "studlab",
       leftlabs = "Study",
       rightcols = c("effect", "ci"),
       rightlabs = c("Sensitivity", "95% CI"),
       xlim = c(0, 1),
       col.diamond = "blue",
       col.diamond.lines = "blue",
       print.I2 = FALSE,
       print.tau2 = FALSE,
       cex = 0.7,
       cex.studlab = 0.6,
       fontsize = 10)

# SPECIFICITY FOREST PLOT
# =======================

# Calculate events and totals for specificity
# Specificity = TN / (TN + FP)
spec_events <- serum_data$TN
spec_total <- serum_data$TN + serum_data$FP

# Create meta-analysis object for specificity
spec_meta <- metaprop(event = spec_events,
                      n = spec_total,
                      studlab = study_labels,
                      sm = "PRAW",  # Use raw proportions
                      fixed = FALSE,
                      random = TRUE,
                      method.ci = "CP",  # Clopper-Pearson CI
                      title = "Specificity Forest Plot")

# Print summary
summary(spec_meta)

# Create specificity forest plot
pdf("specificity_forest.pdf", width = 12, height = 14)  # Optional: save to PDF
forest(spec_meta,
       sortvar = -TE,  # Sort by effect size (descending)
       xlab = "Specificity",
       smlab = "",
       leftcols = "studlab",
       leftlabs = "Study",
       rightcols = c("effect", "ci"),
       rightlabs = c("Specificity", "95% CI"),
       xlim = c(0, 1),
       col.diamond = "red",
       col.diamond.lines = "red",
       print.I2 = FALSE,
       print.tau2 = FALSE,
       cex = 0.7,
       cex.studlab = 0.6,
       fontsize = 10)
dev.off()  # Close PDF

# Also display in R
forest(spec_meta,
       sortvar = -TE,
       xlab = "Specificity",
       smlab = "",
       leftcols = "studlab",
       leftlabs = "Study",
       rightcols = c("effect", "ci"),
       rightlabs = c("Specificity", "95% CI"),
       xlim = c(0, 1),
       col.diamond = "red",
       col.diamond.lines = "red",
       print.I2 = FALSE,
       print.tau2 = FALSE,
       cex = 0.7,
       cex.studlab = 0.6,
       fontsize = 10)

# CREATE COMBINED PLOT (OPTIONAL)
# ================================

# Create both plots side by side
par(mfrow = c(1, 2))

forest(sens_meta,
       sortvar = -TE,
       xlab = "Sensitivity",
       smlab = "",
       leftcols = "studlab",
       leftlabs = "Study",
       rightcols = c("effect", "ci"),
       rightlabs = c("Sens.", "95% CI"),
       xlim = c(0, 1),
       col.diamond = "blue",
       col.diamond.lines = "blue",
       print.I2 = FALSE,
       print.tau2 = FALSE,
       cex = 0.5,
       cex.studlab = 0.4)

forest(spec_meta,
       sortvar = -TE,
       xlab = "Specificity",
       smlab = "",
       leftcols = "studlab",
       leftlabs = "Study",
       rightcols = c("effect", "ci"),
       rightlabs = c("Spec.", "95% CI"),
       xlim = c(0, 1),
       col.diamond = "red",
       col.diamond.lines = "red",
       print.I2 = FALSE,
       print.tau2 = FALSE,
       cex = 0.5,
       cex.studlab = 0.4)

par(mfrow = c(1, 1))  # Reset to single plot

# EXPORT RESULTS TO TABLE
# ========================

# Create a comprehensive results table
results_table <- data.frame(
  Study = study_labels,
  TP = serum_data$TP,
  FP = serum_data$FP,
  FN = serum_data$FN,
  TN = serum_data$TN,
  Sensitivity = round(sens_events / sens_total, 3),
  Sens_Lower_CI = round(sens_meta$lower, 3),
  Sens_Upper_CI = round(sens_meta$upper, 3),
  Specificity = round(spec_events / spec_total, 3),
  Spec_Lower_CI = round(spec_meta$lower, 3),
  Spec_Upper_CI = round(spec_meta$upper, 3)
)

# View the results
View(results_table)

# Save to CSV
write.csv(results_table, "forest_plot_results.csv", row.names = FALSE)

# Print first few rows
head(results_table)
## ggplot workaround-serum- design

library(ggplot2)
library(dplyr)

# Extract data from the bamdit object
# First, check what's in your serum object
str(serum)  # This will show available components

# Typical extraction (adjust based on your actual data structure)
# Method 1: If serum has TP, FP, TN, FN columns
plot_data <- data.frame(
  TP = serum_data$TP,
  FP = serum_data$FP,
  TN = serum_data$TN,
  FN = serum_data$FN,
  design = serum_data$design
)

# Calculate sensitivity and specificity
plot_data <- plot_data %>%
  mutate(
    sensitivity = TP / (TP + FN),
    specificity = TN / (TN + FP),
    FPR = 1 - specificity  # False positive rate
  )

# Create the grouped plot
ggplot(plot_data, aes(x = FPR, y = sensitivity, color = design, shape = design)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  coord_equal() +
  labs(
    x = "False Positive Rate (1 - Specificity)",
    y = "Sensitivity",
    title = "Diagnostic Accuracy by Design Group",
    color = "Design",
    shape = "Design"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  ) +
  scale_color_manual(values = c("0" = "blue", "1" = "red"))



################## graph by antifungal-serum

library(ggplot2)
library(dplyr)

# Extract data from the bamdit object
# First, check what's in your serum object
str(serum)  # This will show available components

# Typical extraction (adjust column names based on your actual data)
plot_data <- data.frame(
  TP = serum$TP,
  FP = serum$FP,
  TN = serum$TN,
  FN = serum$FN,
  antifungal = serum$antifungal  # or serum$antifungal_category, serum$drug, etc.
)

# Calculate sensitivity and specificity
plot_data <- plot_data %>%
  mutate(
    sensitivity = TP / (TP + FN),
    specificity = TN / (TN + FP),
    FPR = 1 - specificity  # False positive rate
  )

# Create the grouped plot by antifungal category
ggplot(plot_data, aes(x = FPR, y = sensitivity, color = antifungal, shape = antifungal)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  coord_equal() +
  labs(
    x = "False Positive Rate (1 - Specificity)",
    y = "Sensitivity",
    title = "Diagnostic Accuracy by Antifungal Category",
    color = "Antifungal",
    shape = "Antifungal"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )


############# Fitting meta-analysis models- no antifungals
library (bamdit)
library (R2jags)

serum_anti0 <- serum %>%
  filter(antifungal == "0")

serum_anti1 <- serum %>%
  filter(antifungal == "1")

serum_anti_b <- serum %>%
  filter(antifungal == 1 | antifungal == 3)

plotdata(serum_anti0, two.by.two = TRUE)

plotdata(serum_anti_b, two.by.two = TRUE)

plotdata(serum_anti1, two.by.two = TRUE)

serum_noantifungal <- metadiag(serum_anti0, two.by.two = TRUE, re = "normal", re.model = "DS",
                     link = "logit", sd.Fisher.rho = 1.7, nr.burnin = 1000,
                     nr.iterations = 10000, nr.chains = 4, r2jags = TRUE)
serum_antifungal <- metadiag(serum_anti0, two.by.two = TRUE, re = "normal", re.model = "DS",
                               link = "logit", sd.Fisher.rho = 1.7, nr.burnin = 1000,
                               nr.iterations = 10000, nr.chains = 4, r2jags = TRUE)

serum_antifungalb <- metadiag(serum_anti_b, two.by.two = TRUE, re = "normal", re.model = "DS",
                             link = "logit", sd.Fisher.rho = 1.7, nr.burnin = 1000,
                             nr.iterations = 10000, nr.chains = 4, r2jags = TRUE)

plot(serum_noantifungal, level = c(0.5, 0.75, 0.95), parametric.smooth = TRUE)
plot(serum_antifungal, level = c(0.5, 0.75, 0.95), parametric.smooth = TRUE)
plot(serum_antifungalb, level = c(0.5, 0.75, 0.95), parametric.smooth = TRUE)
plotsesp(serum_noantifungal)
plotsesp(serum_antifungalb)
plotsesp(serum_antifungal)
bsroc(serum_noantifungal, level = c(0.025, 0.5, 0.975), plot.post.bauc = TRUE,
      fpr.x = seq(0.01, 0.75, 0.01), lower.auc = 0.01, upper.auc = 0.75,
      partial.AUC = FALSE)
bsroc(serum_antifungal, level = c(0.025, 0.5, 0.975), plot.post.bauc = TRUE,
      fpr.x = seq(0.01, 0.75, 0.01), lower.auc = 0.01, upper.auc = 0.75,
      partial.AUC = FALSE)
bsroc(serum_antifungalb, level = c(0.025, 0.5, 0.975), plot.post.bauc = TRUE,
      fpr.x = seq(0.01, 0.75, 0.01), lower.auc = 0.01, upper.auc = 0.75,
      partial.AUC = FALSE)

############# Fitting meta-analysis models- all

library (bamdit)
library (R2jags)



serum_noantifungal <- metadiag(serum_anti0, two.by.two = TRUE, re = "normal", re.model = "DS",
                               link = "logit", sd.Fisher.rho = 1.7, nr.burnin = 1000,
                               nr.iterations = 10000, nr.chains = 4, r2jags = TRUE)
serum_antifungal <- metadiag(serum_anti1, two.by.two = TRUE, re = "normal", re.model = "DS",
                             link = "logit", sd.Fisher.rho = 1.7, nr.burnin = 1000,
                             nr.iterations = 10000, nr.chains = 4, r2jags = TRUE)

plot(serum_noantifungal, level = c(0.5, 0.75, 0.95), parametric.smooth = TRUE)
plot(serum_antifungal, level = c(0.5, 0.75, 0.95), parametric.smooth = TRUE)
plotsesp(serum_noantifungal)
plotsesp(serum_antifungal)
bsroc(serum_noantifungal, level = c(0.025, 0.5, 0.975), plot.post.bauc = TRUE,
      fpr.x = seq(0.01, 0.75, 0.01), lower.auc = 0.01, upper.auc = 0.75,
      partial.AUC = FALSE)
bsroc(serum_antifungal, level = c(0.025, 0.5, 0.975), plot.post.bauc = TRUE,
      fpr.x = seq(0.01, 0.75, 0.01), lower.auc = 0.01, upper.auc = 0.75,
      partial.AUC = FALSE)


