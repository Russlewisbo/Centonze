set.seed(12345)
library (bamdit)
library (ggplot2)


bal <- read_csv("~/Desktop/bal.csv", col_types = cols(TP= col_integer(), 
                                                          FP = col_integer(), FN = col_integer(), 
                                                          TN = col_integer(), year = col_character(), 
                                                          design = col_factor(levels = c("0", "1")), patients = col_integer(), 
                                                          bal = col_factor(levels = c("0", "1")),
                                                          antifungal = col_factor(levels = c("0", "1", "3")), 
                                                          cutoff = col_factor(levels = c("0.5", "1"))))


problems (bal)
# Assuming your data has TP, FP, TN, FN columns
# First, create a mada-compatible data frame
mada_data <- data.frame(
  TP = bal$TP,
  FP = bal$FP,
  TN = bal$TN,
  FN = bal$FN,
  study = "serum$author"  # or serum$author, serum$study_id, etc.
)

# Calculate descriptive statistics
madad_results <- madad(mada_data)

# Create forest plots
forest(madad_results, type = "sens", cex.axis = 0.6,    # X-axis numbers
       cex.lab = 0.7,     # X-axis title
       cex.main = 0.8,    # Main title
       cex = 0.4)         # Study names on Y-axis)  # Sensitivity forest plot
forest(madad_results, type = "spec", cex.axis = 0.6,    # X-axis numbers
       cex.lab = 0.7,     # X-axis title
       cex.main = 0.8,    # Main title
       cex = 0.4)         # Study names on Y-axis)  # Sensitivity forest plot)  # Specificity forest plot


plotdata(, two.by.two = TRUE)                                                                                                                                        

## ggplot workaround-serum- design

library(ggplot2)
library(dplyr)

# Extract data from the bamdit object
# First, check what's in your serum object
str(serum)  # This will show available components

# Typical extraction (adjust based on your actual data structure)
# Method 1: If serum has TP, FP, TN, FN columns
plot_data <- data.frame(
  TP = serum$TP,
  FP = serum$FP,
  TN = serum$TN,
  FN = serum$FN,
  design = serum$design
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
str(bal)  # This will show available components

# Typical extraction (adjust column names based on your actual data)
plot_data <- data.frame(
  TP = bal$TP,
  FP = bal$FP,
  TN = bal$TN,
  FN = bal$FN,
  antifungal = bal$antifungal  # or serum$antifungal_category, serum$drug, etc.
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

bal_anti0 <- bal %>%
  filter(antifungal == "0")

bal_anti1 <- bal %>%
  filter(antifungal == "1")

bal_anti_b <- bal %>%
  filter(antifungal == 1 | antifungal == 3)

plotdata(bal_anti0, two.by.two = TRUE)

plotdata(bal_anti_b, two.by.two = TRUE)

plotdata(bal_anti1, two.by.two = TRUE)

bal_noantifungal <- metadiag(bal_anti0, two.by.two = TRUE, re = "normal", re.model = "DS",
                               link = "logit", sd.Fisher.rho = 1.7, nr.burnin = 1000,
                               nr.iterations = 10000, nr.chains = 4, r2jags = TRUE)
bal_antifungal <- metadiag(bal_anti1, two.by.two = TRUE, re = "normal", re.model = "DS",
                             link = "logit", sd.Fisher.rho = 1.7, nr.burnin = 1000,
                             nr.iterations = 10000, nr.chains = 4, r2jags = TRUE)

bal_antifungalb <- metadiag(bal_anti_b, two.by.two = TRUE, re = "normal", re.model = "DS",
                              link = "logit", sd.Fisher.rho = 1.7, nr.burnin = 1000,
                              nr.iterations = 10000, nr.chains = 4, r2jags = TRUE)

plot(bal_noantifungal, level = c(0.5, 0.75, 0.95), parametric.smooth = TRUE)
plot(bal_antifungal, level = c(0.5, 0.75, 0.95), parametric.smooth = TRUE)
plot(bal_antifungalb, level = c(0.5, 0.75, 0.95), parametric.smooth = TRUE)
plotsesp(bal_noantifungal)
plotsesp(bal_antifungalb)
plotsesp(bal_antifungal)
bsroc(bal_noantifungal, level = c(0.025, 0.5, 0.975), plot.post.bauc = TRUE,
      fpr.x = seq(0.01, 0.75, 0.01), lower.auc = 0.01, upper.auc = 0.75,
      partial.AUC = FALSE)
bsroc(bal_antifungal, level = c(0.025, 0.5, 0.975), plot.post.bauc = TRUE,
      fpr.x = seq(0.01, 0.75, 0.01), lower.auc = 0.01, upper.auc = 0.75,
      partial.AUC = FALSE)
bsroc(bal_antifungalb, level = c(0.025, 0.5, 0.975), plot.post.bauc = TRUE,
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

