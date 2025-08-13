set.seed(12345)
library (bamdit)
library (ggplot2)
serum <- read_csv("~/Desktop/serum.csv", col_types = cols(TP= col_integer(), 
                                         FP = col_integer(), FN = col_integer(), 
                                        TN = col_integer(), year = col_character(), 
                                      design = col_factor(levels = c("0", "1")), patients = col_integer(), 
                                        serum = col_factor(levels = c("0", "1")),
                                      antifungal = col_factor(levels = c("0", "1", "3")), 
                                      cutoff = col_factor(levels = c("0.2", "0.5", "1.5", "1", "0.8"))))


library(mada)


# Assuming your data has TP, FP, TN, FN columns
# First, create a mada-compatible data frame
mada_data <- data.frame(
  TP = serum$TP,
  FP = serum$FP,
  TN = serum$TN,
  FN = serum$FN,
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


plotdata(serum, two.by.two = TRUE)                                                                                                                                        

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


############# Fitting meta-analysis models

serum_anti0 <- serum %>%
  filter(antifungal == "0")

serum_anti1 <- serum %>%
  filter(antifungal == "1")

plotdata(serum_anti0, two.by.two = TRUE)

plotdata(serum_anti1, two.by.two = TRUE)

