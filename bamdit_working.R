library(dplyr)
library(bamdit)
library(ggplot2)
data("glas")
head(glas)


raw <- readxl::read_excel("data.xlsx")
dat0 <- raw %>%
transmute(
study = as.character(author),
tp = as.integer(tp),
     fn = as.integer(n1),
     fp = as.integer(fp),
    tn = as.integer(n2),
    bal = as.integer(bal)
   )
 diag_check <- dat0 %>%
   mutate(
     diseased_total    = tp + fn,
     nondiseased_total = fp + tn,
     reason = dplyr::case_when(
       diseased_total == 0    ~ "no diseased subjects",
       nondiseased_total == 0 ~ "no non-diseased subjects",
     )
   )
good <- diag_check %>% filter(is.na(reason)) %>% select(study, tp, fn, fp, tn)


dat_bamdit <- good %>%
  transmute(
    tp = as.numeric(tp),
    n1 = as.numeric(tp + fn),   # total diseased
    fp = as.numeric(fp),
    n2 = as.numeric(fp + tn),
    # total non-diseased
  ) %>%
  as.data.frame(stringsAsFactors = FALSE)

## 2) Guardrails: no list-cols, all numeric, positive totals
stopifnot(!any(vapply(dat_bamdit, is.list, logical(1))))
stopifnot(all(vapply(dat_bamdit, is.numeric, logical(1))))
dat_bamdit <- subset(dat_bamdit, is.finite(tp) & is.finite(fp) & is.finite(n1) & is.finite(n2))
dat_bamdit <- subset(dat_bamdit, n1 > 0 & n2 > 0 & tp >= 0 & fp >= 0)

## 3) Sanity check the rates we’ll plot
df_plot <- transform(dat_bamdit,
                     TPR = tp / n1,
                     FPR = fp / n2,
                     N   = n1 + n2)
stopifnot(all(is.finite(df_plot$TPR)), all(is.finite(df_plot$FPR)))



## 5) Now Bamdit’s plot (should run cleanly)
plotdata(dat_bamdit)
plotdata(dat_bamdit)

dat_bamdit <- metadiag(dat_bamdit, re = "normal", re.model = "DS",
                     link = "logit", sd.Fisher.rho = 1.7, nr.burnin = 1000,
                    nr.iterations = 10000, nr.chains = 4, r2jags = TRUE)


summary(dat_bamdit, digits = 3)

library("R2jags")
attach.jags(dat_bamdit)
cor(se.pool, sp.pool)

plot(dat_bamdit, level = c(0.5, 0.75, 0.95), parametric.smooth = TRUE)
bsroc(dat_bamdit, level = c(0.025, 0.5, 0.975), plot.post.bauc = TRUE,
       fpr.x = seq(0.01, 0.75, 0.01), lower.auc = 0.01, upper.auc = 0.75,
      partial.AUC = FALSE)
plotsesp(dat_bamdit, binwidth.p = 0.03, CI.level = 0.95)
plotw(m = dat_bamdit)
