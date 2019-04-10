library(tidyverse)
library(rms)
library(cowplot)
library(ROCR)
library(PRROC)

setwd("####/oxygen/")
load("####/cohorts_final.RData")

# Use what we need
for (i in 1:4) {
  temp <- get(paste0("cohort_", i))
  
  temp <- temp %>%
    select(unit_mortality, apache_score, prior_dependency, weight, cpr,
           sex, spell_los, age, is_medical, surgical_classification,
           tw_hyperoxia_13, cumulative_hyperoxia_13, hyperoxia_13, has_copd,
           system_code) %>%
    mutate(resp = factor(
      if_else(system_code == "01", "respiratory", "non-respiratory")))
  
  assign(x = paste0("lrd", i), value = temp)
}
rm(i, temp)
# check how much missingness we have here

misum <- bind_rows(
  lrd1 %>%
  summarise_all(function(x) sum(is.na(x))/length(x)) %>%
    mutate(subset = "day 1"),
  lrd2 %>%
    summarise_all(function(x) sum(is.na(x))/length(x)) %>%
    mutate(subset = "day 3"),
  lrd3 %>%
    summarise_all(function(x) sum(is.na(x))/length(x)) %>%
    mutate(subset = "day 5"),
  lrd4 %>%
    summarise_all(function(x) sum(is.na(x))/length(x)) %>%
    mutate(subset = "day 7")
)

misum <- misum %>%
  select(unit_mortality, apache_score, prior_dependency, weight, cpr,
         sex, spell_los, age, is_medical, tw_hyperoxia_13)

misum <- misum %>%
  mutate_all(.funs = signif, digits = 2)

misum <- misum %>%
  rename(`time weighted mean hyperoxaemia (dose)` = tw_hyperoxia_13,
         `admission type (surgical/medical)` = is_medical,
         `apache II score` = apache_score,
         `prior dependency` = prior_dependency,
         `CPR status` = cpr,
         `unit mortality` = unit_mortality,
         `length of stay` = spell_los)

save(misum, file = "####/missing_data.RData")

for (i in 1:4) {
  temp <- get(paste0("lrd", i))
  
  temp <- temp %>%
    select(unit_mortality, apache_score, prior_dependency, weight, cpr,
           sex, age, is_medical, tw_hyperoxia_13, hyperoxia_13, has_copd,
           system_code, resp) %>%
    na.omit() %>%
    mutate(unit_mortalityL = case_when(
      unit_mortality == "alive" ~ 0L,
      unit_mortality == "deceased" ~ 1L,
      TRUE ~ as.integer(NA)
    ))
  
  assign(x = paste0("lrd", i), value = temp)
}
rm(i, temp)

consort_info <- bind_rows(consort_info, tribble(
  ~stage, ~n,
  "los >= 1 day (no missing)", nrow(lrd1),
  "los >= 3 days (no missing)", nrow(lrd2),
  "los >= 5 days (no missing)", nrow(lrd3),
  "los >= 7 days (no missing)", nrow(lrd4)
))

save(consort_info, file = "####/consort.RData")

## There is a trivial amount of missing data here
## getting 1% back isn't going to change anything

lf <- formula(unit_mortalityL ~ tw_hyperoxia_13 + hyperoxia_13 + resp + apache_score +
                                  is_medical + prior_dependency + rcs(weight) + cpr +
                                  sex + rcs(age))

lf_s <- formula(unit_mortalityL ~ tw_hyperoxia_13 + (hyperoxia_13 * (resp + apache_score +
                              is_medical + prior_dependency + rcs(weight) + cpr +
                              sex + rcs(age))))

lf_i <- formula(unit_mortalityL ~ tw_hyperoxia_13 + (hyperoxia_13 * resp) + apache_score +
                                   is_medical + prior_dependency + rcs(weight) + cpr +
                                 sex + rcs(age))

for(i in 1:4) {
  
  dd <- datadist(get(paste0("lrd", i))); options(datadist = "dd")
  current_data <- get(paste0("lrd", i))

  basic <- lrm(
    lf,
    data = current_data,
    eps = 0.005, maxit = 30, x = TRUE, y = TRUE)

  temp_prob <- Glm(
    formula = lf,
    family = binomial(link = "probit"),
    data = get(paste0("lrd", i)), x = TRUE, y = TRUE)

  temp_clog <- Glm(
    formula = lf,
    family = binomial(link = "cloglog"),
    data = get(paste0("lrd", i)), x = TRUE, y = TRUE)

  assign(x = paste0("prob_", i), value = temp_prob)
  assign(x = paste0("cll_", i), value = temp_clog)
  assign(x = paste0("logb_", i), value = basic)
  
  resp_int <- lrm(
    lf_i, data=current_data,
    eps=0.005, maxit=30)

  assign(x = paste0("logi_", i), value = resp_int)

  current_lr <- lrtest(basic, resp_int)
  assign(paste0("lr", i), current_lr[[1]])
  
  current_aic <- c(AIC(basic), AIC(resp_int))
  assign(paste0("aic", i), current_aic)

  double_con <- contrast(
    resp_int,
         list(hyperoxia_13=0, resp='respiratory'),
         list(hyperoxia_13=1, resp='respiratory'),
         list(hyperoxia_13=0, resp='non-respiratory'),
         list(hyperoxia_13=1, resp='non-respiratory'))

  assign(paste0("dbl_con", i), double_con)
  
  w <- with(current_data, c(sum(resp=='respiratory'), sum(resp=='non-respiratory')))

  single_con <- contrast(
    resp_int,
         list(hyperoxia_13=0, resp=c('respiratory', 'non-respiratory')),
         list(hyperoxia_13=1, resp=c('respiratory', 'non-respiratory')),
         type='average',
         weights=w)
  
  assign(paste0("sgl_con", i), single_con)

  d <- as.data.frame(current_data)
  ## Set all exposure to zero
  
  d$tw_hyperoxia_13 <- 0
  d$hyperoxia_13 <- 0L
  
  p_zero  <- predict(basic, d, type='fitted')

  ## Set all exposure to own
  p_own <- predict(basic, as.data.frame(current_data), type='fitted')

  current_ate <- mean(p_own - p_zero)
  assign(paste0("ate", i), current_ate)

  ## Set cut points
  k  <- c(5,10,15,20,50,100,200,400,800,1600,3200)

  pl <- c(0, quantile(p_zero, 0.99))
  ##
  j  <- pmax(p_zero, p_own) <= pl

  current_plot <- ggfreqScatter(
    x = p_zero[j],
    y = p_own[j],
    cuts = k,
    xlab = 'Predicted mortality with zero exposure',
    ylab = 'Predicted mortality with patient\'s own exposure') +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(0, 0.3)) +
  scale_y_continuous(limits = c(0, 0.3))

  assign(paste0("pplot", i), current_plot)

  xl <- 'Risk Under Patient\'s Own Exposure'
  j <- p_own <= pl

  left_arr <- ggfreqScatter(p_own[j], p_own[j] - p_zero[j], cuts=k, xlab=xl, ylab='ARR with no exposure') +
    scale_x_continuous(limits = c(0, 0.3))

  right_rr <- ggfreqScatter(p_own[j], p_zero[j] / p_own[j], cuts=k, xlab=xl, ylab='RR with no exposure') +
    scale_x_continuous(limits = c(0, 0.3))

  assign(paste0("left_arr", i), left_arr)
  assign(paste0("right_rr", i), right_rr)

  my.valid <- validate(basic, method="boot", B=500)
  my.calib <- calibrate(basic, method="boot", B=500)
  
  assign(x = paste0("log_", i, "_valid"), my.valid)
  assign(x = paste0("log_", i, "_calib"), my.calib)

}
rm(i)

exposures <- tibble(hyperoxia_dose = as.numeric(NULL),
                    prediction = as.numeric(NULL),
                    subset = as.character(NULL))

## Patient Level Exposure Predictions
for(i in 1:4) {
  
  temp_data <- get(paste0("lrd", i))
  temp_model_13 <- get(paste0("logb_", i))
  
  current_exposures <- temp_data %>%
    select(tw_hyperoxia_13) %>%
    rename(hyperoxia_dose = tw_hyperoxia_13) %>%
    mutate(prediction = predict(temp_model_13, type = "fitted"),
           subset = case_when(
             i == 1 ~ "Day 1",
             i == 2 ~ "Day 3",
             i == 3 ~ "Day 5",
             i == 4 ~ "Day 7"
           ))
  
  exposures <- bind_rows(exposures, current_exposures)
  
}

## Scatter plots
exposures %>%
  ggplot(aes(hyperoxia_dose, prediction*100)) +
  geom_point(alpha = 0.5, shape = 1) +
  facet_wrap(~subset, ncol = 2) +
  geom_smooth(method = "loess") +
  scale_y_continuous(limits = c(0, 100)) +
  xlab("time weighted mean hyperoxemia (above 13.3 kPa)") +
  ylab("model predicted probability of hospital mortality (%)")

ggsave(filename = "####/individual_predictions.png",
       width = 10, height = 8, dpi = 300, units = "in")

## Model Fits and Calibrations
model_fits <- tibble(model = as.character(NULL),
                     c_index = as.numeric(NULL),
                     Brier_score = as.numeric(NULL),
                     AIC = as.numeric(NULL))

for(i in 1:4) {
  
  temp <- get(paste0("logb_", i))
  corrected_bier <- get(paste0("log_", i, "_valid"))
  
  additions <- tibble(model = paste0("log_", i),
                      c_index = temp$stats[["C"]],
                      Brier_score = corrected_bier["B","index.corrected"],
                      AIC = AIC(temp))
  
  model_fits <- bind_rows(model_fits, additions)
  
  rm(temp, additions, corrected_bier)
}
rm(i)

model_fits <- model_fits %>%
  mutate(model = paste0("up to day ", seq(from = 1, by = 2, length.out = 4)))

model_fits <- mutate_at(model_fits, vars(c_index:AIC), signif, digits = 3)
names(model_fits) <- c("model", "c-index", "Brier score", "AIC")

write_csv(model_fits, path = "####/model_fits.csv")

## convert the logit link function to probability predictions
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- (odds / (1 + odds))*100
  return(prob)
}

## Average treatment Effect
ate_all <- tibble(group = as.integer(NULL),
                  o2_on = as.numeric(NULL),
                  o2_off = as.numeric(NULL),
                  ate = as.numeric(NULL))

for (i in 1:4) {
  
  current_model <- get(paste0("logb_", i))
  current_data <- get(paste0("lrd", i))
  
  new_data <- current_data %>%
    as.data.frame() # predict doesn't work well with tibbles

  # set O2 to zero for all patients (counterfactual)
  o2_off_raw <- new_data %>%
    mutate(tw_hyperoxia_13 = 0, hyperoxia_13 = 0L)

  # model predicted outcomes
  o2_on_pred <- mutate(new_data, pred = predict(current_model, newdata = new_data, type = "fitted"))
  o2_off_pred <- mutate(new_data, pred = predict(current_model, newdata = o2_off_raw, type = "fitted"))

  
  # average over all patients for ATE
  ate <- bind_cols(
    o2_on_pred %>% select(pred) %>% rename(o2_on = pred),
    o2_off_pred %>% select(pred) %>% rename(o2_off = pred)) %>%
    summarise(o2_on = mean(o2_on),
              o2_off = mean(o2_off)) %>%
    mutate(group = as.integer(i))
  
  # There was some buggy behaviour with dplyr here
  ate$ate <- as.numeric(ate$o2_on) - as.numeric(ate$o2_off)
  ate_all <- bind_rows(ate_all, ate)
  
  ## ROC/PR curves
  
  roc_pred <- prediction(as.numeric(o2_on_pred$pred), as.integer(o2_on_pred$unit_mortality) - 1)
  roc_perf <- performance(roc_pred, measure = "tpr", x.measure = "fpr")
  
  roct <- tibble(fpr = roc_perf@x.values[[1]],
                 tpr = roc_perf@y.values[[1]]) %>%
    mutate(group = i)
  
  roc <- roc.curve(scores.class0 = as.numeric(o2_on_pred$pred[o2_on_pred$unit_mortality == "deceased"]),
                   scores.class1 = as.numeric(o2_on_pred$pred[o2_on_pred$unit_mortality == "alive"]),
                   curve = TRUE)
  
  pr <- pr.curve(scores.class0 = as.numeric(o2_on_pred$pred[o2_on_pred$unit_mortality == "deceased"]),
                 scores.class1 = as.numeric(o2_on_pred$pred[o2_on_pred$unit_mortality == "alive"]),
                 curve = TRUE)
  
  assign(x = paste0("roct", i), value = roct)
  assign(x = paste0("roc", i), value = roc)
  assign(x = paste0("pr", i), value = pr)
  assign(x = paste0("auc", i), value = roc$auc)
  
  rm(ate, new_data, o2_off_raw, o2_off_pred, o2_on_pred)
}
rm(i)

save(ate_all, file = "####/ate.RData")

roc_all <- reduce(list(
    tibble(logx = roc1$curve[,1],
           logy = roc1$curve[,2],
           model = "up to day 1"),
    tibble(logx = roc2$curve[,1],
           logy = roc2$curve[,2],
           model = "up to day 3"),
    tibble(logx = roc3$curve[,1],
           logy = roc3$curve[,2],
           model = "up to day 5"),
    tibble(logx = roc4$curve[,1],
           logy = roc4$curve[,2],
           model = "up to day 7")),
    bind_rows)

pr_all <- reduce(list(
  tibble(logx = pr1$curve[,1],
         logy = pr1$curve[,2],
         model = "up to day 1"),
  tibble(logx = pr2$curve[,1],
         logy = pr2$curve[,2],
         model = "up to day 3"),
  tibble(logx = pr3$curve[,1],
         logy = pr3$curve[,2],
         model = "up to day 5"),
  tibble(logx = pr4$curve[,1],
         logy = pr4$curve[,2],
         model = "up to day 7")),
  bind_rows)

roc_all_plot <- roc_all %>%
  ggplot(aes(x = logx, y = logy, group = model, colour = model)) +
  geom_path() +
  xlab("False Positive Rate") +
  ylab("True Positive Rate (Sensitivity)")

pr_all_plot <- pr_all %>%
  ggplot(aes(x = logx, y = logy, group = model, colour = model)) +
  geom_path() +
  xlab("Recall") +
  ylab("Precision")

plot2by2 <- plot_grid(roc_all_plot, pr_all_plot, ncol = 2)

save_plot("####/prroc.png", plot2by2, ncol = 2,
          nrow = 1, base_aspect_ratio = 1.4, dpi = 300)


for (i in 1:4) {
  for (j in c("left_arr", "right_rr")) {
    ggsave(filename = paste0("./plots/", j, "_", i, ".png"),
           plot = get(paste0(j, i)),
           dpi = 300, width = 8, height = 8, units = "in")
  }
}

for (i in 1:4) {
    ggsave(filename = paste0("./plots/", "cf_", i, ".png"),
           plot = get(paste0("pplot", i)),
           dpi = 300, width = 8, height = 8, units = "in")
  }

## Cleaning up and getting co-efficients for the manuscript
clean_up <- function(mod) {
  
  ## It is quite difficult to extract specific data items from lrm objects
  ## Broom doesn't have a good solution for this yet
  
  these_names <- attr(summary(mod), "dimnames")[[1]]
  these_names <- these_names[1:length(these_names) %% 2 == 1]
  these_names <- rep(these_names, each = 2)
  
  sum_all <- summary(mod)
  
  stb <- tibble(predictor_variables = these_names)
  stb <- bind_cols(stb, tibble(effect = sum_all[,"Effect"],
                               lower.conf = sum_all[,"Lower 0.95"],
                               upper.conf = sum_all[,"Upper 0.95"],
                               type = sum_all[,"Type"])) %>%
    filter(type == 2) %>%
    mutate(predictor_variables = gsub(
      pattern = " [A-Za-z:-]*",
      replacement = "",
      x = predictor_variables))
  
  anova_names <- attr(anova(mod), "dimnames")[[1]]
  
  an_table <- tibble(predictor_variables = anova_names,
                     p_value = anova(mod)[,"P"])
  
  stb <- full_join(stb, an_table, by = "predictor_variables") %>%
    select(-type) %>%
    `[`(1:10,)
  
  return(stb)

}

clean_1 <- clean_up(logb_1)
clean_2 <- clean_up(logb_2)
clean_3 <- clean_up(logb_3)
clean_4 <- clean_up(logb_4)

# Primary
for(i in 1:4) {
  
  if(exists("table_names")) rm(table_names)
  if(exists("temp")) rm(temp)
  
  table_names <- dimnames(summary(get(paste0("logb_", i))))[[1]] %>%
    as.tibble() %>%
    mutate(odds = if_else(row_number() %% 2 == 1, "odd", "even")) %>%
    dplyr::filter(odds == "odd") %>%
    select(value) %>%
    pull() %>%
    rep(each = 2)
  
  temp <- summary(get(paste0("logb_", i))) %>%
    as.data.frame(row.names = FALSE) %>%
    as.tibble() %>%
    mutate(variable = table_names) %>%
    filter(Type == 2) %>%
    rename(low = Low,
           high = High,
           diff = `Diff.`,
           odds_ratio = Effect,
           se = `S.E.`,
           lower_95 = `Lower 0.95`,
           upper_95 = `Upper 0.95`,
           type = Type) %>%
    select(variable, odds_ratio, lower_95, upper_95)
  
  assign(value = temp, x = paste0("log_", i, "_clean"))
  
}

log_prim_clean <- bind_rows(
  log_1_clean %>%
    mutate(subset = "up to 1 day"),
  log_2_clean %>%
    mutate(subset = "up to 3 days"),
  log_3_clean %>%
    mutate(subset = "up to 5 days"),
  log_4_clean %>%
    mutate(subset = "up to 7 days")
)

master_results <- bind_cols(log_prim_clean, bind_rows(select(clean_1, p_value),
                                         select(clean_2, p_value),
                                         select(clean_3, p_value),
                                         select(clean_4, p_value)))

save(master_results, file = "####/master_results.RData")

my_plot <- log_prim_clean %>%
    filter(variable != "age") %>%
    filter(variable != "weight") %>%
    mutate(offset = case_when(
      subset == "up to 1 day" ~ -0.6,
      subset == "up to 3 days" ~ -0.2,
      subset == "up to 5 days" ~ 0.2,
      subset == "up to 7 days" ~ 0.6)) %>%
    mutate(yaxis = rep(rev(seq(from = 1, by = 2, length.out = 8)), 4)) %>%
    mutate(yplace = offset+yaxis) %>%
    ggplot(aes(y = yplace, colour = factor(subset))) +
    geom_point(aes(x = odds_ratio), size = 2) +
    geom_segment(aes(x = lower_95, xend = upper_95, y = yplace, yend = yplace), size = 1) +
    scale_y_discrete(limits = seq(from = 1, by = 2, length.out = 12),
                     labels = rev(c(
                       "dose dependent hyperoxaemia",
                       "dose independent hyperoxaemia",
                       "APACHE II score",
                       "respiratory diagnosis - non-resp:resp",
                       "admission type - surgical:medical",
                       "prior dependency - none:any",
                       "CPR status - no:yes",
                       "sex - male:female"))) +
    scale_x_log10(limits = c(0.4, 4)) +
    annotation_logticks(sides = "b") +
    geom_vline(xintercept = 1, linetype = 2) +
    theme_classic() +
    theme(axis.text.x = element_text(size = rel(1.5)),
          axis.title.x = element_text(size = rel(1.5)),
          axis.text.y = element_text(size = rel(1.5)),
          axis.title.y = element_text(size = rel(1.5)),
          legend.title.align = 0.5) +
    xlab("odds ratio") +
    ylab("independent (predictor) variables") +
    guides(colour = guide_legend(title = "subset"))
  
ggsave(my_plot, filename = "./plots/or.png",
         width = 10, height = 8, dpi = 300, units = "in")

global_calib <- bind_rows(
  log_1_calib[1:50, 1:9] %>% as.tibble() %>% mutate(subset = "day 1 subset"),
  log_2_calib[1:50, 1:9] %>% as.tibble() %>% mutate(subset = "day 3 subset"),
  log_3_calib[1:50, 1:9] %>% as.tibble() %>% mutate(subset = "day 5 subset"),
  log_4_calib[1:50, 1:9] %>% as.tibble() %>% mutate(subset = "day 7 subset")
)

global_calib %>%
  select(predy, calibrated.orig, calibrated.corrected, subset) %>%
  rename(predicted_probability = predy,
         probability_orig = calibrated.orig,
         probability_correct = calibrated.corrected) %>%
  gather(key = prediction_type, value = actual_probability, probability_orig:probability_correct) %>%
  filter(prediction_type != "probability_orig") %>%
  ggplot(aes(x = predicted_probability, y = actual_probability, colour = factor(subset))) +
  geom_smooth(alpha = 0.7, method = "loess") +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(1.5)),
        legend.title.align = 0.5) +
  xlab("predicted probability") +
  ylab("actual probability") +
  guides(colour = guide_legend(title = "subset"))

ggsave(filename = "./plots/calibration.png", width = 10, height = 8, dpi = 300, units = "in")

save.image("####/end_of_experiment.RData")

## :)
