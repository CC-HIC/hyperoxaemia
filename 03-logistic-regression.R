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
    select(unit_mortality, apache_score, prior_dependency, weight,
           sex, spell_los, age, is_medical, surgical_classification,
           tw_hyperoxia_13, cumulative_hyperoxia_13, hyperoxia_13, has_copd,
           system_code, cv, oxy_total) %>%
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
  select(unit_mortality, apache_score, prior_dependency, weight,
         sex, spell_los, age, is_medical, tw_hyperoxia_13, cv)

misum <- misum %>%
  mutate_all(.funs = signif, digits = 2)

misum <- misum %>%
  rename(`time weighted mean hyperoxaemia (dose)` = tw_hyperoxia_13,
         `admission type (surgical/medical)` = is_medical,
         `apache II score` = apache_score,
         `prior dependency` = prior_dependency,
         `ventilated` = cv,
         `unit mortality` = unit_mortality,
         `length of stay` = spell_los)

save(misum, file = "####/missing_data.RData")

for (i in 1:4) {
  temp <- get(paste0("lrd", i))

  temp <- temp %>%
    select(unit_mortality, apache_score, prior_dependency, weight,
           sex, age, is_medical, tw_hyperoxia_13, hyperoxia_13, has_copd,
           system_code, resp, cv, oxy_total) %>%
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

lf <- formula(unit_mortalityL ~ tw_hyperoxia_13 + hyperoxia_13 + resp + rcs(apache_score) + cv +
                                  is_medical + prior_dependency + rcs(weight) +
                                  sex + rcs(age))

# lf_s <- formula(unit_mortalityL ~ tw_hyperoxia_13 + (hyperoxia_13 * (resp + rcs(apache_score) + cv +
#                               is_medical + prior_dependency + rcs(weight) +
#                               sex + rcs(age))))

lf_i <- formula(unit_mortalityL ~ tw_hyperoxia_13 + hyperoxia_13 + resp + rcs(apache_score) +
                                    cv + is_medical + prior_dependency + rcs(weight) +
                                    sex + rcs(age) + hyperoxia_13:resp + hyperoxia_13:cv)


optimal <- formula(unit_mortalityL ~ rcs(oxy_total) + resp + rcs(apache_score) + cv +
                   is_medical + prior_dependency + rcs(weight) +
                   sex + rcs(age))

round_any <- function(x, accuracy, f = round) {
  f(x/accuracy) * accuracy
}

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

  optim_pao2 <- lrm(
    optimal,
    data = current_data,
    eps = 0.005, maxit = 30, x = TRUE, y = TRUE)

  assign(x = paste0("prob_", i), value = temp_prob)
  assign(x = paste0("cll_", i), value = temp_clog)
  assign(x = paste0("logb_", i), value = basic)
  assign(x = paste0("optim_", i), value = optim_pao2)

  # Note the formulation.

  prespec_int <- lrm(
    lf_i, data=current_data,
    eps=0.005, maxit=30)

  current_lr <- lrtest(basic, prespec_int)
  assign(paste0("lr", i), current_lr[[1]])

  prespec_int <- update(prespec_int, x = TRUE, y = TRUE)

  p <- pentrace(prespec_int,
                list(simple = 0,
                     interaction = c(1, 5, 10, 15, 20, 30, 40, 100, 200, 500, 1000,
                                     5000, 10000, 20000, 50000, 100000)),
                maxit = 30)

  assign(x = paste0("penalty_", i), value = p$penalty["interaction"][1,1])

  # Use best penalty

  prespec_int_p <- update(prespec_int,
              penalty = list(
                simple = 0,
                interaction = p$penalty["interaction"][1,1]
              ), maxit = 30, eps = 0.005)

  assign(x = paste0("logip_", i), value = prespec_int_p)
  assign(x = paste0("logi_", i), value = prespec_int)
  current_aic <- c(AIC(basic), AIC(prespec_int), AIC(prespec_int_p))
  names(current_aic) <- c("basic", "interaction", "pen_interaction")

  assign(paste0("aic", i), current_aic)

  double_con_resp <- contrast(
    prespec_int_p,
         list(hyperoxia_13=0, resp='respiratory'),
         list(hyperoxia_13=1, resp='respiratory'),
         list(hyperoxia_13=0, resp='non-respiratory'),
         list(hyperoxia_13=1, resp='non-respiratory'))

  assign(paste0("dbl_con_resp", i), double_con_resp)

  double_con_vent <- contrast(
    prespec_int_p,
    list(hyperoxia_13=0, cv=1),
    list(hyperoxia_13=1, cv=1),
    list(hyperoxia_13=0, cv=0),
    list(hyperoxia_13=1, cv=0))

  assign(paste0("dbl_con_vent", i), double_con_vent)

  ia <- grep('\\*', names(coef(prespec_int_p)))
  x <- coef(prespec_int_p)[-c(1, ia)]
  y <- coef(prespec_int)[-c(1, ia)]
  r <- range(c(x, y))

  png(filename = paste0("./plots/penalty_main_", i, ".png"))
  plot(x, y,
       xlab = "Unpenalized Coefficients",
       ylab = "Penalized Coefficients",
       main = paste("Main Effects: Model", i),
       xlim = r, ylim = r)
  abline(a = 0, b = 1)
  dev.off()

  x <- coef(prespec_int_p)[ia]
  y <- coef(prespec_int)[ia]
  r <- range(c(x, y))

  png(filename = paste0("./plots/penalty_int_", i, ".png"))
  plot(x, y,
       xlab = "Unpenalized Coefficients",
       ylab = "Penalized Coefficients",
       main = paste("Interaction Effects: Model", i),
       xlim = r, ylim = r)
  abline(a = 0, b = 1)
  dev.off()

  d <- as.data.frame(current_data)

  d$tw_hyperoxia_13 <- 0
  d$hyperoxia_13 <- 0L

  p_zero  <- predict(basic, d, type='fitted')
  p_own <- predict(basic, as.data.frame(current_data), type='fitted')

  oxygen_exposure <- current_data %>%
    select(tw_hyperoxia_13, hyperoxia_13) %>%
    mutate(arr = p_own - p_zero) %>%
    mutate(arr = round_any(arr, 0.005),
           p_own = round_any(p_own, 0.005)) %>%
    group_by(arr, p_own) %>%
    tally() %>%
    mutate(bins = cut(n, breaks = c(0, 5, 10, 15, 20, 50, 100, 200, 400, 800, 1600),
                      labels = c(
                        "(  0,   5]",
                        "(  5,  10]",
                        "( 10,  15]",
                        "( 15,  20]",
                        "( 20,  50]",
                        "( 50, 100]",
                        "(100, 200]",
                        "(200, 400]",
                        "(400, 800]",
                        "(800,1600]"
                      )))

  current_plot <- ggplot(data = oxygen_exposure,
         mapping = aes(x = p_own, y = arr, colour = bins)) +
    geom_point() +
    scale_colour_viridis_d(drop = FALSE, guide = FALSE) +
    geom_abline(intercept = 0, slope = 0, linetype = 2) +
    xlim(c(0, 0.3)) +
    ylim(c(-0.06, 0.1)) +
    ylab(label = "") +
    xlab(label = "") +
    theme(axis.text = element_text(size = rel(1.5)))

  assign(paste0("arr_plot", i), current_plot)

  current_ate <- mean(p_own - p_zero)
  assign(paste0("ate", i), current_ate)

  my.valid <- validate(prespec_int_p, method="boot", B=500)
  my.calib <- calibrate(prespec_int_p, method="boot", B=500)

  assign(x = paste0("log_", i, "_valid"), my.valid)
  assign(x = paste0("log_", i, "_calib"), my.calib)

}
rm(i)

## Remeber to remove guide when creating the grid plot

plot2by2 <- plot_grid(arr_plot1, arr_plot2,
                      arr_plot3, arr_plot4, ncol = 2)

save_plot("./plots/arr.png", plot2by2, ncol = 2,
          nrow = 2, base_aspect_ratio = 1.3, dpi = 300)

ggsave(arr_plot1, filename = "./plots/arr_guide.png")

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

## Model Fits and Calibrations
model_fits <- tibble(model = as.character(NULL),
                     c_index = as.numeric(NULL),
                     Brier_score = as.numeric(NULL),
                     AIC = as.numeric(NULL))

for(i in 1:4) {

  temp <- get(paste0("logb_", i))
  corrected_bier <- get(paste0("log_", i, "_valid"))

  additions <- tibble(model = paste0("logb_", i),
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

write_csv(model_fits, path = "./manuscript/model_fits.csv")

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

  # plots
  # temp <- new_data %>%
  #   mutate(pred = predict(get(paste0("log_", i)), newdata = new_data, type = "response")) %>%
  #   ggplot(aes(x = tw_hyperoxia_13, y = pred)) + geom_point() +
  #   geom_smooth(method = "loess") +
  #   xlab(label = "time weighted mean hyperoxia") +
  #   ylab(label = "predicted probability of death") +
  #   ylim(c(0, 0.1)) +
  #   xlim(c(0, 5))
  #
  # assign(x = paste0("indv_", i), value = temp)

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

save(ate_all, file = "./manuscript/ate.RData")

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

save_plot("./plots/prroc.png", plot2by2, ncol = 2,
          nrow = 1, base_aspect_ratio = 1.4, dpi = 300)

clean_up <- function(model, ddata){
  mod <- model
  dd <- datadist(ddata); options(datadist = "dd")

  sum_mod <- broom::tidy(summary(mod))
  anova_mod <- broom::tidy(anova(mod))

  val_names <- sum_mod %>%
    filter(Type == 1) %>%
    select(.rownames) %>%
    pull() %>%
    str_extract("\\w*")

  clean_base <- sum_mod %>%
    filter(Type == 2) %>%
    mutate(.rownames = val_names) %>%
    select(.rownames, Effect, Lower.0.95, Upper.0.95) %>%
    rename(variable = .rownames,
           odds_ratio = Effect,
           lower_95 = Lower.0.95,
           upper_95 = Upper.0.95)

  result <- anova_mod %>%
    mutate(variable = case_when(
      grepl("[*]", .rownames) ~ .rownames,
      grepl("Nonlinear", .rownames) ~ "x",
      grepl("All Interactions", .rownames) ~ "x",
      grepl("TOTAL", .rownames) ~ "x",
      TRUE ~ str_extract(.rownames, "\\w*"))) %>%
    filter(variable != "x") %>%
    select(-.rownames) %>%
    full_join(clean_base, .)

  return(result)
}

clean_1 <- clean_up(logb_1, lrd1)
clean_2 <- clean_up(logb_2, lrd2)
clean_3 <- clean_up(logb_3, lrd3)
clean_4 <- clean_up(logb_4, lrd4)

master_results <- bind_rows(
  clean_1 %>%
    mutate(subset = "up to 1 day"),
  clean_2 %>%
    mutate(subset = "up to 3 days"),
  clean_3 %>%
    mutate(subset = "up to 5 days"),
  clean_4 %>%
    mutate(subset = "up to 7 days")
)

save(master_results, file = "./manuscript/master_results.RData")

my_plot <- master_results %>%
    filter(variable != "age") %>%
    filter(variable != "weight") %>%
    filter(variable != "apache_score") %>%
    filter(!grepl("[*]", variable)) %>%
    mutate(offset = case_when(
      subset == "up to 1 day" ~ -0.6,
      subset == "up to 3 days" ~ -0.2,
      subset == "up to 5 days" ~ 0.2,
      subset == "up to 7 days" ~ 0.6)) %>%
    mutate(yaxis = rep(rev(seq(from = 1, by = 2, length.out = 7)), 4)) %>%
    mutate(yplace = offset+yaxis) %>%
    ggplot(aes(y = yplace, colour = factor(subset))) +
    geom_point(aes(x = odds_ratio), size = 2) +
    geom_segment(aes(x = lower_95, xend = upper_95, y = yplace, yend = yplace), size = 1) +
    scale_y_discrete(limits = seq(from = 1, by = 2, length.out = 12),
                     labels = rev(c(
                       "dose dependent hyperoxemia",
                       "dose independent hyperoxemia",
                       "continuous ventilation",
                       "respiratory diagnosis - non-resp:resp",
                       "admission type - surgical:medical",
                       "prior dependency - none:any",
                       "sex - male:female"))) +
    scale_x_log10(limits = c(0.4, 5)) +
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
  geom_smooth(method = "loess", se = FALSE) +
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

save(lr1, lr2, lr3, lr4, file = "./manuscript/likelihoods.RData")
save(dbl_con_resp1, dbl_con_resp2, dbl_con_resp3, dbl_con_resp4,
     dbl_con_vent1, dbl_con_vent2, dbl_con_vent3, dbl_con_vent4,
     file = "./manuscript/contrasts.RData")
save(aic1, aic2, aic3, aic4, file = "./manuscript/aics.RData")

save.image("./data/end_of_experiment.RData")

## :)
