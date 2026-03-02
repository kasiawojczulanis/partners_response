# Script for the study published in XXX

# Title: "Limited capacity for parental compensation under experimental handicap in a small pelagic seabird, the Little Auk" 
# Authors: Katarzyna Wojczulanis-Jakubas1*, Pauline Bodson1, Jerome Fort2, David Gremillet3,4, Andrea Grunst5, Melissa Grunst5, Kristin Piening1, Martyna Syposz1, Dariusz Jakubas1

# Study Question # 2 
# Do behavioural effects persist after the burden is removed? 
  
# Libs ----

library(tidyverse)
library(ggdist)
library(ggalluvial)

library(brms)
library(bayesplot)
library(bayestestR)
library(tidybayes)

library(emmeans)

library(patchwork)

# ---- Read data ----------------------------------

nfeeds_chr2_df <- readRDS("./Data/nfeeds_chr2_df.rds")




# ---- Sample size ----------------------------------

# examine sample size
nfeeds_chr2_df %>%
  distinct(season, nest, nest_device_stat) %>%
  count(season, nest_device_stat)



# ---- Reshape data for modelling -------------------------------
# top = higher-feeding partner
# bottom = lower-feeding partner

# controls
partners_control <- nfeeds_chr2_df %>%
  filter(nest_device_stat == "control") %>%
  select(-nest_device_stat, -partner_device_stat) %>%
  pivot_wider(
    names_from = sx,
    values_from = nfeeds,
    names_prefix = "nfeeds_"
  ) %>%
  mutate(
    nfeeds_f = replace_na(nfeeds_f, 0),
    top    = pmax(nfeeds_f, nfeeds_m),
    bottom = pmin(nfeeds_f, nfeeds_m),
    group  = "control"
  )

# deployment
partners_device <- nfeeds_chr2_df %>%
  filter(nest_device_stat %in% c("acc", "gps")) %>%
  select(-nest_device_stat, -sx) %>%
  pivot_wider(
    names_from = partner_device_stat,
    values_from = nfeeds,
    names_prefix = "nfeeds_"
  ) %>%
  mutate(
    top    = pmax(nfeeds_deploy, nfeeds_partner),
    bottom = pmin(nfeeds_deploy, nfeeds_partner),
    group  = "deployment"
  )


# combine control-deployment data sets
top_bottom_chr2 <- bind_rows(
  partners_control  %>% select(season, nest, top, bottom, group),
  partners_device %>% select(season, nest, top, bottom, group)
  
)


# ---- Analysis #1 -------------------------------
# Comparison of provisioning rates of higher‑ and lower‑feeding parents between experimental and control groups.

# prepare data for the model
top_bottom_all_long <- top_bottom_chr2 %>%
  pivot_longer(names_to = "parent", values_to = "nfeeds", cols = -c(season, nest, group)) 


# priors ----
priors <- c(
  set_prior("normal(2, 0.3)", class = "Intercept"),
  set_prior("normal(0, 0.2)", class = "b"),
  set_prior("exponential(5)", class = "sd", group = "season"),
  set_prior("exponential(5)", class = "sd", group = "nest")
)

model_prior_nfeeds <- brm(
  formula = nfeeds ~ parent * group + (1 | season) + (1|nest),
  data = top_bottom_all_long,
  family = poisson(link = "log"),
  prior = priors,
  sample_prior = "only",
  seed = 1313,
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.995, max_treedepth = 12), refresh = 0
)


# priors check
pp_check(model_prior_nfeeds, fun = "dens_overlay") +
  xlim(-1, 50) +
  labs(
    x = "Number of feedings per 48h",
    y = "Density",
    title = "Prior predictive distribution\nfor mid chick rearing (devices removed)",
    subtitle = "Model: Poisson with log link"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(face = "bold"),
    legend.position = "top",
    legend.title = element_blank()
  )

# ggsave("./Plots/Q2a_priors_nfeeds.png", width = 8, height = 9, dpi = 300)

 

# model ----
model_nfeeds <- brm(
  formula = nfeeds ~ parent * group + (1 | season) + (1 | nest),
  data = top_bottom_all_long,
  family = poisson(link = "log"),
  prior = priors,
  seed = 1313,
  chains = 4, cores = 4,iter = 4000, warmup = 2000,     
  control = list(adapt_delta = 0.995, max_treedepth = 12), refresh = 0  # better exploration
)


# model diagnostics

# posterior check
pp_check(model_nfeeds, type = "dens_overlay")


# overdispersion check
pp_check(model_nfeeds, type = "stat", stat = "sd")
pp_check(model_nfeeds, type = "stat", stat = "var")

# model check
plot(model_nfeeds)

mcmc_plot(model_nfeeds, type = "acf")
mcmc_plot(model_nfeeds, type = "pairs")

summary(model_nfeeds)$warnings

nuts_params(model_nfeeds) %>%  subset(Parameter == "divergent__") %>% 
  group_by(Value) %>% 
  summarise(n = n())

nuts_params(model_nfeeds) %>%  subset(Parameter == "treedepth__") %>% 
  group_by(Value) %>% 
  summarise(n = n())


# model summary
describe_posterior(model_nfeeds)
bayes_R2(model_nfeeds)
summary(model_nfeeds)



# plotting posterior predictions -----

# data prep
newdata <- expand.grid(
  parent = unique(top_bottom_all_long$parent),
  group  = unique(top_bottom_all_long$group),
  season = NA
)

pred <- posterior_epred(model_nfeeds, allow_new_levels = TRUE, newdata = newdata)

# means, quantiles
pred_summary <- apply(pred, 2, function(x) {
  c(mean = mean(x),
    lwr = quantile(x, 0.025),
    upr = quantile(x, 0.975))
})

pred_df <- cbind(newdata, t(pred_summary))
pred_df <- pred_df %>%
  dplyr::mutate(row_id = dplyr::row_number())


# pred posterior vals
pred_long <- pred %>%
  as.data.frame() %>%
  mutate(draw = row_number()) %>%   # numer próbki posterioru
  pivot_longer(
    cols = starts_with("V"),        # V1, V2, V3, V4...
    names_to = "row_id",
    values_to = "value"
  ) %>%
  mutate(row_id = as.numeric(gsub("V", "", row_id))) %>%  # V1 → 1
  left_join(pred_df, by = "row_id") %>% 
  mutate(parent = if_else(parent == "top", "higher-feeding", "lower-feeding")) # for plotting/ms adjust labels



# plot
posteriors1_plot <- ggplot(pred_long, aes(x = parent, y = value, fill = group)) +
  stat_halfeye(alpha = 0.6, position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = c(
    "control"    = "steelblue",  
    "deployment" = "tomato4"  
  )) +
  labs(x = "Parent type",
       y = "Posterior predicted number of feedings (per 48h)",
       title = "Carry-over effects: number of feedings in mid chick rearing (after device removal)",
       subtitle = "posterior distributions and medians with 95 CI",
       fill = "Experimental group") +
  
  theme_bw()


# ggsave(posteriors1_plot, filename = "./Plots/Q2a_posteriors_nfeeds.png", width = 8, height = 9, dpi = 300)

# contrasts -----

em <- emmeans(model_nfeeds, ~ group * parent, type = "response")

contr <- pairs(em)
contr_df <- as.data.frame(contr) %>%
  rename(
    lower = lower.HPD,
    upper = upper.HPD
  )

# delete cross-level contrasts
contr_df <- contr_df %>%
  filter(
    !(grepl("bottom", contrast) &
        grepl("top", contrast) &
        grepl("control", contrast) &
        grepl("deployment", contrast))
  )

# contrasts order
desired_order <- c(
  "control bottom / control top",
  "deployment bottom / deployment top",
  "control bottom / deployment bottom",
  "control top / deployment top"
)

all_contrasts <- unique(contr_df$contrast)
desired_order <- desired_order[desired_order %in% all_contrasts]

contr_df$contrast <- factor(contr_df$contrast, levels = desired_order)

# relabel for plotting
contr_df$contrast <- recode(
  contr_df$contrast,
  "control bottom / deployment bottom" = "control vs experimental, lower-feeding parents",
  "control bottom / control top"       = "control; lower-feeding vs higher-feeding parents",
  "deployment bottom / deployment top" = "experimental; lower-feeding vs higher-feeding parents",
  "control top / deployment top"       = "control vs experimental, higher-feeding parents"
)


# plotting contrasts - forest plot
plot1_contrasts <- ggplot(contr_df, aes(
  x = contrast,
  y = ratio,
  ymin = lower,
  ymax = upper
)) +
  geom_pointrange(size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Parent effect within each group",
    x = "",
    y = "Ratio (lower-feeding vs higher-feeding parent)"
  )

ggsave(plot1_contrasts, filename = "Q2a_contrasts_nfeeds.png", width = 10, height = 9, dpi = 300)

# combine posteriors and contrasts plots (for ms)

combined_plot <- posteriors1_plot + plot1_contrasts +
  plot_layout(ncol = 2)
# ggsave(combined_plot, filename = "./Plots/Q2a_res_plots.png", width = 14, height = 9, dpi = 300)





# ---- Analysis #2 -------------------------------
# Analysis  of the relationship between provisioning rates of the parners 


# priors ----

# same priors as the data does not change
priors <- c(
  set_prior("normal(2, 0.3)", class = "Intercept"),
  set_prior("normal(0, 0.2)", class = "b"),
  set_prior("exponential(5)", class = "sd", group = "season")
)

model_prior_topbottom <- brm(
  formula = top ~ bottom * group + (1 | season),
  data = top_bottom_chr2,
  family = poisson(link = "log"),
  prior = priors,
  seed = 1313,
  sample_prior = "only",
  chains = 4, cores = 4, iter = 4000,
  refresh = 0
)

# priors check
pp_check(model_prior_topbottom, fun = "dens_overlay") +
  xlim(-1, 50) +
  labs(
    x = "Number of feedings per 48h",
    y = "Density",
    title = "Prior predictive distribution\nfor mid chick rearing (devices removed)",
    subtitle = "Model: Poisson with log link"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(face = "bold"),
    legend.position = "top",
    legend.title = element_blank()
  )

# ggsave("./Plots/Q2b_priors_feeds_relation.png", width = 7, height = 9, dpi = 300)


# model----

model_topbottom <- brm(
  formula = top ~ bottom * group + (1 | season),
  data = top_bottom_chr2,
  family = poisson(link = "log"),
  prior = priors,
  seed = 1313,
  chains = 4, cores = 4,iter = 4000, warmup = 2000,     
  control = list(adapt_delta = 0.995, max_treedepth = 12), refresh = 0  # better exploration
)



# model diagnostics

# posterior check
pp_check(model_topbottom, type = "dens_overlay")

# overdispersion check
pp_check(model_topbottom, type = "stat", stat = "sd")
pp_check(model_topbottom, type = "stat", stat = "var")

# model check
plot(model_topbottom)

mcmc_plot(model_topbottom, type = "acf")
mcmc_plot(model_topbottom, type = "pairs")

summary(model_topbottom)$warnings

nuts_params(model_topbottom) %>%  subset(Parameter == "divergent__") %>% 
  group_by(Value) %>% 
  summarise(n = n())

nuts_params(model_topbottom) %>%  subset(Parameter == "treedepth__") %>% 
  group_by(Value) %>% 
  summarise(n = n())

loo(model_topbottom)



# model summary
describe_posterior(model_topbottom)
bayes_R2(model_topbottom)

summary(model_topbottom)


# plotting posterior predictions ----

# prep data
newdata <- expand.grid(
  bottom = seq(min(top_bottom_chr2$bottom), max(top_bottom_chr2$bottom), length.out = 100),
  group = c("control", "deployment")
)

pred <- posterior_epred(
  model_topbottom,
  newdata = newdata,
  re_formula = NA   
) %>% 
  as.data.frame()

pred_summary <- newdata %>%
  mutate(
    pred_mean  = apply(pred, 2, mean),
    pred_lower = apply(pred, 2, quantile, probs = 0.025),
    pred_upper = apply(pred, 2, quantile, probs = 0.975)
  )

# plot
posteriors2_plot <- ggplot() +
  geom_jitter(
    data = top_bottom_chr2,
    aes(x = bottom, y = top, color = group),
    alpha = 0.5, size = 2
  ) +
  geom_ribbon(
    data = pred_summary,
    aes(x = bottom, ymin = pred_lower, ymax = pred_upper, fill = group),
    alpha = 0.2, color = NA
  ) +
  geom_line(
    data = pred_summary,
    aes(x = bottom, y = pred_mean, color = group),
    linewidth = 1.2
  ) +
  scale_color_manual(
    values = c(
      "control"    = "steelblue",
      "deployment" = "tomato4"
    )
  ) +
  scale_fill_manual(
    values = c(
      "control"    = "steelblue",
      "deployment" = "tomato4"
    )
  ) +
  labs(
    title = "Model predictions for feeding rate of higher-feeding parent",
    x = "Feeding rate of lower-feeding parent count",
    y = "Predicted feeding rate of higher-feeding parent",
    color = "Experimental group",
    fill = "Experimental group"
  ) +
  theme_minimal(base_size = 14)


# ggsave(posteriors2_plot, filename = "./Plots/Q2b_posteriors_feeds_relation.png", width = 9, height = 9, dpi = 300)


# ---- Analysis #3 -------------------------------
# Analysis of parents feeding position in the early (burden deployment) and mid chick rearing (burden removal)

# data
all_nfeeds <- readRDS("./Data/nfeeds_chr12_df.rds")

# data adjustment
top_bottom_temp <- all_nfeeds %>%
  group_by(session, season, nest, sx) %>%
  summarise(total_feeds = sum(nfeeds, na.rm = TRUE), .groups = "drop_last") %>%
  arrange(desc(total_feeds)) %>%
  mutate(rank = row_number()) %>%   # 1 = top feeder, 2 = bottom feeder
  mutate(feeder_type = if_else(rank == 1, "top", "bottom")) %>%
  ungroup() 

top_bottom_df <- all_nfeeds %>%
  left_join(top_bottom_temp, by = c("session", "season", "nest", "sx")) %>%
  mutate(group = if_else(nest_device_stat == "control", "control", "deployment")) 

df_wide <- top_bottom_df %>%
  select(season, session, nest_device_stat, nest, sx, group, feeder_type) %>%
  distinct() %>%
  pivot_wider(
    names_from  = session,
    values_from = feeder_type,
    names_prefix = "session_"
  ) %>% 
  rename(
    session_cr1 = `session_chick rearing1`,
    session_cr2 = `session_chick rearing2`
  )


# alluvial plot
df_wide <- df_wide %>%
  mutate(group = factor(group, levels = c("control", "deployment")))

allu_df <- df_wide %>% 
  filter(!is.na(group)) %>%
  filter(session_cr1 == "top") %>%   # single (top) parent from the nest (the other would be redundant)
  # relabelling for ms
  mutate(session_cr1 = if_else(session_cr1 == "top", "higher-feeding", NA_character_)) %>% 
  mutate(session_cr2 = if_else(session_cr2 == "top", "higher-feeding", 
                               if_else(session_cr2 == "bottom", "lower-feeding", NA_character_)))
  
  
allu_plot <-   ggplot(allu_df, aes(axis1 = session_cr1,
           axis2 = session_cr2,
           y = 1,
           fill = group,
           order = as.numeric(group))) +
  geom_alluvium(alpha = 0.7, na.rm = TRUE) +
  geom_stratum(width = 0.2, na.rm = TRUE) +
  geom_text(stat = "stratum",            # <- use geom_text instead of geom_label
            aes(label = after_stat(stratum)),
            angle = 90,        # rotate text
            na.rm = TRUE) +
  scale_x_discrete(limits = c("early (burden on) ", "mid (burden off)")) +
  scale_fill_manual(values = c(
    "control"    = "steelblue",
    "deployment" = "tomato4"
  ), na.translate = FALSE) +
  labs(
    title = "Change in the feeder status between early and mid chick rearing",
    x = "Chick rearing phase",
    y = "Number of individuals",
    fill = ""
  ) +
  theme_minimal()

# ggsave(allu_plot, filename = "./Plots/Q2c_allu_plot.png", width = 7, height = 9, dpi = 300)


# data prep for the model
df_model <- df_wide %>%
  mutate(
    stayed_top = case_when(
      session_cr1 == "top" & session_cr2 == "top"    ~ 1,
      session_cr1 == "top" & session_cr2 != "top"    ~ 0,
      TRUE                                           ~ NA_real_
    )
  ) %>%
  filter(!is.na(stayed_top))


df_model <- df_model %>%
  mutate(
    group  = factor(group),     # control / deployment
    nest   = factor(nest),
    season = factor(season)
  )


# bayes model

priors <- c(
  set_prior("normal(0, 1.5)", class = "Intercept"),
  set_prior("normal(0, 1)", class = "b"),               
  set_prior("student_t(3, 0, 2.5)", class = "sd")       
)

fit_bayes <- brm(
  formula = stayed_top ~ group + (1 | season),
  data    = df_model,
  family  = bernoulli(link = "logit"),
  prior   = priors,
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  seed    = 1234,
  control = list(adapt_delta = 0.995, max_treedepth = 12)
)


# model summary
summary(fit_bayes)
fixef(fit_bayes)
posterior_summary(fit_bayes)


# diagnostics
plot(fit_bayes)
neff_ratio(fit_bayes)
pp_check(fit_bayes)


ce_group <- conditional_effects(
  fit_bayes,
  effects = "group",
  method  = "posterior_epred" # means a posteriori on probabilistic scale
)

ce_data <- as.data.frame(ce_group$group)

# plot poterior props 
pp1 <- ggplot(ce_data, aes(x = group, y = estimate__, fill = group)) +
  geom_col(alpha = 0.7) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__), width = 0.2, ) +
  scale_fill_manual(
    values = c(
      "control"    = "steelblue",
      "deployment" = "tomato4"
    ),
    na.translate = FALSE
  ) +
  ylim(0, 1) +
  theme_minimal() +
  labs(
    x = "Group",
    y = "Posterior mean P(stayed higher-feeding parent)",
    title = "Probability of staying higher-feeding parent\n(control vs deployment, Bayesian model)"
  )


# contrasts

# generated quantities – epred for each group:
nd <- expand.grid(
  group  = factor(c("control", "deployment"), levels = levels(df_model$group)),
  season = NA,
  nest   = NA
)

# posterior predictions on prob scale
ep <- posterior_epred(fit_bayes, allow_new_levels = TRUE, newdata = nd)

# col 1 = control, 2 = deployment
diff_prob <- ep[, 2] - ep[, 1]

# plot diff distribution
pp2 <- tibble(diff_prob = as.numeric(diff_prob)) %>%
  ggplot(aes(x = diff_prob)) +
  geom_histogram(bins = 40, fill = "grey40", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Deployment - Control (P(stayed higher-feeding parent))",
    y = "Posterior density (approx.)",
    title = "Poster distribution of the difference between the treatment groups\nin probability of staying higher-feeding parent"
  )


combined_plot <- pp1 + pp2

# ggsave(
#   filename ="./Plots/Q2c_posterior_props_diff_combined.png",
#   plot = combined_plot,
#   width = 14,
#   height = 5,
#   dpi = 300
# )

