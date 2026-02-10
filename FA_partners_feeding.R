# PARTNERS 

# Libs ----

library(tidyverse)

library(brms)
library(bayesplot)
library(bayestestR)
library(emmeans)

library(tidybayes)
library(ggdist)

library(patchwork)


# ---- Read data ----------------------------------

nfeeds_chr1_df <-   readRDS("C:/Users/DELL/Dropbox/Working files/Projects/GPSeffect/GPS_effect/nfeeds_chr1_df.rds")


# examine sample size
nfeeds_chr1_df %>%
  distinct(season, nest, nest_device_stat) %>%
  count(season, nest_device_stat)


# examine sex of the deployed ind
nfeeds_chr1_df %>%
  filter(nest_device_stat %in% c("acc", "gps"),
         partner_device_stat == "deploy") %>% 
  group_by(sx, nest_device_stat) %>% 
  summarise(n = n())


# examine sex of the top feeders
paired <- nfeeds_chr1_df %>%
  select(season, nest, nest_device_stat, sx, nfeeds) %>%
  pivot_wider(
    names_from = sx,
    values_from = nfeeds,
    names_prefix = "feeds_"
  )

top_bottom_sex <- paired %>%
  mutate(
    top_feeder_sex = case_when(
      feeds_f > feeds_m ~ "female",
      feeds_m > feeds_f ~ "male",
      feeds_f == feeds_m ~ "same"
    ),
    bottom_feeder_sex = case_when(
      feeds_f < feeds_m ~ "female",
      feeds_m < feeds_f ~ "male",
      feeds_f == feeds_m ~ "same"
    )
  )

top_bottom_sex %>%
  group_by(nest_device_stat) %>% 
  count(top_feeder_sex) 
  

female_top <- 20
male_top   <- 18

binom.test(female_top, female_top + male_top, p = 0.5)

# ---- Reshape data -------------------------------

# controls
partners_control <- nfeeds_chr1_df %>%
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
partners_device <- nfeeds_chr1_df %>%
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
top_bottom_chr1 <- bind_rows(
  partners_control  %>% select(season, nest, top, bottom, group),
  partners_device %>% select(season, nest, top, bottom, group)
  
)


# ---- Analysis #1 -------------------------------

# prepare data for the model
top_bottom_all_long <- top_bottom_chr1 %>%
  pivot_longer(names_to = "parent", values_to = "nfeeds", cols = -c(season, nest, group)) 


# priors ----

priors <- c(
  set_prior("normal(1.8, 0.3)", class = "Intercept"),
  set_prior("normal(0, 0.5)", class = "b"),
  set_prior("exponential(5)", class = "sd", group = "season"),
  set_prior("exponential(5)", class = "sd", group = "nest")
)



# priors check
model_prior_nfeeds <- brm(
  formula = nfeeds ~ parent * group + (1 | season) + (1 | nest),
  data = top_bottom_all_long,
  family = poisson(link = "log"),
  prior = priors,
  sample_prior = "only",
  seed = 1313,
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.995, max_treedepth = 12), refresh = 0
)


pp_check(model_prior_nfeeds, fun = "dens_overlay") +
  xlim(-1, 50) +
  labs(
    x = "Number of feedings per 48h",
    y = "Density",
    title = "Prior predictive distribution",
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

# ggsave("prior_feeds_plot.png", width = 7, height = 9, dpi = 300)


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
  mutate(draw = row_number()) %>%  
  pivot_longer(
    cols = starts_with("V"),        
    names_to = "row_id",
    values_to = "value"
  ) %>%
  mutate(row_id = as.numeric(gsub("V", "", row_id))) %>%
  left_join(pred_df, by = "row_id")


# plot
ggplot(pred_long, aes(x = parent, y = value, fill = group)) +
  stat_halfeye(alpha = 0.6, position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = c(
    "control"    = "steelblue",  
    "deployment" = "tomato4"  
  )) +
  labs(x = "Parent type",
     y = "Posterior predicted number of feedings (per 48h)",
     title = "Posterior distributions and medians with 95 CI",
     fill = "Experimental group") +
  
    theme_bw()


# ggsave("posterior_feeds_plot.png", width = 7, height = 9, dpi = 300)

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

# plotting contrasts - forest plot
ggplot(contr_df, aes(
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
    y = "Ratio (bottom / top)"
  )

# ggsave("contrasts_feeds_plot.png", width = 7, height = 9, dpi = 300)





# ---- Analysis #2 -------------------------------


# priors ----

# same priors as the data does not change
priors <- c(
  set_prior("normal(1.8, 0.3)", class = "Intercept"),
  set_prior("normal(0, 0.5)", class = "b"),
  set_prior("exponential(5)", class = "sd", group = "season")
)

model_prior_topbottom <- brm(
  formula = top ~ bottom * group + (1 | season),
  data = top_bottom_chr1,
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
    x = "Feeding rate (n feeds/48h) of top feeders",
    y = "Density",
    title = "Prior predictive distribution",
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

ggsave("prior_feeds_relation_plot.png", width = 7, height = 9, dpi = 300)



# model----

model_topbottom <- brm(
  formula = top ~ bottom * group + (1 | season),
  data = top_bottom_chr1,
  family = poisson(link = "log"),
  prior = priors,
  seed = 1313,
  chains = 4, cores = 4,iter = 4000, warmup = 2000,     
  control = list(adapt_delta = 0.995, max_treedepth = 12), refresh = 0  # better exploration
)


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


# plotting posterior preditions ----

# prep data
newdata <- expand.grid(
  bottom = seq(min(top_bottom_chr1$bottom), max(top_bottom_chr1$bottom), length.out = 100),
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
ggplot() +
  geom_jitter(
    data = top_bottom_chr1,
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
    title = "Model predictions for bottom feeding rate",
    x = "Bottom partner feeding count",
    y = "Predicted top partner feeding count",
    color = "Experimental group",
    fill = "Experimental group"
  ) +
  theme_minimal(base_size = 14)


# posterior slopes
posterior <- as_draws_df(model_topbottom)

control_slope     <- posterior$b_bottom
deployment_slope  <- posterior$b_bottom + posterior$`b_bottom:groupdeployment`

slopes_df <- tibble(
  slope = c(control_slope, deployment_slope),
  group = rep(c("control", "deployment"), each = length(control_slope))
)

p1 <- ggplot(slopes_df, aes(x = group, y = slope, fill = group)) +
  stat_halfeye(alpha = 0.6) +
  scale_fill_manual(values = c("control" = "steelblue", "deployment" = "tomato4")) +
  labs(
    x = "",
    y = "Posterior slope (effect of bottom on top)",
    title = "Posterior distribution of slopes by group"
  ) +
  theme_bw()

# ggsave("posterior_slopes_plot.png", width = 7, height = 9, dpi = 300)



# posterior slopes - diff
post <- as_draws_df(model_topbottom)

control_slope     <- post$b_bottom
deployment_slope  <- post$b_bottom + post$`b_bottom:groupdeployment`

slope_diff <- control_slope - deployment_slope
df_diff <- tibble(slope_diff = slope_diff)

p2 <- ggplot(df_diff, aes(x = slope_diff)) +
  geom_histogram(bins = 40, fill = "grey40", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", size = 1) +
  labs(
    title = "Posterior distribution of the difference in slopes",
    subtitle = "Difference = slope_control − slope_deployment",
    x = "Difference in slope",
    y = "Posterior density"
  ) +
  theme_bw(base_size = 14)


# ggsave("posterior_slopes_diff_plot.png", width = 7, height = 9, dpi = 300)


combined_plot <- p1 + p2

# ggsave(
#   filename ="posterior_slopes_combined.png",
#   plot = combined_plot,
#   width = 12,     
#   height = 5,     
#   dpi = 300
# )


