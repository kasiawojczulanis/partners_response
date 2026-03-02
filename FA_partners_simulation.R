# Script for the study published in XXX

# Title: "Limited capacity for parental compensation under experimental handicap in a small pelagic seabird, the Little Auk" 
# Authors: Katarzyna Wojczulanis-Jakubas1*, Pauline Bodson1, Jerome Fort2, David Gremillet3,4, Andrea Grunst5, Melissa Grunst5, Kristin Piening1, Martyna Syposz1, Dariusz Jakubas1

# Study Question # 3 
# Simulation of the foraging trips number

# Data preparation ----

# Libs ----
library(data.table)
library(tidyverse)
library(patchwork)

# Real data -----
ff_all <- readRDS("./Data/activity_controls_chr1.rds") %>%
  filter(behav_blocks == "foraging_flight")

# short and long artifacts - viz

# plot short
ff25 <- ff_all %>% 
  filter(tot_dur_min <= quantile(tot_dur_min, 0.25)) %>%
  ggplot(aes(x = tot_dur_min)) +
  geom_histogram(
    bins = 30,
    fill = "grey",
    color = "white",
    alpha = 0.85
  ) +
  geom_vline(
    xintercept = 20,
    color = "tomato",
    linewidth = 1.1,
    linetype = "dashed"
  ) +
  scale_x_continuous(
    breaks = seq(0, max(ff_all$tot_dur_min, na.rm = TRUE), by = 5),
    expand = c(0.01, 0.01)
  ) +
  labs(
    x = "Duration (minutes)",
    y = "Count",
    title = "Distribution of short flights durations",
    subtitle = "below 25th percentile of the whole distribution") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray30"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "bold"),
    plot.caption = element_text(size = 10, color = "gray40")
  )


# plot long
ff75 <- ff_all %>% 
  filter(tot_dur_min >= quantile(tot_dur_min, 0.75)) %>%
  ggplot(aes(x = tot_dur_min/60)) +
  geom_histogram(
    bins = 30,
    fill = "grey",
    color = "white",
    alpha = 0.85
  ) +
  geom_vline(
    xintercept = 24,
    color = "tomato",
    linewidth = 1.1,
    linetype = "dashed"
  ) +
  scale_x_continuous(
    breaks = seq(0, max(ff_all$tot_dur_min, na.rm = TRUE), by = 5),
    expand = c(0.01, 0.01)
  ) +
  labs(
    x = "Duration (hours)",
    y = "Count",
    title = "Distribution of long flights durations",
    subtitle = "above 75th percentile of the whole distribution") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray30"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "bold"),
    plot.caption = element_text(size = 10, color = "gray40")
  )

# combined short and long
ff_25_75 <- ff25 + ff75

# ggsave(
#   filename ="./Plots/ff_artifacts.png",
#   plot = ff_25_75,
#   width = 12,
#   height = 5,
#   dpi = 300
# )

# filter out the artifacts

ff_dur <- ff_all %>% 
  filter(behav_blocks =="foraging_flight") %>% 
  filter(tot_dur_min > 20) %>% # short artifacts
  filter(tot_dur_min <= 24*60) # long artifacts

# n flights excluded
nrow(ff_all) - nrow(ff_dur)


# foraging trips distribution - viz

# seasons pooled
p1 <- ggplot(ff_dur, aes(x = log(tot_dur_min))) +
  geom_density(fill = "#4C72B0", alpha = 0.7, color = NA) +
  labs(
    title = "Years pooled",
    x = "log(duration in minutes)",
    y = "density"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# by season
p2 <- ggplot(ff_dur, aes(x = log(tot_dur_min))) +
  geom_density(fill = "#55A868", alpha = 0.7, color = NA) +
  facet_wrap(~season) +
  labs(
    title = "By year",
    x = "log(duration in minutes)",
    y = "density"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )


# combined
ff_dist <- p1 + p2 + 
  plot_annotation(
    title = "Density distributions of foraging trips duration",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 14, color = "gray30")
    )
  )

# ggsave(
#   filename ="ff_distributions.png",
#   plot = ff_dist,
#   width = 12,
#   height = 5,
#   dpi = 300
# )


# other actitivies distribution - viz

other_dt <- readRDS("./Data/activity_controls_chr1.rds") 
other_dt <- other_dt %>% 
  filter(behav_blocks %in% c("colony", "nest", "latency"))

ggplot(other_dt, aes(x = log(tot_dur_min))) +
  facet_grid(season ~ behav_blocks) +
  geom_density(fill = "#4C72B0", alpha = 0.7, color = NA) +
  labs(
    x = "log(duration in minutes)",
    y = "density"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# ggsave(
#   filename ="./Plots/other_distributions.png",
#   plot = ff_dist,
#   width = 12,
#   height = 5,
#   dpi = 300
# )




# Simulation ---- 

sim_dt <- readRDS("./Data/activity_controls_chr1.rds") 


# Prepare the vals ----
get_pool <- function(df, state) {
  df %>%
    filter(behav_blocks == state) %>%
    pull(tot_dur_min)
}

pool_forage_dt <- get_pool(sim_dt, "foraging_flight")
pool_latency_dt <- get_pool(sim_dt, "latency")
pool_nest_dt    <- get_pool(sim_dt, "nest")
pool_colony_dt  <- get_pool(sim_dt, "colony")

# option 1 
# sampling from the original vector
set.seed(1313)
N <- 10000

pool_forage_bootstrap <- sample(pool_forage_dt, size = N, replace = TRUE)
pool_latency_bootstrap <- sample(pool_latency_dt, size = N, replace = TRUE)
pool_nest_bootstrap <- sample(pool_nest_dt, size = N, replace = TRUE)
pool_colony_bootstrap <- sample(pool_colony_dt, size = N, replace = TRUE)


# option 2
# Kernel density–based random generator (smooth non‑parametric)
set.seed(1313)
N <- 10000

kde_sample <- function(n, dens) {
  keep <- dens$x >= 0
  sample(dens$x[keep], size = n, replace = TRUE, prob = dens$y[keep])
}

pool_forage_kde <- kde_sample(N, density(pool_forage_dt))
pool_latency_kde <- kde_sample(N, density(pool_latency_dt))
pool_nest_kde    <- kde_sample(N, density(pool_nest_dt))
pool_colony_kde  <- kde_sample(N, density(pool_colony_dt))


# Visualize generated distributions ---

df_all <- tibble(
  forage = c(pool_forage_bootstrap, pool_forage_kde),
  latency = c(pool_latency_bootstrap, pool_latency_kde),
  nest = c(pool_nest_bootstrap, pool_nest_kde),
  colony = c(pool_colony_bootstrap, pool_colony_kde),
  method = rep(c("Bootstrap", "KDE"),
               each = length(pool_forage_bootstrap))
) %>%
  pivot_longer(cols = c(forage, latency, nest, colony),
               names_to = "variable",
               values_to = "value")

boot_kde <- ggplot(df_all, aes(x = value, color = method)) +
  geom_density(size = 1.1) +
  facet_wrap(~ variable, scales = "free") +
  theme_minimal() +
  labs(
    x = "Value",
    y = "Density",
    title = "Bootstrap vs KDE – Density curves for all variables"
  )


# ggsave(
#   filename ="./Plots/bootstrap_vs_kde.png",
#   plot = boot_kde,
#   width = 12,
#   height = 5,
#   dpi = 300
# )


# Simulation 1 - full time budget ---

# data ----

# bootstrap dt - activate for boostrap sampling
# pool_forage <- pool_forage_bootstrap
# pool_latency <- pool_latency_bootstrap
# pool_colony <- pool_colony_bootstrap
# pool_nest   <- pool_nest_bootstrap

# kde dt - activate for kde sampling
pool_forage <- pool_forage_kde
pool_latency <- pool_latency_kde
pool_colony <- pool_colony_kde
pool_nest   <- pool_nest_kde


# funs and params ---

T_max <- 48 * 60 # time windonw (48h) in minutes

simulate_one_day <- function(T_max,
                             pool_forage,
                             pool_latency,
                             pool_nest,
                             pool_colony) {
  t <- 0
  n_flights <- 0
  
  events <- data.frame(
    event = character(),
    duration = numeric(),
    t_start = numeric(),
    t_end = numeric(),
    stringsAsFactors = FALSE
  )
  
  while (TRUE) {
    
    forage_time <- sample(pool_forage, 1, replace = TRUE)
    if (t + forage_time > T_max) break
    events <- rbind(events, data.frame(
      event = "forage",
      duration = forage_time,
      t_start = t,
      t_end = t + forage_time
      
    ))
    
    t <- t + forage_time
    n_flights <- n_flights + 1

    latency_time <- sample(pool_latency, 1, replace = TRUE)
    if (t + latency_time > T_max) break
    events <- rbind(events, data.frame(
      event = "latency",
      duration = latency_time,
      t_start = t,
      t_end = t + latency_time
    ))
    t <- t + latency_time
    
    nest_time <- sample(pool_nest, 1, replace = TRUE)
    if (t + nest_time > T_max) break
    events <- rbind(events, data.frame(
      event = "nest",
      duration = nest_time,
      t_start = t,
      t_end = t + nest_time
    ))
    t <- t + nest_time
    
    colony_time <- sample(pool_colony, 1, replace = TRUE)
    if (t + colony_time > T_max) break
    events <- rbind(events, data.frame(
      event = "colony",
      duration = colony_time,
      t_start = t,
      t_end = t + colony_time
    ))
    t <- t + colony_time
    
  }
  
  list(
    n_flights = n_flights,
    total_time = t,
    events = events
  )
}

# n runs simulation - n flights
n_sim <- 1000

results <- replicate(
  n_sim,
  simulate_one_day(T_max, pool_forage, pool_latency, pool_nest, pool_colony)$n_flights
)


df_results <- data.frame(n_flights = results)

# Simulation results - viz ----

# plot params
min_f <- min(df_results$n_flights)
max_f <- max(df_results$n_flights)
mean_f <- mean(df_results$n_flights)
q05   <- quantile(df_results$n_flights, 0.05)
q95   <- quantile(df_results$n_flights, 0.95)

# plot
ps1 <- ggplot(df_results, aes(x = n_flights)) +
  geom_histogram(binwidth = 1, fill = "grey", color = "black", alpha = 0.5) +
  geom_vline(xintercept = min_f, color = "gray", linetype = "dotted", size = 1) +
  geom_vline(xintercept = max_f, color = "gray", linetype = "dotted", size = 1) +
  geom_vline(xintercept = q05, color = "tomato4", linetype = "dashed", size = 1) +
  geom_vline(xintercept = q95, color = "tomato4", linetype = "dashed", size = 1) +
  geom_vline(xintercept = mean_f, color = "tomato", linetype = "solid", size = 1) +
  scale_x_continuous(breaks = seq(min(df_results$n_flights),
                                  max(df_results$n_flights),
                                  by = 1), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(
    title = "Simulated over 1000 iter, full time-budget (trip/latency/nest/colony)",
    x = "Number of foraging trips/feeds",
    y = "Frequency"
  ) + 
  theme_classic()




# Simulation 2 - restricted time budget  ---
# with latency, nest and colony lasting of minimum vals <= Q1

# data ----

# bootstrap dt - activate for boostrap sampling
# pool_latency <- pool_latency_bootstrap
# pool_nest   <- pool_nest_bootstrap
# pool_colony <- pool_colony_bootstrap

# kde dt - activate for kde sampling
pool_latency <- pool_latency_kde
pool_colony <- pool_colony_kde
pool_nest   <- pool_nest_kde


pools <- list(
  forage  = pool_forage,
  latency = pool_latency,
  nest    = pool_nest,
  colony  = pool_colony
)

pool_25 <- lapply(pools, function(x) {
  q25 <- quantile(x, 0.25)
  x[x <= q25]
})


pool_latency25 <- pool_25$latency
pool_nest25 <- pool_25$nest
pool_colony25 <- pool_25$colony

# n runs simulation ----
set.seed(1313)
n_sim <- 1000

results_flightsonly <- replicate(
  n_sim,
  simulate_one_day(T_max, pool_forage, pool_latency25, pool_nest25, pool_colony25)$n_flights
)

df_results_flightonly <- data.frame(n_flights = results_flightsonly)

# Simulation results - viz ----

# plot params
min_f <- min(df_results_flightonly$n_flights)
max_f <- max(df_results_flightonly$n_flights)
mean_f <- mean(df_results_flightonly$n_flights)
q05   <- quantile(df_results_flightonly$n_flights, 0.05)
q95   <- quantile(df_results_flightonly$n_flights, 0.95)

# plot
ps2 <- ggplot(df_results_flightonly, aes(x = n_flights)) +
  geom_histogram(binwidth = 1, fill = "grey", color = "black", alpha = 0.5) +
  geom_vline(xintercept = min_f, color = "gray", linetype = "dotted", size = 1) +
  geom_vline(xintercept = max_f, color = "gray", linetype = "dotted", size = 1) +
  geom_vline(xintercept = q05, color = "tomato4", linetype = "dashed", size = 1) +
  geom_vline(xintercept = q95, color = "tomato4", linetype = "dashed", size = 1) +
  geom_vline(xintercept = mean_f, color = "tomato", linetype = "solid", size = 1) +
  scale_x_continuous(breaks = seq(min_f, max_f,by = 1), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(
    title = "Simulated over 1000 iter, restricted time-budget (q25 of latency/nest/colony)",
    x = "Number of foraging trips/feeds",
    y = "Frequency"
  ) + 
  theme_classic()

combined_plot <- ps1 + ps2

# ggsave(
#   filename ="./Plots/Q3_simulation_bootstrap.png",
#   plot = combined_plot,
#   width = 14,
#   height = 5,
#   dpi = 300
# )

# ggsave(
#   filename ="./Plots/Q3_simulation_kde.png",
#   plot = combined_plot,
#   width = 14,
#   height = 5,
#   dpi = 300
# )
