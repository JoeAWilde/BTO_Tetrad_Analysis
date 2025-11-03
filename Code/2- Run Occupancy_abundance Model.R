#-----------Intro-------------#

# Author: Joe A. Wilde
# Creation: 20/11/2024

#----------------------------#

#Load required packages ####
library(tidyverse)
library(tidylog)
library(cmdstanr)
library(bayesplot)
library(ggExtra)
library(progress)
library(ggdist)
library(sf)
library(terra)
library(tidyterra)
library(ggspatial)
library(HDInterval)
library(ggpubr)
library(ggnewscale)
library(ggblend)
library(kableExtra)
options(kableExtra.latex.load_packages = FALSE)

#Read in the data ####
df <- read.table("Outputs/script_1/ModelData.csv", sep = ",", header = T, fill = T) %>%
  na.omit() %>%
  pivot_longer(cols = c(pheasant_breed, pheasant_wint), names_to = "time_period", values_to = "abund") %>%
  filter(abund != "M") %>%
  mutate(new_abund = 0)

#Convert the categorical abundance measures to continuous ####
for(i in 1:nrow(df)) {
  if(i == 1) {
    pb <- progress_bar$new(total = nrow(df), format = "[:bar] :percent eta::eta", clear = F)
    pb$tick(0)
  }
  if(df$County[i] == "Cornwall") {
    df$new_abund[i] <- case_when(df$abund[i] == "0" ~ 0,
                                 df$abund[i] == "0.11-1" ~ 0.55,
                                 df$abund[i] == "0.2-1" ~ 0.6,
                                 df$abund[i] == "1" ~ 1,
                                 df$abund[i] == "1.33-1.5" ~ 1.415,
                                 df$abund[i] == "1.25-2" ~ 1.625,
                                 df$abund[i] == "1.6-3" ~ 2.3,
                                 df$abund[i] == "2" ~ 2,
                                 df$abund[i] == "2.25-3" ~ 2.625,
                                 df$abund[i] == "3" ~ 3,
                                 df$abund[i] == "3.25-5" ~ 4.125,
                                 df$abund[i] == "3.5-6" ~ 4.75,
                                 df$abund[i] == "4" ~ 4,
                                 df$abund[i] == "5-32" ~ 18.5,
                                 df$abund[i] == "5.33-35" ~ 20.165, 
                                 df$abund[i] == "7-85" ~ 46)
  } else if(df$County[i] == "Devon") {
    df$new_abund[i] <- case_when(df$abund[i] == "0" ~ 0,
                                 df$abund[i] == "1" ~ 0.5,
                                 df$abund[i] == "2" ~ 1.5,
                                 df$abund[i] == "3" ~ 3.5,
                                 df$abund[i] == "4" ~ 7.5,
                                 df$abund[i] == "5" ~ 15,
                                 df$abund[i] == "6" ~ 35,
                                 df$abund[i] == "7" ~ 50)
  } else if(df$County[i] == "Berkshire") {
    df$new_abund[i] <- case_when(df$abund[i] == "0" ~ (0*7.3),
                                 df$abund[i] == "1" ~ (0.1*7.3),
                                 df$abund[i] == "2" ~ (0.25*7.3),
                                 df$abund[i] == "3" ~ (0.75*7.3),
                                 df$abund[i] == "4" ~ (1.25*7.3),
                                 df$abund[i] == "5" ~ (1.75*7.3),
                                 df$abund[i] == "6" ~ (2.25*7.3),
                                 df$abund[i] == "7" ~ (2.5*7.3))
  } else {
    df$new_abund[i] <- as.numeric(df$abund[i])
  }
  
  pb$tick()
}

#Tidy the dataframe to be used in the analysis ####
post_df <- df %>%
  mutate(
    new_tet = as.integer(factor(X, levels = unique(X))), 
    new_county = as.integer(factor(County, levels = unique(County)))
    time_period = factor(time_period, levels = c("pheasant_wint", "pheasant_breed"))
  ) %>%
  group_by(County, time_period) %>% 
  mutate(
    N_birds = sum(new_abund), 
    rel_abund = new_abund / sum(N_birds)
  ) %>%
  ungroup() %>%
  group_by(County) %>%
  mutate(
    N_tetrads = length(unique(new_tet))
  ) %>%
  ungroup()

#Create a dataframe of predictors
options(na.action='na.pass')
preds <- data.frame(model.matrix(~ County + time_period + broad_wood + conif_wood + arable + imp_grass + 
                                   semnat_grass + mountain + saltwater + freshwater +
                                   coastal + built_up + habitat_div_shannon + PA_suit,
                                 data = post_df)) %>%
  select(-X.Intercept.) %>%
  scale()

#Save the centre and scale of the predictors to back-transform the coefficients later ####
pred_mus <- attr(preds, "scaled:center")
pred_sds <- attr(preds, "scaled:scale")

#Create the two- and three-way interactions ####
preds <- preds %>%
  data.frame() %>%
  mutate(
    berk_dummy = if_else(CountyCornwall < 0 & CountyDevon < 0 & CountyHertfordshire < 0, 0.5, -0.5),
    `pheasant_wint:PA_suit` = time_periodpheasant_breed * PA_suit, 
    `CountyBerkshire:time_periodpheasant_wint:PA_suit` = berk_dummy * (time_periodpheasant_breed*-1) * PA_suit,
    `CountyCornwall:time_periodpheasant_wint:PA_suit` = CountyCornwall * (time_periodpheasant_breed*-1) * PA_suit,  
    `CountyDevon:time_periodpheasant_wint:PA_suit` =   CountyDevon * (time_periodpheasant_breed*-1) * PA_suit,
    `CountyHertfordshire:time_periodpheasant_wint:PA_suit` = CountyHertfordshire * (time_periodpheasant_breed*-1) * PA_suit,
    `CountyBerkshire:time_periodpheasant_breed:PA_suit` = berk_dummy * time_periodpheasant_breed * PA_suit,
    `CountyCornwall:time_periodpheasant_breed:PA_suit` =   CountyCornwall * time_periodpheasant_breed * PA_suit,
    `CountyDevon:time_periodpheasant_breed:PA_suit` =  CountyDevon * time_periodpheasant_breed * PA_suit,    
    `CountyHertfordshire:time_periodpheasant_breed:PA_suit` = CountyHertfordshire * time_periodpheasant_breed * PA_suit,
    saltwater = if_else(is.na(saltwater), 0, saltwater),
    coastal = if_else(is.na(coastal), 0, coastal), 
    log_ntetrads = log(post_df$N_tetrads)
  ) %>%
  select(-berk_dummy)

#Create a list of data for the model ####
stan_data <- list(
  N = nrow(post_df), 
  N_tets = length(unique(post_df$new_tet)), 
  tets = post_df$new_tet,
  N_counties = length(unique(post_df$new_county)),
  county = post_df$new_county,
  X = ncol(preds), 
  preds = preds, 
  occ = ifelse(post_df$rel_abund> 0, 1, 0), 
  abund = post_df$rel_abund
)

#Compile the model ####
model <- cmdstan_model("Code/Stan script/occ_abund_model.stan")

#Sample from the model ####
fit <- model$sample(
  data = stan_data,
  iter_warmup = 1000,
  iter_sampling = 1000,
  chains = 4,
  parallel_chains = 4
)

#Save the model object ####
fit$save_object("Outputs/script_2/occ_abund_run.rds")

#Load model object
fit <- readRDS("Outputs/script_2/occ_abund_run.rds")

#Create a list of parameters of interest ####
pars <- c("a_occ_mu", "sd_occ_tets", "a_abund_mu", "sd_abund_tets", "beta", "p_det_eta", "phi_county")

#Visualise a trace plot of model run ####
mcmc_trace(fit$draws(pars))

#Output summary table for model run ####
fit$print(pars, max_rows = 34)

#Visualise posterior distributions of fixed effect coefficients ####
mcmc_areas(fit$draws("beta")) + 
  scale_y_discrete(labels = names(preds))

#Visualise posterior predictive check ####
pp_check(post_df$rel_abund , fit$draws("occ_sim", format = "matrix"), 
         ppc_dens_overlay)


#Extract highest density interval ####
all_draws <- fit$draws(pars, format = "df")
names(all_draws)[5:29] <- names(preds)
hdi_alldraws <- HDInterval::hdi(all_draws)

#Make plots of effect of PA ####
#Create a dataframe to simulate data from ####
sim_df <- data.frame(
  PA_suit_scaled = rep_len(seq(min(preds$PA_suit), max(preds$PA_suit), length.out = 50), length.out = 100),
  time_period_scaled = rep(c(-1, 1), each = 50)) %>%
  mutate(
    PA_suit = (PA_suit_scaled * pred_sds[16]) + pred_mus[16], 
    time_period = if_else(time_period_scaled == -1, "Winter", "Breeding")
  )

#Extract parameter draws ####
int_draws <- c(fit$draws("a_abund_mu", format = "matrix"))
pa_beta <- c(fit$draws("beta[16]", format = "matrix"))
time_beta <- c(fit$draws("beta[4]", format = "matrix"))

#Create a matrix to simulate response data into ####
y_sim <- matrix(nrow = length(int_draws), ncol = nrow(sim_df))

#Simulate response data ####
for(c in 1:nrow(sim_df)) {
  if(c == 1) {
    pb <- progress_bar$new(total = nrow(sim_df), 
                           format = "[:bar] :percent eta::eta", 
                           clear = F)
    pb$tick(0)
  }
  for(r in 1:length(int_draws)) {
    y_sim[r, c] <- boot::inv.logit(int_draws[r] + pa_beta[r] * sim_df$PA_suit_scaled[c] + 
                                    time_beta[r] * sim_df$time_period_scaled[c])
  }
  pb$tick()
}

#Extract the mean and upper and lower HDL of posterior prediction ####
sim_df$y_mu <- colMeans(y_sim)
sim_df$y_low <- HDInterval::hdi(y_sim)[1, ]
sim_df$y_upp <- HDInterval::hdi(y_sim)[2, ]

#Convert time period variable to clean text ####
post_df$time_period <- ifelse(post_df$time_period == "pheasant_wint", "Winter", "Breeding")

#Plot the predcicted effect of protected area coverage ####
p1 <- ggplot() + 
  # geom_point(data = post_df, aes(x = PA_suit, y = rel_abund / N_Birds, shape = factor(County), colour = time_period),
  #            alpha = 0.6, size = 5) + 
  geom_ribbon(data = sim_df, aes(x = PA_suit, ymin = y_low, ymax = y_upp, group = time_period, fill = time_period), alpha = 0.33)  |> blend("multiply") + 
  geom_smooth(data = sim_df, aes(x = PA_suit, y = y_mu, colour = time_period), se = F, linewidth = 3) + 
  scale_color_manual(name = "Time period", values = c("orange2", "purple3")) +
  scale_fill_manual(name = "Time period", values = c("orange2", "purple3")) +
  scale_shape_manual(name = "County", values = c(4, 0, 2, 1)) + 
  scale_y_continuous(name = "Relative pheasant abundance") + #, limits = c(0, 0.01)) + 
  scale_x_continuous(name = "Proportion of protected area coverage") + 
  theme_classic(base_size = 45) + 
  theme(legend.key.width = unit(2, "cm")) + 
  coord_cartesian(ylim = c(0, 0.00125))

p1
ggsave(p1, filename = "Outputs/plots/MS_final_PA_effect.png", units = "px", height = 4320, width = 7890)

#Make plot of effect of bw ####
#Create dataframe to simulate from ####
bw_sim_df <- data.frame(
  broad_wood_scaled = rep_len(seq(min(preds$broad_wood), max(preds$broad_wood), length.out = 50), length.out = 100),
  time_period_scaled = rep(c(-1, 1), each = 50)) %>%
  mutate(
    broad_wood = (broad_wood_scaled * pred_sds[5]) + pred_mus[5], 
    time_period = if_else(time_period_scaled == -1, "Winter", "Breeding")
  )
#Extract parameter draws ####
bw_beta <- c(fit$draws("beta[5]", format = "matrix"))

#Create a matrix to simulate response data into ####
bw_y_sim <- matrix(ncol = nrow(bw_sim_df), nrow = length(bw_beta))

#Simulate response data ####
for(r in 1:length(bw_beta)) {
  for(c in 1:nrow(bw_sim_df)) {
    bw_y_sim[r, c] <- boot::inv.logit(
      int_draws[r] + 
        bw_beta[r] * bw_sim_df$broad_wood_scaled[c] + 
        time_beta[r] * bw_sim_df$time_period_scaled[c]
    )
  }
}

#Extract simulated mean and upper and lower HDL ####
bw_sim_df$y_mu <- colMeans(bw_y_sim)
bw_sim_df$y_low <- HDInterval::hdi(bw_y_sim)[1, ]
bw_sim_df$y_upp <- HDInterval::hdi(bw_y_sim)[2, ]

#Plot the effect of broadlead woodland ####
p2 <- ggplot() + 
  # geom_point(data = post_df, aes(x = broad_wood/100, y = rel_abund / N_Birds, shape = factor(County), colour = time_period),
  #            alpha = 0.4, size = 5) + 
  geom_ribbon(data = bw_sim_df, aes(x = broad_wood/100, ymin = y_low, ymax = y_upp, group = time_period, fill = time_period), alpha = 0.33) |> blend("multiply") + 
  geom_smooth(data = bw_sim_df, aes(x = broad_wood/100, y = y_mu, colour = time_period), se = F, linewidth = 3) + 
  scale_color_manual(name = "Time period", values = c("orange2", "purple3")) +
  scale_fill_manual(name = "Time period", values = c("orange2", "purple3")) +
  scale_shape_manual(name = "County", values = c(4, 0, 2, 1)) + 
  scale_y_continuous(name = "Relative pheasant abundance") + #, limits = c(0, 0.01)) + 
  scale_x_continuous(name = "Proportion of broadleaf woodland coverage") + 
  theme_classic(base_size = 45) +
  theme(legend.key.width = unit(2, "cm")) + 
  coord_cartesian(ylim = c(0, 0.002))
p2
ggsave(p2, filename = "Outputs/plots/MS_final_bw_effect.png", units = "px", height = 4320, width = 7890)


#Create table of model summary
m1tidy <- fit$summary(pars)

names(m1tidy)<-c("Parameter", 
                 "Mean", "Median", "SD", "Mad", 
                 "5% quantile", "95% quantile", "R-hat", "Bulk", "Tail")

m1hdis <- hdi(fit$draws(pars, format = "df")[, 1:34], credMass=0.95)

m1tidy$Upper <- m1hdis[2,]
m1tidy$Lower <- m1hdis[1,]


m1tidy$Group <- c(
  rep("Pr(Occupied)", 2), 
  rep("Abundance intercept", 2), 
  rep("Abundance covariates", 25), 
  "Pr(Detection)", 
  rep("Family-specific parameter (Beta)", 4)
)

m1tidy$Parameter <- c(
  "Population mean", 
  "Population SD", 
  "Population mean", 
  "Population SD", 
  "County: Berkshire → Cornwall", 
  "County: Berkshire → Devon", 
  "County: Berkshire → Hertfordshire", 
  "Time period: winter → breeding", 
  "Habitat: Broadleaf woodland", 
  "Habitat: Coniferous woodland", 
  "Habitat: Arable", 
  "Habitat: Improved grassland", 
  "Habitat: Semi-natural grassland", 
  "Habitat: Mountain, heath and bog", 
  "Habitat: Saltwater", 
  "Habitat: Freshwater", 
  "Habitat: Coastal", 
  "Habitat: Built-up areas and gardens", 
  "Shannon habitat diversity index", 
  "Protected area (PA) coverage", 
  "PA coverage × Time period: winter", 
  "PA coverage × Time period: winter × County: Berkshire", 
  "PA coverage × Time period: winter × County: Cornwall", 
  "PA coverage × Time period: winter × County: Devon", 
  "PA coverage × Time period: winter × County: Hertfordshire", 
  "PA coverage × Time period: breeding × County: Berkshire", 
  "PA coverage × Time period: breeding × County: Cornwall", 
  "PA coverage × Time period: breeding × County: Devon", 
  "PA coverage × Time period: breeding × County: Hertfordshire", 
  "Pr(detection)", 
  "Phi (precision): Devon", 
  "Phi (precision): Hertfordshire", 
  "Phi (precision): Berkshire", 
  "Phi (precision): Cornwall"
)

p_det_draws <- fit$draws("p_det_eta", format = "df") %>%
  pull(p_det_eta) %>%
  boot::inv.logit(.)

m1tidy[30, 2] <- mean(p_det_draws)
m1tidy[30, 4] <- sd(p_det_draws)
m1tidy[30, 12] <- hdi(p_det_draws)[1]
m1tidy[30, 11] <- hdi(p_det_draws)[2]

m1tidy[,2:(length(m1tidy)-1)] <- round(m1tidy[2:(length(m1tidy)-1)], digits = 2)
m1tidy$`R-hat` <- round(m1tidy$`R-hat`, digits = 1)
m1tidy[,9:10] <- round(m1tidy[,9:10], digits = 0)
m1tidy[,11:12] <- round(m1tidy[, 11:12], digits = 2)

m1clean <- m1tidy %>%
  select(-c(`5% quantile`, `95% quantile`, Median, Mad, Group)) %>%
  select(Parameter, Mean, SD, Lower, Upper, `R-hat`, Bulk, Tail) %>%
  .[c(30, 1:29, 31:34), ] %>%
  data.frame()

table1 <- knitr::kable(m1clean, row.names = FALSE) %>%
  kableExtra::kable_classic(full_width = F) %>%
  pack_rows("Pr(occupied) (logit)",start_row = 2, end_row = 3, 
            label_row_css = "background-color: white;") %>%
  pack_rows("Abundance (logit)", start_row = 4, end_row = 30, 
            label_row_css = "background-color: white;") %>%
  pack_rows("Intercept",start_row = 4, end_row = 5, 
            label_row_css = "background-color: white;") %>%
  pack_rows("Covariates",start_row = 6, end_row = 30, 
            label_row_css = "background-color: white;") %>%
  pack_rows("Two-way interactions",start_row = 22, end_row = 22, 
            label_row_css = "background-color: white;") %>%
  pack_rows("Three-way interactions",start_row = 23, end_row = 30, 
            label_row_css = "background-color: white;") %>%
  pack_rows("Family-specific parameters",start_row = 31, end_row = 34, 
            label_row_css = "background-color: white;") %>%
  row_spec(0,background = "white") %>%
  add_header_above(c(" " = 1, "Posterior values" = 2, "95% Highest Density Limits" = 2," " = 1, "Effective sample size" = 2)) %>%
  cat(., file = "Outputs/tables/model_output_tables.html")


#extract residuals and plot spatially
y_pred_abund <- fit$draws("occ_sim", format = "df") %>%
  data.frame() %>%
  select(-c(.chain, .draw, .iteration)) %>%
  pivot_longer(cols = everything()) %>%
  # filter(value != 0) %>%
  group_by(name) %>%
  summarise(mean_value = mean(value), .groups = "drop") %>%
  mutate(
    idx = substr(name, 9, 999), 
    idx = as.integer(substr(idx, 1, (nchar(idx)-1)))
  ) %>%
  select(-name) %>%
  arrange(idx)

head(y_pred_abund)

post_df$abund_pred <- y_pred_abund$mean_value

abund_resid_df <- post_df %>%
  mutate(
    abund_resid = rel_abund - abund_pred
  )

cols <- RColorBrewer::brewer.pal(5, "Dark2") %>%
  .[4:5]

#Devon resid plot
wint_dev_df <- abund_resid_df %>% 
  filter(time_period == "Winter" & County == "Devon")

wint_dev_plot <- ggplot() + 
  geom_point(data = wint_dev_df, aes(x = X, y = Y, colour = abund_resid), size = 1) + 
  scale_colour_gradientn(name = "Relative abundance residuals \n(grey = not present/detected)", 
    colors = c("low" = cols[1], "high" = cols[2])) + 
  ggtitle("Devon - winter") + 
  theme_bw(base_size = 5)

ggsave(wint_dev_plot, path = "Outputs/plots/MS_final_plots/residual_plots/", filename = "devon_winter_resid.png", 
  units = "px", height = 990, width = 1280)

breed_dev_df <- abund_resid_df %>% 
  filter(time_period == "Breeding" & County == "Devon")

breed_dev_plot <- ggplot() + 
  geom_point(data = breed_dev_df, aes(x = X, y = Y, colour = abund_resid), size = 1) + 
  scale_colour_gradientn(name = "Relative abundance residuals \n(grey = not present/detected)", 
    colors = c("low" = cols[1], "high" = cols[2])) + 
  ggtitle("Devon - breeding") + 
  theme_bw(base_size = 5)

ggsave(breed_dev_plot, path = "Outputs/plots/MS_final_plots/residual_plots/", filename = "devon_breed_resid.png", 
  units = "px", height = 990, width = 1280)

#Cornwall resid plot
wint_cw_df <- abund_resid_df %>% 
  filter(time_period == "Winter" & County == "Cornwall")

wint_cw_plot <- ggplot() + 
  geom_point(data = wint_cw_df, aes(x = X, y = Y, colour = abund_resid), size = 1) + 
  scale_colour_gradientn(name = "Relative abundance residuals \n(grey = not present/detected)", 
    colors = c("low" = cols[1], "high" = cols[2])) + 
  ggtitle("Cornwall - winter") + 
  theme_bw(base_size = 5)

ggsave(wint_cw_plot, path = "Outputs/plots/MS_final_plots/residual_plots/", filename = "cornwall_winter_resid.png", 
  units = "px", height = 990, width = 1280)

breed_cw_df <- abund_resid_df %>% 
  filter(time_period == "Breeding" & County == "Cornwall")

breed_cw_plot <- ggplot() + 
  geom_point(data = breed_cw_df, aes(x = X, y = Y, colour = abund_resid), size = 1) + 
  scale_colour_gradientn(name = "Relative abundance residuals \n(grey = not present/detected)", 
    colors = c("low" = cols[1], "high" = cols[2])) + 
  ggtitle("Cornwall - breeding") + 
  theme_bw(base_size = 5)

ggsave(breed_cw_plot, path = "Outputs/plots/MS_final_plots/residual_plots/", filename = "cornwall_breed_resid.png", 
  units = "px", height = 990, width = 1280)


#Berkshire resid plot
wint_bk_df <- abund_resid_df %>% 
  filter(time_period == "Winter" & County == "Berkshire")

wint_bk_plot <- ggplot() + 
  geom_point(data = wint_bk_df, aes(x = X, y = Y, colour = abund_resid), size = 2) + 
  scale_colour_gradientn(name = "Relative abundance residuals \n(grey = not present/detected)", 
    colors = c("low" = cols[1], "high" = cols[2])) + 
  ggtitle("Berkshire - winter") + 
  theme_bw(base_size = 3)

ggsave(wint_bk_plot, path = "Outputs/plots/MS_final_plots/residual_plots/", filename = "Berkshire_winter_resid.png", 
  units = "px", height = 420, width = 1280)

breed_bk_df <- abund_resid_df %>% 
  filter(time_period == "Breeding" & County == "Berkshire")

breed_bk_plot <- ggplot() + 
  geom_point(data = breed_bk_df, aes(x = X, y = Y, colour = abund_resid), size = 2) + 
  scale_colour_gradientn(name = "Relative abundance residuals \n(grey = not present/detected)", 
    colors = c("low" = cols[1], "high" = cols[2])) + 
  ggtitle("Berkshire - breeding") + 
  theme_bw(base_size = 3)

ggsave(breed_bk_plot, path = "Outputs/plots/MS_final_plots/residual_plots/", filename = "Berkshire_breed_resid.png", 
  units = "px", height = 420, width = 1280)

#Hertfordshire resid plot
wint_hf_df <- abund_resid_df %>% 
  filter(time_period == "Winter" & County == "Hertfordshire")

wint_hf_plot <- ggplot() + 
  geom_point(data = wint_hf_df, aes(x = X, y = Y, colour = abund_resid), size = 2) + 
  scale_colour_gradientn(name = "Relative abundance residuals \n(grey = not present/detected)", 
    colors = c("low" = cols[1], "high" = cols[2])) + 
  ggtitle("Hertfordshire - winter") + 
  theme_bw(base_size = 5)

ggsave(wint_hf_plot, path = "Outputs/plots/MS_final_plots/residual_plots/", filename = "Hertfordshire_winter_resid.png", 
  units = "px", height = 780, width = 1280)

breed_hf_df <- abund_resid_df %>% 
  filter(time_period == "Breeding" & County == "Hertfordshire")

breed_hf_plot <- ggplot() + 
  geom_point(data = breed_hf_df, aes(x = X, y = Y, colour = abund_resid), size = 2) + 
  scale_colour_gradientn(name = "Relative abundance residuals \n(grey = not present/detected)", 
    colors = c("low" = cols[1], "high" = cols[2])) + 
  ggtitle("Hertfordshire - breeding") + 
  theme_bw(base_size = 5)

ggsave(breed_hf_plot, path = "Outputs/plots/MS_final_plots/residual_plots/", filename = "Hertfordshire_breed_resid.png", 
  units = "px", height = 780, width = 1280)
