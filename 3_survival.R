# LOAD REQUIRED PACKAGES ####
library(lubridate)
library(MASS)
library(bbmle)
library(lsmeans)
library(ggplot2)

# LOAD DATA ####
# ensure PIT tag numbers do not present in scientific notation
options(scipen=999)

# load comprehensive mark data set (all fish tagged at RST with predicted length at marine entry)
### individual fish are identified by PIT tag number
dat_rst <-read.csv("dat_rst_mar.csv") 

# load comprehensive mark-recapture data set (only tagged fish that were recaptured in the estuary)
### individual fish are identified by PIT tag number
### columns preceded by "tag_" indicate initial mark data
### columns preceded by "rc_" indicate recapture data
dat_rc <- read.csv("dat_rc.csv")

# load comprehensive marine survival data set (all fish tagged at RST with adult survival data from RFID array detections)
### 0 = died, 1 = redetected as adult 
dat_surv <- read.csv("dat_surv.csv")

# CHECK: SIZE-SELECTIVE MORTALITY BETWEEN FRESHWATER EXIT AND MARINE ENTRY ####
### Kolmogorov-Smirnov tests to check for annual differences in the freshwater exit size distribution for 1) all individuals tagged at the RST and 2) the subset of individuals recaptured in the estuary 
### no significant difference (p > 0.05) inidicates a lack of size-selective mortality in the estuary

# create annual RST tagging data frames
dat_rst_17 <- subset(dat_rst, tag_year==2017)
dat_rst_18 <- subset(dat_rst, tag_year==2018)
dat_rst_19 <- subset(dat_rst, tag_year==2019)

# create annual Reach 1 recapture data frames
dat_rc <- subset(dat_rc, tag_loc == "smolt trap" & rc_loc_1 == 1) # fish tagged at RST and recaptured in Reach 1 (estuary mouth / marine entry)
dat_rc$tag_year <- as.factor(dat_rc$tag_year)
dat_rc$tag_year <- droplevels(dat_rc$tag_year)

dat_rc_17 <- subset(dat_rc, tag_year==2017) # Reach 1 recaptures in 2017
dat_rc_18 <- subset(dat_rc, tag_year==2018) # Reach 1 recaptures in 2018
dat_rc_19 <- subset(dat_rc, tag_year==2019) # Reach 1 recaptures in 2019

# Kolmogorov-Smirnov similarity tests - distribution of tag lengths at RST tagging vs Reach 1 recaptures
ks.test(dat_rst_17$tag_length, dat_rc_17$tag_length) # 2017; p = 0.4883 (distributions not significantly different)
ks.test(dat_rst_18$tag_length, dat_rc_18$tag_length) # 2018; p = 0.4131 (distributions not significantly different)
ks.test(dat_rst_19$tag_length, dat_rc_19$tag_length) # 2019; p = 0.8555 (distributions not significantly different)

# PREPARE MARINE SURVIVAL DATA FOR SURVIVAL MODELS ####
# keep only columns relevant for analyses
dat_surv <- dat_surv[,c("PIT_number", "tag_year", "doy", "length", "surv", "recap_year")]
colnames(dat_surv) <- c("PIT_number", "tag_year", "tag_doy", "length_fw", "surv", "rc_year")

# check number of tagged spawners returning each year
table(dat_surv$tag_year, dat_surv$rc_year)

# append marine entry sizes to marine survival data
### create columns for predicted length at marine entry (length_mar) and standard deviation of growth estimate (sd_growth)
dat_surv$length_mar <- NA
dat_surv$sd_growth <- NA

# append length at marine entry and growth SD values to survival data frame
### for each individual in the survival data set
for(i in 1:nrow(dat_surv)) {
  ### identify tag year 
  tag_year <- dat_surv$tag_year[i]
  
  ### identify length at freshwater exit (tag length)
  tag_length <- dat_surv$length_fw[i]
  
  ### identify fish in RST data set that matches tag year and tag length
  id <- which(dat_rst$tag_year == tag_year & dat_rst$tag_length == tag_length)
  
  ### for that fish, assign predicted length at marine entry
  dat_surv$length_mar[i] <- unique(dat_rst[id,]$pred_length)
  
  ### for that fish, assign growth SD
  dat_surv$sd_growth[i] <- unique(dat_rst[id,]$pred_growth_sd)
}

# restrict survival data to fish tagged in 2017 or 2018
### these are smolt years with complete adult return data available
dat_surv <- dat_surv[dat_surv$tag_year == 2017 | dat_surv$tag_year == 2018,]

# create two data frames for marine survival models
### both data frames include all fish tagged at the RST in 2017 and 2018
### one data frame is for observed survival model (estuary-present model) -- length is observed length at freshwater exit; estuary growth is implicit in survival data because fish transit the estuary en route to the ocean
### one data frame is for predicted survival model (estuary-absent model) -- length is modeled length at marine entry; no estuary growth between freshwater exit and marine entry b/c these points are coincident in space and time
surv_obs <- data.frame("tag_year" = as.factor(dat_surv$tag_year),
                      "length" = dat_surv$length_fw,
                      "surv" = dat_surv$surv,
                      "estuary" = 1,
                      "sd_growth" = NA,
                      "weight" = NA,
                      "weight_sc" = 1 # all fish have equal weight in the model
                      )

surv_pred <- data.frame("tag_year" = as.factor(dat_surv$tag_year),
                        "length" = dat_surv$length_mar, 
                        "surv" = dat_surv$surv,
                        "estuary" = 0,
                        "sd_growth" = dat_surv$sd_growth,#standard deviation associated with predicted size at marine entry
                        "weight" = NA, 
                        "weight_sc" = NA
                        )

# for estuary-absent model, propagate error from residence and growth model predictions to weight individual fish
### each individual's weight is the inverse of the variance associated with its predicted size at marine entry -- this propagates error associated with both estuary residence and growth models, such that fish with more variance about their predicted length receive less weight in model fitting
### weights are scaled to 1 for consistency with observed survival data
surv_pred$weight <- 1/(surv_pred$sd_growth^2) 
surv_pred$weight_sc <- surv_pred$weight/max(surv_pred$weight) 

# combine observed and predicted data sets into one data frame
dat <- data.frame(rbind(surv_obs, surv_pred))
dat$surv <- as.integer(dat$surv)
dat$estuary <- as.factor(dat$estuary)
dat$tag_year <- as.factor(dat$tag_year)

# RUN CANDIDATE MARINE SURVIVAL MODELS
# RESPONSE VARAIABLE: observed survival (surv)
# CANDIDATE PREDICTORS: length (length), tag year (tag_year), and estuary (1 = estuary present, 0 = estuary absent)

### model family is quasibinomial rather than binomial to accomodate weights argument; statistical handling for quasibinomial distribution is the same as binomial (see https://github.com/alan-turing-institute/PosteriorBootstrap/issues/16)

mod_null <- glm(surv ~ 1, data = dat, family = quasibinomial("logit"), weights = weight_sc);summary(mod_null) # null model - constant survival rate
mod_l <- glm(surv ~ length, data = dat, family = quasibinomial("logit"), weights = weight_sc);summary(mod_l) # effect of length
mod_l_y <- glm(surv ~ length + tag_year, data = dat, family = quasibinomial("logit"), weights = weight_sc);summary(mod_l_y) # effect of length + tag year
mod_l_y_est <- glm(surv ~ length + tag_year + estuary, data = dat, family = quasibinomial("logit"), weights = weight_sc);summary(mod_l_y_est) # effect of length + tag year + estuary
mod_l_y_est_int <- glm(surv ~ length + tag_year*estuary, data = dat, family = quasibinomial("logit"), weights = weight_sc);summary(mod_l_y_est_int)

# model comparison: forwards fit
### AICc cannot be calculated for quasibinomial distribution, so used forwards fit to determine most parsimonous model
summary(mod_null)
summary(mod_l) # length significant, p < 0.001
summary(mod_l_y) # year significant, p < 0.001
summary(mod_l_y_est) # estuary significant, p = 0.0257
summary(mod_l_y_est_int) # estuary:year interaction term not significant

# model comparison: ANOVA
### same result as forwards fitting above
mod_anova <- anova(mod_null, mod_l, mod_l_y, mod_l_y_est, mod_l_y_est_int, test="F")
mod_anova

### BEST FIT SURVIVAL MODEL: surv ~ length + tag year + estuary
summary(mod_l_y_est)

# PAIRWISE ANNUAL SURVIVAL COMPARISONS ####
#  tukey test to determine whether annual survival rates are significantly different from one another
### 2017 differs from 2018
mod_surv_lsm <- lsmeans(mod_l_y_est, ~tag_year)
contrast(mod_surv_lsm, "tukey")

# PREDICT MARINE SURVIVAL PROBABILITY FOR EACH FISH TAGGED AT THE RST #### 
# standardized lenght inputs for model predictions
length_vals <- seq(min(dat_surv$length_fw[dat_surv$surv == 1]), max(dat_surv$length_fw[dat_surv$surv == 1]), 0.1)

# categorical estuary and year inputs for model predictions
est_0 <- as.factor(rep(0, length(length_vals)))
est_1 <- as.factor(rep(1, length(length_vals)))
year_17 <- as.factor(rep(2017, length(length_vals)))
year_18 <- as.factor(rep(2018, length(length_vals)))

# predict survival probabilities; create data frames with CIs
### 2017 estuary-present
surv_obs_17 <- predict(mod_l_y_est, list(length = length_vals, estuary = est_1, tag_year = year_17), se.fit = T, type = "response")
surv_obs_17 <- as.data.frame(surv_obs_17)
surv_obs_17$lower <- surv_obs_17$fit - (1.96*surv_obs_17$se.fit) 
surv_obs_17$upper <- surv_obs_17$fit + (1.96*surv_obs_17$se.fit) 

### 2018 estuary-present
surv_obs_18 <- predict(mod_l_y_est, list(length = length_vals, estuary = est_1, tag_year = year_18), se.fit = T, type = "response")
surv_obs_18 <- as.data.frame(surv_obs_18)
surv_obs_18$lower <- surv_obs_18$fit - (1.96*surv_obs_18$se.fit) 
surv_obs_18$upper <- surv_obs_18$fit + (1.96*surv_obs_18$se.fit) 

### 2017 estuary-absent
surv_pred_17 <- predict(mod_l_y_est, list(length = length_vals, estuary = est_0, tag_year = year_17), se.fit = T, type = "response")
surv_pred_17 <- as.data.frame(surv_pred_17)
surv_pred_17$lower <- surv_pred_17$fit - (1.96*surv_pred_17$se.fit) 
surv_pred_17$upper <- surv_pred_17$fit + (1.96*surv_pred_17$se.fit) 

### 2018 estuary-absent
surv_pred_18 <- predict(mod_l_y_est, list(length = length_vals, estuary = est_0, tag_year = year_18), se.fit = T, type = "response")
surv_pred_18 <- as.data.frame(surv_pred_18)
surv_pred_18$lower <- surv_pred_18$fit - (1.96*surv_pred_18$se.fit) 
surv_pred_18$upper <- surv_pred_18$fit + (1.96*surv_pred_18$se.fit) 

# SAVE DATA, MODELS, AND PREDICTIONS ####
saveRDS(dat, "dat_surv_obs_pred.rds")
saveRDS(mod_l_y_est, "mod_surv.rds")
saveRDS(surv_obs_17, "surv_obs_17.rds")
saveRDS(surv_obs_18, "surv_obs_18.rds")
saveRDS(surv_pred_17, "surv_pred_17.rds")
saveRDS(surv_pred_18, "surv_pred_18.rds")

# TABLE 1: MARINE SURVIVAL MODEL EFFECT SIZES, CONFIDENCE INTERVALS, P-VALUES ####
# run model without intercept (for interpretability)
mod_l_y_est <- glm(surv ~ 0 + length + tag_year + estuary, data = dat, family = quasibinomial("logit"), weights = weight_sc)

# determine effect sizes, p-values
summary(mod_l_y_est)

# determine confidence intervals
round(confint(mod_l_y_est), 3)

# FIGURE 6: MARINE SURVIVAL ####
# create survival data frames for easy plotting
### estuary-present model
p_surv_obs <- rbind(surv_obs_17, surv_obs_18)
p_surv_obs$length <- rep(length_vals,2)
p_surv_obs$year <- c(as.character(year_17), as.character(year_18))

### estuary-absent model
p_surv_pred <- rbind(surv_pred_17, surv_pred_18)
p_surv_pred$length <- rep(length_vals,2)
p_surv_pred$year <- c(as.character(year_17), as.character(year_18))

# establish plot inputs
### points
pt_size <- 1.5
pt_co <- "grey50"

### lines
line_size <- 1
obs_co <- "#023a45"
pred_co <- "#2baa8e"

### axes
ax_title <- 10
ax_text <- 8

### quantile points
quantile(dat_rst$tag_length, c(0.1, 0.5, 0.9))
q0.1 <- 78
q0.5 <- 93
q0.9 <- 112

# draw figure 6: marine survival plot
fig_6 <- ggplot(p_surv_obs) + 
  facet_grid(. ~ year) + 
  coord_cartesian(xlim = c(60, 142), ylim = c(0, 0.16)) +
  scale_x_continuous(name = "Fork length at at freshwater exit (mm)", breaks = seq(60, 140, 20)) +
  scale_y_continuous(
  
    # Features of the primary y-axis
    name = "Pr(survival)", breaks = seq(0, 0.16, 0.04),
    
    # Features of the secondary y-axis
    sec.axis = sec_axis( trans=~.*6.25, name="Absolute survival", breaks = c(0,1))
  ) +
  
  # estuary-absent model predictions
  geom_line(data = p_surv_pred, aes(x = length, y = fit, color = pred_co), lwd=line_size) +
  geom_ribbon(data = p_surv_pred, aes(ymin = lower, ymax = upper, x = length), fill = adjustcolor(pred_co, 0.4)) +
  
  # add points for size quantiles
  geom_point(data = p_surv_pred[p_surv_pred$length==q0.1,], aes(x = length, y = fit), shape = 21, size = 2.5, color = pred_co, fill = pred_co) +
  geom_point(data = p_surv_pred[p_surv_pred$length==q0.5,], aes(x = length, y = fit), shape = 23, size = 2.5, color = pred_co, fill = pred_co) +
  geom_point(data = p_surv_pred[p_surv_pred$length==q0.9,], aes(x = length, y = fit), shape = 22, size = 2.5, color = pred_co, fill = pred_co) +
  
  # estuary-present model predictions
  geom_line(data = p_surv_obs, aes(x = length, y = fit, color = obs_co), lwd=line_size) +
  geom_ribbon(data = p_surv_obs, aes(ymin = lower, ymax = upper, x = length), fill = adjustcolor(obs_co, 0.4)) +
  
  # add points for size quantiles
  geom_point(data = p_surv_obs[p_surv_obs$length==q0.1,], aes(x = length, y = fit), shape = 21, size = 2.5, color = obs_co, fill = obs_co) +
  geom_point(data = p_surv_obs[p_surv_obs$length==q0.5,], aes(x = length, y = fit), shape = 23, size = 2.5, color = obs_co, fill = obs_co) +
  geom_point(data = p_surv_obs[p_surv_obs$length==q0.9,], aes(x = length, y = fit), shape = 22, size = 2.5, color = obs_co, fill = obs_co) +
  
  # add secondary axis points indicating absolute individual survival (0/1) 
  geom_point(data = dat_surv[dat_surv$surv==0,], aes(x = length_fw, y = 0, alpha = 0.6), position = position_jitter(w = 1, h = 0), col = pt_co, shape = 22, size = 2.5, show.legend = F) +
  geom_point(data = dat_surv[dat_surv$surv==1,], aes(x = length_fw, y = 0.16, alpha = 0.6), position = position_jitter(w = 1, h = 0), col = pt_co, shape = 22, size = 2.5, show.legend = F) +
  scale_color_identity(name = "Model",
                       breaks = c("#023a45", "#2baa8e"), 
                       labels = c("Marine survival w/ estuary rearing (obs)", 
                                  "Marine survival w/o estuary rearing (pred)"), 
                       guide = "legend") +
  
  # map theme
  theme_classic() +
  theme(axis.title = element_text(size = ax_title, face = "bold"), 
        axis.text = element_text(size = ax_text, color = "black"),
        axis.line = element_line(color = NA),
        strip.background = element_rect(fill = "grey70"),
        strip.text = element_text(size = ax_title, face = "bold"),
        panel.border = element_rect(fill = NA, size = 1, color = "black"),
        legend.title = element_text(size = ax_text, color = "black"),
        legend.text = element_text(size = ax_text-0.65, color = "black"),
        legend.position = c(0.25, 0.8))

fig_6

# save marine survival plot
ggsave("fig_6.pdf", dpi = 300, height = 4.5, width = 6, units = "in")

