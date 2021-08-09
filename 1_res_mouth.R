# LOAD REQUIRED PACKAGES ####
library(bbmle) 
library(MASS) 
library(ggplot2)

# LOAD DATA ####
# ensure PIT tag numbers do not present in scientific notation
options(scipen=999)

# load comprehensive mark-recapture data set
### individual fish are identified by PIT tag number
### columns preceded by "tag_" indicate initial mark data
### columns preceded by "rc_" indicate recapture data
dat_rc <- readRDS("dat_rc.rds")
rownames(dat_rc) <- c()

# create Reach 1 mark-recapture data set
### keep only fish tagged at the RST and recaptured in Reach 1
### Reach 1 is a proxy for marine entry
dat <- dat_rc[which(dat_rc$tag_loc == "smolt trap" & dat_rc$rc_loc_1 == "1"),]

### RESPONSE VARAIABLE: observed estuary residence (days)
### CANDIDATE PREDICTORS: length at tagging (tag_length), tag year (tag_year), tag date (tag_doy), and daily growth rate in the estuary (growth_d)

# RUN CANDIDATE LINEAR MODELS ####
### one predictor
mod_lm_l <- glm(days ~ tag_length, data=dat);summary(mod_lm_l) 
mod_lm_y <- glm(days ~ tag_year, data=dat);summary(mod_lm_y) 
mod_lm_doy <- glm(days ~ tag_doy, data=dat);summary(mod_lm_doy) 
mod_lm_gd <- glm(days ~ growth_d, data=dat);summary(mod_lm_gd) 

### two predictors
mod_lm_l_y <- glm(days ~ tag_length + tag_year, data=dat);summary(mod_lm_l_y) 
mod_lm_l_doy <- glm(days ~ tag_length + tag_doy, data=dat);summary(mod_lm_l_doy) 
mod_lm_l_gd <- glm(days ~ tag_length + growth_d, data=dat);summary(mod_lm_l_gd) 
mod_lm_y_doy <- glm(days ~ tag_year + tag_doy, data=dat);summary(mod_lm_y_doy) 
mod_lm_y_gd <- glm(days ~ tag_year + growth_d, data=dat);summary(mod_lm_y_gd) 
mod_lm_doy_gd <- glm(days ~ tag_doy + growth_d, data=dat);summary(mod_lm_doy_gd) 

#### three predictors
mod_lm_l_y_doy <- glm(days ~ tag_length + tag_year + tag_doy, data=dat);summary(mod_lm_l_y_doy) 
mod_lm_l_y_gd <- glm(days ~ tag_length + tag_year + growth_d, data=dat);summary(mod_lm_l_y_gd) 
mod_lm_l_doy_gd <- glm(days ~ tag_length + tag_doy + growth_d, data=dat);summary(mod_lm_l_doy_gd) 
mod_lm_y_doy_gd <- glm(days ~ tag_year + tag_doy + growth_d, data=dat);summary(mod_lm_y_doy_gd) 

### four predictors
mod_lm_all <- glm(days ~ tag_length + tag_year + tag_doy + growth_d, data=dat);summary(mod_lm_all) 

# determine best-fit linear model
### aggregate candidate models into list
mod_list <- list("mod_lm_l" = mod_lm_l, # effect of length at tagging
                 "mod_lm_y" = mod_lm_y, # effect of tag year
                 "mod_lm_doy" = mod_lm_doy, # effect of tag date
                 "mod_lm_gd" = mod_lm_gd, # effect of daily growth rate
                 "mod_lm_l_y" = mod_lm_l_y, # effect of length at tagging + tag year
                 "mod_lm_l_doy" = mod_lm_l_doy, # effect of length at tagging + tag date
                 "mod_lm_l_gd" = mod_lm_l_gd, # effect of length at tagging + daily growth rate
                 "mod_lm_y_doy" = mod_lm_y_doy, # effect of tag year + tag date
                 "mod_lm_y_gd" = mod_lm_y_gd, # effect of tag year + daily growth rate
                 "mod_lm_doy_gd" = mod_lm_doy_gd, # effect of tag date + daily growth rate 
                 "mod_lm_l_y_doy" = mod_lm_l_y_doy, # effect of length at tagging + tag year + tag date
                 "mod_lm_l_y_gd" = mod_lm_l_y_gd, # effect of length at tagging + tag year + daily growth rate
                 "mod_lm_l_doy_gd" = mod_lm_l_doy_gd, # effect of length at tagging + tag date + daily growth rate
                 "mod_lm_y_doy_gd" = mod_lm_y_doy_gd, # effect of tag year + tag date + daily growth rate
                 "mod_lm_all" = mod_lm_all # effect of length at tagging + tag year + tag date + daily growth rate
                 )

### compare models using AICc
AICctab(mod_list, delta = TRUE, weights = TRUE, base = TRUE) 

### BEST FIT LINEAR MODEL: days ~ tag_length
summary(mod_lm_l)
              
# RUN CANDIDATE GENERALIZED LINEAR MODELS - POISSON DISTRIBUTION ####
### one predictor
mod_pois_l <- glm(days ~ tag_length, data=dat, family = poisson(link="log"));summary(mod_pois_l) 
mod_pois_y <- glm(days ~ tag_year, data=dat, family = poisson(link="log"));summary(mod_pois_y) 
mod_pois_doy <- glm(days ~ tag_doy, data=dat, family = poisson(link="log"));summary(mod_pois_doy) 
mod_pois_gd <- glm(days ~ growth_d, data=dat, family = poisson(link="log"));summary(mod_pois_gd) 
### two predictors
mod_pois_l_y <- glm(days ~ tag_length + tag_year, data=dat, family = poisson(link="log"));summary(mod_pois_l_y) 
mod_pois_l_doy <- glm(days ~ tag_length + tag_doy, data=dat, family = poisson(link="log"));summary(mod_pois_l_doy) 
mod_pois_l_gd <- glm(days ~ tag_length + growth_d, data=dat, family = poisson(link="log"));summary(mod_pois_l_gd) 
mod_pois_y_doy <- glm(days ~ tag_year + tag_doy, data=dat, family = poisson(link="log"));summary(mod_pois_y_doy) 
mod_pois_y_gd <- glm(days ~ tag_year + growth_d, data=dat, family = poisson(link="log"));summary(mod_pois_y_gd) 
mod_pois_doy_gd <- glm(days ~ tag_doy + growth_d, data=dat, family = poisson(link="log"));summary(mod_pois_doy_gd) 

### three predictors
mod_pois_l_y_doy <- glm(days ~ tag_length + tag_year + tag_doy, data=dat, family = poisson(link="log"));summary(mod_pois_l_y_doy) 
mod_pois_l_y_gd <- glm(days ~ tag_length + tag_year + growth_d, data=dat, family = poisson(link="log"));summary(mod_pois_l_y_gd) 
mod_pois_l_doy_gd <- glm(days ~ tag_length + tag_doy + growth_d, data=dat, family = poisson(link="log"));summary(mod_pois_l_doy_gd) 
mod_pois_y_doy_gd <- glm(days ~ tag_year + tag_doy + growth_d, data=dat, family = poisson(link="log"));summary(mod_pois_y_doy_gd) 

### four predictors
mod_pois_all <- glm(days ~ tag_length + tag_year + tag_doy + growth_d, data=dat, family = poisson(link="log"));summary(mod_pois_all) 

# determine best-fit poisson-distributed GLM
### aggregate candidate models into list
mod_list <- list("mod_pois_l" = mod_pois_l, # effect of length at tagging
                 "mod_pois_y" = mod_pois_y, # effect of tag year
                 "mod_pois_doy" = mod_pois_doy, # effect of tag date
                 "mod_pois_gd" = mod_pois_gd, # effect of daily growth rate
                 "mod_pois_l_y" = mod_pois_l_y, # effect of length at tagging + tag year
                 "mod_pois_l_doy" = mod_pois_l_doy, # effect of length at tagging + tag date
                 "mod_pois_l_gd" = mod_pois_l_gd, # effect of length at tagging + daily growth rate
                 "mod_pois_y_doy" = mod_pois_y_doy, # effect of tag year + tag date
                 "mod_pois_y_gd" = mod_pois_y_gd, # effect of tag year + daily growth rate
                 "mod_pois_doy_gd" = mod_pois_doy_gd, # effect of tag date + daily growth rate 
                 "mod_pois_l_y_doy" = mod_pois_l_y_doy, # effect of length at tagging + tag year + tag date
                 "mod_pois_l_y_gd" = mod_pois_l_y_gd, # effect of length at tagging + tag year + daily growth rate
                 "mod_pois_l_doy_gd" = mod_pois_l_doy_gd, # effect of length at tagging + tag date + daily growth rate
                 "mod_pois_y_doy_gd" = mod_pois_y_doy_gd, # effect of length at tagging + tag year + tag date + daily growth rate
                 "mod_pois_all" = mod_pois_all # effect of length at tagging + tag year + tag date + daily growth rate
                 )

### compare models using AICc
AICctab(mod_list, delta = TRUE, weights = TRUE, base = TRUE) 

### BEST FIT POISSON GENERALIZED LINEAR MODEL: days ~ tag_length + tag_year
summary(mod_pois_l_y)

# RUN CANDIDATE GENERALIZED LINEAR MODELS - NEGATIVE BINOMIAL DISTRIBUTION ####
### one predictor
mod_nb_l <- glm.nb(days ~ tag_length, data=dat);summary(mod_nb_l) 
mod_nb_y <- glm.nb(days ~ tag_year, data=dat);summary(mod_nb_y) 
mod_nb_doy <- glm.nb(days ~ tag_doy, data=dat);summary(mod_nb_doy) 
mod_nb_gd <- glm.nb(days ~ growth_d, data=dat);summary(mod_nb_gd) 

# two predictors
mod_nb_l_y <- glm.nb(days ~ tag_length + tag_year, data=dat);summary(mod_nb_l_y)
mod_nb_l_doy <- glm.nb(days ~ tag_length + tag_doy, data=dat);summary(mod_nb_l_doy)
mod_nb_l_gd <- glm.nb(days ~ tag_length + growth_d, data=dat);summary(mod_nb_l_gd)
mod_nb_y_doy <- glm.nb(days ~ tag_year + tag_doy, data=dat);summary(mod_nb_y_doy)
mod_nb_y_gd <- glm.nb(days ~ tag_year + growth_d, data=dat);summary(mod_nb_y_gd)
mod_nb_doy_gd <- glm.nb(days ~ tag_doy + growth_d, data=dat);summary(mod_nb_doy_gd)

# three predictors
mod_nb_l_y_doy <- glm.nb(days ~ tag_length + tag_year + tag_doy, data=dat);summary(mod_nb_l_y_doy) 
mod_nb_l_y_gd <- glm.nb(days ~ tag_length + tag_year + growth_d, data=dat);summary(mod_nb_l_y_gd) 
mod_nb_l_doy_gd <- glm.nb(days ~ tag_length + tag_doy + growth_d, data=dat);summary(mod_nb_l_doy_gd) 
mod_nb_y_doy_gd <- glm.nb(days ~ tag_year + tag_doy + growth_d, data=dat);summary(mod_nb_y_doy_gd) 

# four predictors
mod_nb_all <- glm.nb(days ~ tag_length + tag_year + tag_doy + growth_d, data=dat);summary(mod_nb_all) 

# determine best-fit negative binomial-distributed GLM
### aggregate candidate models into list
mod_list <- list("mod_nb_l" = mod_nb_l, # effect of length at tagging
                 "mod_nb_y" = mod_nb_y, # effect of tag year
                 "mod_nb_doy" = mod_nb_doy, # effect of tag date
                 "mod_nb_gd" = mod_nb_gd, # effect of daily growth rate
                 "mod_nb_l_y" = mod_nb_l_y, # effect of length at tagging + tag year
                 "mod_nb_l_doy" = mod_nb_l_doy, # effect of length at tagging + tag date
                 "mod_nb_l_gd" = mod_nb_l_gd, # effect of length at tagging + daily growth rate
                 "mod_nb_y_doy" = mod_nb_y_doy, # effect of tag year + tag date
                 "mod_nb_y_gd" = mod_nb_y_gd, # effect of tag year + daily growth rate
                 "mod_nb_doy_gd" = mod_nb_doy_gd, # effect of tag date + daily growth rate 
                 "mod_nb_l_y_doy" = mod_nb_l_y_doy, # effect of length at tagging + tag year + tag date
                 "mod_nb_l_y_gd" = mod_nb_l_y_gd, # effect of length at tagging + tag year + daily growth rate
                 "mod_nb_l_doy_gd" = mod_nb_l_doy_gd, # effect of length at tagging + tag date + daily growth rate
                 "mod_nb_y_doy_gd" = mod_nb_y_doy_gd, # effect of tag year + tag date + daily growth rate
                 "mod_nb_all" = mod_nb_all # effect of length at tagging + tag year + tag date + daily growth rate
                 )

### compare models using AICc
AICctab(mod_list, delta = TRUE, weights = TRUE, base = TRUE) 

### BEST FIT NEGATIVE BINOMIAL GENERALIZED LINEAR MODEL: days ~ tag_length
summary(mod_nb_l) 

# COMPARE CANDIDATE MODELS - ALL DISTRIBUTIONS ####
# aggregate all candidate models into list
mod_list <- list("mod_lm_l" = mod_lm_l, 
                 "mod_lm_y" = mod_lm_y,
                 "mod_lm_doy" = mod_lm_doy,
                 "mod_lm_gd" = mod_lm_gd,
                 "mod_lm_l_y" = mod_lm_l_y,
                 "mod_lm_l_doy" = mod_lm_l_doy,
                 "mod_lm_l_gd" = mod_lm_l_gd,
                 "mod_lm_y_doy" = mod_lm_y_doy,
                 "mod_lm_y_gd" = mod_lm_y_gd,
                 "mod_lm_doy_gd" = mod_lm_doy_gd,
                 "mod_lm_l_y_doy" = mod_lm_l_y_doy,
                 "mod_lm_l_y_gd" = mod_lm_l_y_gd,
                 "mod_lm_l_doy_gd" = mod_lm_l_doy_gd,
                 "mod_lm_y_doy_gd" = mod_lm_y_doy_gd,
                 "mod_lm_all" = mod_lm_all,
                 "mod_pois_l" = mod_pois_l, 
                 "mod_pois_y" = mod_pois_y,
                 "mod_pois_doy" = mod_pois_doy,
                 "mod_pois_gd" = mod_pois_gd,
                 "mod_pois_l_y" = mod_pois_l_y,
                 "mod_pois_l_doy" = mod_pois_l_doy,
                 "mod_pois_l_gd" = mod_pois_l_gd,
                 "mod_pois_y_doy" = mod_pois_y_doy,
                 "mod_pois_y_gd" = mod_pois_y_gd,
                 "mod_pois_doy_gd" = mod_pois_doy_gd,
                 "mod_pois_l_y_doy" = mod_pois_l_y_doy,
                 "mod_pois_l_y_gd" = mod_pois_l_y_gd,
                 "mod_pois_l_doy_gd" = mod_pois_l_doy_gd,
                 "mod_pois_y_doy_gd" = mod_pois_y_doy_gd,
                 "mod_pois_all" = mod_pois_all,
                 "mod_nb_l" = mod_nb_l, 
                 "mod_nb_y" = mod_nb_y,
                 "mod_nb_doy" = mod_nb_doy,
                 "mod_nb_gd" = mod_nb_gd,
                 "mod_nb_l_y" = mod_nb_l_y,
                 "mod_nb_l_doy" = mod_nb_l_doy,
                 "mod_nb_l_gd" = mod_nb_l_gd,
                 "mod_nb_y_doy" = mod_nb_y_doy,
                 "mod_nb_y_gd" = mod_nb_y_gd,
                 "mod_nb_doy_gd" = mod_nb_doy_gd,
                 "mod_nb_l_y_doy" = mod_nb_l_y_doy,
                 "mod_nb_l_y_gd" = mod_nb_l_y_gd,
                 "mod_nb_l_doy_gd" = mod_nb_l_doy_gd,
                 "mod_nb_y_doy_gd" = mod_nb_y_doy_gd,
                 "mod_nb_all" = mod_nb_all)

### compare models using AICc
AICctab(mod_list, delta = TRUE, weights = TRUE, base = TRUE) 

### BEST FIT MODEL IS NEGATIVE BINOMIAL GLM: days ~ tag_length
summary(mod_nb_l)
### negative binomial models receive vast majority of weight
### dAICc of next best model is >2

# PREDICT ESTUARY RESIDENCE VALUES ####
# standardized length inputs for model predictions
length_vals <- seq(min(dat$tag_length), max(dat$tag_length), 1)

# predict estuary residence values from best-fit model
days_nb_l <- predict.glm(mod_nb_l, se.fit=T, list(tag_length=length_vals), interval=c("confidence"), level=0.95)

# predict upper and lower 95% confidence intervals
days_nb_l$lwr <- days_nb_l$fit - (1.96*days_nb_l$se.fit)
days_nb_l$upr <- days_nb_l$fit + (1.96*days_nb_l$se.fit)

# create data frame with predicted estuary residence days for fish of a given length at tagging
### exponentiate model predictions to translate from log space into real space
days_nb_l <- as.data.frame(days_nb_l) # create data frame from model predictions
days_nb_l$tag_length <- length_vals # add length values used to predict estuary residence
days_nb_l$fit <- exp(days_nb_l$fit) # exponentiate model fit
days_nb_l$lwr <- exp(days_nb_l$lwr) # exponentiate lower 95% CI
days_nb_l$upr <- exp(days_nb_l$upr) # exponentiate upper 95% CI

# SAVE BEST-FIT MODEL AND MODEL PREDICTIONS ####
saveRDS(mod_nb_l, "mod_res_glm.rds")
saveRDS(days_nb_l, "days_res_glm.rds")

# FIGURE 3: ESTUARY RESIDENCE ####
# establish plot inputs
### points
pt_size <- 3
pt_co <- "grey"

### lines
line_size <- 1
line_co <- "#023a45"

### axes
ax_title <- 10
ax_text <- 8

# draw figure 3: estuary residence plot
fig_3 <- ggplot() + 
  geom_point(data = dat, aes(x = tag_length, y = days), col = adjustcolor(pt_co, 0.6), size = pt_size) +
  labs(x = "Fork length at freshwater exit (mm)", y = "Days in estuary") +
  # coord_cartesian(ylim = c(0, 50)) +
  scale_x_continuous(breaks=seq(70, 120, 10)) +
  geom_line(data = days_nb_l, aes(x = tag_length, y = fit), col = line_co, lwd = line_size) +
  geom_ribbon(data = days_nb_l, aes(x = tag_length, ymin = lwr, ymax = upr), fill = adjustcolor(line_co, 0.4), inherit.aes = F) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = ax_title), 
        axis.text = element_text(color = "black", size = ax_text))

fig_3

# save estuary residence plot
ggsave("fig_3.pdf", dpi = 300, height = 3, width = 3, units = "in")
