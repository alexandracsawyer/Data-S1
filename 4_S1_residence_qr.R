# LOAD REQUIRED PACKAGES ####
library(quantreg)
library(Qtools)

# LOAD DATA ####
# ensure PIT tag numbers do not present in scientific notation
options(scipen=999)

# load comprehensive mark-recapture data set
### individual fish are identified by PIT tag number
### columns preceded by "tag_" indicate initial mark data
### columns preceded by "rc_" indicate recapture data
dat_rc <- read.csv("dat_rc.csv")
rownames(dat_rc) <- c()

# restrict mark-recapture data set to fish tagged at the RST (freshwater exit)
dat_rc <- subset(dat_rc, tag_loc == "smolt trap")

# RUN CANDIDATE QUANTILE REGRESSION MODELS ####
### tau argument identifies data quantile
### tested models for multiple tau levels 0.75-1.00, representing upper quantile of estuary residence values
mod_qr <- rq(formula = log(days) ~ tag_length, tau = seq(0.75, 1, 0.05), data = dat_rc)
mod_qr

# COMPARE CANDIDATE MODELS (TAU 0.75-1.0)
### create data frame to store model comparison values
mod_compare <- data.frame("tau" = c("0.75", "0.8", "0.85", "0.9", "0.95", "1.0"))

### calculate  and store AIC for all models
### tau = 0.75 has lowest AIC
mod_compare$AIC <- AIC.rq(mod_qr)

### calculate AICc from AIC for all models
### tau = 0.75 has lowest AICc
n <- nrow(dat_rc) # sample size
K <- 1 # number of free parameters in the model
mod_compare$AICc <- mod_compare$AIC * (n/(n-K-1))

### tau goodness of fit test: a large test statistic (small p value) is evidence of lack of fit
### tau = 0.75 has small test statistic and large p value
gof <- GOFTest(mod_qr, alpha=0.05, B = 1000)
mod_compare$gof_ts <- gof[[1]]$Tn
mod_compare$gof_p <- gof[[1]]$p.value

mod_compare

### BEST FIT QUANTILE REGRESSION MODEL: tau = 0.75
### rerun standalone 0.75 QR model
mod_qr75 <- rq(formula = log(days) ~ tag_length, tau = 0.75, data = dat_rc) 
mod_qr75

# PREDICT ESTUARY RESIDENCE VALUES ####
# standardized length inputs for model predictions
length_vals <- seq(min(dat_rc$tag_length), max(dat_rc$tag_length), 1)

# predict estuary residence values from best-fit model
days_qr75 <- predict.rq(mod_qr75, list(tag_length=length_vals), interval = c("confidence"), level = 0.95)

# create data frame with predicted estuary residence days for fish of a given length at tagging
### exponentiate model predictions to translate from log space into real space
days_qr75 <- as.data.frame(days_qr75) # create data frame from model predictions
days_qr75$fit <- exp(days_qr75$fit) # exponentiate model fit
days_qr75$lower <- exp(days_qr75$lower) # exponentiate lower 95% CI
days_qr75$higher <- exp(days_qr75$higher) # exponentiate upper 95% CI
days_qr75$tag_length <- length_vals # append length values to data frame

# SAVE BEST-FIT MODEL AND MODEL PREDICTIONS ####
saveRDS(mod_qr75, "mod_res_qr.rds")
write.csv(days_qr75, "days_res_qr.csv")

# FIGURE S1 - PANEL B: ESTUARY RESIDENCE (QUANTILE REGRESSION) ####
# load recapture data
dat_rc <- read.csv("dat_rc.csv")
dat_rc <- subset(dat_rc, tag_loc == "smolt trap") # keep only fish tagged at the RST

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
fig_s1_b <- ggplot() + 
  geom_point(data = dat_rc, aes(x = tag_length, y = days), col = adjustcolor(pt_co, 0.6), size = pt_size) +
  labs(x = "Fork length at freshwater exit (mm)", y = "Days in estuary") +
  coord_cartesian(ylim = c(0, 50)) +
  scale_x_continuous(breaks=seq(70, 120, 10)) +
  geom_line(data = days_qr75, aes(x = tag_length, y = fit), col = line_co, lwd = line_size) +
  geom_ribbon(data = days_qr75, aes(x = tag_length, ymin = lower, ymax = higher), fill = adjustcolor(line_co, 0.4), inherit.aes = F) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = ax_title), 
        axis.text = element_text(color = "black", size = ax_text))

fig_s1_b
