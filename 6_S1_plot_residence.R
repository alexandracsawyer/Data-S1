# LOAD REQUIRED PACKAGES ####
library(ggplot2)
library(patchwork)

# LOAD DATA ####
dat_rst <- readRDS("dat_rst.rds")
dat_rc <- readRDS("dat_rc.rds")

# load negative binomial GLM model and predictions
mod_res_glm <- readRDS("mod_res_glm.rds")
days_res_glm <- readRDS("days_res_glm.rds")

# load 75% quantile regression model and predictions
mod_res_qr <- readRDS("mod_res_qr.rds")
days_res_qr <- readRDS("days_res_qr.rds")

# load CJS model
mod_res_cjs <- readRDS("mod_res_cjs.rds")
days_res_cjs <- readRDS("days_res_cjs.rds")

# TABLE S1: RESIDENCE MODEL PARAMETER ESTIMATES, RESIDENCE VALUES, CONFIDENCE INTERVALS ####
#establish quantiles
quantile(dat_rst$length, c(0.1, 0.5, 0.9))

q0.1 <- 78
q0.5 <- 93
q0.9 <- 112

# extract parameter estimates
round(mod_res_glm$coefficients, 3) # negative binomial GLM
round(mod_res_qr$coefficients, 3) # 75% quantile regression
round(mod_res_cjs$results$beta$Phi, 3) # CJS mark-recapture model

# extract 1) estuary residence estimates 2) lower 95% confidence bound, and 3) upper 95% confidence bound by size quantile (0.1, 0.5, 0.9)
### negative binomial GLM
round(days_res_glm$fit[days_res_glm$tag_length==q0.1 | days_res_glm$tag_length==q0.5 | days_res_glm$tag_length==q0.9], 1)
round(days_res_glm$lwr[days_res_glm$tag_length==q0.1 | days_res_glm$tag_length==q0.5 | days_res_glm$tag_length==q0.9], 1)
round(days_res_glm$upr[days_res_glm$tag_length==q0.1 | days_res_glm$tag_length==q0.5 | days_res_glm$tag_length==q0.9], 1)

### 75% quantile regression
round(days_res_qr$fit[days_res_qr$tag_length==q0.1 | days_res_qr$tag_length==q0.5 | days_res_qr$tag_length==q0.9], 1)
round(days_res_qr$lower[days_res_qr$tag_length==q0.1 | days_res_qr$tag_length==q0.5 | days_res_qr$tag_length==q0.9], 1)
round(days_res_qr$higher[days_res_qr$tag_length==q0.1 | days_res_qr$tag_length==q0.5 | days_res_qr$tag_length==q0.9], 1)

### CJS mark-recapture model
### no CIs for CJS model; error did not converge well due to low recapture rate
round(days_res_cjs$days[days_res_cjs$tag_length==q0.1 | days_res_cjs$tag_length==q0.5 | days_res_cjs$tag_length==q0.9], 1)

# FIGURE S1: COMPARISON OF ESTUARY RESIDENCE ESTIMATES ####
# subset data sets for 1) reach 1 and 2) reach 2-5 recaptures
dat_R1 <- subset(dat_rc, tag_loc == "smolt trap" & rc_loc_1 == 1)
dat_all <- subset(dat_rc, tag_loc == "smolt trap" & (rc_loc_1 == 2 | 
                                                       rc_loc_1 == 3 | 
                                                       rc_loc_1 == 4 | 
                                                       rc_loc_1 == 5))

# create standardized length values based on subset data sets (for easy plotting)
length_vals_R1 <- seq(min(dat_R1$tag_length), max(dat_R1$tag_length), 1)
length_vals_all <- seq(min(dat_all$tag_length), max(dat_all$tag_length), 1) 
days_res_cjs <- days_res_cjs[which(days_res_cjs$tag_length %in% length_vals_all),]

# establish plot inputs
### points
pt_size <- 3
pt_co <- "grey"

### lines
line_size <- 1
line_co <- "#023a45"

### axes & labels
ax_title <- 8
ax_text <- 8
sub_title <- 2.75


# draw base plot: Reach 1 recaptures for negative binomial GLM
p <- ggplot() + 
  geom_point(data = dat_all, aes(x = tag_length, y = days), col = adjustcolor(pt_co, 0.6), fill = adjustcolor("white", 0.0), size = pt_size, shape = 21) +
  geom_point(data = dat_R1, aes(x = tag_length, y = days), col = adjustcolor(pt_co, 0.6), size = pt_size) +
  coord_cartesian(xlim =c(65, 125), ylim = c(0, 50)) +
  scale_x_continuous(breaks=seq(70, 120, 10)) +
  theme_classic() +
  theme(axis.title = element_blank())

# draw base plot: all recaptures for quantile regression and CJS mark-recapture models
p2 <- ggplot() + 
  geom_point(data = dat_R1, aes(x = tag_length, y = days), col = adjustcolor(pt_co, 0.6), size = pt_size) +
  geom_point(data = dat_all, aes(x = tag_length, y = days), col = adjustcolor(pt_co, 0.6), size = pt_size) +
  coord_cartesian(xlim =c(65, 125), ylim = c(0, 50)) +
  scale_x_continuous(breaks=seq(70, 120, 10)) +
  theme_classic() + 
  theme(axis.title = element_blank())

# draw negative binonial GLM residence plot
plot_glm <- p + 
  geom_line(data = days_res_glm, aes(x = length_vals_R1, y = fit), col = line_co, lwd = line_size) +
  geom_ribbon(data = days_res_glm, aes(x = length_vals_R1, ymin = lwr, ymax = upr), fill = adjustcolor(line_co, 0.4), inherit.aes = F) + 
  labs(x = "Length at freshwater exit (mm)", y = "Days in estuary") +
  annotate("text", x = 93, y = 50, label = "a) Negative binomial GLM", size = sub_title) +
  theme(axis.title.y = element_text(face = "bold", size = ax_title), 
        axis.title.x = element_blank(),
        axis.text = element_text(size = ax_text, color = "black"),
        plot.title = element_text(size = ax_text, color = "black"))

# 75% quantile regression residence plot
plot_qr <- p2 + 
  geom_line(data = days_res_qr, aes(x = length_vals_all, y = fit), col = line_co, lwd = line_size) +
  geom_ribbon(data = days_res_qr, aes(x = length_vals_all, ymin = lower, ymax = higher), fill = adjustcolor(line_co, 0.4), inherit.aes = F) +
  labs(x = "Length at freshwater exit (mm)", y = "Days in estuary") +
  annotate("text", x = 94, y = 50, label = "b) 75% quantile regression", size = sub_title) +
  theme(axis.title.x = element_text(face = "bold", size = ax_title), 
        axis.title.y = element_blank(),
        axis.text = element_text(size = ax_text, color = "black"), 
        plot.title = element_text(size = ax_text, color = "black"))

# cjs residence plot
plot_cjs <- p2 +
  geom_line(data = days_res_cjs, aes(x = length_vals_all, y = days_cjs), col = line_co, lwd = line_size) +
  labs(x = "Length at freshwater exit (mm)", y = "Days in estuary") +
  # ggtitle("c) CJS mark-recapture") +
  annotate("text", x = 90, y = 50, label = "c) CJS mark-recapture", size = sub_title) +
  theme(axis.title = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = ax_text, color = "black"),
        plot.title = element_text(size = ax_text, color = "black"))

# plot with all models
plot_glm + plot_qr + plot_cjs

# save comparative estuary residence plot
ggsave("fig_s1.pdf", dpi = 300, height = 4.5, width = 6, units = "in")

