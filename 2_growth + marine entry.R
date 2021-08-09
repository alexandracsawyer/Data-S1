# LOAD REQUIRED PACKAGES ####
library(bbmle)
library(lsmeans)
library(ggplot2)
library(patchwork)

# LOAD DATA ####
# ensure PIT tag numbers do not present in scientific notation
options(scipen=999)

# load comprehensive mark-recapture data set (only tagged fish that were recaptured in the estuary)
### individual fish are identified by PIT tag number
### columns preceded by "tag_" indicate initial mark data
### columns preceded by "rc_" indicate recapture data
dat_rc <- readRDS("1_data/processed/dat_rc.rds")
rownames(dat_rc) <- c()

# load comprehensive mark data set (all fish that were tagged at the RST)
### individual fish are identified by PIT tag number
dat_rst <- readRDS("1_data/processed/dat_rst.rds")
rownames(dat_rst) <- c()

# LOAD ESTUARY RESIDENCE MODEL ####
mod_res <- readRDS("2_analyses/0_final/mod_res_glm.rds")

# PREPARE DATA FOR GROWTH RATE MODELS ####
# mark-recapture data set 
### keep only fish tagged at the RST
dat_rc <- subset(dat_rc, tag_loc == "smolt trap")

### identify growth rate anomalies
### remove fish whose growth exceeds 2mm/day -- clear measurement error
plot(dat_rc$growth_d)
dat_rc <- dat_rc[-which(dat_rc$growth_d > 2),] 

### drop unused tag year factors
dat_rc$tag_year <- droplevels(dat_rc$tag_year) 

### keep only columns relevant for analyses
dat_rc <- dat_rc[,c("PIT_number", "tag_year", "tag_doy", "tag_length", "rc_doy", "rc_length", "days", "growth", "growth_d")] 

# RST data set
### drop unused tag year factors
dat_rst$tag_year <- droplevels(dat_rst$tag_year)

### keep only columns relevant for analyses
### adjust row and column names for consistency with small data set
dat_rst <- dat_rst[,c("PIT_number", "tag_year", "doy", "length")]
colnames(dat_rst) <- c("PIT_number", "tag_year", "tag_doy", "tag_length")

# RUN CANDIDATE GROWTH RATE MODELS ####
# null model - constant growth rate
mod_growth_null <- lm(growth_d ~ 1, data = dat_rc); summary(mod_growth_null)

# growth rate model with effect of tag year
mod_growth_y <- lm(growth_d ~ 0 + factor(tag_year), data = dat_rc); summary(mod_growth_y)

# determine best-fit annual growth rate model
### aggregate candidate models into list
mod_list <- list("mod_growth_null" = mod_growth_null, # null model - constant growth rate
                 "mod_growth_y" = mod_growth_y # effect of tag year
                 )

### compare models using AICc
AICctab(mod_list, delta = TRUE, weights = TRUE, base = TRUE)

### BEST-FIT GROWTH RATE MODEL: growth ~ tag_year
summary(mod_growth_y)

# PAIRWISE ANNUAL GROWTH RATE COMPARISONS ####
#  tukey test to determine whether annual growth rates are significantly different from one another
### 2017 differs from 2018; p = 0.0264
### 2017 does not differ from 2019; p = 0.9577
### 2018 differs from 2019; p = 0.0032
mod_growth_lsm <- lsmeans(mod_growth_y, ~factor(tag_year))
contrast(mod_growth_lsm, "tukey")

# PREDICT ESTUARY RESIDENCE FOR ALL FISH TAGGED AT THE RST ####
# predict estuary residence for each fish
pred_days <- predict.glm(mod_res, list(tag_length = dat_rst$tag_length), se.fit = T, interval = "confidence", level = 0.95) 
pred_days <- as.data.frame(pred_days)

# append predicted estuary residence to RST data set
### exponentiate model predictions to translate from log space into real space
dat_rst$pred_days <- exp(pred_days$fit)
dat_rst$pred_days_sd <- exp(pred_days$se.fit)

# PREDICT ESTUARY GROWTH FOR ALL FISH TAGGED AT THE RST ####
# standardized years inputs for model predictions
year_vals <- as.factor(c(2017, 2018, 2019))

# predict annual growth rate for each fish 
pred_growth_d <- predict.lm(mod_growth_y, list(tag_year = dat_rst$tag_year), se.fit = T, interval = "confidence", level = 0.95)
pred_growth_d <- as.data.frame(pred_growth_d)

# calculate standard deviation in annual mean growth rate
### append tag year to data frame
### N is number of fish recaptured in each sampling year, for use in calculating SD from SE
pred_growth_d$tag_year <- dat_rst$tag_year
N <- table(dat_rc$tag_year)
pred_growth_d$sd <- pred_growth_d$se.fit * sqrt(N[match(pred_growth_d$tag_year, names(N))])

# append predicted estuary growth to RST data set
### predicted growth is predicted daily growth * predicted estuary residence
### predicted growth SD incorporates uncertainy associated with both predicted annual growth rates and predicted estuary residence
dat_rst$pred_growth <- pred_growth_d$fit.fit * dat_rst$pred_days
dat_rst$pred_growth_sd <- sqrt((pred_growth_d$sd/pred_growth_d$fit.fit)^2 + (dat_rst$pred_days_sd/dat_rst$pred_days)^2)

# PREDICT LENGTH AT MARINE ENTRY FOR ALL FISH TAGGED AT THE RST
# calculate predicted length at marine entry
### add predicted estuary growth to length at RST tagging
dat_rst$pred_length <- dat_rst$tag_length + dat_rst$pred_growth

# SAVE DATA ####
saveRDS(mod_growth_y, "2_analyses/0_final/mod_growth_d.rds")
saveRDS(dat_rst, "1_data/processed/dat_rst_mar.rds")

# FIGURE 4: ANNUAL ESTUARY GROWTH RATES ####
# create growth data frame for easy plotting
dat_growth <- data.frame("year" = c(2017, 2018, 2019))
dat_growth$growth_d <- mod_growth_y$coefficients
dat_growth$lower <- confint(mod_growth_y)[1:3]
dat_growth$upper <- confint(mod_growth_y)[4:6]

# make years numeric for plotting
dat_rc$tag_year <- as.numeric(as.character(dat_rc$tag_year))

# establish plot inputs
### points
pt_size <- 1.5
pt_co <- "grey"

# lines
line_size <- 1
fw_co <- "#2baa8e"
mar_co <- "#023a45"

# axes
ax_title <- 10
ax_text <- 8

# draw figure 4: estuary growth plot
fig_4 <- ggplot() +
  geom_jitter(data = dat_rc, aes(x = tag_year, y = growth_d), width = 0.15, color = adjustcolor(pt_co, 0.6), size = pt_size) +
  geom_point(data = dat_growth, aes(x = year, y = growth_d), col = mar_co, fill = adjustcolor("white", 0.0), size = pt_size, ) +
  geom_pointrange(data = dat_growth, aes(x = year, y= growth_d, ymin = lower, ymax = upper), col = mar_co, lwd = line_size) +
  labs(x = "Year", y = "Growth rate (mm/day)") +
  coord_cartesian(xlim =c(2016.5, 2019.5), ylim = c(0, 2.05)) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = ax_title), 
        axis.text = element_text(color = "black", size = ax_text))

fig_4

# save estuary growth plot
ggsave("3_plots/0_final/fig_4.pdf", dpi = 300, height = 3, width = 3, units = "in")

# FIGURE 5: SIZE AT MARINE ENTRY ####
# create annual rst data sets
dat_2017 <- dat_rst[dat_rst$tag_year=="2017",]
dat_2018 <- dat_rst[dat_rst$tag_year=="2018",]
dat_2019 <- dat_rst[dat_rst$tag_year=="2019",]

# determine 10%, 50%, and 90% size quantiles (in mm) at freshwater exit
quantile(dat_rst$tag_length, c(0.1, 0.5, 0.9))
q0.1 <- 78
q0.5 <- 93
q0.9 <- 112

# calculate annual predicted length at marine entry by size quantile 
end0.1 <- round(unique(dat_rst$pred_length[dat_rst$tag_length==q0.1]))
end0.5 <- round(unique(dat_rst$pred_length[dat_rst$tag_length==q0.5]))
end0.9 <- round(unique(dat_rst$pred_length[dat_rst$tag_length==q0.9]))

# calculate annual % growth by size quantile
growth0.1 <- round(unique(dat_rst$pred_growth[dat_rst$tag_length==q0.1])/q0.1*100, 0)
growth0.5 <- round(unique(dat_rst$pred_growth[dat_rst$tag_length==q0.5])/q0.5*100, 0)
growth0.9 <- round(unique(dat_rst$pred_growth[dat_rst$tag_length==q0.9])/q0.9*100, 0)

# create size quantile data frame for easy plotting
## create data frame
dat_size <- data.frame("year" = c(rep(2017, 3), rep(2018, 3), rep(2019, 3)), # year of tagging
                       "q" = rep(c(0.1, 0.5, 0.9), 3), # size quantile
                       "start" = NA, # start point for arrows (length at freshwater exit)
                       "end" = NA, # end point for arrows (length at marine entry)
                       "growth" = NA, # percent growth
                       "text" = NA, # text label for percent growth
                       "text_x" = NA, # plot location for arrow label
                       "text_y" = NA) # plot location for arrow label

### start size is observed size at freshwater exit
dat_size$start <- rep(c(q0.1, q0.5, q0.9), 3)

### end size is predicted size at marine entry
dat_size$end[dat_size$q == 0.1] <- end0.1
dat_size$end[dat_size$q == 0.5] <- end0.5
dat_size$end[dat_size$q == 0.9] <- end0.9

### growth is annual % growth
dat_size$growth[dat_size$q == 0.1] <- growth0.1
dat_size$growth[dat_size$q == 0.5] <- growth0.5
dat_size$growth[dat_size$q == 0.9] <- growth0.9

### add growth label information
dat_size$text <- paste0("+", dat_size$growth, "%")
dat_size$text_x <- dat_size$end-((dat_size$end-dat_size$start)/2)
dat_size$text_y <- rep(c(0.08, 0.09, 0.09), 3)

# draw base plot
### points
pt_size <- 1.5
pt_co <- "black"

### lines
line_size <- 1.25/2
arrow_size <- 1.25/2
fw_co <- "#2baa8e"
mar_co <- "#023a45"

## axes and labels
ax_title <- 10
ax_text <- 8
sub_title <- 2.75
height <- 0.09

# draw base plot
p <- ggplot() +
  geom_blank() +
  labs(x = "Fork Length (mm)", y = "Density") +
  coord_cartesian(xlim =c(60, 150), ylim = c(0,0.1)) +
  scale_y_continuous(breaks = seq(0, 0.1, 0.05)) +
  scale_x_continuous(breaks = seq(60, 150, 20)) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = ax_title), 
        axis.text = element_text(color = "black", size = ax_text))

# draw figure 5: size at marine entry
### 2017 plot
p_2017 <- p + 
  geom_density(data = dat_2017, aes(x=tag_length), col = fw_co, lwd = line_size, fill = adjustcolor(fw_co, 0.5)) +
  geom_density(data = dat_2017, aes(x=pred_length), col = mar_co, lwd = line_size, fill = adjustcolor(mar_co, 0.5)) +
  geom_point(data = dat_size[dat_size$year == 2017,], aes(x = start, y = text_y), size = pt_size, pch=c(21, 23, 22)) +
  geom_point(data = dat_size[dat_size$year == 2017,], aes(x = end, y = text_y), size = pt_size, pch=c(21, 23, 22), fill= pt_co) +
  geom_segment(data = dat_size[dat_size$year == 2017,], aes(x = start + 1, xend = end - 1, y = text_y, yend = text_y), size = arrow_size, arrow = arrow(length = unit(0.1, "cm"), type = "closed")) +
  annotate("text", x = dat_size$text_x[dat_size$year==2017], y = dat_size$text_y[dat_size$year==2017] + 0.009, label = dat_size$text[dat_size$year==2017], size = sub_title) +
  annotate("text", x = 60, y = height + 0.01, label = "a) 2017", size = sub_title) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

### 2018 plot
p_2018 <- p + 
  geom_density(data = dat_2018, aes(x=tag_length), col = fw_co, lwd = line_size, fill = adjustcolor(fw_co, 0.5)) +
  geom_density(data = dat_2018, aes(x=pred_length), col = mar_co, lwd = line_size, fill = adjustcolor(mar_co, 0.5)) +
  geom_point(data = dat_size[dat_size$year == 2018,], aes(x = start, y = text_y), size = pt_size, pch=c(21, 23, 22)) +
  geom_point(data = dat_size[dat_size$year == 2018,], aes(x = end, y = text_y), size = pt_size, pch=c(21, 23, 22), fill = pt_co) +
  geom_segment(data = dat_size[dat_size$year == 2018,], aes(x = start + 1, xend = end - 1, y = text_y, yend = text_y), size = arrow_size, arrow = arrow(length = unit(0.1, "cm"), type = "closed")) +
  annotate("text", x = dat_size$text_x[dat_size$year==2018], y = dat_size$text_y[dat_size$year==2018] + 0.009, label = dat_size$text[dat_size$year==2018], size = sub_title) +
  annotate("text", x = 60, y = height + 0.01, label = "b) 2018", size = sub_title) +
  theme(axis.title.x = element_blank())

### 2019 pot
p_2019 <- p + 
  geom_density(data = dat_2019, aes(x=tag_length), col = fw_co, lwd = line_size, fill = adjustcolor(fw_co, 0.5)) +
  geom_density(data = dat_2019, aes(x=pred_length), col = mar_co, lwd = line_size, fill = adjustcolor(mar_co, 0.5)) +
  geom_point(data = dat_size[dat_size$year == 2019,], aes(x = start, y = text_y), size = pt_size, pch=c(21, 23, 22)) +
  geom_point(data = dat_size[dat_size$year == 2019,], aes(x = end, y = text_y), size = pt_size, pch=c(21, 23, 22), fill=c("black")) +
  geom_segment(data = dat_size[dat_size$year == 2019,], aes(x = start + 1, xend = end - 1, y = text_y, yend = text_y), size = arrow_size, arrow = arrow(length = unit(0.1, "cm"), type = "closed")) +
  annotate("text", x = dat_size$text_x[dat_size$year==2019], y = dat_size$text_y[dat_size$year==2019] + 0.009, label = dat_size$text[dat_size$year==2019], size = sub_title) +
  annotate("text", x = 60, y = height + 0.01, label = "c) 2019", size = sub_title) +
  theme(axis.title.y = element_blank())

### plot for all years
fig_5 <- p_2017 / p_2018 / p_2019
fig_5

# save size at marine entry plot
ggsave("3_plots/0_final/fig_5.pdf", dpi = 300, height = 4.5, width = 6, units = "in")

