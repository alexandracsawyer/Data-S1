# LOAD REQUIRED PACKAGES ####
library(lubridate)
library(marked)
library(ggplot2)

# LOAD DATA ####
# ensure PIT tag numbers do not present in scientific notation
options(scipen=999)

# load comprehensive mark-recapture data set (only tagged fish that were recaptured in the estuary)
### individual fish are identified by PIT tag number
### columns preceded by "tag_" indicate initial mark data
### columns preceded by "rc_" indicate recapture data
dat_rc <- read.csv("dat_rc.csv")
rownames(dat_rc) <- c()
dat_rc$tag_year <- as.factor(dat_rc$tag_year)

# load comprehensive initial mark data set (all fish that were tagged at the RST)
### individual fish are identified by PIT tag number
dat_rst <- read.csv("dat_rst.csv")
rownames(dat_rst) <- c()
dat_rst$tag_year <- as.factor(dat_rst$tag_year)

# load day of year (DOY) data - this has adjusted DOY for different pooled data scenarios
doy <- read.csv("date_adjust.csv")

# load estuary sampling data - this has dates of estuary sampling for each year
sets <- read.csv("estuary_sets.csv")
sets$doy_adj3 <- doy$doy_adj3[match(sets$doy, doy$doy)] # adjust dates to reflect pooled sampling occasions
sets_2017 <- subset(sets, year==2017)
sets_2018 <- subset(sets, year==2018)
sets_2019 <- subset(sets, year==2019)

# PREPARE DATA FOR INDIVIDUAL ENCOUNTER HISTORY MATRIX ####
# mark-recapture data set 
dat_rc <- subset(dat_rc, tag_loc == "smolt trap") # keep only fish tagged at the RST
dat_rc <- dat_rc[,c("PIT_number", "tag_year", "tag_doy", "tag_length", "rc_year", "rc_doy", "rc_length", "rc_loc_1", "rc_loc_2")] # keep only columns relevant for analyses

# RST data set
dat_rst <- dat_rst[,c("PIT_number", "tag_year", "doy", "location", "length")] # keep only columns relevant for analyses
colnames(dat_rst) <- c("PIT_number", "tag_year", "tag_doy", "tag_loc", "tag_length") # adjust row and column names for consistency with mark-recapture data set

# merge recapture data into RST data set
dat_rst$rc_year <- dat_rc$rc_year[match(dat_rst$PIT_number, dat_rc$PIT_number)] # append recapture year
dat_rst$rc_doy <- dat_rc$rc_doy[match(dat_rst$PIT_number, dat_rc$PIT_number)] # append recapture doy
dat_rst$rc_length <- dat_rc$rc_length[match(dat_rst$PIT_number, dat_rc$PIT_number)] # append recapture length
dat_rst$rc_loc_1 <- dat_rc$rc_loc_1[match(dat_rst$PIT_number, dat_rc$PIT_number)] # append recapture location 1 (reach)
dat_rst$rc_loc_2 <- dat_rc$rc_loc_2[match(dat_rst$PIT_number, dat_rc$PIT_number)] # append recapture location 2 (site)

# adjust tag and recapture dates to match pooled dates (pooled for every 3 days)
### pooling sampling occasions enables model convergence
### otherwise too many sampling dates with too few recaptures
dat_rst$tag_doy <- doy$doy_adj3[match(dat_rst$tag_doy, doy$doy)]
dat_rst$rc_doy <- doy$doy_adj3[match(dat_rst$rc_doy, doy$doy)]

# CREATE INDIVIDUAL ENCOUNTER HISTORY MATRIX ####
# create empty data frame for individual encounter history matrix (IEH)
IEH <- data.frame("PIT_number" = dat_rst$PIT_number,
                  "tag_year" = dat_rst$tag_year,
                  "tag_doy" = dat_rst$tag_doy,
                  "tag_length" = dat_rst$tag_length,
                  "tag_loc" = dat_rst$tag_loc,
                  "rc_year" = dat_rst$rc_year,
                  "rc_doy" = dat_rst$rc_doy,
                  "rc_length" = dat_rst$rc_length,
                  "rc_loc" = dat_rst$rc_loc_2)

# check for correct number of recaps, confirm no duplicate PIT numbers
length(dat_rst$PIT[!is.na(dat_rst$rc_loc_2)==T]) # 64 recaps
setdiff(dat_rst$PIT[!is.na(dat_rst$rc_loc_2)==T], dat_rc$PIT_number) # same PIT numbers as recapture data set
which(duplicated(IEH$PIT)==T) # no duplicate PIT numbers

# create sampling occasion columns for each unique date RST and/or estuary were sampled 
rst_dates <- sort(unique(dat_rst$tag_doy)) # dates where fish were tagged at the RST
est_dates <- sort(unique(doy$doy_adj3)) # dates where where estuary seining occurred (for recaptures)
cols <- sort(unique(c(rst_dates, est_dates))) # create columns with unique sampling dates
IEH <- cbind(IEH, setNames(lapply(cols, function(x) x=NA), cols)) # add sampling occasion columns to IEH 

# run for-loop to establish encounter histories
for(i in 1:nrow(IEH)) {
  
  # determine tag and recap dates  of ith fish
  # if the PIT tagged fish matches an initial tag...
    if(IEH$PIT_number[i] == dat_rst$PIT_number[i]) {
      
      tag_date <- as.character(dat_rst$tag_doy[i]) # id the tag date
      rc_date <- ifelse(!is.na(dat_rst$rc_doy[i])==T, as.character(dat_rst$rc_doy[i]), NA) # id the recap date, NA if no recap
      
      # for each sampling occasion
      for(j in 10:ncol(IEH)) {
        # assign a 1 when tag_date = column name (this indicates initial tag)
        if(tag_date == colnames(IEH)[j]) {
          IEH[i,j] <- 1
        }
        # assign a 1 when rc_date is not NA and rc_date = column name (this indicates estuary recapture)
        else if(!is.na(rc_date) == T & rc_date == colnames(IEH)[j]) {
          IEH[i,j] <- 1
        }
        # otherwise assign a 0
        else IEH[i,j] <- 0
    }
    }
}

# manually adjust IEH for single fish that was recaptured twice in the estuary
IEH[IEH$PIT_number==226001029713,"157"] <- 1 # this fish was captured for a second time on doy 156 (adjusted doy = 157)
IEH[IEH$PIT_number==226001029713,] # confirm that IEH has been updated

# format IEH for CJS model run
IEH$ch <- do.call(paste0, IEH[,10:35]) # paste all encounter 1/0s to form ch column for mark-recapture model run
IEH <- IEH[,c(1, 36, 2:35)] # save new IEH data frame with full encounter histories
IEH$tag_year <- as.factor(IEH$tag_year)

# SAVE INDIVIDUAL ENCOUNTER HISTORY MATRIX ####
saveRDS(IEH, file="mod_res_cjs_ieh.rds")

# PREPARE DATA FOR CJS MODEL RUN ####
# model inputs
days <- as.numeric(colnames(IEH[,11:36])) # sampling occasions
begin <- days[1] # date of first sampling occasion
int <- as.vector(diff(days)) # interval between sampling occasions

# process data for 'marked' CJS model run
pr <- process.data(IEH, model="CJS", groups=c("tag_year"), begin.time=begin, time.intervals = int, accumulate = F)

# create model design matrices for Phi (apparent survival probability) and p (detection probability)
design.Phi <- list(static=c("tag_length", "tag_year")) # design matrix for Phi, including individual parameters for tag length and tag year
design.p <- list(static=c("tag_length", "tag_year")) # design matrix for p, including individual parameters for tag length and tag year
design.parameters <- list(Phi=design.Phi, p=design.p)
ddl <- make.design.data(pr, parameters=design.parameters)

# fix p (probability of recapture) at 0 for dates when no estuary sampling occurred
ddl_p_fix <- ddl$p[,c("tag_year", "time", "fix")] # create subset matrix with only columns of interest

### fix p (probability of recapture) at 0 for dates when only RST sampling occurred
### if p for a given sampling occasion is NA (estimable) and the estuary was not sampled on that day, fix p at 0
ddl_p_fix$fix <- ifelse(is.na(ddl_p_fix$fix == T) & !(ddl_p_fix$time %in% sets$doy_adj3), 0, ddl_p_fix$fix)

### fix p (probability of recapture) at 0 for annually-specific dates when estuary sampling did not occur
### if p for a given sampling occasion is NA (estimable) and the estuary was not sampled on that day in that year, fix p at 0
ddl_p_fix$fix <- ifelse(is.na(ddl_p_fix$fix == T) & ddl_p_fix$tag_year==2017 & !(ddl_p_fix$time %in% sets_2017$doy_adj3), 0, 
                        ifelse(is.na(ddl_p_fix$fix == T) & ddl_p_fix$tag_year==2018 & !(ddl_p_fix$time %in% sets_2018$doy_adj3), 0, 
                               ifelse(is.na(ddl_p_fix$fix == T) & ddl_p_fix$tag_year==2019 & !(ddl_p_fix$time %in% sets_2019$doy_adj3), 0, ddl_p_fix$fix)))

### replace fixed values in design matrix for p (probability of detection) with updated values that reflect annually-specific estuary sampling occasions
ddl$p$fix <- ddl_p_fix$fix
head(ddl$p)

# RUN CJS MARK-RECAPTURE MODEL ####
# establish model inputs
Phi.length <- list(formula=~tag_length) # apparent survival probability ~ length at tagging
p.dot <- list(formula=~1) # constant probability of detection

# run cjs mark-recapture model
mod_length_dot <- crm(pr, ddl, model.parameters=list(Phi=Phi.length,p=p.dot))
mod_length_dot <- cjs.hessian(mod_length_dot)
mod_length_dot

# PREDICT ESTUARY RESIDENCE VALUES ####
### phi derivation to predict estuary residence days: days = -1/log(survival probability)

# standardized length inputs for model predictions
length_vals <- seq(min(IEH$tag_length[is.na(IEH$rc_length)==F]), max(IEH$tag_length[is.na(IEH$rc_length)==F]), 1)

# predict length values in logit space, transform to real (probability) space
logit <- mod_length_dot$results$beta$Phi[1]+(mod_length_dot$results$beta$Phi[2]*length_vals)
odds <- exp(logit)
prob <- odds/(1+odds)
days_cjs <- -1/log(prob) # this is estuary residence in days

# create data frame with predicted estuary residence days for fish of a given length at tagging
### exponentiate model predictions to translate from log space into real space
days_cjs <- as.data.frame(days_cjs) # create data frame from model predictions
days_cjs$tag_length <- length_vals # add length values used to predict estuary residence
days_cjs

# SAVE BEST-FIT MODEL AND MODEL PREDICTIONS ####
saveRDS(mod_length_dot, file="mod_res_cjs.rds")
write.csv(days_cjs, "days_res_cjs.csv")

# FIGURE S1 - PANEL C: ESTUARY RESIDENCE (CJS MARK-RECAPTURE) ####
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
fig_s1_c <- ggplot() + 
  geom_point(data = dat_rc, aes(x = tag_length, y = days), col = adjustcolor(pt_co, 0.6), size = pt_size) +
  labs(x = "Fork length at freshwater exit (mm)", y = "Days in estuary") +
  coord_cartesian(ylim = c(0, 50)) +
  scale_x_continuous(breaks=seq(70, 120, 10)) +
  geom_line(data = days_cjs, aes(x = tag_length, y = days_cjs), col = line_co, lwd = line_size) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = ax_title), 
        axis.text = element_text(color = "black", size = ax_text))

fig_s1_c

