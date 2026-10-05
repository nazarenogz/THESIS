# Behavioral state occupancy of ungulates in the Northern Apennines: 
# a multi-state camera trap approach to the foraging-vigilance trade-off
# Thesis of Nazareno Gimenez Zapiola

# script: Charlotte Vanderlocht

# IMPORTANT: 
# 0 = camera operational + species not detected 
# 1 = camera operational + Moving 
# 2 = camera operational + Vigilant 
# 3 = camera operational + Foraging 
# 4 = camera operational + Other 
# NA = camera unavailable / videos not reviewed 
# 
# Camera effort is retained as an observation-level covariate. 
# This means that "no detection" is NOT assigned when the camera 
# was unavailable.

rm(list = ls())
gc()


# Packages ----------------------------------------------------------------
library(tidyr)
library(lubridate)
library(tidyverse)
library(unmarked)
library(corrplot)
library(RColorBrewer)
library(DHARMa)
library(sf)
library(sjPlot)
library(htmltools)
library(sjmisc)
library(readxl)
library(mclogit)
library(forestplot)

# Functions ---------------------------------------------------------------

# Function to combine sets of X Julian days
combine_julian_days <- function(data, 
                                step) {
  new_ncol <- ceiling(ncol(data) / step) # Number of new columns will be the original number of columns divided by the step
  new_data <- matrix(NA, 
                     nrow = nrow(data), 
                     ncol = new_ncol) # Initialize the new matrix to store the results
  
  for (i in 1:new_ncol) { # Define the column indices for the current X-day block
    start_col <- (i - 1) * step + 1
    end_col <- min(i * step, 
                   ncol(data))
    new_data[, i] <- apply(data[, start_col:end_col, drop = FALSE], 1, function(x) { # Apply the logic for combining the Julian days
      if (any(is.na(x))) {
        return(NA)
      } else if (any(x == 1)) {
        return(1)
      } else {
        return(0)
      }
    })
  }
  
  return(new_data)
}

# Function to make species matrix at specific temporal resolution
make_species_matrix <- function(species_name, 
                                data, 
                                effort_data, 
                                Sites, 
                                time_values, 
                                effort_threshold) { 
  # Species observations 
  species_data <- data %>%
    filter(taxon == species_name) %>% 
    group_by(CameraTrapID, period) %>% 
    summarise(
      n_obs = n(), # all detections of the species, classified or not
      n_LOCO = sum(main_behav == "LOCO", na.rm = TRUE ),
      n_VIG = sum(main_behav == "VIG", na.rm = TRUE ),
      n_FOR = sum(main_behav == "FOR", na.rm = TRUE ),
      n_classified = n_LOCO + n_VIG + n_FOR,
      .groups = "drop") %>%
    mutate(
      main_behav = case_when(
        # Detected, but no behaviour in our 3 categories = NA (NOT Absent)
        n_classified == 0 ~ NA_character_,
        n_LOCO >= n_VIG & n_LOCO >= n_FOR ~ "Moving",
        n_VIG > n_LOCO & n_VIG >= n_FOR ~ "Vigilant",
        n_FOR > n_VIG & n_FOR > n_LOCO ~ "Foraging",
        TRUE ~ NA_character_)
      ) %>% 
    dplyr::select(CameraTrapID, period, n_obs, main_behav)

  # Complete camera x period matrix
  species_data <- expand_grid(CameraTrapID = Sites, 
                              period = time_values) %>% 
    left_join(species_data, by = c("CameraTrapID", "period")) 
  
  # Add camera effort 
  species_data <- species_data %>% 
    left_join(effort_data, by = c("CameraTrapID", "period")) 
  
  # effort not ok                                   = NA
  # effort ok + no detections (n_obs == 0)           = Absent
  # effort ok + detections, behaviour not classified = NA
  # effort ok + detections, behaviour classified     = Moving / Vigilant / Foraging
  species_data <- species_data %>%
    mutate(main_behav = case_when(
      # 1) Camera insufficiently sampled
      is.na(effort) ~ NA_character_,
      effort < effort_threshold ~ NA_character_,
      # 2) Camera sufficiently sampled and species not detected
      is.na(n_obs) ~ "Absent",
      # 3) Species detected: classified state, or NA if unclassified
      TRUE ~ main_behav
    )
    )
  
  # Wide behavioural matrix 
  species_wide <- species_data %>% 
    dplyr::select(CameraTrapID, period, main_behav) %>% 
    pivot_wider(names_from = period, 
                values_from = main_behav ) %>% 
    as.data.frame() 
  rownames(species_wide) <- species_wide$CameraTrapID 
  species_wide$CameraTrapID <- NULL 
  
  # Wide effort matrix 
  effort_wide <- species_data %>% 
    dplyr::select(CameraTrapID, period, effort) %>% 
    pivot_wider(names_from = period, values_from = effort) %>% 
    as.data.frame() 
  rownames(effort_wide) <- effort_wide$CameraTrapID 
  effort_wide$CameraTrapID <- NULL 
  
  # Convert to matrices 
  species_matrix <- as.matrix(species_wide) 
  effort_matrix <- as.matrix(effort_wide) 
  
  # Return both 
  return(list(observations = species_matrix, 
                effort = effort_matrix, 
                long = species_data) 
  ) 
  }

# Function to convert behavioural states to numeric
# 0 = Absent 
# 1 = Moving
# 2 = Vigilant 
# 3 = Foraging 
# NA = camera or behaviour unavailable 
make_numeric_states <- function(x) { 
  x_num <- matrix( 
    NA_integer_, 
    nrow = nrow(x), 
    ncol = ncol(x), 
    dimnames = dimnames(x) 
  ) 
  x_num[x == "Absent"] <- 0 
  x_num[x == "Moving"] <- 1 
  x_num[x == "Vigilant"] <- 2 
  x_num[x == "Foraging"] <- 3
  return(x_num) 
  }


# Parameters --------------------------------------------------------------

## Temporal resolution -----------------------------------------------------
study_start <- as.Date("2024-03-01") 
study_end <- as.Date("2024-06-21")

time_unit_days <- 3
period_range <- seq_len(ceiling(as.integer(study_end - study_start + 1) / time_unit_days))

effort_threshold <- 0.3 # camera-trap operational at least 30% of unit
# important: effort is also retained as a continuous observation-level covariate!

## Behavioural states -----------------------------------------------------

# 0 = Absent / non-detection 
# 1 = Moving # REMOVE?
# 2 = Vigilant 
# 3 = Foraging 
# 4 = Other
# NA = No data 

# Import data ---------------------------------------------------------------
setwd("D:/TIROCINIO/TESI/DATA ANALYSIS")
RawData <- read_excel("Primavera_2024_(01-03-24_21-03_24).xlsx")
str(RawData)

covs <- read_excel("COVS_SITES_TESI.xlsx") %>%
  rename(CameraTrapID = FT) %>%
  mutate(dens_road = as.numeric(dens_road))

# Clean data --------------------------------------------------------------
data_clean <- RawData %>%
  mutate(
    date_character_in = NA,
    date_character_in = case_when(
      substring(DateMDY, 3, 3) == "/" ~ as.character(mdy(DateMDY)),
      str_length(DateMDY) == 5 ~ as.character(as.Date(as.numeric(DateMDY), origin = "1899-12-30"))
    )) %>%
  mutate(
    date = NA,
    date = case_when(
      substring(DateMDY, 3, 3) == "/" ~ ymd(date_character_in),
      str_length(DateMDY) == 5 ~ ydm(date_character_in)
    )) %>%
  rowwise() %>%
  mutate(
    date_character = as.character(date),
    timestamp_character = paste0(date_character, " ", substring(as.character(Time), 12, 19)),
    timestamp = ymd_hms(timestamp_character),
    main_behav = Bdisplay1
  ) %>%
  ungroup()

is.timepoint(data_clean$timestamp)
is.timepoint(data_clean$date)
range(data_clean$date)
range(month(data_clean$date))

# These have a date but no timestamp
which(is.na(data_clean$timestamp))
data_clean[c(1039,1197,2606,4184,4816,5502,7007,7210,7214),]

range(data_clean$date)

# Species
table(data_clean$Species)

data_clean_tax <- data_clean %>%
  mutate(
    taxon = NA,
    taxon = case_when(
      Species == "Fallow_deer" ~ "Dama dama",
      Species == "Red_deer" ~ "Cervus elaphus",
      Species == "Roe_deer" ~ "Capreolus capreolus",
      Species == "Wild_boar" | Species == "wild_boar" ~ "Sus scrofa",
      Species == "Hare" ~ "Lepus europaeus",
      Species == "Red_fox" ~ "Vulpes vulpes",
      Species == "Porcupine" | Species == "porcupine" ~ "Hystrix cristata",
      Species == "Wolf" | Species == "wolf" ~ "Canis lupus",
      Species == "Badger" | Species == "badger" ~ "Meles meles",
      str_contains(Species, c("Human", "HUman", "human"), logic = "or") ~ "Homo sapiens"#,
      # str_contains(Species, c("car", "Car"), logic = "or") ~ "Vehicle",
      # str_contains(Species, c("Martes", "marten"), logic = "or") ~ "Martes sp."
      ))
table(data_clean_tax$taxon)
table(data_clean_tax$taxon, data_clean_tax$main_behav)

# Count number of videos of wolves and humans per cameratrap
str(data_clean_tax)
data_clean_tax_reddeer_pred_humans <- data_clean_tax %>%
  group_by(CameraTrapID) %>%
  mutate(
    n_video_Homo_sapiens = sum(`N. video`[taxon == "Homo sapiens"], na.rm = TRUE),
    n_video_Canis_lupus = sum(`N. video`[taxon == "Canis lupus"], na.rm = TRUE)
  ) %>%
  ungroup() %>%
  dplyr::select(CameraTrapID, n_video_Homo_sapiens, n_video_Canis_lupus) %>%
  distinct()
table(data_clean_tax_reddeer_pred_humans$CameraTrapID, data_clean_tax_reddeer_pred_humans$n_video_Canis_lupus)
table(data_clean_tax_reddeer_pred_humans$CameraTrapID, data_clean_tax_reddeer_pred_humans$n_video_Homo_sapiens)

# Focus on ungulates with behaviours
ungul <- data_clean_tax %>%
  filter(
    taxon %in% c("Cervus elaphus", "Capreolus capreolus", "Dama dama", "Sus scrofa")
  ) %>%
  mutate(
    jday = yday(date),
    period = floor(as.integer(date - study_start) / time_unit_days) + 1
  )
range(ungul$jday)
range(ungul$period)

# Camera-trap outages -----------------------------------------------------

# Based on message from Umberto
camera_outages <- tribble(
  ~CameraTrapID, ~start, ~end, ~reason,
  "MM07", "2024-03-15", "2024-04-04", "video_not_reviewed",
  "MM08", "2024-04-05", "2024-05-03", "malfunction",
  "MM08", "2024-05-20", "2024-06-18", "malfunction",
  "MM09", "2024-05-11", "2024-05-30", "dead_battery",
  "MM10", "2024-03-06", "2024-04-05", "malfunction",
  "MM10", "2024-05-03", "2024-05-30", "malfunction",
  "MM11", "2024-05-16", "2024-06-10", "dead_battery",
  "MM12", "2024-05-18", "2024-06-10", "dead_battery",
  "MM13", "2024-05-26", "2024-06-10", "dead_battery",
  "MM15", "2024-03-06", "2024-04-26", "malfunction",
  "MM16", "2024-05-12", "2024-05-22", "malfunction",
  "MM18", "2024-03-08", "2024-04-26", "video_not_reviewed",
  "MM21", "2024-03-06", "2024-04-05", "malfunction",
  "MM22", "2024-03-06", "2024-04-05", "malfunction",
  "MM22", "2024-04-21", "2024-05-03", "dead_battery",
  "MM22", "2024-05-20", "2024-05-30", "dead_battery",
  "MM25", "2024-04-24", "2024-05-22", "video_not_reviewed",
  "MM27", "2024-03-15", "2024-04-11", "malfunction",
  "MM27", "2024-04-30", "2024-06-14", "malfunction",
  "MM28", "2024-04-11", "2024-05-27", "stolen",
  "MM29", "2024-04-11", "2024-04-26", "video_not_reviewed",
  "MM30", "2024-04-11", "2024-04-26", "video_not_reviewed",
  "MM31", "2024-03-15", "2024-04-04", "video_not_reviewed",
  "MM32", "2024-03-15", "2024-04-04", "video_not_reviewed",
  "MM33", "2024-03-15", "2024-04-04", "video_not_reviewed",
  "MM33", "2024-05-30", "2024-06-07", "dead_battery",
  "MM34", "2024-03-23", "2024-04-24", "malfunction",
  "MM35", "2024-03-07", "2024-04-24", "malfunction",
  "MM37", "2024-03-07", "2024-04-26", "malfunction",
  "MM38", "2024-04-12", "2024-05-27", "malfunction"
) %>%
  mutate(
    start = as.Date(start),
    end   = as.Date(end)
  )
camera_outages

# Daily effort -------------------------------------------------------------

Sites <- unique(data_clean_tax$CameraTrapID)  # Consider each site-per-year as a "site" in occupancy
Sites <- Sites[!is.na(Sites)]

# Template: all TRUE
all_camera_days <- 
  expand_grid( 
    CameraTrapID = Sites, 
    date = seq(study_start, study_end, by = "day" )
  ) %>% 
  mutate(operational = TRUE, 
         outage_reason = NA_character_ )

# Fill with outages: FALSE
for(i in seq_len(nrow(camera_outages))) { 
  all_camera_days$operational[
    all_camera_days$CameraTrapID == camera_outages$CameraTrapID[i] & 
      all_camera_days$date >= camera_outages$start[i] & 
      all_camera_days$date <= camera_outages$end[i] 
  ] <- FALSE 
  all_camera_days$outage_reason[ 
    all_camera_days$CameraTrapID == camera_outages$CameraTrapID[i] & 
      all_camera_days$date >= camera_outages$start[i] & 
      all_camera_days$date <= camera_outages$end[i] 
  ] <- camera_outages$reason[i] 
  }

# Check
table(all_camera_days$operational, useNA = "ifany") 
table(all_camera_days$outage_reason, useNA = "ifany")
# View(all_camera_days)

all_camera_days <- all_camera_days %>% 
  mutate( 
    jday = yday(date), 
    period = floor(as.integer(date - study_start) / time_unit_days) + 1,, 
    week_start = floor_date(date, unit = "week", week_start = 1), 
    month = month(date) 
    )


# Period effort -----------------------------------------------------------

period_effort <- all_camera_days %>% 
  filter(period %in% period_range) %>% 
  group_by(CameraTrapID, period) %>% 
  summarise(n_days = n(), 
            n_operational = sum(operational), 
            n_unavailable = sum(!operational), 
            effort = n_operational / n_days, 
            .groups = "drop" ) 
# Check
summary(period_effort$effort) 
table(period_effort$effort == 0) 
table(period_effort$effort == 1)
# View(period_effort)


# Make species matrices ---------------------------------------------------

## Red deer --------------------------------------------------------------- 
RedDeer_data <- make_species_matrix(species_name = "Cervus elaphus", 
                                    data = ungul, 
                                    effort_data = period_effort, 
                                    Sites = Sites, 
                                    time_values = period_range, 
                                    effort_threshold = effort_threshold) 
RedDeer <- RedDeer_data$observations 
RedDeer_effort <- RedDeer_data$effort 

## Roe deer ---------------------------------------------------------------- 
RoeDeer_data <- make_species_matrix(species_name = "Capreolus capreolus", 
                                    data = ungul, 
                                    effort_data = period_effort, 
                                    Sites = Sites, 
                                    time_values = period_range, 
                                    effort_threshold = effort_threshold) 
RoeDeer <- RoeDeer_data$observations 
RoeDeer_effort <- RoeDeer_data$effort 

## Fallow deer ------------------------------------------------------------- 
FallowDeer_data <- make_species_matrix(species_name = "Dama dama", 
                                       data = ungul, 
                                       effort_data = period_effort, 
                                       Sites = Sites, 
                                       time_values = period_range, 
                                       effort_threshold = effort_threshold) 
FallowDeer <- FallowDeer_data$observations 
FallowDeer_effort <- FallowDeer_data$effort 

## Wild boar --------------------------------------------------------------- 
WildBoar_data <- make_species_matrix(species_name = "Sus scrofa", 
                                     data = ungul, 
                                     effort_data = period_effort, 
                                     Sites = Sites, 
                                     time_values = period_range, 
                                     effort_threshold = effort_threshold) 
WildBoar <- WildBoar_data$observations 
WildBoar_effort <- WildBoar_data$effort

# Various checks ------------------------------------------------------------------

## Check dimensions --------------------------------------------------------
rownames(RedDeer) == rownames(RoeDeer) 
rownames(RedDeer) == rownames(FallowDeer) 
rownames(RedDeer) == rownames(WildBoar)
dim(RedDeer) 
dim(RedDeer_effort)
dim(RoeDeer) 
dim(RoeDeer_effort)
dim(FallowDeer) 
dim(FallowDeer_effort)
dim(WildBoar) 
dim(WildBoar_effort)

## Check observations ------------------------------------------------------
table(RedDeer, useNA = "ifany") 
table(RoeDeer, useNA = "ifany") 
table(FallowDeer, useNA = "ifany") 
table(WildBoar, useNA = "ifany")


## Check effort ------------------------------------------------------------
summary(as.vector(RedDeer_effort)) 
table(RedDeer_effort == 0, useNA = "ifany" ) 
table(RedDeer_effort < effort_threshold, useNA = "ifany")

# Site covariates -------------------------------------------------------------------

# Add empty rows for camera traps MM09 MM11 MM16
to_add <- covs[(1:3),]
to_add[,] <- NA
to_add$CameraTrapID <- c("MM09", "MM11", "MM16")

all_covs <- rbind(covs[-c(28:30),], to_add)

all_sites <- all_covs %>%
  filter(CameraTrapID %in% rownames(RedDeer)) %>%
  arrange(match(CameraTrapID, rownames(RedDeer))) %>%
  left_join(data_clean_tax_reddeer_pred_humans, 
            by = "CameraTrapID")

# Check site order
table(rownames(RedDeer) == all_sites$CameraTrapID) # order is perfect
table(rownames(RoeDeer) == all_sites$CameraTrapID) # order is perfect
table(rownames(FallowDeer) == all_sites$CameraTrapID) # order is perfect
table(rownames(WildBoar) == all_sites$CameraTrapID) # order is perfect


# Remove site covariate columns that are missing --------------------------
all_sites_model <- all_sites %>% 
  dplyr::select(-CameraTrapID) 
all_sites_model <- all_sites_model[
  ,
  sapply(all_sites_model, function(x) 
    !all(is.na(x)) ), 
  drop = FALSE 
  ]

# Check site covariates
summary(all_sites_model)


# Make numeric matrices ---------------------------------------------------
RedDeer_num <- make_numeric_states(RedDeer) 
table(RedDeer_num, useNA = "ifany") # 0,1,2,3 NA

RoeDeer_num <- make_numeric_states(RoeDeer) 
table(RoeDeer_num, useNA = "ifany") # 0,1,2,3, NA

FallowDeer_num <- make_numeric_states(FallowDeer) 
table(FallowDeer_num, useNA = "ifany") # 0, NA

WildBoar_num <- make_numeric_states(WildBoar)
table(WildBoar_num, useNA = "ifany") # 0,1,2,3 NA
# Too little "2" --> NA
# WildBoar_num[WildBoar_num == 2] <- NA_integer_
# table(WildBoar_num, useNA = "ifany") # 0,1,3 NA

# Observation covariates --------------------------------------------------
RedDeer_obsCovs <- list(effort = RedDeer_effort) 
RoeDeer_obsCovs <- list(effort = RoeDeer_effort) 
FallowDeer_obsCovs <- list(effort = FallowDeer_effort) 
WildBoar_obsCovs <- list(effort = WildBoar_effort)


# umf: unmarked frames ----------------------------------------------------
umf_reddeer <- unmarkedFrameOccuMS(y = RedDeer_num, 
                                   siteCovs = all_sites_model, 
                                   obsCovs = RedDeer_obsCovs) 

umf_roedeer <- unmarkedFrameOccuMS(y = RoeDeer_num, 
                                   siteCovs = all_sites_model, 
                                   obsCovs = RoeDeer_obsCovs) 

# umf_fallowdeer <- unmarkedFrameOccuMS(y = FallowDeer_num, 
#                                       siteCovs = all_sites_model, 
#                                       obsCovs = FallowDeer_obsCovs)
# Less than 3 states.. doesn't work.

umf_wildboar <- unmarkedFrameOccuMS(y = WildBoar_num, 
                                    siteCovs = all_sites_model, 
                                    obsCovs = WildBoar_obsCovs)

# Check
summary(umf_reddeer) 
summary(umf_roedeer) 
# summary(umf_fallowdeer) # nope
summary(umf_wildboar)


# Model formulas ----------------------------------------------------------
# TBD with Nazareno

detection_formula <- "~ scale(effort)" # + scale(veg_h)" # + scale(dens_trees)" #  
occupancy_formula <- "~ scale(wood_biomass) + scale(n_video_Homo_sapiens) + scale(n_video_Canis_lupus)" # + scale(dens_road)" #+ scale(slope)"

# For red deer, roe deer, wild boar at periodly scale: 
# we have
# 0 = Absent 
# 1 = Moving 
# 2 = Vigilant 
# 3 = Foraging 
# NA

S <- 4
n_psi <- S-1  # occupancy parameters, for each comparison with baseline
n_det <- S*(S-1)/2 # detection parameters
n_psi 
n_det

# Simplify detection formulas. 
# Let's assume that WHEN an animal is detected, the behavioural classification is correct.
detformulas_simple <- c(
  detection_formula,  # p11 ~ effort
  "~1",               # p12
  "~1",               # p13
  detection_formula,  # p22
  "~1",               # p23
  detection_formula   # p33
)

detformulas_supersimple <- c(
  detection_formula,  # p11 ~ effort
  "~1",               # p12
  "~1",               # p13
  "~1",               # p22
  "~1",               # p23
  "~1"               # p33
)


# RED DEER ----------------------------------------------------------------

## Null model
model_null_reddeer <- occuMS( 
  detformulas = rep("~1", n_det), 
  psiformulas = rep("~1", n_psi), 
  data = umf_reddeer 
  ) 
summary(model_null_reddeer)

## Full model
model_all_reddeer <- occuMS( 
  detformulas = rep(detection_formula, n_det), #detformulas_supersimple, # rep("~1", n_det), #  
  psiformulas = rep(occupancy_formula, n_psi), 
  data = umf_reddeer, 
  parameterization = "multinomial",
  control = list(maxit = 5000)
) 
summary(model_all_reddeer)

## Plot

plot(model_all_reddeer) # Predicted values vs Residuals

# Prepare dataset for forest plot
mod_summary <- summary(model_all_reddeer)
df_state <- mod_summary$state
df_state$submodel <- "State (Occupancy)"
df_state$variable <- rownames(df_state)
df_det <- mod_summary$det
df_det$submodel <- "Detection"
df_det$variable <- rownames(df_det)
plot_data <- rbind(df_state, df_det)
colnames(plot_data) <- c("Estimate", "SE", "z", "p_value", "Submodel", "Variable")
plot_data$Lower <- plot_data$Estimate - (1.96 * plot_data$SE)
plot_data$Upper <- plot_data$Estimate + (1.96 * plot_data$SE)
table_text <- cbind(
  c("Variable", plot_data$Variable),
  c("Estimate", round(plot_data$Estimate, 2))
)
# You can further change the names here, to make it more readable! 
# "psi[1]" into "Occupancy - Moving", for example
# removing "scale()" for example
# etc.

# Forest plot option 1 'classical'
plot_data <- plot_data %>%
  mutate(Variable = factor(Variable, levels = rev(unique(Variable))),    
         Submodel = factor(Submodel, levels = c("State (Occupancy)", "Detection")))

ggplot(plot_data,
       aes(x = Estimate, y = Variable)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") + # Reference line = no effect
  geom_errorbarh(aes(xmin = Lower, xmax = Upper),
                 height = 0.15,
                 linewidth = 0.7) + # 95% CI
  geom_point(size = 3) + # Estimate
  facet_grid(Submodel ~ ., scales = "free_y", space = "free_y") +  # Separate State and Detection
  labs(x = "Coefficient", y = NULL, title = "Red deer") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 13),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        axis.title.x = element_text(size = 13),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()
  )

# ROE DEER ----------------------------------------------------------------

## Null model
model_null_roedeer <- occuMS(
  detformulas = rep("~1", n_det), 
  psiformulas = rep("~1", n_psi), 
  data = umf_roedeer ) 
summary(model_null_roedeer)

## Full model
model_all_roedeer <- occuMS( 
  detformulas = detformulas_supersimple, # detformulas_simple, #rep(detection_formula, n_det), 
  psiformulas = rep(occupancy_formula, n_psi), 
  data = umf_roedeer, 
  parameterization = "multinomial",
  control = list(maxit = 5000)
) 
summary(model_all_roedeer)

# ROE DEER FOREST PLOT -----------------------------------------------------

mod_summary_roedeer <- summary(model_all_roedeer)
df_state_roedeer <- mod_summary_roedeer$state
df_state_roedeer$submodel <- "State (Occupancy)"
df_state_roedeer$variable <- rownames(df_state_roedeer)
df_det_roedeer <- mod_summary_roedeer$det
df_det_roedeer$submodel <- "Detection"
df_det_roedeer$variable <- rownames(df_det_roedeer)
plot_data_roedeer <- rbind(df_state_roedeer, df_det_roedeer)
colnames(plot_data_roedeer) <- c("Estimate", "SE", "z", "p_value", "Submodel", "Variable")
plot_data_roedeer$Lower <- plot_data_roedeer$Estimate - (1.96 * plot_data_roedeer$SE)
plot_data_roedeer$Upper <- plot_data_roedeer$Estimate + (1.96 * plot_data_roedeer$SE)

plot_data_roedeer <- plot_data_roedeer %>%
  mutate(Variable = factor(Variable, levels = rev(unique(Variable))),    
         Submodel = factor(Submodel, levels = c("State (Occupancy)", "Detection")))

ggplot(plot_data_roedeer,
       aes(x = Estimate, y = Variable)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper),
                 height = 0.15,
                 linewidth = 0.7) +
  geom_point(size = 3) +
  facet_grid(Submodel ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Coefficient", y = NULL, title = "Roe deer") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 13),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        axis.title.x = element_text(size = 13),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()
  )

# WILD BOAR ----------------------------------------------------------------

## Null model
model_null_wildboar <- occuMS(
  detformulas = rep("~1", n_det), 
  psiformulas = rep("~1", n_psi), 
  data = umf_wildboar ) 
summary(model_null_wildboar)

## Full model
model_all_wildboar <- occuMS( 
  detformulas = detformulas_supersimple, #rep(detection_formula, n_det), 
  psiformulas = rep(occupancy_formula, n_psi), 
  data = umf_wildboar, 
  parameterization = "multinomial",
  control = list(maxit = 5000)
) 
summary(model_all_wildboar)

# WILD BOAR FOREST PLOT -----------------------------------------------------

mod_summary_wildboar <- summary(model_all_wildboar)
df_state_wildboar <- mod_summary_wildboar$state
df_state_wildboar$submodel <- "State (Occupancy)"
df_state_wildboar$variable <- rownames(df_state_wildboar)
df_det_wildboar <- mod_summary_wildboar$det
df_det_wildboar$submodel <- "Detection"
df_det_wildboar$variable <- rownames(df_det_wildboar)
plot_data_wildboar <- rbind(df_state_wildboar, df_det_wildboar)
colnames(plot_data_wildboar) <- c("Estimate", "SE", "z", "p_value", "Submodel", "Variable")
plot_data_wildboar$Lower <- plot_data_wildboar$Estimate - (1.96 * plot_data_wildboar$SE)
plot_data_wildboar$Upper <- plot_data_wildboar$Estimate + (1.96 * plot_data_wildboar$SE)

plot_data_wildboar <- plot_data_wildboar %>%
  mutate(Variable = factor(Variable, levels = rev(unique(Variable))),    
         Submodel = factor(Submodel, levels = c("State (Occupancy)", "Detection")))

ggplot(plot_data_wildboar,
       aes(x = Estimate, y = Variable)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper),
                 height = 0.15,
                 linewidth = 0.7) +
  geom_point(size = 3) +
  facet_grid(Submodel ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Coefficient", y = NULL, title = "Wild boar") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 13),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        axis.title.x = element_text(size = 13),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()
  )

# SIMPLE MULTINOMIAL ------------------------------------------------------

ungul_glmm <- as.data.frame(ungul) 
ungul_glmm$main_behav[ungul_glmm$main_behav == "FOR"] <- "Foraging"
ungul_glmm$main_behav[ungul_glmm$main_behav == "LOCO"] <- "Moving"
ungul_glmm$main_behav[ungul_glmm$main_behav == "VIG"] <- "Vigilant"
table(ungul_glmm$main_behav)

ungul_glmm <- ungul_glmm %>%
  mutate(
    main_behav_f = relevel(as.factor(main_behav), ref = "Foraging"),
    # Video-level covariates
    log_group = log(N_tot),             # group size (log: very skewed, 1-15)
    young = as.integer(N_ageYoung > 0)  # at least one young in the video
  ) %>%
  left_join(all_sites, by = "CameraTrapID")

table(ungul_glmm$main_behav_f)


## Settings ----------------------------------------------------------------

# Species to model (label = taxon)
mn_species <- c("Red deer"  = "Cervus elaphus",
                "Roe deer"  = "Capreolus capreolus",
                "Wild boar" = "Sus scrofa")

# Covariates (name = axis label). All are site-level.
# grasslands removed: > 0 at only 3 of 27 sites (MM23, MM38, MM29),
# so its effect rested almost entirely on MM29.
# Other habitat covariates screened (each added one at a time to these 6):
# dens_trees, veg_h, veg_h_rug, sh_uso, dist_forest, mixed, dist_water.
# None significant after Holm correction (7 covariates x 3 species) --> not included.
mn_covs <- c(n_video_Homo_sapiens = "Number of local human sightings",
             n_video_Canis_lupus  = "Number of local wolf sightings",
             wood_biomass         = "Woody biomass (Mg/ha)",
             dens_road            = "Road density (km/km²)",
             dist_urb             = "Distance to urban areas (m)",
             slope                = "Slope (%)") # max 117.7 --> percent, not degrees

# Video-level covariates (also screened: climate bio03/06/08/09, dens_pop, nl,
# dist_road, diel period, julian day, camera temperature --> not included).
# Group size and presence of young: strongest effects of all candidates
# (wild boar group size, red deer young; both p < 0.002) and a priori
# expected to affect vigilance (group-size effect, maternal vigilance).
# log_group is continuous (scaled, plotted); young is binary (not scaled).
mn_video_covs <- c(log_group = "Group size (individuals)")
mn_controls <- c(young = "Young present")

mn_covs_all <- c(mn_covs, mn_video_covs)

mn_out_dir <- "D:/TIROCINIO/TESI/DATA ANALYSIS/multinomial_output"
dir.create(mn_out_dir, showWarnings = FALSE)

behav_colours <- c(Foraging = "#287C78", Moving = "#C78332", Vigilant = "#625A91")
behav_fills   <- c(Foraging = "#9FC5C2", Moving = "#E1C18A", Vigilant = "#AAA3CA")


## Collinearity (site level) -----------------------------------------------
# Covariates are measured per camera trap, so correlations are computed
# across sites (not across videos, which would weight sites by n videos).

site_covs_mn <- all_sites %>%
  dplyr::select(CameraTrapID, all_of(names(mn_covs))) %>%
  filter(complete.cases(.))

cor_sites <- cor(site_covs_mn[, names(mn_covs)])
round(cor_sites, 2)

vif_sites <- sapply(names(mn_covs), function(v) {
  1 / (1 - summary(lm(reformulate(setdiff(names(mn_covs), v), response = v),
                      data = site_covs_mn))$r.squared)
})
round(vif_sites, 2) # all < 2 --> no collinearity problem

png(file.path(mn_out_dir, "Collinearity_site_covariates.png"),
    width = 2400, height = 2400, res = 300)
corrplot(cor_sites, method = "number", type = "upper",
         tl.col = "black", tl.cex = 0.8, number.cex = 0.8)
dev.off()


## Functions ---------------------------------------------------------------

# Multinomial mixed model, random intercept per camera trap.
# estimator = "REML": with the default "ML" the camera-trap variance
# collapses to ~0 (false convergence), which makes it a plain multinomial
# model treating every video as independent --> SEs about 2x too small.
fit_mn <- function(data, covs, controls = NULL) {
  f <- reformulate(c(paste0("scale(", covs, ")"), controls), response = "main_behav_f")
  eval(bquote(
    mblogit(formula = .(f),
            data = data,
            random = ~1|CameraTrapID,
            estimator = "REML")
  ))
}

# Joint Wald test per covariate (both equations together, df = 2).
# Replaces anova()-based pLRTs: for mixed mblogit models the deviance is a
# quasi-likelihood, and anova() itself warns that the LRT is unreliable.
wald_tests_mn <- function(model, covs, controls = NULL) {
  b <- coef(model)
  V <- vcov(model)
  keys <- c(setNames(paste0("scale(", covs, ")"), covs),
            setNames(paste0("~", controls), controls))
  bind_rows(lapply(names(keys), function(v) {
    idx <- grep(keys[[v]], names(b), fixed = TRUE)
    chi2 <- as.numeric(t(b[idx]) %*% solve(V[idx, idx]) %*% b[idx])
    tibble(covariate = v,
           chi2 = chi2,
           df = length(idx),
           p_value = pchisq(chi2, df = length(idx), lower.tail = FALSE))
  }))
}

# Predicted probabilities over the observed range of one covariate,
# all other covariates at their mean (young = proportion of videos with
# young), random effect = 0.
predict_mn <- function(model, data, cov, n = 50) {
  at <- list(seq(min(data[[cov]]), max(data[[cov]]), length.out = n))
  names(at) <- cov
  p <- emmeans::emmip(model,
                      as.formula(paste("main_behav_f ~", cov)),
                      at = at,
                      data = data,
                      plotit = FALSE,
                      CIs = TRUE)
  tibble(covariate = cov,
         x = p[[cov]],
         behaviour = p$main_behav_f,
         prob = p$yvar,
         lower = pmax(p$LCL, 0), # delta-method CIs can cross 0/1
         upper = pmin(p$UCL, 1))
}

plot_mn <- function(pred, title, xlab = NULL) {
  g <- ggplot(pred, aes(x = x, y = prob)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = behaviour), alpha = 0.3) +
    geom_line(aes(colour = behaviour), linewidth = 1) +
    scale_colour_manual(name = "Main behaviour", values = behav_colours) +
    scale_fill_manual(name = "Main behaviour", values = behav_fills) +
    labs(x = xlab, y = "Predicted probability", title = title) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw() +
    theme(axis.title = element_text(face = "plain", size = 14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          title = element_text(face = "bold", size = 16),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14, face = "bold"),
          panel.grid.minor = element_line(colour = "grey93"),
          panel.grid.major = element_line(colour = "grey93"),
          strip.background = element_rect(colour = "white", fill = "white"),
          strip.text = element_text(face = "bold", size = 12))
  if (is.null(xlab)) {
    # Panel labels below each x axis, so they work as x-axis titles
    g <- g + facet_wrap(~ covariate_label, scales = "free_x", ncol = 3,
                        strip.position = "bottom") +
      labs(x = NULL) +
      theme(strip.placement = "outside",
            strip.text = element_text(face = "plain", size = 14))
  }
  g
}

# Predicted probabilities without / with young (other covariates at mean)
predict_young_mn <- function(model, data) {
  p <- emmeans::emmip(model,
                      main_behav_f ~ young,
                      at = list(young = c(0, 1)),
                      data = data,
                      plotit = FALSE,
                      CIs = TRUE)
  tibble(young = factor(p$young, levels = c(0, 1), labels = c("No young", "Young present")),
         behaviour = p$main_behav_f,
         prob = p$yvar,
         lower = pmax(p$LCL, 0),
         upper = pmin(p$UCL, 1))
}

plot_young_mn <- function(pred, title) {
  ggplot(pred, aes(x = young, y = prob, colour = behaviour)) +
    geom_pointrange(aes(ymin = lower, ymax = upper),
                    position = position_dodge(width = 0.4), size = 0.8, linewidth = 1) +
    scale_colour_manual(name = "Main behaviour", values = behav_colours) +
    labs(x = "Presence of young", y = "Predicted probability", title = title) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(face = "plain", size = 14),
          title = element_text(face = "bold", size = 16),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14, face = "bold"),
          panel.grid.minor = element_blank())
}


## Fit, test and plot all species ------------------------------------------

mn_models <- list()
mn_tests <- list()
mn_preds <- list()
mn_preds_young <- list()

for (sp_label in names(mn_species)) {

  # Data: only videos with all covariates (MM09, MM11, MM16 have none)
  data_sp <- ungul_glmm %>%
    filter(taxon == mn_species[[sp_label]]) %>%
    filter(if_all(all_of(c(names(mn_covs_all), names(mn_controls))), ~ !is.na(.x)))
  print(table(data_sp$main_behav_f))

  # Model
  model_sp <- fit_mn(data_sp, names(mn_covs_all), names(mn_controls))
  print(summary(model_sp))
  print(performance::check_collinearity(model_sp))
  mn_models[[sp_label]] <- model_sp

  # Tests
  mn_tests[[sp_label]] <- wald_tests_mn(model_sp, names(mn_covs_all), names(mn_controls)) %>%
    mutate(species = sp_label, .before = 1)

  # Predicted probabilities for every continuous covariate
  pred_sp <- bind_rows(lapply(names(mn_covs_all), function(v)
    predict_mn(model_sp, data_sp, v))) %>%
    mutate(x = ifelse(covariate == "log_group", exp(x), x), # back to n individuals
           x = ifelse(covariate == "dens_road", x * 1000, x), # m/m2 --> km/km2
           species = sp_label,
           covariate_label = factor(mn_covs_all[covariate], levels = mn_covs_all))
  mn_preds[[sp_label]] <- pred_sp

  # One panel with all covariates
  file_sp <- gsub(" ", "", sp_label)
  png(file.path(mn_out_dir, paste0(file_sp, "_AllCovariates.png")),
      width = 4200, height = 3600, res = 300)
  print(plot_mn(pred_sp, title = sp_label))
  dev.off()

  # One figure per covariate
  for (v in names(mn_covs_all)) {
    png(file.path(mn_out_dir, paste0(file_sp, "_", v, ".png")),
        width = 2400, height = 1800, res = 300)
    print(plot_mn(filter(pred_sp, covariate == v),
                  title = sp_label,
                  xlab = mn_covs_all[[v]]))
    dev.off()
  }

  # Presence of young
  pred_young_sp <- predict_young_mn(model_sp, data_sp) %>%
    mutate(species = sp_label, .before = 1)
  mn_preds_young[[sp_label]] <- pred_young_sp
  png(file.path(mn_out_dir, paste0(file_sp, "_young.png")),
      width = 2400, height = 1800, res = 300)
  print(plot_young_mn(pred_young_sp, title = sp_label))
  dev.off()
}


## Results tables ----------------------------------------------------------

# Joint Wald tests (df = 2) per species and covariate
mn_tests_all <- bind_rows(mn_tests)
mn_tests_all
write.csv(mn_tests_all, file.path(mn_out_dir, "Wald_tests.csv"), row.names = FALSE)

# Coefficients
mn_coefs_all <- bind_rows(lapply(names(mn_models), function(s) {
  cf <- summary(mn_models[[s]])$coefficients
  tibble(species = s, term = rownames(cf),
         estimate = cf[, 1], se = cf[, 2], z = cf[, 3], p_value = cf[, 4])
}))
write.csv(mn_coefs_all, file.path(mn_out_dir, "Coefficients.csv"), row.names = FALSE)

# Predicted probability at the minimum and maximum of each covariate
# (replaces the "Average percentages" filters)
mn_minmax <- bind_rows(mn_preds) %>%
  group_by(species, covariate, behaviour) %>%
  filter(x == min(x) | x == max(x)) %>%
  mutate(at = ifelse(x == min(x), "min", "max")) %>%
  ungroup() %>%
  dplyr::select(species, covariate, behaviour, at, x, prob, lower, upper)
mn_minmax
write.csv(mn_minmax, file.path(mn_out_dir, "Predicted_min_max.csv"), row.names = FALSE)

# Predicted probability without / with young
mn_young_all <- bind_rows(mn_preds_young)
mn_young_all
write.csv(mn_young_all, file.path(mn_out_dir, "Predicted_young.csv"), row.names = FALSE)
