library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(readr)
library(lubridate)
library(zoo)
library(ggplot2)
library(patchwork)
library(data.table)
library(sf)
library(ggforce)

# Define reusable path
fig_path <- "C:/Users/angalletti/Desktop/DOUBLECHECK"

# Define drought color palette
drought_colors <- c(
  "Rainfall deficit"  = "#0072B2",
  "Rain to snow"      = "#56B4E9",
  "Winter recession"  = "#CC79A7",
  "Cold snow season"  = "#009E73",
  "Warm snow season"  = "#E69F00",
  "Hot and dry"          = "#DC143C",
  "Composite"         = "#FFD700"
)

##############################################################################
# RUN SWITCHES (set TRUE/FALSE to control what executes)
##############################################################################
RUN_PREPROCESSING <- TRUE
RUN_DECISION_TREE <- TRUE
RUN_CLASSIFICATION <- TRUE
RUN_FIGURES_CORE <- TRUE
RUN_FIGURES_EXTRA <- FALSE

# If TRUE, use cached MC result instead of running MC
USE_MC_CACHE <- TRUE





##############################################################################
# 01 PREPROCESSING (read inputs, aggregate, build anomalies)
##############################################################################
if (RUN_PREPROCESSING) {
##### Read basins properties #####
basins <- read_delim("basin_specs/basins.csv", delim = ";", 
                     locale = locale(decimal_mark = ".")) %>% 
  slice(1:88)


# Read study_basins.txt
study_basins <- read_csv("basin_specs/study_basins.txt", 
                         col_types = cols(
                           name = col_character(),
                           ID = col_integer(),
                           class = col_character(),
                           contrib = col_character(),
                           nodelink = col_character()
                         ))
study_basins <- study_basins %>%
  # Clean the contrib column: remove curly braces, split on commas, convert to integer vector
  mutate(contrib = str_remove_all(contrib, "[{}]"),
         contrib = str_split(contrib, ",") %>% map(~ as.integer(trimws(.x)))) %>%
  rowwise() %>%
  mutate(
    total_area = sum(basins$AREA[basins$ID %in% contrib], na.rm = TRUE),
    weighted_elev = sum(basins$AREA[basins$ID %in% contrib] *
                          basins$AVGQUOTE[basins$ID %in% contrib], na.rm = TRUE) / total_area
  ) %>%
  ungroup()

##### Read hydrometeo inputs #####

daily_avg_data <- read_csv("OUT_CONC/combined_out_average_snow_rain_we_snowmelt_sca_et_daily.txt")
temp_data <- read_csv("daily_data/temperature.txt")
flow_data <- read_csv("OUT_CONC/combined_out_flow_daily.txt")
##### aggregate at the study basin scale #####

# Create membership table from study_basins (handles nested basins)
study_membership <- study_basins %>%
  select(name, contrib) %>% 
  unnest(contrib) %>%          # Each row is a subbasin from the study basin
  rename(subbasin_id = contrib)

# Join daily_avg_data with study_membership and then with basins (to get area)
daily_join <- daily_avg_data %>%
  rename(subbasin_id = idba) %>%  # Rename for consistency
  inner_join(study_membership, by = "subbasin_id") %>%
  inner_join(basins %>% select(ID, AREA), by = c("subbasin_id" = "ID"))

# Group by study basin and time, and compute weighted averages for each variable
daily_avg_agg <- daily_join %>%
  group_by(name, time) %>%
  summarize(
    rain = sum(val_rain * AREA, na.rm = TRUE) / sum(AREA, na.rm = TRUE),
    snow = sum(val_snow * AREA, na.rm = TRUE) / sum(AREA, na.rm = TRUE),
    we   = sum(value_we * AREA, na.rm = TRUE) / sum(AREA, na.rm = TRUE),
    melt = sum(val_melt * AREA, na.rm = TRUE) / sum(AREA, na.rm = TRUE),
    et   = sum(value_et * AREA, na.rm = TRUE) / sum(AREA, na.rm = TRUE),
    pet  = sum(val_pet * AREA, na.rm = TRUE) / sum(AREA, na.rm = TRUE),
    .groups = "drop"
  )
daily_avg_agg <- daily_avg_agg %>% 
  mutate(time = as.Date(time))

# Loop over each study basin to compute its hourly adjusted temperature
temp_daily_list <- lapply(1:nrow(study_basins), function(i) {
  # Extract study basin info (including its weighted elevation and name)
  basin <- study_basins[i, ]
  
  # Copy the hourly temperature data
  temp_hourly <- temp_data
  
  # Compute the adjusted temperature for this basin using the lapse rate:
  # temp_adjusted = int + slo * weighted_elev
  temp_hourly$temp_adjusted <- temp_hourly$int + temp_hourly$slo * basin$weighted_elev
  
  # Append the study basin name to each record
  temp_hourly$name <- basin$name
  
  # Create a date column from the hourly timestamp
  temp_hourly$date <- as.Date(temp_hourly$time)
  
  # Aggregate hourly data to daily averages
  daily_temp <- aggregate(temp_adjusted ~ name + date, data = temp_hourly, FUN = mean, na.rm = TRUE)
  
  return(daily_temp)
})

# Combine all the results into one data frame
temp_daily <- do.call(rbind, temp_daily_list)

temp_daily <- temp_daily %>% mutate(date = as.Date(date))

daily_avg_agg <- left_join(daily_avg_agg, temp_daily %>% 
                             rename(time = date) %>% 
                             select(name, time, temp_adjusted),
                           by = c("name", "time"))


# Prepare the flow data: rename idno to nodelink and convert time to Date.
flow_data_agg <- flow_data %>%
  rename(nodelink = idno) %>%
  mutate(time = as.Date(time)) %>%
  # Convert study_basins nodelink to numeric before the join.
  inner_join(study_basins %>% mutate(nodelink = as.numeric(nodelink)) %>% select(name, nodelink),
             by = "nodelink") %>%
  select(name, time, value) %>%
  rename(daily_flow = value)

# Merge with daily_avg_agg
daily_avg_agg <- left_join(daily_avg_agg, flow_data_agg, by = c("name", "time"))


# Remove Feb 29 redistributing it on Feb 28 Mar 1
# Step 1: Create helper columns: year and date_str (month-day)
daily_avg_agg <- daily_avg_agg %>%
  mutate(year = format(time, "%Y"),
         date_str = format(time, "%m-%d"))

# Step 2: Extract Feb 29 rows for each study basin
feb29_df <- daily_avg_agg %>%
  filter(date_str == "02-29") %>%
  select(name, year, daily_flow, temp_adjusted, we, rain, snow, melt, et, pet) %>%
  # Rename each variable to have a "_feb29" suffix
  rename_with(~ paste0(.x, "_feb29"), -c(name, year))

# Step 3: Join the feb29 data back to the main dataset by name and year
daily_avg_agg <- daily_avg_agg %>%
  left_join(feb29_df, by = c("name", "year"))

# Step 4: Adjust values for Feb 28 and Mar 1
daily_avg_agg <- daily_avg_agg %>%
  mutate(
    # Group A variables: daily_flow, temp_adjusted, we are averaged with feb 29
    daily_flow = case_when(
      date_str %in% c("02-28", "03-01") & !is.na(daily_flow_feb29) ~ (daily_flow + daily_flow_feb29) / 2,
      TRUE ~ daily_flow
    ),
    temp_adjusted = case_when(
      date_str %in% c("02-28", "03-01") & !is.na(temp_adjusted_feb29) ~ (temp_adjusted + temp_adjusted_feb29) / 2,
      TRUE ~ temp_adjusted
    ),
    we = case_when(
      date_str %in% c("02-28", "03-01") & !is.na(we_feb29) ~ (we + we_feb29) / 2,
      TRUE ~ we
    ),
    # Group B variables: rain, snow, melt, et, pet: value of feb 29 is split between the two days to preserve volumes
    rain = case_when(
      date_str %in% c("02-28", "03-01") & !is.na(rain_feb29) ~ rain + (rain_feb29 / 2),
      TRUE ~ rain
    ),
    snow = case_when(
      date_str %in% c("02-28", "03-01") & !is.na(snow_feb29) ~ snow + (snow_feb29 / 2),
      TRUE ~ snow
    ),
    melt = case_when(
      date_str %in% c("02-28", "03-01") & !is.na(melt_feb29) ~ melt + (melt_feb29 / 2),
      TRUE ~ melt
    ),
    et = case_when(
      date_str %in% c("02-28", "03-01") & !is.na(et_feb29) ~ et + (et_feb29 / 2),
      TRUE ~ et
    ),
    pet = case_when(
      date_str %in% c("02-28", "03-01") & !is.na(pet_feb29) ~ pet + (pet_feb29 / 2),
      TRUE ~ pet
    )
  )

# Step 5: Remove Feb 29 rows and helper columns
daily_avg_agg <- daily_avg_agg %>%
  filter(date_str != "02-29") %>%
  select(-ends_with("_feb29"), -year, -date_str)

daily_avg_agg <- daily_avg_agg %>%
  # Compute total precipitation as rain + snow
  mutate(precip = rain + snow) %>%
  group_by(name) %>%
  arrange(time) %>%
  ungroup()


# Create a new data frame with smoothed daily variables
daily_smooth_agg <- daily_avg_agg %>%
  group_by(name) %>%
  arrange(time) %>%
  mutate(
    # Smooth the daily variables using a 30-day centered moving average (partial = TRUE to avoid NA at edges)
    rain_smoothed    = rollapply(rain, width = 30, FUN = mean, align = "center", partial = TRUE),
    snow_smoothed    = rollapply(snow, width = 30, FUN = mean, align = "center", partial = TRUE),
    melt_smoothed    = rollapply(melt, width = 30, FUN = mean, align = "center", partial = TRUE),
    precip_smoothed  = rollapply(precip, width = 30, FUN = mean, align = "center", partial = TRUE),
    swe_smoothed     = rollapply(we, width = 30, FUN = mean, align = "center", partial = TRUE),
    et_smoothed      = rollapply(et, width = 30, FUN = mean, align = "center", partial = TRUE),
    pet_smoothed     = rollapply(pet, width = 30, FUN = mean, align = "center", partial = TRUE),
    flow_smoothed    = rollapply(daily_flow, width = 30, FUN = mean, align = "center", partial = TRUE),
    temp_smoothed    = rollapply(temp_adjusted, width = 30, FUN = mean, align = "center", partial = TRUE)
  ) %>%
  ungroup() %>%
  select(name, time, rain_smoothed, snow_smoothed, melt_smoothed, 
         precip_smoothed, swe_smoothed, et_smoothed, pet_smoothed, 
         flow_smoothed, temp_smoothed)

# For cumulative melt, one approach is to compute it from the smoothed daily melt.
# Since cumulative melt resets on October 1, we need a water-year indicator and a water-day count.
daily_smooth_agg <- daily_smooth_agg %>%
  group_by(name) %>%
  arrange(time) %>%
 
  mutate(
    melt_1_smoothed = if_else(
      is.na(lag(swe_smoothed)), 
      0, 
      pmax(0, lag(swe_smoothed) - swe_smoothed)
    )
  ) %>%
  mutate(
    water_year  = if_else(month(time) >= 10, year(time) + 1, year(time)),
    water_start = if_else(month(time) >= 10,
                          as.Date(paste0(year(time), "-10-01")),
                          as.Date(paste0(year(time) - 1, "-10-01"))),
    water_day   = as.integer(time - water_start) + 1
  ) %>%
  group_by(name, water_year) %>%
  arrange(time) %>%
  mutate(
    cum_melt_smoothed = cumsum(melt_smoothed),
    cum_melt_1_smoothed = cumsum(melt_1_smoothed)
  ) %>%
  ungroup()

# This table groups by basin and by calendar day (mm-dd), and computes reference values.
# Add the new reference for cumulative melt calculated from the calendar-year cumulative melt (cum_melt_1_smoothed).
ref_daily <- daily_smooth_agg %>%
  mutate(day = format(time, "%m-%d")) %>%
  group_by(name, day) %>%
  summarize(
    ref_rain       = quantile(rain_smoothed, probs = 0.20, na.rm = TRUE),
    ref_snow       = quantile(snow_smoothed, probs = 0.20, na.rm = TRUE),
    ref_precip     = quantile(precip_smoothed, probs = 0.20, na.rm = TRUE),
    ref_flow       = quantile(flow_smoothed, probs = 0.20, na.rm = TRUE),
    ref_swe        = mean(swe_smoothed, na.rm = TRUE),
    ref_temp       = mean(temp_smoothed, na.rm = TRUE),
    ref_cum_melt   = mean(cum_melt_smoothed, na.rm = TRUE),  # water-year based
    ref_et         = mean(et_smoothed, na.rm = TRUE),
    ref_pet        = mean(pet_smoothed, na.rm = TRUE),
    ref_cum_melt_1 = mean(cum_melt_1_smoothed, na.rm = TRUE), # calendar-year based
    .groups = "drop"
  ) %>%
  select(name, day, starts_with("ref_"))

# Create anomalies by joining daily_avg_agg with ref_daily on 'name' and calendar day.
anomalies <- daily_smooth_agg %>%
  # Create a "day" column in "mm-dd" format from the time column
  mutate(day = format(time, "%m-%d")) %>%
  left_join(ref_daily %>% 
              select(name, day, ref_rain, ref_snow, ref_precip, ref_flow, 
                     ref_swe, ref_temp, ref_cum_melt, ref_et, ref_pet, ref_cum_melt_1),
            by = c("name", "day")) %>%
  # Compute anomalies: observed minus reference
  mutate(
    anomaly_rain       = rain_smoothed - ref_rain,
    anomaly_snow       = snow_smoothed - ref_snow,
    anomaly_precip     = precip_smoothed - ref_precip,
    anomaly_flow       = flow_smoothed - ref_flow,
    anomaly_swe        = swe_smoothed - ref_swe,
    anomaly_temp       = temp_smoothed - ref_temp,
    anomaly_cum_melt   = cum_melt_smoothed - ref_cum_melt,
    anomaly_cum_melt_1 = cum_melt_1_smoothed - ref_cum_melt_1,
    anomaly_et         = et_smoothed - ref_et,
    anomaly_pet        = pet_smoothed - ref_pet
  ) %>%
  # Keep only the identifier columns and anomaly columns, plus the flow one for drought characteristics calc
  select(name, time, flow_smoothed, starts_with("anomaly_"))
}


##############################################################################
# 02 DECISION TREE + WRAPPERS (functions only)
##############################################################################
if (RUN_DECISION_TREE) {




compute_deficits <- function(anomalies, params) {
  
  # 1. Flow deficits (using anomaly_flow and daily_flow)
  flow_deficits <- anomalies %>%
    filter(anomaly_flow < 0) %>%
    arrange(name, time) %>%
    group_by(name) %>%
    # Identify gaps: if the previous date is not exactly one day before, start a new group.
    mutate(gap = if_else(is.na(lag(time)) | 
                           ((lag(time) + 1 != time) & !(format(lag(time), "%m-%d") == "02-28" & format(time, "%m-%d") == "03-01")),
                         1, 0),
           group = cumsum(gap)) %>%
    group_by(name, group) %>%
    summarize(
      date_start         = min(time),
      date_end           = max(time),
      duration           = n(),
      total_deficit_flow = sum(anomaly_flow, na.rm = TRUE),
      avg_flow           = mean(flow_smoothed, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(duration >= params$runoff_duration_threshold) %>%
    mutate(
      event_id      = row_number(),
      avg_intensity = total_deficit_flow / duration
    )
  
  # 2. Precipitation deficits (using anomaly_precip)
  precip_deficits <- anomalies %>%
    filter(anomaly_precip < 0) %>%
    arrange(name, time) %>%
    group_by(name) %>%
    mutate(gap = if_else(is.na(lag(time)) | 
                           ((lag(time) + 1 != time) & !(format(lag(time), "%m-%d") == "02-28" & format(time, "%m-%d") == "03-01")),
                         1, 0),
           group = cumsum(gap)) %>%
    group_by(name, group) %>%
    summarize(
      date_start            = min(time),
      date_end              = max(time),
      duration              = n(),
      total_deficit_precip  = sum(anomaly_precip, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(duration >= params$rainfall_duration_threshold)
  
  # 3. Cumulative melt deficits (using anomaly_cum_melt)
  melt_deficits <- anomalies %>%
    filter(anomaly_cum_melt < 0) %>%
    arrange(name, time) %>%
    group_by(name) %>%
    mutate(gap = if_else(is.na(lag(time)) | 
                           ((lag(time) + 1 != time) & !(format(lag(time), "%m-%d") == "02-28" & format(time, "%m-%d") == "03-01")),
                         1, 0),
           group = cumsum(gap)) %>%
    group_by(name, group) %>%
    summarize(
      date_start         = min(time),
      date_end           = max(time),
      duration           = n(),
      total_deficit_melt = sum(anomaly_cum_melt, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(duration > params$melt_duration_threshold)
  
  # 4. SWE anomalies (storage deficits: below-normal snowpack)
  swe_anomalies <- anomalies %>%
    filter(anomaly_swe < 0) %>%
    arrange(name, time) %>%
    group_by(name) %>%
    mutate(
      gap   = if_else(is.na(lag(time)) | lag(time) + 1 != time, 1, 0),
      group = cumsum(gap)
    ) %>%
    group_by(name, group) %>%
    summarise(
      date_start        = min(time),
      date_end          = max(time),
      duration          = n(),
      mean_anom_SWE     = mean(anomaly_swe, na.rm = TRUE),
      min_anom_SWE      = min(anomaly_swe, na.rm = TRUE),
      total_deficit_SWE = sum(anomaly_swe, na.rm = TRUE),  # cumulative missing snow
      .groups = "drop"
    ) %>%
    filter(duration >= params$melt_duration_threshold)
  
  # 5. Positive temperature anomalies (hot spells)
  temperature_positive_anomalies <- anomalies %>%
    filter(anomaly_temp > 0) %>%   # <-- use anomaly_temp here
    arrange(name, time) %>%
    group_by(name) %>%
    mutate(
      gap   = if_else(is.na(lag(time)) | lag(time) + 1 != time, 1, 0),
      group = cumsum(gap)
    ) %>%
    group_by(name, group) %>%
    summarise(
      date_start   = min(time),
      date_end     = max(time),
      duration     = n(),
      mean_anom_T  = mean(anomaly_temp, na.rm = TRUE),
      max_anom_T   = max(anomaly_temp, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(
      duration >= 5,
      mean_anom_T > params$high_T_anom
    )
  
  # 6. Negative temperature anomalies (cold spells)
  temperature_negative_anomalies <- anomalies %>%
    filter(anomaly_temp < 0) %>%   # <-- and here
    arrange(name, time) %>%
    group_by(name) %>%
    mutate(
      gap   = if_else(is.na(lag(time)) | lag(time) + 1 != time, 1, 0),
      group = cumsum(gap)
    ) %>%
    group_by(name, group) %>%
    summarise(
      date_start   = min(time),
      date_end     = max(time),
      duration     = n(),
      mean_anom_T  = mean(anomaly_temp, na.rm = TRUE),
      min_anom_T   = min(anomaly_temp, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(
      duration >= 5,
      mean_anom_T < -params$high_T_anom
    )
  
  # Return a list containing all anomaly types
  return(list(
    flow_deficits                  = flow_deficits,
    precip_deficits                = precip_deficits,
    melt_deficits                  = melt_deficits,
    swe_anomalies                  = swe_anomalies,
    temperature_positive_anomalies = temperature_positive_anomalies,
    temperature_negative_anomalies = temperature_negative_anomalies
  ))
}

# Row‐wise classification wrapper: apply classification to each row of the deficits.
classification_row_wrapper <- function(deficit_data, params) {
  # For each row in the deficit data, apply classify_drought_events.
  classifications <- pmap_chr(deficit_data, function(...) {
    row <- as_tibble(list(...))
    classify_drought_events(row, params)
  })
  
  deficit_data <- ungroup(deficit_data)
  mutate(deficit_data, classification = classifications)
}

# Classification function for a single drought event.
classify_drought_events <- function(event, params) {
  if (params$verbose)
    message(paste("Checking event", event$event_id, "basin", event$name ,": from", event$date_start, "to", event$date_end))
  
  if (precipitation_deficit_MC(event, params)) {
    if (params$verbose)
      message("Precipitation deficit detected, checking deficit duration similarity...")
    if (precipitation_duration_matches_runoff_MC(event, params)) {
      if (params$verbose)
        message("Deficit duration matches, classifying as Rainfall deficit.\n")
      return("Rainfall deficit")
    } else {
      if (params$verbose)
        message("Deficit duration does not match, checking if drought in autumn.")
      if (drought_extends_to_snow_season_MC(event, params)) {
        if (params$verbose)
          message("Drought extends to snow season, classifying as Rain to snow.\n")
        return("Rain to snow")
      } else {
        if (params$verbose)
          message("Drought not continued because of snow accumulation start, classifying as Composite.\n")
        return("Composite")
      }
    }
  } else {
    if (params$verbose)
      message("No precipitation deficit detected, checking melt...")
    if (insufficient_melt_MC(event, params)) {
      if (params$verbose)
        message("Melt insufficient, checking if early melt or delayed...")
      if (swe_melted_MC(event, params)) {
        if (params$verbose)
          message("SWE already melted, classifying as Warm snow season.\n")
        return("Warm snow season")
      } else {
        if (params$verbose)
          message("SWE still present, checking drought start season...")
        if (drought_start_in_spring_MC(event, params)) {
          if (params$verbose)
            message("Cold snow season spring drought.\n")
          return("Cold snow season")
        } else {
          if (params$verbose)
            message("Checking if drought occurred in winter...")
          if (drought_in_winter_MC(event, params)) {
            if (params$verbose)
              message("Winter recession drought.\n")
            return("Winter recession")
          } else {
            if (params$verbose)
              message("Checking ET conditions...")
            if (temperature_high_anom_MC(event, params)) {
              if (params$verbose)
                message("High T detected, summer Hot and dry.\n")
              return("Hot and dry")
            } else {
              if (params$verbose)
                message("Unrelated to melt, classifying as Composite.\n")
              return("Composite")
            }
          }
        }
      }
    } else {
      if (params$verbose)
        message("Melt not an issue, checking ET...")
      if (temperature_high_anom_MC(event, params)) {
        if (params$verbose)
          message("High T detected, classifying as Hot and dry.\n")
        return("Hot and dry")
      } else {
        if (params$verbose)
          message("No temperature-related issues, classifying as Composite.\n")
        return("Composite")
      }
    }
  }
}

##### MC revised decision functions #####
precipitation_deficit_MC <- function(event, params) {
  event_start <- event$date_start
  event_end <- event$date_end
  
  overlapping_deficit <- deficits$precip_deficits %>%
    filter(name == event$name,
           date_end >= (event_start - params$P_lookback),
           date_start <= event_start)
  
  if (params$verbose) {
    message("Found ", nrow(overlapping_deficit), " overlapping precipitation deficit events.")
  }
  
  return(nrow(overlapping_deficit) > 0)
}

precipitation_duration_matches_runoff_MC <- function(event, params) {
  event_start <- event$date_start
  
  overlapping_deficit <- deficits$precip_deficits %>%
    filter(name == event$name,
           date_end >= (event_start - params$P_lookback),
           date_start <= event_start)
  
  preceding_rainfall_deficit <- overlapping_deficit %>%
    arrange(desc(date_end)) %>%
    slice(1)
  
  if (nrow(preceding_rainfall_deficit) == 0) {
    message("WARNING: No preceding rainfall deficit event found for comparison.")
    return(FALSE)
  }
  
  rainfall_deficit_duration <- preceding_rainfall_deficit$duration
  drought_duration <- event$duration
  
  is_similar <- rainfall_deficit_duration >= (params$Q_deficit_ratio * drought_duration)
  
  if (params$verbose) {
    if (is_similar) {
      message("Durations are similar: rainfall deficit duration = ", rainfall_deficit_duration, 
              ", drought duration = ", drought_duration)
    } else {
      message("Durations are NOT similar: rainfall deficit duration = ", rainfall_deficit_duration, 
              ", drought duration = ", drought_duration)
    }
  }
  
  return(is_similar)
}

drought_in_autumn_MC <- function(event, params) {
  event_start <- event$date_start
  event_end <- event$date_end
  
  event_dates <- seq.Date(event_start, event_end, by = "day")
  event_months <- unique(format(event_dates, "%m"))
  
  return(any(event_months %in% c("09", "10", "11", "12")))
}

drought_extends_to_snow_season_MC <- function(event, params) {
  event_end_month <- as.numeric(format(event$date_end, "%m"))
  snow_season_months <- c(11, 12, 1, 2, 3)
  
  return(event_end_month %in% snow_season_months)
}

insufficient_melt_MC <- function(event, params) {
  # Create a sequence of dates for the event period.
  event_dates <- seq(event$date_start, event$date_end, by = "day")
  # Remove Feb 29 if present (since your datasets have been cleansed of that date).
  event_dates <- event_dates[format(event_dates, "%m-%d") != "02-29"]
  
  event_day_str <- format(event_dates, "%m-%d")
  
  # Filter the global daily_smooth_agg for this event's basin and dates, and pull the cum_melt_smoothed values.
  melt_values <- daily_smooth_agg %>% 
    filter(name == event$name, time %in% event_dates) %>% 
    arrange(time) %>%
    pull(cum_melt_1_smoothed)
  
  # From ref_daily, filter for the same basin and for the corresponding calendar days (in "mm-dd" format).
  ref_series <- sapply(event_day_str, function(d) {
    ref_value <- ref_daily %>% 
      filter(name == event$name, day == d) %>% 
      pull(ref_cum_melt_1)
    if (length(ref_value) == 0) NA else ref_value
  })
  
  # Error checking: ensure we have values for each event day.
  if (length(melt_values) == 0 && length(ref_series) == 0) {
    stop("Both melt values and reference values are missing for the specified event duration.")
  } else if (length(melt_values) == 0) {
    stop("Melt values are missing for the specified event duration.")
  } else if (length(ref_series) == 0) {
    stop("Melt reference values are missing for the specified event duration.")
  }
  
  if (length(melt_values) != length(event_dates) || length(ref_series) != length(event_dates)) {
    stop("Mismatch in the number of event dates and the corresponding melt or reference values.")
  }
  
  # Calculate the number of days where the observed melt is below the reference.
  deficit_days <- sum(melt_values < ref_series)
  
  if (params$verbose) {
    message("Deficit days: ", deficit_days, " out of ", length(event_dates), " days.")
  }
  
  # Return TRUE if the fraction of deficit days meets or exceeds the threshold.
  return(deficit_days / length(event_dates) >= params$Melt_contrib)
}

# insufficient_melt_MC <- function(event, params) { # only considers melt in march-july
#   
#   event_dates <- seq(event$date_start, event$date_end, by = "day")
#   event_dates <- event_dates[format(event_dates, "%m-%d") != "02-29"]
#   
#   # Finestra rilevante: 01-mar .. 31-lug (valida per ogni anno)
#   md <- format(event_dates, "%m-%d")
#   is_relevant <- md >= "03-01" & md <= "07-31"
#   
#   # Se non ci sono giorni rilevanti nell'evento, decidi cosa fare:
#   # (qui: ritorno FALSE + optional messaggio)
#   if (!any(is_relevant)) {
#     if (isTRUE(params$verbose)) message("No relevant melt days (Mar 1–Jul 31) within the event window.")
#     return(FALSE)
#   }
#   
#   melt_values <- daily_smooth_agg %>%
#     filter(name == event$name, time %in% event_dates) %>%
#     arrange(time) %>%
#     pull(cum_melt_1_smoothed)
#   
#   ref_series <- sapply(md, function(d) {
#     v <- ref_daily %>%
#       filter(name == event$name, day == d) %>%
#       pull(ref_cum_melt_1)
#     if (length(v) == 0) NA else v
#   })
#   
#   if (length(melt_values) == 0 && length(ref_series) == 0) stop("Both melt values and reference values are missing.")
#   if (length(melt_values) == 0) stop("Melt values are missing.")
#   if (length(ref_series) == 0) stop("Melt reference values are missing.")
#   if (length(melt_values) != length(event_dates) || length(ref_series) != length(event_dates)) {
#     stop("Mismatch in the number of event dates and the corresponding melt or reference values.")
#   }
#   
#   # Considera SOLO i giorni rilevanti
#   melt_rel <- melt_values[is_relevant]
#   ref_rel  <- ref_series[is_relevant]
#   
#   # Se vuoi essere prudente anche con NA:
#   ok <- !is.na(melt_rel) & !is.na(ref_rel)
#   melt_rel <- melt_rel[ok]
#   ref_rel  <- ref_rel[ok]
#   
#   if (length(melt_rel) == 0) stop("No valid melt/reference pairs in the relevant window (after NA removal).")
#   
#   deficit_days <- sum(melt_rel < ref_rel)
#   
#   if (isTRUE(params$verbose)) {
#     message("Relevant deficit days: ", deficit_days, " out of ", length(melt_rel), " relevant days (Mar 1–Jul 31).")
#   }
#   
#   return(deficit_days / length(melt_rel) >= params$Melt_contrib)
# }



swe_melted_MC <- function(event, params) {
  event_start <- event$date_start
  
  # Old version: pulling the value exactly 20 days ago from daily_smooth_agg
  # swe_20_days_ago <- daily_smooth_agg %>% 
  #   filter(name == event$name, time == (event_start - 20)) %>% 
  #   pull(swe_smoothed)
  # if (params$verbose) {
  #   message("OLD: SWE at event start: ", swe_start_old, "; 20 days ago: ", swe_20_days_ago)
  # }
  # decision <- swe_start_old < params$SWE_melted_ratio * swe_20_days_ago
  
  # New version: use daily_smooth_agg and get the maximum SWE from the past 20 days
  swe_start <- daily_smooth_agg %>%
    filter(name == event$name, time == event_start) %>%
    pull(swe_smoothed)
  
  # Filter the SWE values for the 20 days prior to event_start (excluding event_start)
  swe_past <- daily_smooth_agg %>%
    filter(name == event$name, time >= (event_start - 60), time <= event_start) %>%
    pull(swe_smoothed)
  
  max_swe_past <- if (length(swe_past) > 0) max(swe_past, na.rm = TRUE) else NA
  
  if (params$verbose) {
    message("NEW: SWE at event start: ", swe_start, "; Maximum SWE in past 60 days: ", max_swe_past)
  }
  
  # Return TRUE if the current SWE is less than the threshold (SWE_melted_ratio * max_swe_past)
  return(swe_start < params$SWE_melted_ratio * max_swe_past)
}


drought_start_in_spring_MC <- function(event, params) {
  event_month <- as.numeric(format(event$date_start, "%m"))
  return(event_month %in% c(3, 4, 5))
}

temperature_high_MC <- function(event, params) {
  event_start <- event$date_start
  event_end <- event$date_end
  
  temp_above_threshold <- daily_smooth_agg %>%
    filter(name == event$name, time >= event_start, time <= event_end) %>% 
    filter(temp_smoothed > params$high_T) %>%
    nrow()
  
  total_days <- as.numeric(event_end - event_start) + 1
  
  if (params$verbose) {
    message("Temperature above threshold days: ", temp_above_threshold, " of ", total_days)
  }
  
  return(temp_above_threshold > (0.5 * total_days))
}

temperature_high_anom_MC <- function(event, params) {
  event_start <- event$date_start
  event_end <- event$date_end
  
  temp_above_threshold <- daily_smooth_agg %>%
    filter(name == event$name, time >= event_start, time <= event_end) %>%
    # Create a calendar-day identifier to join with the reference table
    mutate(day = format(time, "%m-%d")) %>%
    left_join(ref_daily %>% 
                filter(name == event$name) %>% 
                select(day, ref_temp), 
              by = "day") %>%
    # Compute the temperature anomaly for each day.
    mutate(temp_anomaly = temp_smoothed - ref_temp) %>%
    summarise(num_above = sum(temp_anomaly > params$high_T_anom, na.rm = TRUE)) %>%
    pull(num_above)
  
  
  total_days <- as.numeric(event_end - event_start) + 1
  
  if (params$verbose) {
    message("Temperature anomaly above threshold days: ", temp_above_threshold, " of ", total_days)
  }
  
  return(temp_above_threshold > (0.5 * total_days))
}


drought_in_summer_MC <- function(event, params) {
  event_start <- event$date_start
  event_end <- event$date_end
  
  event_months <- seq(event_start, event_end, by = "day") %>% 
    format("%m") %>% as.numeric()
  
  return(any(event_months %in% c(7, 8, 9)))
}

drought_in_winter_MC <- function(event, params) {
  event_start <- event$date_start
  event_end <- event$date_end
  
  event_months <- seq(event_start, event_end, by = "day") %>% 
    format("%m") %>% as.numeric()
  
  return(any(event_months %in% c(12, 1, 2)))
}



# Full workflow: compute deficits (from pre-computed anomalies) then classify events.
full_classification_workflow <- function(params, anomalies) {
  # Compute deficits (for flow, precipitation, and melt)
  deficits <<- compute_deficits(anomalies, params) # WARNING: deficits created globally, check here for debugging
  
  # Classify the flow deficits (you can later expand to classify precip or melt deficits)
  classified_events <- classification_row_wrapper(deficits$flow_deficits, params)
  
  return(classified_events)
}


# Monte Carlo helper: sample parameters from a parameter list.
sample_parameters <- function(param_list) {
  sampled_params <- lapply(param_list, function(x) {
    if (length(x) > 1) sample(x, 1) else x
  })
  return(sampled_params)
}

# Monte Carlo simulation function using the full workflow.
monte_carlo_simulation <- function(n_iter = 1000, param_list, anomalies) {
  results <- vector("list", n_iter)
  progress_step <- n_iter / 10
  
  for (i in 1:n_iter) {
    # Sample a new parameter set
    sampled_params <- sample_parameters(param_list)
    
    # Run the full workflow with the sampled parameters
    classified_events <- full_classification_workflow(sampled_params, anomalies)
    total_events <- nrow(classified_events)
    
    if (total_events == 0) {
      counts <- data.frame(classification = character(), count = numeric(), stringsAsFactors = FALSE)
      proportions <- data.frame(classification = character(), proportion = numeric(), stringsAsFactors = FALSE)
    } else {
      counts <- classified_events %>%
        group_by(classification) %>%
        summarise(count = n(), .groups = "drop")
      
      proportions <- counts %>%
        mutate(proportion = count / total_events)
    }
    
    # Save full output including events
    results[[i]] <- list(
      params      = sampled_params,
      counts      = counts,
      proportions = proportions,
      events      = classified_events  # <--- NEW FIELD ADDED HERE
    )
    
    if (i %% progress_step == 0) {
      cat(sprintf("Completed %d%% of iterations (%d/%d)\n", (i / n_iter) * 100, i, n_iter))
    }
  }
  
  return(results)
}


# Unified wrapper for direct analysis and Monte Carlo simulation.
# Unified wrapper: runs direct analysis (with fixed parameters) or Monte Carlo simulation (with parameter ranges)
unified_classification_wrapper <- function(mode = "direct", anomalies, fixed_params = NULL, param_list = NULL, n_iter = 1000) {
  if (mode == "direct") {
    if (is.null(fixed_params)) stop("For direct mode, please provide fixed_params.")
    classified_events <- full_classification_workflow(fixed_params, anomalies)
    return(classified_events)
  } else if (mode == "MC") {
    if (is.null(param_list)) stop("For MC mode, please provide param_list.")
    mc_results <- monte_carlo_simulation(n_iter = n_iter, param_list = param_list, anomalies = anomalies)
    return(mc_results)
  } else {
    stop("Invalid mode: choose 'direct' or 'MC'.")
  }
}



}


##############################################################################
# 03 CLASSIFICATION RUN (direct + MC/cache)
##############################################################################
if (RUN_CLASSIFICATION) {

# Fixed parameters for direct analysis
fixed_params <- list(
  start_date                 = as.Date("1992-10-01"),
  end_date                   = as.Date("2023-09-30"),
  runoff_duration_threshold  = 20,
  rainfall_duration_threshold= 20,
  melt_duration_threshold    = 20,
  P_lookback                 = 20,
  Q_deficit_ratio            = 0.6,
  Melt_contrib               = 0.8,
  SWE_melted_ratio           = 0.6,
  high_T                     = 10,
  high_T_anom                = 0.5,
  verbose                    = TRUE
)

# Parameter list with ranges for Monte Carlo simulation
param_list <- list(
  start_date                 = as.Date("1992-10-01"),
  end_date                   = as.Date("2023-09-30"),
  runoff_duration_threshold  = 20,        # range from 3 to 40
  rainfall_duration_threshold= 20,       # p_deficit range or value
  melt_duration_threshold    = 15,       # cum_melt range or value
  P_lookback                 = 15:30,       # range from 15 to 30
  Q_deficit_ratio            = seq(0.5, 0.9, by = 0.05),
  Melt_contrib               = seq(0.6, 0.9, by = 0.05),
  SWE_melted_ratio           = seq(0.5, 0.8, by = 0.05),
  high_T                     = 10,        # range from 7 to 15
  high_T_anom                = seq(0.5, 2, by = 0.1),
  verbose                    = FALSE        # typically keep verbose off in MC mode
)

# Example computation of deficits for checking variable attributes
deficits <- compute_deficits(anomalies,fixed_params)

# In direct mode, if verbose = TRUE, the two lines below log the output in an external file
log_con <- file("mylog.txt", open = "wt")
sink(log_con, type = "message")
# Call the unified classification wrapper for direct mode
result_direct <- unified_classification_wrapper(mode = "direct", anomalies = anomalies, fixed_params = fixed_params)
sink(type = "message")
save_result_dir <- function(result_direct, dir_path = "dir_outputs", prefix = "result_direct")  {
  # Ensure output directory exists
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
  
  # Create filename with timestamp
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  file_path <- file.path(dir_path, paste0(prefix, "_", timestamp, ".rds"))
  
  # Save to disk
  saveRDS(result_direct, file = file_path)
  message("Saved successfully to: ", file_path)
}
save_result_dir(result_direct)
# Call the unified classification wrapper for Monte Carlo simulation
# (Here, n_iter is set to 1000; adjust as needed.)
result_MC <- unified_classification_wrapper(mode = "MC", anomalies = anomalies, param_list = param_list, n_iter = 1000)
# OR!!!!!!!!!!! Load results
result_MC <- readRDS("MC_outputs/result_MC_20250720_214641.rds")

save_result_MC <- function(result_MC, dir_path = "MC_outputs", prefix = "result_MC") {
  # Ensure output directory exists
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
  
  # Create filename with timestamp
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  file_path <- file.path(dir_path, paste0(prefix, "_", timestamp, ".rds"))
  
  # Save to disk
  saveRDS(result_MC, file = file_path)
  message("Saved successfully to: ", file_path)
}

save_result_MC(result_MC)



##### MC GRID PLOT #####
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(patchwork)  # for plot_layout and wrap_plots

# Subset the first 10 iterations to inspect
result_MC_subset <- result_MC[1:10]
str(result_MC_subset)
params <- param_list
# For testing, we can create a smaller dataframe from a subset of MC results:
mc_df <- map_dfr(seq_along(result_MC), function(i) {
  res <- result_MC[[i]]
  # Convert sampled parameters (a list of scalars) into a one-row tibble,
  # then add the iteration number.
  param_df <- as_tibble(res$params) %>% mutate(iteration = i)
  # Get the counts data and add iteration info.
  count_df <- res$counts %>% mutate(iteration = i)
  # Join the parameter values with the counts data.
  left_join(count_df, param_df, by = "iteration")
})

# Alternatively, if you're using the full result_MC, substitute `result_MC_subset` with `result_MC`.
# For ease of plotting, you can use either raw counts or proportions.
# Let's suppose you want to switch between plotting 'count' or a new variable 'proportion'.
# For now, we assume your counts dataframe has a column 'count' from the MC simulation.

# Extract the names of parameters that vary (using the previously provided method)
original_varying_params <- names(param_list)[sapply(param_list, function(x) length(x) > 1)]
print(original_varying_params)

# Pivot the data to long format for the originally varying parameters.
mc_long <- mc_df %>%
  pivot_longer(
    cols = all_of(original_varying_params),
    names_to = "parameter",
    values_to = "param_value"
  )
# If you also want to be able to switch between plotting count or proportion:
# Suppose the counts data frame in each iteration includes the column 'count',
# and you want to be able to choose the y-axis variable.
y_axis_var <- "count"  # change to "proportion" if you want to plot proportions

# Optional: If your classifications in mc_df are stored in a column named "classification",
# you might need to ensure they're factors with the desired order.
desired_order <- c("Rainfall deficit", "Rain to snow", "Composite", 
                   "Warm snow season", "Cold snow season", "Winter recession", 
                   "Hot and dry")
mc_long <- mc_long %>% 
  mutate(classification = factor(classification, levels = desired_order))

# Create an empty list to store composite plots (one per varying parameter)
plot_list <- list()

for (i in seq_along(original_varying_params)) {
  var <- original_varying_params[i]
  
  # Subset mc_long for the current parameter
  mc_long_subset <- mc_long %>% filter(parameter == var)
  
  # Create a boxplot faceted by drought classification.
  # Adjust the theme and facet layout as needed.
  p_box <- ggplot(mc_long_subset, aes(x = factor(param_value), y = .data[[y_axis_var]])) +
    geom_boxplot() +
    facet_wrap(~ classification, scales = "free_x", nrow = 1) +
    labs(x = "", y = y_axis_var) +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA),
          strip.text = element_text(size = 10),
          plot.margin = margin(0, 1, 0, 1),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  # Optionally, add a narrow label plot on the right with the parameter name.
  p_label <- ggplot() +
    annotate("text", x = 1, y = 1, label = var, angle = 90, hjust = 0.5, size = 3.5) +
    theme_void() +
    theme(plot.margin = margin(0, 1, 0, 1))
  
  # Combine boxplot and label plot horizontally.
  combined <- p_box + p_label + plot_layout(widths = c(1, 0.05))
  plot_list[[var]] <- combined
}

# Stack all composite plots vertically (one per varying parameter)
final_plot <- wrap_plots(plot_list, ncol = 1) +
  plot_layout(heights = rep(1, length(plot_list)))

# Print final composite plot
print(final_plot)

# Optionally, save the plot:
png("plots_sensitivity/6param_multi_1000.png", width = 9000, height = 6500, res = 600)
print(final_plot)
dev.off()


}


##############################################################################
# 04 FIGURES
##############################################################################
if (RUN_FIGURES_CORE || RUN_FIGURES_EXTRA) {
##### Event count per duration per catchment ######
# Define the range of duration thresholds to test
duration_values <- 5: 40

# Initialize list to collect results
duration_results <- list()

for (dur in duration_values) {
  # Update threshold in fixed parameter set
  current_params <- fixed_params
  current_params$runoff_duration_threshold <- dur
  
  # Compute deficits
  deficits <- compute_deficits(anomalies, current_params)
  
  # Extract flow deficits
  flow_def <- deficits$flow_deficits
  
  # Count events per catchment
  count_df <- flow_def %>%
    group_by(name) %>%
    summarise(n_events = n(), .groups = "drop") %>%
    mutate(threshold = dur)
  
  duration_results[[as.character(dur)]] <- count_df
}

# Bind all results into one data frame
duration_summary <- bind_rows(duration_results)
# Ensure names are in the desired order
ordered_names <- study_basins$name
# Convert to factor with explicit levels in the desired order
duration_summary$name <- factor(duration_summary$name, levels = ordered_names)
summary_stats <- duration_summary %>%
  group_by(threshold) %>%
  summarise(
    mean_events = mean(n_events, na.rm = TRUE),
    lower = quantile(n_events, 0.10, na.rm = TRUE),
    upper = quantile(n_events, 0.90, na.rm = TRUE),
    .groups = "drop"
  )

ggplot() +
  # Ribbon (10th–90th percentile)
  geom_ribbon(
    data = summary_stats,
    aes(x = threshold, ymin = lower, ymax = upper, fill = "10th–90th percentile"),
    alpha = 0.5, show.legend = TRUE
  ) +
  
  # Individual catchments
  geom_line(
    data = duration_summary,
    aes(x = threshold, y = n_events, group = name, color = "Individual catchments"),
    size = 0.1, alpha = 0.35, show.legend = TRUE
  ) +
  
  # Mean line
  geom_line(
    data = summary_stats,
    aes(x = threshold, y = mean_events, color = "Mean"),
    size = 1.2, show.legend = TRUE
  ) +
  
  # Manual legends
  scale_color_manual(
    name = NULL,
    values = c("Individual catchments" = "grey40", "Mean" = "black"),
    guide = guide_legend(
      order = 1,
      override.aes = list(
        fill = NA,       # removes square
        shape = NA,      # removes markers
        size = c(0.4, 1.2),
        alpha = c(0.6, 1)
      )
    )
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("10th–90th percentile" = "grey80"),
    guide = guide_legend(
      order = 2,
      override.aes = list(
        color = NA,      # removes line
        fill = "grey80",
        size = NA
      )
    )
  ) +
  
  # Labels + Theme
  labs(
    x = "Duration threshold (days)",
    y = "Number of events"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.97, 0.97),
    legend.justification = c(1, 1),
    legend.box.background = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.spacing.y = unit(-10, "pt")
  )

# Save plot using ggsave
if (RUN_FIGURES_CORE) {
ggsave(
  filename = file.path(fig_path, "duration_sensitivity.png"),
  width = 3.3, height = 3, units = "in", dpi = 600, scale = 1.5
)
}


##### variability of type-parameter-event count #####

# Step 1: Compute SD and IQR of drought type counts for each param value
sensitivity_summary <- mc_long %>%
  group_by(parameter, param_value, classification) %>%
  summarise(
    sd_count = sd(count, na.rm = TRUE),
    iqr_count = IQR(count, na.rm = TRUE),
    .groups = "drop"
  )

# Step 2: Plot: bars = SD, points = IQR
ggplot(sensitivity_summary, aes(x = factor(param_value))) +
  geom_col(aes(y = sd_count, fill = parameter), position = "dodge", width = 0.8, show.legend = FALSE) +
  geom_point(aes(y = iqr_count), color = "black", shape = 21, size = 2, stroke = 0.3) +
  facet_grid(classification ~ parameter, scales = "free_x", switch = "y") +
  labs(
    x = "Parameter value",
    y = "Spread of drought type count (SD = bar, IQR = point)"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 9),
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA)
  )



##### event-based classification stability #####

# Create a single dataframe from all events across MC iterations
all_events_df <- map2_dfr(
  .x = result_MC,
  .y = seq_along(result_MC),
  .f = function(res, i) {
    res$events %>%
      mutate(
        iteration = i,
        event_uid = paste(name, group, sep = "_")
      )
  }
)

# Check how many times each event_uid appears
check_uid_counts <- all_events_df %>%
  count(event_uid) %>%
  filter(n != 1000)

# View unexpected cases
print(check_uid_counts)

# Summarize frequencies per event_uid and classification
classification_summary <- all_events_df %>%
  count(event_uid, classification) %>%
  group_by(event_uid) %>%
  mutate(
    total = sum(n),
    prop = n / total
  ) %>%
  arrange(event_uid, desc(prop))

# Get dominant and secondary classification per event
classification_ranks <- classification_summary %>%
  group_by(event_uid) %>%
  arrange(desc(prop)) %>%
  slice(1:2) %>%
  mutate(rank = row_number()) %>%
  select(event_uid, rank, classification, prop) %>%
  ungroup() %>%
  pivot_wider(
    names_from = rank,
    values_from = c(classification, prop),
    names_glue = "{if_else(rank == 1, 'dominant', 'secondary')}_{.value}"
  )



# Attach metadata (like basin name, group, start date, etc.)
event_metadata <- all_events_df %>%
  select(event_uid, name, group, date_start, date_end, duration, avg_flow, avg_intensity) %>%
  distinct()

# Final dataframe
event_stability_df <- classification_ranks %>%
  left_join(event_metadata, by = "event_uid")

# Extract basin name from event_uid
classification_ranks <- classification_ranks %>%
  mutate(basin = sub("_[0-9]+$", "", event_uid))

# Define thresholds
stability_thresholds <- c(0.9, 0.8, 0.7, 0.6)

# Function to compute basin stability for a given threshold
compute_basin_stability <- function(threshold) {
  classification_ranks %>%
    group_by(basin) %>%
    summarise(
      total_events = n(),
      stable_events = sum(dominant_prop >= threshold, na.rm = TRUE),
      stability_pct = stable_events / total_events * 100
    ) %>%
    arrange(desc(stability_pct))
}

# Apply the function to all thresholds
basin_stability_list <- set_names(
  map(stability_thresholds, compute_basin_stability),
  paste0("≥", stability_thresholds * 100, "%")
)

# Combine all into one long table with an added column for threshold
basin_stability_all <- bind_rows(
  map2(basin_stability_list, names(basin_stability_list), ~ mutate(.x, threshold = .y))
)



##### PLOTTING WATER YEARS DROUGHT VARIABLES #####

# Function to compute water year boundaries for a given water year label.
compute_water_year_bounds <- function(wy) {
  # For water year wy, the period runs from October 1 of (wy - 1) to September 30 of wy.
  water_year_start <- as.Date(paste0(wy - 1, "-10-01"))
  water_year_end   <- as.Date(paste0(wy, "-09-30"))
  list(start = water_year_start, end = water_year_end)
}

# Function to compute the intersection (cropped to water year) of a drought event.
# Input: an event (with event$date_start and event$date_end) and water year boundaries.
intersect_event_with_wy <- function(event, wy_bounds) {
  # Intersection start: later of the event start and water year start.
  int_start <- max(event$date_start, wy_bounds$start)
  # Intersection end: earlier of event end and water year end.
  int_end <- min(event$date_end, wy_bounds$end)
  if (int_start <= int_end) {
    return(list(start = int_start, end = int_end))
  } else {
    return(NULL)
  }
}

# Main plotting framework function.
# Input: catchment name.
plot_by_catchment <- function(catchment_name) {
  # Filter daily_smooth_agg for the selected catchment.
  catch_data <- daily_smooth_agg %>% filter(name == catchment_name)
  catch_deficits <<- result_direct %>% filter(name == catchment_name)
  
  # Determine the overall time span from the data.
  overall_start <- min(catch_data$time, na.rm = TRUE)
  overall_end   <- max(catch_data$time, na.rm = TRUE)
  
  # Determine the water years covered.
  # We'll compute the water year as: if month(time)>=10, then water_year = year(time)+1, else water_year = year(time)
  catch_data <- catch_data %>%
    mutate(water_year = if_else(month(time) >= 10, year(time) + 1, year(time)),
           day = format(time, "%m-%d"))  # add day for joining references
  
  water_years <- sort(unique(catch_data$water_year))
  
  # Create output directory for this catchment if it doesn't exist.
  out_dir <- file.path("subbasin_plots", catchment_name)
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # For each water year, produce a plot.
  for (wy in water_years) {
    bounds <- compute_water_year_bounds(wy)
    
    # Subset data for this water year.
    wy_data <- catch_data %>% filter(time >= bounds$start, time <= bounds$end)
    
    if(nrow(wy_data) == 0) next  # skip if no data
    
    # For reference values, we assume ref_daily has columns: name, day, ref_rain, ref_cum_melt, etc.
    wy_data <- wy_data %>%
      left_join(ref_daily %>% filter(name == catchment_name), by = c("name", "day"))
    
    # For flow deficits, subset the global deficits$flow_deficits for this catchment that intersect this water year.
    # Assume deficits is global.
    wy_flow_deficits <- deficits$flow_deficits %>%
      filter(name == catchment_name) %>%
      # Only keep those that intersect the water year.
      filter(date_end >= bounds$start & date_start <= bounds$end) %>%
      rowwise() %>%
      mutate(
        # Compute intersection boundaries.
        int_start = max(date_start, bounds$start),
        int_end   = min(date_end, bounds$end)
      ) %>%
      ungroup()
    
    # Now, create individual plots for each variable.
    # For each plot, the x-axis is time (or water day).
    
    # Precipitation: observed (precip_smoothed) vs reference (ref_precip)
    p_precip <- ggplot(wy_data, aes(x = time)) +
      geom_line(aes(y = precip_smoothed, color = "Obs")) +
      geom_line(aes(y = ref_precip, color = "Ref"), linetype = "dashed") +
      labs(title = paste(catchment_name, "Water Year", wy, "- Precipitation"),
           y = "Precipitation (smoothed)") +
      theme_minimal()
    
    # Cumulative snowmelt: observed (cum_melt_smoothed) vs reference (ref_cum_melt)
    p_melt <- ggplot(wy_data, aes(x = time)) +
      geom_line(aes(y = cum_melt_1_smoothed, color = "Obs")) +
      geom_line(aes(y = ref_cum_melt_1, color = "Ref"), linetype = "dashed") +
      labs(title = paste(catchment_name, "Water Year", wy, "- Cumulative Snowmelt"),
           y = "Cum. Melt (smoothed)") +
      theme_minimal()
    
    # Snow water equivalent: observed (swe_smoothed) vs reference (ref_swe)
    p_swe <- ggplot(wy_data, aes(x = time)) +
      geom_line(aes(y = swe_smoothed, color = "Obs")) +
      geom_line(aes(y = ref_swe, color = "Ref"), linetype = "dashed") +
      labs(title = paste(catchment_name, "Water Year", wy, "- SWE"),
           y = "SWE (smoothed)") +
      theme_minimal()
    
    # Temperature: observed (temp_smoothed) vs reference (ref_temp)
    p_temp <- ggplot(wy_data, aes(x = time)) +
      geom_line(aes(y = temp_smoothed, color = "Obs")) +
      geom_line(aes(y = ref_temp, color = "Ref"), linetype = "dashed") +
      labs(title = paste(catchment_name, "Water Year", wy, "- Temperature"),
           y = "Temp (smoothed)") +
      theme_minimal()
    
    # Flow: observed (flow_smoothed) vs reference (ref_flow), with drought ribbons.
    p_flow <- ggplot(wy_data, aes(x = time)) +
      geom_line(aes(y = flow_smoothed, color = "Obs")) +
      geom_line(aes(y = ref_flow, color = "Ref"), linetype = "dashed") +
      labs(title = paste(catchment_name, "Water Year", wy, "- Flow"),
           y = "Flow (smoothed)") +
      theme_minimal()
    
    # Add red ribbons to highlight drought events for flow.
    # Loop over each drought event in wy_flow_deficits and add a geom_rect.
    if(nrow(wy_flow_deficits) > 0){
      for(j in 1:nrow(wy_flow_deficits)){
        event <- wy_flow_deficits[j, ]
        p_flow <- p_flow +
          annotate("rect",
                   xmin = event$int_start,
                   xmax = event$int_end,
                   ymin = -Inf, ymax = Inf,
                   alpha = 0.2, fill = "red")
      }
    }
    
    # Combine all five plots into one composite (e.g., 5 rows)
    composite <- (p_precip / p_melt / p_swe / p_temp / p_flow) +
      plot_annotation(title = paste("Hydrological Variables for", catchment_name, "Water Year", wy))
    
    # Save the composite plot in the subdirectory.
    outfile <- file.path(out_dir, paste0("plot_", wy, ".png"))
    ggsave(outfile, composite, width = 10, height = 15, dpi = 300)
    
    message("Saved plot for water year ", wy, " to ", outfile)
  }
}

# Example call for a catchment (adjust catchment name as needed)
plot_by_catchment("Adige_Bronzolo")

for (name in study_basins$name){
  plot_by_catchment(name)
}

##### Summary statistics per catchment or up/down catchment class #####

summarize_drought_statistics <- function(classified_events, catchment = NULL, catchment_class = NULL) {
  # If a specific catchment is provided, filter by its name.
  if (!is.null(catchment)) {
    classified_events <- classified_events %>% filter(name %in% catchment)
  }
  
  # If a catchment class is provided, first extract the catchment names corresponding to that class from study_basins.
  if (!is.null(catchment_class)) {
    catchment_names <- study_basins %>% 
      filter(class %in% catchment_class) %>% 
      pull(name)
    
    classified_events <- classified_events %>% filter(name %in% catchment_names)
  }
  
  total_events <- nrow(classified_events)
  
  # Summarize statistics for each drought type.
  summary_df <- classified_events %>%
    group_by(classification) %>%
    summarize(
      count = n(),
      percent = if(total_events > 0) 100 * n() / total_events else NA,
      avg_duration = mean(duration, na.rm = TRUE),
      avg_intensity = mean(avg_intensity, na.rm = TRUE),
      avg_severity = mean(abs(total_deficit_flow), na.rm = TRUE),
      .groups = "drop"
    )
  
  # Optionally, if a single catchment is provided, compute average relative flow.
  # Relative flow = event's average flow divided by the overall average flow at that station.
  if (!is.null(catchment) && length(catchment) == 1) {
    overall_flow <- daily_smooth_agg %>% 
      filter(name == catchment) %>% 
      summarize(overall = mean(flow_smoothed, na.rm = TRUE)) %>% 
      pull(overall)
    
    rel_flow_df <- classified_events %>%
      mutate(relative_flow = avg_flow / overall_flow) %>%
      group_by(classification) %>%
      summarize(avg_relative_flow = mean(relative_flow, na.rm = TRUE), .groups = "drop")
    
    summary_df <- left_join(summary_df, rel_flow_df, by = "classification")
  }
  
  return(summary_df)
}
# The instruction below can be customized to summarize (then plot) statistics for:
# a single catchment, specifying "catchment = "Adige_Bronzolo" " or
# a class of catchments (upstream-downstream), specifying catchment_class = c("up","down") or catchment_class = "up"
stats_class <- summarize_drought_statistics(classified_events = result_direct, catchment_class = c("up","down"))

# Pie chart for the "up" class (stats_class)
p_up <- ggplot(stats_class, aes(x = "", y = percent, fill = classification)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = drought_colors, breaks = names(drought_colors), drop = FALSE) +
  labs(title = "", fill = "Drought Type") +
  theme_void() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5)) +
  guides(fill = guide_legend(nrow = 7, byrow = TRUE, title.position = "top"))

print(p_up)


##### METRICS BARPLOTS #####
# Assume result_direct is your classified events data frame, with at least these columns:
# name, duration, avg_intensity, total_deficit_flow, avg_flow, classification

# Compute overall average flow for each catchment from daily_smooth_agg.
overall_flow <- daily_smooth_agg %>%
  group_by(name) %>%
  summarise(overall_flow = mean(flow_smoothed, na.rm = TRUE), .groups = "drop")

# Add additional metrics to your events data: 
# - Relative flow: avg_flow / overall_flow
# - Abs intensity: absolute value of avg_intensity
# - Severity: absolute value of total_deficit_flow
events_df <- result_direct %>%
  left_join(overall_flow, by = "name") %>%
  mutate(relative_flow = avg_flow / overall_flow,
         abs_intensity = abs(avg_intensity),
         severity = abs(total_deficit_flow),
         norm_intensity = abs_intensity / overall_flow,
         norm_severity  = severity      / overall_flow)

# Define the metrics to plot; here, we use:
# duration, abs_intensity, severity, and relative_flow.
metric_cols <- c("duration", "abs_intensity", "severity", "relative_flow")


# Pre-squish the values for each metric in events_long
events_long <- events_long %>%
  group_by(metric) %>%
  mutate(
    # Compute the 5th and 95th percentiles for the current metric
    q_low = quantile(value, 0.05, na.rm = TRUE),
    q_high = quantile(value, 0.95, na.rm = TRUE),
    # Squish the value: values below q_low are set to q_low, those above q_high become q_high.
    value_squished = squish(value, range = c(q_low, q_high))
  ) %>%
  ungroup()


# Assume events_df has columns:
#   classification, duration, abs_intensity, severity, relative_flow
order_levels <- names(drought_colors)

events_df      <- dplyr::mutate(events_df,      classification = factor(classification, levels = order_levels))
stats_class <- dplyr::mutate(stats_class, classification = factor(classification, levels = order_levels))


theme_set(theme_minimal(base_size = 6))

# 1) Duration
p1 <- ggplot(events_df, aes(x = classification, y = duration, fill = classification)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.2) +
  coord_cartesian(ylim = c(0, 250)) +
  scale_fill_manual(values = drought_colors, breaks = order_levels) +
  scale_x_discrete(breaks = NULL) + # remove this line to restore x axis labels
  labs(y = "Duration (days)", x = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = rel(0.8)),
        axis.text.y = element_text(size = rel(0.7)),
        axis.text.x = element_text(angle = 45, hjust = 1)  
  )

# 2) Absolute Intensity
p2 <- ggplot(events_df, aes(x = classification, y = norm_intensity, fill = classification)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.2) +
  coord_cartesian(ylim = c(0, 0.3)) +
  scale_fill_manual(values = drought_colors, breaks = order_levels) +
  scale_x_discrete(breaks = NULL) + # remove this line to restore x axis labels
  labs(y = "Normalized Intensity (-)", x = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = rel(0.8)),
        axis.text.y = element_text(size = rel(0.7)),
        axis.text.x = element_text(angle = 45, hjust = 1) 
  )

# 3) Severity
p3 <- ggplot(events_df, aes(x = classification, y = norm_severity, fill = classification)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.2) +
  coord_cartesian(ylim = c(0, 30)) +
  scale_fill_manual(values = drought_colors, breaks = order_levels) +
  scale_x_discrete(breaks = NULL) + # remove this line to restore x axis labels
  labs(y = "Normalized Severity (s)", x = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = rel(0.8)),
        axis.text.y = element_text(size = rel(0.7)),
        axis.text.x = element_text(angle = 45, hjust = 1) 
  )

# 4) Relative Flow
p4 <- ggplot(events_df, aes(x = classification, y = relative_flow, fill = classification)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.2) +
  coord_cartesian(ylim = c(0, 2)) +
  scale_fill_manual(values = drought_colors, breaks = order_levels) +
  scale_x_discrete(breaks = NULL) + # remove this line to restore x axis labels
  labs(y = "Relative Flow (-)", x = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = rel(0.8)),
        axis.text.y = element_text(size = rel(0.7)),
        axis.text.x = element_text(angle = 45, hjust = 1)  
        )

# Combine into a 2×2 grid
final_plot <- (p1 | p2) / (p3 | p4) 

print(final_plot)

p_up <- {
  total <- sum(stats_class$percent, na.rm = TRUE)  # use 1 if percent in 0–1
  
  donut <- ggplot(stats_class, aes(x = 2, y = percent, fill = classification)) +
    geom_col(width = 1, color = "white", linewidth = 0.6) +            # colored slices + radial separators
    geom_col(data = data.frame(percent = total),
             aes(x = 2, y = percent),
             inherit.aes = FALSE, fill = NA, color = "black",          # inner & outer ring border
             linewidth = 0.3, width = 1) +
    coord_polar(theta = "y") +
    xlim(-0.2, 2.5) +                                                    # first number ↑ => bigger hole
    scale_fill_manual(values = drought_colors,
                      breaks = names(drought_colors), drop = FALSE) +
    labs(title = "", fill = "Drought Type") +
    theme_void() +
    theme(
      plot.margin       = margin(t = -300, r = 0, b = 60, l = 0),
      legend.position   = "bottom",
      legend.box.margin = margin(t = -100, r = 0, b = -100, l = 0),
      legend.text       = element_text(size = 8),
      legend.title      = element_text(size = 9),
      legend.key.height = unit(0.3, "lines"),  # ⬅ smaller vertical padding
      legend.key.width  = unit(0.5, "lines")
    ) +
    guides(fill = guide_legend(nrow = 14, byrow = TRUE, title.position = "top"))
  
  (donut) / guide_area() +
    plot_layout(heights = c(0.75, 0.25), guides = "collect")
}

# Side-by-side with 70/30 width split
final_pie <- (final_plot | p_up) + plot_layout(widths = c(0.7, 0.3))
print(final_pie)

if (RUN_FIGURES_CORE) {
png(file.path(fig_path, "unica.png"), width = 6.85, height = 3.6, units = "in", res = 600)
print(final_pie)       # or print(p) if p is a ggplot

}
dev.off()

##### PDF kernel density of events per type
library(dplyr)
library(ggplot2)

set.seed(1)

df <- events_df %>%
  mutate(
    classification = factor(classification, levels = order_levels),
    date_mid = as.Date(date_start) + floor((as.Date(date_end) - as.Date(date_start)) / 2),
    md = format(date_mid, "%m-%d")
  ) %>%
  filter(!is.na(date_mid), md != "02-29") %>%
  mutate(doy = as.integer(strftime(date_mid, "%j")))

bin_width <- 7
x_jitter  <- 0.9   # smaller jitter -> tighter piles (days)

df_stack <- df %>%
  mutate(week_bin = floor((doy - 1) / bin_width) * bin_width + 1) %>%
  group_by(week_bin) %>%
  arrange(classification, .by_group = TRUE) %>%     # pile by type order within each bin
  mutate(
    stack_y = row_number(),
    x_plot  = week_bin + runif(n(), -x_jitter, x_jitter)
  ) %>%
  ungroup()

p <- ggplot(df_stack, aes(x = x_plot, y = stack_y, color = classification)) +
  geom_point(size = 5, alpha = 0.85) +            # bigger dots
  scale_color_manual(values = drought_colors, breaks = order_levels, drop = FALSE) +
  scale_x_continuous(
    limits = c(1, 365),
    breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 365),
    labels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","")
  ) +
  labs(x = "Day of year (weekly bins; event midpoint)", y = "Stack height (events per week)", color = "Drought Type") +
  theme_minimal(base_size = 8) +
  theme(
    legend.position = "right",
    panel.grid.minor.y = element_blank()
  )

print(p)
#RUG PLOT
library(dplyr)
library(ggplot2)

df <- events_df %>%
  mutate(
    classification = factor(classification, levels = order_levels),
    date_mid = as.Date(date_start) + floor((as.Date(date_end) - as.Date(date_start)) / 2),
    md = format(date_mid, "%m-%d")
  ) %>%
  filter(!is.na(date_mid), md != "02-29") %>%
  mutate(doy = as.integer(strftime(date_mid, "%j")))

bin_width <- 30

rug_df <- df %>%
  mutate(week_bin = floor((doy - 1) / bin_width) * bin_width + 1) %>%
  group_by(week_bin) %>%
  arrange(classification, .by_group = TRUE) %>%     # stack in type order
  mutate(stack_y = row_number()) %>%
  ungroup()

# Each event becomes a short horizontal segment ("rug") centered on the bin
seg_width <- 26.0  # days; keep <= bin_width (7) for a clean weekly rug

rug_df <- rug_df %>%
  mutate(
    x0 = week_bin - seg_width / 2,
    x1 = week_bin + seg_width / 2,
    y0 = stack_y,
    y1 = stack_y
  )

p <- ggplot(rug_df) +
  geom_segment(aes(x = x0, xend = x1, y = y0, yend = y1, color = classification),
               linewidth = 2.9, lineend = "square") +
  scale_color_manual(values = drought_colors, breaks = order_levels, drop = FALSE) +
  scale_x_continuous(
    limits = c(1, 365),
    breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 365),
    labels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","")
  ) +
  labs(x = "Day of year (weekly bins; event midpoint)", y = "Stack height (events per week)", color = "Drought Type") +
  theme_minimal(base_size = 8) +
  theme(
    legend.position = "right",
    panel.grid.minor.y = element_blank()
  )

print(p)

library(dplyr)
library(ggplot2)

df <- events_df %>%
  mutate(
    classification = factor(classification, levels = order_levels),
    date_mid = as.Date(date_start) + floor((as.Date(date_end) - as.Date(date_start)) / 2),
    md = format(date_mid, "%m-%d")
  ) %>%
  filter(!is.na(date_mid), md != "02-29") %>%
  mutate(
    doy   = as.integer(strftime(date_mid, "%j")),
    month = as.integer(format(date_mid, "%m"))
  )

# Month boundaries in a non-leap year (doy start/end)
month_starts <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
month_ends   <- c(31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365)

month_df <- data.frame(
  month = 1:12,
  m_start = month_starts,
  m_end   = month_ends
) %>%
  mutate(
    m_center = (m_start + m_end) / 2,
    m_width  = (m_end - m_start + 1)
  )

rug_df <- df %>%
  left_join(month_df, by = "month") %>%
  group_by(month) %>%
  arrange(classification, .by_group = TRUE) %>%
  mutate(stack_y = row_number()) %>%
  ungroup() %>%
  mutate(
    seg_width = pmax(m_width - 2, 10),      # leave 1 day margin each side; min width safeguard
    x0 = m_center - seg_width/2,
    x1 = m_center + seg_width/2,
    y0 = stack_y,
    y1 = stack_y
  )

p <- ggplot(rug_df) +
  geom_segment(
    aes(x = x0, xend = x1, y = y0, yend = y1, color = classification),
    linewidth = 2.9, lineend = "square"
  ) +
  scale_color_manual(values = drought_colors, breaks = order_levels, drop = FALSE) +
  scale_x_continuous(
    limits = c(1, 365),
    breaks = month_starts,
    labels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  ) +
  labs(x = "Month (event midpoint)", y = "Stack height (events per month)", color = "Drought Type") +
  theme_minimal(base_size = 8) +
  theme(
    legend.position = "right",
    panel.grid.minor.y = element_blank()
  )

print(p)


# ggsave(file.path(fig_path, "events_doy_weekly_stack.png"), p, width = 6.85, height = 3.0, dpi = 600)




##### flow duration curves if ever  needed #####
plot_fdc_normalized <- function(basin_name) {
  # Filter for the selected basin
  basin_data <- daily_avg_agg[daily_avg_agg$name == basin_name, ]
  
  # Clean and normalize flow data
  flows <- basin_data$daily_flow
  flows <- flows[!is.na(flows)]
  mean_flow <- mean(flows)
  flows_norm <- flows / mean_flow
  
  # Sort and calculate exceedance
  flows_sorted <- sort(flows_norm, decreasing = TRUE)
  exceedance <- 100 * (rank(-flows_sorted, ties.method = "first") / length(flows_sorted))
  
  # Prepare data frame
  fdc_df <- data.frame(
    exceedance = exceedance,
    flow_norm = flows_sorted
  )
  
  # Plot with ggplot2
  ggplot(fdc_df, aes(x = exceedance, y = flow_norm)) +
    geom_line(color = "steelblue", size = 1) +
    scale_y_log10(limits = c(0.01, 10), breaks = c(0.01, 0.1, 1, 10)) +
    labs(
      title = paste("Flow Duration Curve -", basin_name),
      x = "Exceedance Probability (%)",
      y = "Normalized Flow (log scale)"
    ) +
    theme_minimal()
}


plot_fdc_normalized("Aurino_Cadipietra")

##### Conditional probability heatmaps #####
compute_conditional_probabilities <- function(
    anomaly_events, drought_events, anomaly_type,
    lead_times = c(30,60,90,120)) {
  
  results <- list()
  
  for (t in lead_times) {
    
    # assign median month of the event
    anomaly_events_mod <- anomaly_events %>%
      mutate(
        mid_date = date_start + floor(duration / 2),
        month    = lubridate::month(mid_date)
      )
    
    # ---- BINARY VERSION (default) ----
    # for each anomaly: did at least one drought start within [date_end, date_end + t] ?
    anomaly_events_mod <- anomaly_events_mod %>%
      rowwise() %>%
      mutate(
        drought_followed = any(
          drought_events$name == name &
            drought_events$date_start >= date_start &
            drought_events$date_start <= date_end + t
        )
      ) %>%
      ungroup()
    
    # ---- ALTERNATIVE: COUNT TOTAL DROUGHTS ----
    # anomaly_events_mod <- anomaly_events_mod %>%
    #   rowwise() %>%
    #   mutate(
    #     drought_count = sum(
    #       drought_events$name == name &
    #       drought_events$date_start >= date_end &
    #       drought_events$date_start <= date_end + t
    #     )
    #   ) %>%
    #   ungroup()
    # # Note: then you would compute prob = mean(drought_count) instead of mean(drought_followed)
    
    # aggregate pooled across all basins by month
    summary <- anomaly_events_mod %>%
      group_by(month) %>%
      summarise(
        total_anoms = n(),
        followed    = sum(drought_followed, na.rm = TRUE),  # binary version
        prob        = ifelse(total_anoms > 0, followed / total_anoms, NA_real_),
        # prob        = ifelse(total_anoms > 0, mean(drought_count), NA_real_),  # total-count version
        .groups = "drop"
      ) %>%
      mutate(
        type = anomaly_type,
        lead = t
      )
    
    results[[as.character(t)]] <- summary
  }
  
  bind_rows(results)
}

all_probs <- bind_rows(
  compute_conditional_probabilities(deficits$precip_deficits,  deficits$flow_deficits, "P"),
  compute_conditional_probabilities(deficits$melt_deficits,    deficits$flow_deficits, "MELT"),
  compute_conditional_probabilities(deficits$swe_anomalies,    deficits$flow_deficits, "SWE"),
  compute_conditional_probabilities(deficits$temperature_positive_anomalies, deficits$flow_deficits, "Thigh"),
  compute_conditional_probabilities(deficits$temperature_negative_anomalies, deficits$flow_deficits, "Tlow")
)


heatmap_data <- all_probs %>%
  unite("row_id", type, lead, sep = "_") %>%  # e.g. P_7, P_14 ...
  select(row_id, month, prob) %>%
  pivot_wider(names_from = month, values_from = prob)

heatmap_long <- heatmap_data %>%
  pivot_longer(cols = `1`:`12`, names_to = "month", values_to = "prob") %>%
  mutate(
    month = as.integer(month),
    row_id = factor(row_id, levels = unique(row_id))  # keep order
  )

# derive group from row_id (text before first "_"); adjust if your rule differs
# Define row order
row_levels <- unlist(
  lapply(c("P_low","MELT_low","SWE_low","SWE_high","T_low","T_high"),
         function(drv) paste(drv, c(7,14,30,60,90), sep = "_"))
)
sep_lines <- seq(5, length(row_levels)-1, by = 5) + 0.5

ggplot(heatmap_long, aes(x = month, y = row_id, fill = prob)) +
  geom_tile(color = "white", width = 1, height = 0.98) +
  geom_hline(yintercept = sep_lines, color = "black", size = 0.3) +
  scale_fill_gradient(low = "white", high = "red", na.value = "grey90") +
  scale_x_continuous(breaks = 1:12, labels = month.abb, expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Month", y = "Anomaly × Lead time", fill = "P(drought)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(heatmap_long, aes(x = month, y = row_id, fill = prob)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "green",
    mid = "yellow",
    high = "red",
    midpoint = 0.4,
    limits = c(0.2, 0.8),
    oob = scales::squish
  ) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(x = "Month", y = "Anomaly × Lead time", fill = "P(drought)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##### second take on conditional probability - quantiles monitoring #####
monthly_anomalies <- anomalies %>%
  mutate(year  = year(time),
         month = month(time)) %>%
  group_by(name, year, month) %>%
  summarise(
    # Flux: Precip (monthly sum of anomalies)
    precip_month = sum(anomaly_precip, na.rm = TRUE),
    
    # Flux: Melt (delta of cumulative anomaly = last - first of month)
    melt_month   = last(anomaly_cum_melt) - first(anomaly_cum_melt),
    
    # Storage: SWE (monthly mean anomaly)
    swe_month    = mean(anomaly_swe, na.rm = TRUE),
    
    # State: Temperature (monthly mean anomaly)
    temp_month   = mean(anomaly_temp, na.rm = TRUE),
    
    .groups = "drop"
  )

monthly_anomalies_q <- monthly_anomalies %>%
  group_by(name, month) %>%
  mutate(
    q_precip = ecdf(precip_month)(precip_month),
    q_melt   = ecdf(melt_month)(melt_month),
    q_swe    = ecdf(swe_month)(swe_month),
    q_temp   = ecdf(temp_month)(temp_month)
  ) %>%
  ungroup()

monthly_anomalies_flags <- monthly_anomalies_q %>%
  mutate(
    P_low    = q_precip < 0.2,
    MELT_low = q_melt   < 0.2,
    SWE_low  = q_swe    < 0.2,
    SWE_high = q_swe    > 0.8,
    T_low    = q_temp   < 0.2,
    T_high   = q_temp   > 0.8
  )

compute_conditional_probs_monthly <- function(flagged_months, drought_events,
                                              flag_var, anomaly_type,
                                              lead_times = c(7,14,30,60,90)) {
  
  results <- list()
  
  for (t in lead_times) {
    
    # Select months flagged for this anomaly type
    anoms <- flagged_months %>%
      filter(!!sym(flag_var)) %>%
      mutate(
        # Use last day of month as anchor
        month_end = as.Date(paste0(year, "-", month, "-", days_in_month(month)))
      )
    
    # For each anomaly month: check if a drought starts within [month_end, month_end + t]
    anoms <- anoms %>%
      rowwise() %>%
      mutate(
        drought_followed = any(
          drought_events$name == name &
            drought_events$date_start >= month_end &
            drought_events$date_start <= month_end + t
        )
      ) %>%
      ungroup()
    
    # Summarise pooled across basins by month
    summary <- anoms %>%
      group_by(month) %>%
      summarise(
        total_anoms = n(),
        followed    = sum(drought_followed, na.rm = TRUE),
        prob        = ifelse(total_anoms > 0, followed / total_anoms, NA_real_),
        .groups = "drop"
      ) %>%
      mutate(type = anomaly_type, lead = t)
    
    results[[as.character(t)]] <- summary
  }
  
  bind_rows(results)
}

all_probs <- bind_rows(
  compute_conditional_probs_monthly(monthly_anomalies_flags, deficits$flow_deficits, "P_low",    "P_low"),
  compute_conditional_probs_monthly(monthly_anomalies_flags, deficits$flow_deficits, "MELT_low", "MELT_low"),
  compute_conditional_probs_monthly(monthly_anomalies_flags, deficits$flow_deficits, "SWE_low",  "SWE_low"),
  compute_conditional_probs_monthly(monthly_anomalies_flags, deficits$flow_deficits, "SWE_high", "SWE_high"),
  compute_conditional_probs_monthly(monthly_anomalies_flags, deficits$flow_deficits, "T_low",    "T_low"),
  compute_conditional_probs_monthly(monthly_anomalies_flags, deficits$flow_deficits, "T_high",   "T_high")
)

heatmap_data <- all_probs %>%
  unite("row_id", type, lead, sep = "_") %>%   # e.g. "P_low_7", "SWE_low_30"
  select(row_id, month, prob) %>%
  pivot_wider(names_from = month, values_from = prob) %>%
  arrange(row_id)


# Reshape to long format
heatmap_long <- all_probs %>%
  mutate(row_id = paste(type, lead, sep = "_")) %>%
  select(row_id, type, lead, month, prob) %>%
  pivot_longer(cols = prob, names_to = "var", values_to = "prob") %>%
  mutate(
    month = as.integer(month),
    # Order rows: type blocks, lags increasing
    type = factor(type, levels = c("P_low","MELT_low","SWE_low","SWE_high","T_low","T_high")),
    lead = factor(lead, levels = c(7,14,30,60,90)),
    row_id = paste(type, lead, sep = "_"),
    row_id = factor(row_id, levels = unique(row_id))
  )

# Define row order
row_levels <- unlist(
  lapply(c("P_low","MELT_low","SWE_low","SWE_high","T_low","T_high"),
         function(drv) paste(drv, c(7,14,30,60,90), sep = "_"))
)
# Separator positions (after each block of 5 rows)
sep_lines <- seq(5, length(row_levels)-1, by = 5) + 0.5

# Plot
p <- ggplot(heatmap_long, aes(x = month, y = row_id, fill = prob)) +
  geom_tile(color = "white", width = 1, height = 0.98) +
  geom_hline(yintercept = sep_lines, color = "black", linewidth = 0.3) +
  scale_fill_gradient(low = "white", high = "red", na.value = "grey90") +
  scale_x_continuous(breaks = 1:12, labels = month.abb, expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Month", y = "Anomaly × Lead time", fill = "P(drought)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save output as PDF
if (RUN_FIGURES_CORE) {
pdf(file.path(fig_path, "cond_prob_heatmap.pdf"), 
    width = 9, height = 5)  # dimensions in inches
print(p)

}
dev.off()


##### per type trigger onset analysis #####
compute_conditional_probs_monthly_types <- function(flagged_months, drought_events,
                                                    flag_var, anomaly_type,
                                                    lead_times = c(7,14,30,60,90)) {
  results <- list()
  
  for (t in lead_times) {
    anoms <- flagged_months %>%
      filter(!!sym(flag_var)) %>%
      mutate(month_end = as.Date(paste0(year, "-", month, "-", days_in_month(month))),
             id = paste(name, year, month))
    
    linkages <- anoms %>%
      rowwise() %>%
      mutate(
        drought_classes = list(
          drought_events %>%
            filter(.data$name == name,
                   date_start >= month_end,
                   date_start <= month_end + t) %>%
            pull(classification)
        )
      ) %>%
      ungroup()
    
    expanded <- linkages %>%
      select(id, month, drought_classes) %>%
      tidyr::unnest_longer(drought_classes, values_to = "drought_class", keep_empty = TRUE) %>%
      mutate(drought_class = ifelse(is.na(drought_class), "none", drought_class))
    
    # P(any drought) per month
    any_summary <- expanded %>%
      group_by(month, id) %>%
      summarise(any_drought = any(drought_class != "none"), .groups = "drop") %>%
      group_by(month) %>%
      summarise(p_any = mean(any_drought), .groups = "drop")
    
    # Shares of drought types
    all_classes <- unique(drought_events$classification)
    
    type_summary <- expanded %>%
      filter(drought_class != "none") %>%
      distinct(id, month, drought_class) %>%
      count(month, drought_class) %>%
      group_by(month) %>%
      mutate(share = n / sum(n)) %>%
      ungroup() %>%
      complete(month, drought_class = all_classes, fill = list(share = 0))
    
    out <- type_summary %>%
      left_join(any_summary, by = "month") %>%
      mutate(type = anomaly_type, lead = t)
    
    results[[as.character(t)]] <- out
  }
  
  bind_rows(results)
}

# Anomaly–drought probabilities
all_probs_types <- bind_rows(
  compute_conditional_probs_monthly_types(monthly_anomalies_flags, result_direct, "P_low",    "P_low"),
  compute_conditional_probs_monthly_types(monthly_anomalies_flags, result_direct, "MELT_low", "MELT_low"),
  compute_conditional_probs_monthly_types(monthly_anomalies_flags, result_direct, "SWE_low",  "SWE_low"),
  compute_conditional_probs_monthly_types(monthly_anomalies_flags, result_direct, "SWE_high", "SWE_high"),
  compute_conditional_probs_monthly_types(monthly_anomalies_flags, result_direct, "T_low",    "T_low"),
  compute_conditional_probs_monthly_types(monthly_anomalies_flags, result_direct, "T_high",   "T_high")
) %>%
  mutate(
    share = ifelse(is.na(share), 0, share),
    p_any = ifelse(is.na(p_any), 0, p_any)
  )

# Define row order
row_levels <- unlist(
  lapply(c("P_low","MELT_low","SWE_low","SWE_high","T_low","T_high"),
         function(drv) paste(drv, c(7,14,30,60,90), sep = "_"))
)

# Plotting data
heatmap_plotdata <- all_probs_types %>%
  group_by(type, lead, month) %>%
  arrange(drought_class) %>%
  mutate(
    xmin = month - 0.5 + c(0, head(cumsum(share), -1)),
    xmax = month - 0.5 + cumsum(share),
    row_id = paste(type, lead, sep = "_")
  ) %>%
  ungroup() %>%
  mutate(row_id = factor(row_id, levels = row_levels),
         row_num = as.numeric(row_id))


# Separator positions (after each block of 5 rows)
sep_lines <- seq(5, length(row_levels)-1, by = 5) + 0.5

p <- ggplot(heatmap_plotdata) +
  geom_rect(aes(xmin = xmin, xmax = xmax,
                ymin = row_num - 0.5, ymax = row_num + 0.5,
                fill = drought_class)) +
  geom_hline(yintercept = sep_lines, color = "black", size = 0.3) +
  scale_x_continuous(breaks = 1:12, labels = month.abb, expand = c(0, 0)) +  # NEW
  scale_y_continuous(breaks = 1:length(row_levels), labels = row_levels, expand = c(0, 0)) +
  scale_fill_manual(values = c(
    "Hot and dry"="#DC143C","Warm snow season"="#E69F00","Rain to snow"="#56B4E9",
    "Rainfall deficit"="#0072B2","Cold snow season"="#009E73",
    "Winter recession"="#CC79A7","Composite"="#FFD700"
  ),
  breaks = c("Rainfall deficit","Rain to snow","Winter recession",
             "Cold snow season","Warm snow season","Hot and dry","Composite")) +
  labs(x = "Month", y = "Anomaly × Lead time", fill = "Drought type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8))

# Save output as PDF
if (RUN_FIGURES_CORE) {
pdf(file.path(fig_path, "shares_heatmap.pdf"), 
    width = 9, height = 5)  # dimensions in inches
print(p)

}
dev.off()

##### normalized line_ribbon MC plots + rotated #####
# Load results
result_MC <- readRDS("MC_outputs/result_MC_20250720_214641.rds")
original_varying_params <- names(param_list)[sapply(param_list, function(x) length(x) > 1)]
print(original_varying_params)
recode_class <- function(x) {
  dplyr::recode(trimws(x),
                "Heatwave" = "Hot and dry",
                "heatwave" = "Hot and dry")
}

result_MC <- lapply(result_MC, function(run) {
  if (!is.null(run$counts) && "classification" %in% names(run$counts)) {
    run$counts <- run$counts %>% mutate(classification = recode_class(classification))
  }
  if (!is.null(run$proportions) && "classification" %in% names(run$proportions)) {
    run$proportions <- run$proportions %>% mutate(classification = recode_class(classification))
  }
  if (!is.null(run$events) && "classification" %in% names(run$events)) {
    run$events <- run$events %>% mutate(classification = recode_class(classification))
  }
  run
})



mc_df <- map_dfr(seq_along(result_MC), function(i) {
  res <- result_MC[[i]]
  # Convert sampled parameters (a list of scalars) into a one-row tibble,
  # then add the iteration number.
  param_df <- as_tibble(res$params) %>% mutate(iteration = i)
  # Get the counts data and add iteration info.
  count_df <- res$counts %>% mutate(iteration = i)
  # Join the parameter values with the counts data.
  left_join(count_df, param_df, by = "iteration")
})
# Pivot the data to long format for the originally varying parameters.
mc_long <- mc_df %>%
  pivot_longer(
    cols = all_of(original_varying_params),
    names_to = "parameter",
    values_to = "param_value"
  )
# If you also want to be able to switch between plotting count or proportion:
# Suppose the counts data frame in each iteration includes the column 'count',
# and you want to be able to choose the y-axis variable.
y_axis_var <- "count"  # change to "proportion" if you want to plot proportions

# Optional: If your classifications in mc_df are stored in a column named "classification",
# you might need to ensure they're factors with the desired order.
desired_order <- c("Rainfall deficit", "Rain to snow", "Composite", 
                   "Warm snow season", "Cold snow season", "Winter recession", 
                   "Hot and dry")
mc_long <- mc_long %>% 
  mutate(classification = factor(classification, levels = desired_order))




# Define drought color palette
drought_colors <- c(
  "Rainfall deficit"  = "#0072B2",
  "Rain to snow"      = "#56B4E9",
  "Winter recession"  = "#CC79A7",
  "Cold snow season"  = "#009E73",
  "Warm snow season"  = "#E69F00",
  "Hot and dry"          = "#DC143C",
  "Composite"         = "#FFD700"
)

# Step 1 — Summarize medians and IQR from mc_long
mc_summary <- mc_long %>%
  group_by(classification, parameter, param_value) %>%
  summarize(
    Q25    = quantile(.data[[y_axis_var]], 0.25, na.rm = TRUE),
    median = median(.data[[y_axis_var]], na.rm = TRUE),
    Q75    = quantile(.data[[y_axis_var]], 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# Step 2 — Compute normalization baseline: median of per-param medians for each class
class_baselines <- mc_summary %>%
  group_by(classification) %>%
  summarize(median_baseline = median(median, na.rm = TRUE), .groups = "drop")

# Step 3 — Compute true per-class median event count across MC runs
library(purrr)
counts_df <- map_dfr(seq_along(result_MC), function(i) {
  result_MC[[i]]$counts %>%
    mutate(run_id = i)
})

class_median_n <- counts_df %>%
  group_by(classification) %>%
  summarize(median_n = median(count), .groups = "drop") %>%
  mutate(label = paste0("median events = ", median_n))

# Step 4 — Normalize and join labels
mc_summary_normalized <- mc_summary %>%
  left_join(class_baselines, by = "classification") %>%
  left_join(class_median_n, by = "classification") %>%
  mutate(
    norm_median = median / median_baseline,
    norm_Q25    = Q25    / median_baseline,
    norm_Q75    = Q75    / median_baseline
  )

# Desired drought type order
desired_order <- c(
  "Rainfall deficit",
  "Rain to snow",
  "Composite",
  "Warm snow season",
  "Cold snow season",
  "Winter recession",
  "Hot and dry"
)

mc_summary_normalized <- mc_summary_normalized %>%
  mutate(classification = factor(classification, levels = desired_order))

# Step 5 — Plot: one ribbon+line per parameter
plot_list <- list()
for (var in original_varying_params) {
  
  df <- mc_summary_normalized %>% filter(parameter == var)
  
  p_ribbon <- ggplot(df, aes(x = param_value, group = classification,
                             fill = classification, color = classification)) +
    geom_ribbon(aes(ymin = norm_Q25, ymax = norm_Q75), alpha = 0.2, color = NA) +
    geom_line(aes(y = norm_median), size = 1) +
    geom_text(data = df %>% distinct(classification, label),
              aes(x = Inf, y = Inf, label = label),
              hjust = 1.1, vjust = 1.2, size = 3, inherit.aes = FALSE) +
    facet_wrap(~ classification, scales = "free_x", nrow = 1) +
    scale_color_manual(values = drought_colors) +
    scale_fill_manual(values = drought_colors) +
    labs(x = "", y = "Relative frequency") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      strip.text = element_text(size = 10),
      plot.margin = margin(0, 1, 0, 1),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "none"
    )
  
  # Label panel with parameter name
  p_label <- ggplot() +
    annotate("text", x = 1, y = 1, label = var, angle = 90, hjust = 0.5, size = 3.5) +
    theme_void() +
    theme(plot.margin = margin(0, 1, 0, 1))
  
  plot_list[[var]] <- p_ribbon + p_label + plot_layout(widths = c(1, 0.05))
}

# Step 6 — Stack plots vertically
final_plot <- wrap_plots(plot_list, ncol = 1) +
  plot_layout(heights = rep(1, length(plot_list)))

# Step 7 — Save output
png(file.path(fig_path, "5_param_lines_ribbons_norm_annotated_1000_a.png"), width = 11.7, height = 8.3, units = "in", res = 600)
print(final_plot)
dev.off()

final_plot <- ggplot(mc_summary_normalized, aes(x = param_value, group = classification,
                                  fill = classification, color = classification)) +
  geom_ribbon(aes(ymin = norm_Q25, ymax = norm_Q75), alpha = 0.2, color = NA) +
  geom_line(aes(y = norm_median), linewidth = 1) +
  geom_text(data = mc_summary_normalized %>% distinct(classification, label),
            aes(x = Inf, y = Inf, label = label),
            hjust = 1.03, vjust = 1.2, size = 3, inherit.aes = FALSE) +
  facet_grid(classification ~ parameter, scales = "free") +
  scale_color_manual(values = drought_colors) +
  scale_fill_manual(values = drought_colors) +
  labs(x = NULL, y = "Relative frequency") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    strip.text = element_text(size = 9),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 7),
    legend.position = "none"
  )

if (RUN_FIGURES_CORE) {
png(file.path(fig_path, "rotated_7x5_mc_plot_a.png"), width = 6.85, height = 9.5, units = "in", res = 600)
print(final_plot)

}
dev.off()




##### Flow regime during each drought type #####
# --- Libraries
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(ggplot2)

# --- Parameters
station <- "Adige_Bronzolo"
out_dir <- "plots_q_deficits"   # create if missing
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Helpers: non-leap DOY axis (month ticks/labels)
month_starts <- yday(seq(as.Date("2001-01-01"), as.Date("2001-12-01"), by = "1 month"))
month_labs   <- month.abb

# --- Climatology (black line) for the station
clim <- ref_daily %>%
  filter(name == station) %>%
  mutate(doy = yday(as.Date(paste0("2001-", day)))) %>%  # non-leap reference year
  distinct(doy, ref_flow)

# --- Full daily flow series for the station (+ DOY and ref_flow join)
flow_st <- flow_data_agg %>%
  filter(name == station) %>%
  mutate(year = year(time),
         doy  = yday(time),
         day  = format(time, "%m-%d")) %>%
  left_join(ref_daily %>%
              filter(name == station) %>%
              select(day, ref_flow),
            by = "day")

# --- Drought types available for the station
types <- result_direct %>%
  filter(name == station) %>%
  distinct(classification) %>%
  pull(classification)

# --- Loop over drought types
for (target_typ in types) {
  
  # Select events of target type for the station
  ev_sel <- result_direct %>%
    filter(name == station,
           str_to_lower(classification) == str_to_lower(target_typ)) %>%
    mutate(y_start = year(date_start),
           y_end   = year(date_end))
  
  # Skip if none (defensive)
  if (nrow(ev_sel) == 0) next
  
  # Expand each event to the years it spans (1 or 2), then attach daily flow
  ev_year_flow <- ev_sel %>%
    rowwise() %>%
    mutate(years = list(unique(c(y_start, y_end)))) %>%
    unnest(years) %>%
    rename(ev_year = years) %>%
    ungroup() %>%
    left_join(flow_st, by = c("name", "ev_year" = "year"))
  
  # Ribbon segments: only deficit days within the event window and where flow < climatology
  ribbons <- ev_year_flow %>%
    filter(time >= date_start, time <= date_end) %>%
    filter(!is.na(daily_flow), !is.na(ref_flow), daily_flow <= ref_flow) %>%
    mutate(ymin = daily_flow,
           ymax = ref_flow,
           grp  = interaction(event_id, ev_year, drop = TRUE))
  
  # Red lines: full-year daily flow for all years involved in those events
  lines_red <- ev_year_flow %>%
    filter(!is.na(daily_flow)) %>%
    mutate(grp = interaction(event_id, ev_year, drop = TRUE))
  
  # Build plot
  p <- ggplot() +
    # Climatology (black)
    geom_line(data = clim, aes(doy, ref_flow), color = "black", linewidth = 0.7) +
    # Deficit ribbons (between climatology and event flow, only over deficit dates)
    geom_ribbon(data = ribbons, aes(doy, ymin = ymin, ymax = ymax, group = grp),
                fill = "grey", alpha = 0.35) +
    # Event flow lines (red, full year)
    #geom_line(data = lines_red, aes(doy, daily_flow, group = grp),
              #color = "red", alpha = 0.3, linewidth = 0.5) +
    # Axes and limits
    scale_x_continuous(" ", breaks = month_starts, labels = month_labs, limits = c(1, 366)) +
    scale_y_continuous(limits = c(0, 125)) +
    # Labels and theme
    ylab(expression(paste("Portata (m"^3," s"^-1,")"))) +
    ggtitle(paste0('Typical flow deficit during "', target_typ, '" droughts — ', station)) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.title.x = element_text(margin = margin(t = 6))
    )
  
  # Save
  file_safe <- gsub("[^[:alnum:]_]+", "_", target_typ)
  ggsave(filename = file.path(out_dir, paste0("Q_deficit_", file_safe, ".png")),
         plot = p, width = 8, height = 5, dpi = 300)
}

}
