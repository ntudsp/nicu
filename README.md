# README
Bhan, Lam

# Replication code for “Do neonates hear what we measure? Assessing neonatal ward soundscapes at the neonates’ ears”

The GitHub repository contains the code to replicate the analysis,
figures and tables for the paper titled: “Do neonates hear what we
measure? Assessing neonatal ward soundscapes at the neonates’ ears”.

The data that support the findings of this study are openly available in
NTU research  
data repository DR-NTU (Data) at <https://doi.org/10.21979/N9/8GHNGX>

The subheadings in this repository follows the headings in the paper
(after the \[Data Loading\] section) for consistency.

The following figures are produced by this replication code:

## Initialisation

## Data Loading

Download the dataset from Dataverse if it does not exist. Note that this
code repository is configure to ignore the <kbd>data</kbd> folder during
git commits due to the large size (3.4 GB) of the dataset <kbd>csv</kbd>
file.

    [1] "The dataset exists"

Load the “timeSeries1sec.csv” dataset file containing the measurement
data

## Data Preparation

### Adding sub location data due to bed position changes

The bed positions in the wards were shifted according to the schedule
below.

From NICU:

1.  Start to 2022-03-17 11:20: NICU-A
2.  2022-03-17 11:20 to 2022-03-21 12:45: NICU-B
3.  2022-03-21 12:45 onwards: NICU-A

From HD:

1.  Start to 2022-04-03 09:00: HD-A
2.  2022-04-03 09:00 to 2022-04-03 10:00: HD-B
3.  2022-04-03 10:00 onwards: HD-A

Changeover dates to mark out:

1.  start: 20220303 16:00:00 /
2.  change: 20220318 16:00:00 to 17:00:00
3.  change: 20220324 16:00:00 to 17:00:00

``` r
timeseries.dt <- timeseries.dt |>
        #account for transition period of \pm 15 mins
        dplyr::mutate(
                sub_location = case_when(
                        #trim start for both
                        datetime < ymd_hms(
                                "2022-03-03 17:00:00", tz = "Singapore"
                        ) ~ "transition",
                        #NICU-A
                        datetime < ymd_hms(
                                "2022-03-17 11:05:00", tz = "Singapore"
                                ) & location == "NICU" ~ "NICU-A",
                        #30 min gap; NICU-B
                        datetime > ymd_hms(
                                "2022-03-17 11:35:00", tz = "Singapore"
                                ) &
                                datetime < ymd_hms(
                                "2022-03-21 12:45:00", tz = "Singapore"
                                ) & 
                                location == "NICU" ~ "NICU-B",
                        #30 min gap; NICU-A and trim end for nicu
                        datetime > ymd_hms(
                                "2022-03-21 13:15:00", tz = "Singapore"
                                ) & 
                                datetime < ymd_hms(
                                "2022-03-24 21:00:00", tz = "Singapore"
                                ) &
                                location == "NICU" ~ "NICU-A",
                        #HD-A
                        datetime < ymd_hms(
                                "2022-04-03 08:45:00", tz = "Singapore"
                                ) & location == "HD" ~ "HD-A",
                        #30 min gap; HD-B
                        datetime > ymd_hms(
                                "2022-04-03 09:15:00", tz = "Singapore"
                                ) &
                                datetime < ymd_hms(
                                "2022-04-03 09:45:00", tz = "Singapore"
                                ) & location == "HD" ~ "HD-B",
                        #30 min gap: HD-A and trim end for HD
                        datetime > ymd_hms(
                                "2022-04-03 10:15:00", tz = "Singapore"
                                ) & 
                                datetime < ymd_hms(
                                "2022-04-13 15:00:00", tz = "Singapore"
                                ) &
                                location == "HD" ~ "HD-A",
                       .default = "transition"
                )
        ) |>
        dplyr::filter(!sub_location == "transition")
```

### Summarising metrics by time and computation of summary statistics

``` r
#metrics to be summarised
params.art.df <- data.frame(
        aweight=c("L[AS]","L[ASmax]","L[AS10]","L[AS50]","L[AS90]"),
        cweight=c("L[CS]","L[CSmax]","L[CS10]","L[CS50]","L[CS90]"),
        tuhms=c("T","T[max]","T[10]","T[50]","T[90]")
        )

# Function to calculate metrics
calculate_metrics <- function(data, prefix) {
  data |> summarise(
          !!sym(prefix) := if(prefix %in% c("LAS", "LCS")) meandB(score_1min) else mean(score_1min),
          !!sym(paste0(prefix, "max")) := max(score_1min),
          !!sym(paste0(prefix, "10")) := quantile(score_1min, 0.90),
          !!sym(paste0(prefix, "50")) := quantile(score_1min, 0.50),
          !!sym(paste0(prefix, "90")) := quantile(score_1min, 0.10)
    )
}

# Process and combine all metrics
timeseries_1hour_metrics <- timeseries.dt |>
        dplyr::mutate(hour = floor_date(datetime,unit = "hour"))
        
timeseries_1hour_metrics <-
        bind_cols(
                #a-weighted
                timeseries_1hour_metrics |> 
                        dplyr::filter(acoUnit == "dBA") |>
                        group_by(sub_location, microphone, hour) |>
                        group_modify(~ calculate_metrics(., "LAS")) |>
                        ungroup(),
                #c-weighted
                timeseries_1hour_metrics |>
                        dplyr::filter(acoUnit == "dBC") |>
                        group_by(sub_location, microphone, hour) |>
                        group_modify(~ calculate_metrics(., "LCS")) |>
                        ungroup() |>
                        dplyr::select(c("LCS":"LCS90")),
                #tonality
                timeseries_1hour_metrics |>
                        dplyr::filter(acoUnit == "tuHMS") |>
                        group_by(sub_location, microphone, hour) |>
                        group_modify(~ calculate_metrics(., "T")) |>
                        ungroup() |>
                        dplyr::select(c("T":"T90"))
                ) |>
        dplyr::mutate(
                microphone = as.factor(microphone),
                sub_location = as.factor(sub_location),
                `LCS-LAS` = LCS - LAS,
                `LAS10-LAS90` = LAS10 - LAS90,
                location = ifelse(
                        grepl("NICU",sub_location),
                        "NICU",
                        "HD"
                        ), .before="sub_location"
                )
```
