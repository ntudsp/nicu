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
(after the [Data Preparation](#dataprep) section) for consistency.

The following figures are produced by this replication code:

- [`Figure 4`](#fig_4) in [4.3 Occurence rates](#sec4.3)

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

## Acoustic variation within-between wards

### 4.3 Occurrence rates

#### `Figure 4`: Occurrence rate of $\textit{OR}^h_\text{SNR}(5)$ and $\textit{OR}_{T}^h(0.4)$ averaged over the same daily 1-h period throughout the entire measurement campaign. A-weighted decibel metrics and tonality metrics averaged by hour of the day across the entire measurement duration at NICU-A, NICU-B and HD-A measurement points.

``` r
SNR = 5 #SNR threshold for dBA occurence rate

#compute occurrence rate
OR_dBA_tuHMS <- timeseries.dt |>
        #filter only binaural and dBA and tuHMS metrics
        dplyr::filter(
                acoUnit %in% c("dBA","tuHMS") &
                        microphone %in% c("binL","binR") #&
                        # datetime < ymd_hms("2022-03-05 23:00:00",
                        #                tz = "Singapore") #for testing
                        ) |> 
        #prepare data frame
        pivot_wider(
                names_from = microphone,
                values_from = score_1min,
                values_fn = ~ mean(.x, na.rm = TRUE) #to handle duplicates
        ) |>
        #aggregate both binL and binR into binLR
        dplyr::mutate(
                binLR = ifelse(
                        acoUnit == "dBA",
                        10*log10((10^(binL/10)+10^(binR/10))/2),
                        (binL+binR)/2
                )
        ) |>
        #drop binL and binR columns
        dplyr::select(!c(binL,binR)) |>
        #filter only required sub locations
        dplyr::filter(sub_location %in% c("NICU-A","NICU-B","HD-A")) |>
        #factorise sub_location, create hour and minute column
        dplyr::mutate(
                sub_location = as.factor(sub_location),
                minute = lubridate::floor_date(datetime, unit = "minute"),
                hour = lubridate::floor_date(datetime, unit = "hour")
                ) |>
        #aggregate by minute
        group_by(sub_location, acoUnit, minute) |>
        dplyr::mutate(
                # binLR_1h_mean = ifelse(
                #         acoUnit == "dBA",
                #         meandB(binLR_1min_mean),
                #         mean(binLR_1min_mean)
                # 50% exceedance level
                binLR_1min_max = max(binLR)
        ) |>
        # summarise_by_time(
        #         .date_var = datetime,
        #         .by       = "minute",
        #         binLR_1min_max=max(binLR),
        #         binLR_1min_mean=ifelse(
        #                 acoUnit == "dBA",
        #                 meandB(binLR),
        #                 mean(binLR)
        #                 )
        #         ) |>
        ungroup() |>
        distinct() |>
        #find 1h average as the noise floor
        dplyr::group_by(sub_location,acoUnit,hour) |>
        dplyr::mutate(
                # binLR_1h_mean = ifelse(
                #         acoUnit == "dBA",
                #         meandB(binLR_1min_mean),
                #         mean(binLR_1min_mean)
                # 50% exceedance level
                binLR_50 = quantile(binLR, 0.5)
        ) |>
        ungroup() |>
        distinct () |>
        summarise(
                .by = c(sub_location,acoUnit,hour),
                count = n(), #number of events in an hour
                OR = ifelse(
                        acoUnit == "dBA", #if dBA
                        #find the number of events > SNR of 5 dBA
                        sum(
                                binLR_1min_max > (binLR_50 + SNR),
                                na.rm = TRUE
                                )*(100/count), #1min max SPL > 1h SPL
                        sum(binLR_1min_max > 0.4)*(100/count) #tonality > 0.4
                )
        ) |>
        distinct() |>
        #add column to aggregate hour and convert to 0 to 23 range
        dplyr::mutate(hour_only = as.numeric(hour(hour)))

or_plot <- 
        ggplot(
                OR_dBA_tuHMS, 
                aes(x = as.numeric(hour_only), #convert to numeric for plot
                    y = OR, 
                    color = sub_location, 
                    fill = sub_location)
                )  + 
        facet_wrap(
                ~acoUnit, 
                nrow = 2,
                scales = "free",
                labeller = labeller(
                        acoUnit = as_labeller(
                                c("dBA" = "italic(L)[AS]", "tuHMS" = "italic(T)"),
                                label_parsed)
                        )
                ) + 
        stat_summary(
                fun = mean,
                geom = "line",
                linewidth = 1
        ) +
        stat_summary(
                fun = mean,
                geom = "ribbon",
                alpha = .3,
                linetype = 0,
                fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),
                fun.min = function(x) mean(x) - sd(x) / sqrt(length(x))
        ) +
        scale_color_paletteer_d("Redmonder::qPBI", direction = -1) +
        scale_fill_paletteer_d("Redmonder::qPBI", direction = -1) +
        xlim(0,24) + ylim(0,100) +
        ylab("% of time over threshold") +
        xlab("Hour") +
        scale_x_continuous(
                breaks = seq(0, 24, by = 6),
                minor_breaks = seq(0, 24, 2)
                ) +
        ggthemes::theme_hc() + 
        theme(panel.grid.minor.y = element_line(color = 1,
                                                size = 0.1,
                                                linetype = 3)) +
        theme(panel.grid.major.x = element_line(color = 1,
                                                size = 0.1,
                                                linetype = 1)) +
        theme(panel.grid.minor.x = element_line(color = 1,
                                                size = 0.1,
                                                linetype = 3))

or_plot
```

![](README_files/figure-commonmark/occurence-1.png)

``` r
ggsave(path = "output/",
       plot = or_plot,
       filename = "or_dBA_tuHMS.pdf",
       device = "pdf",
       units = "px",
       width = 1600,
       height = 1500,
       #scale = 3.5,
       scale = 2.0,
       dpi = 600) 
```
