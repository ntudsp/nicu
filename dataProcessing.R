#This is a script containing functions to extract the data processed by 
#HEADAcoustics ARTEMIS
#Author: Bhan Lam (NTU)
#Date: 26 Apr 2023

# Extraction of Level vs. Time data
# Arguments:
#       filename: full file path and filename
#       col_names: list(character()) | list of column names
#       col_types: character() | string of characters indicating data types (see:read_table)
#       origin_datetime: POSIXct() | start datetime of measurement first session/run
#       location: character() | location name
#       acousticUnit: character() | acoustic parameter units
artemisTimeSeriesExtract<-function(filename,origin_datetime,
                                   location,acoUnit){
        
        #If not L(A) or L(C)
        if (origin_datetime==FALSE) {
                #skip 19 rows instead of 20
                skiprows = 19;
        } else {
                skiprows = 20;
        }
        
        #update colnames based on location
        if (location=="NICU") {
                colNames<-c("time","binL","binR","146AEIn","146AEOut")
        } else {
                colNames<-c("time","binL","binR","146AE")
        }
        
        #retrieve recording start time from "(A) Slow"
        #read time series data header to extract start time
        start.datetime<-read_table(
                file = paste0(paste0(fileName[1:3],collapse =""),
                              ".Level vs. Time (A) Slow.csv"),
                skip = 6,n_max = 1,col_names = c("name","ymd","hms","tz"),
                col_types = "cccc")
        start.datetime<-ymd_hms(paste0(start.datetime$ymd, start.datetime$hms),
                                tz = "Singapore")
        
        #read csv file        
        timeseries.df<-read_table(file = paste0(filename,collapse = ""), 
                   skip = skiprows, 
                   col_names = colNames,
                   col_types = paste0(rep("d",length(colNames)),collapse=""),
                   skip_empty_rows = TRUE) %>%
                #convert to datetime and set origin time
                dplyr::mutate(time=as.duration(time),
                              datetime=start.datetime+time) %>%
                pivot_longer(cols=colNames[2:length(colNames)],
                             names_to = "microphone", 
                             values_to = "score") %>%
                dplyr::mutate(microphone=as.factor(microphone),
                              location=as.factor(location),
                              acoUnit=acoUnit)
        
        return(timeseries.df)
}


oldartemisTimeSeriesExtract<-function(filename,origin_datetime,
                                   location,acoUnit){
        
        #read csv file        
        timeseries.df<-read_table(file = filename, 
                                  skip = 19, 
                                  col_names = c("time","binL","binR","146AEIn","146AEOut"),
                                  col_types = "ddddd",
                                  skip_empty_rows = TRUE) %>%
                #convert to datetime and set origin time
                dplyr::mutate(time=as.duration(time),
                              datetime=as_datetime(time,origin=origindatetime),
                              .after="time") %>%
                pivot_longer(cols=c("binL","binR","146AEIn","146AEOut"),
                             names_to = "microphone", 
                             values_to = "score") %>%
                dplyr::mutate(microphone=as.factor(microphone),
                              location=as.factor(location),
                              acoUnit=acoUnit)
        
        return(timeseries.df)
}


#       location: character() | location name
#       acousticParam: character() | first acoustic parameter in artemis

artemisSingleValExtract<-function(filename,location,origin_datetime,
                                  acousticParam){
        
        #If not L(A) or L(C)
        if (origin_datetime==FALSE) {
                #skip 7 instead of 8 rows as there is no "recording date" row
                skiprows <- 7
                #if location is "HD" there are only 3 mics
                if (location == "HD") {nmax <- 4} else {nmax <- 3}
        } else {
                skiprows <- 8
                #if location is "HD" there are only 3 mics
                if (location == "HD") {nmax <- 4} else {nmax <- 3}
        }
        
        
        #retrieve recording start time from "(A) Slow"
        #read time series data header to extract start time
        start.datetime<-read_table(
                file = paste0(paste0(fileName[1:3],collapse =""),
                              ".Level vs. Time (A) Slow.csv"),
                skip = 6,n_max = 1,col_names = c("name","ymd","hms","tz"),
                col_types = "cccc")
        start.datetime<-ymd_hms(paste0(start.datetime$ymd, start.datetime$hms),
                                tz = "Singapore")
        
        df.singleValue<-read_delim(file = paste0(filename,collapse = ""),
                                   delim = ";",
                                   skip = skiprows,
                                   n_max = nmax,
                                   col_names = FALSE,
                                   show_col_types = FALSE) %>%
                
                #separate single value data in first col with "," delimiter
                tidyr::separate_wider_delim(X1, ",",
                                            names=c("microphone",NA,NA,
                                                    acousticParam),
                                            too_few = "align_start") %>%
                #separate numeric data
                tidyr::separate_wider_delim(acousticParam:X8,"=",
                                            names_sep = ":") %>%
                #replace tuHMS units if exists
                dplyr::mutate(!! names(.)[3] := gsub(" tuHMS","",.[[3]])) %>%
                #separate microphone data
                tidyr::separate_wider_delim(microphone,":",
                                            names=c(NA,"microphone")) %>%
                #separate SPL data
                tidyr::separate_wider_delim(c(`X6:2`,`X7:2`,`X8:2`)," ",
                                            names_sep=":") %>%
                #remove redundant cols
                dplyr::select(!c(`X6:2:2`,`X7:2:2`,`X7:2:3`,`X7:2:4`,
                                 `X8:2:2`,`X8:2:3`,`X8:2:4`)) %>%
                #tidy columns with names
                pivot_wider(names_from = paste0(acousticParam,":1"), 
                            values_from = paste0(acousticParam,":2")) %>%
                pivot_wider(names_from = `X2:1`, values_from = `X2:2`) %>%
                pivot_wider(names_from = `X3:1`, values_from = `X3:2`) %>%
                pivot_wider(names_from = `X4:1`, values_from = `X4:2`) %>%
                pivot_wider(names_from = `X5:1`, values_from = `X5:2`) %>%
                pivot_wider(names_from = `X6:1`, values_from = `X6:2:1`) %>%
                pivot_wider(names_from = `X7:1`, values_from = `X7:2:1`) %>%
                pivot_wider(names_from = `X8:1`, values_from = `X8:2:1`) %>%
                #time column
                dplyr::mutate(datetime=start.datetime,.before=microphone,
                              location=as.factor(location)) %>%
                `colnames<-`(gsub(" ","",colnames(.))) %>%
                #clean_names() %>% #clean up col names
                #change to numeric
                dplyr::mutate(across(!c(datetime,location,microphone),
                                     as.numeric)) %>%
                #convert microphone names
                dplyr::mutate(microphone=case_when(
                        grepl("3128104_L",microphone) ~ as.factor("bin_L"),
                        grepl("3128104_R",microphone) ~ as.factor("bin_R"),
                        grepl("445966",microphone) ~ as.factor("146AEIn"),
                        grepl("445974",microphone) ~ as.factor("146AEOut"),
                        grepl("3149650_L",microphone) ~ as.factor("bin_L"),
                        grepl("3149650_R",microphone) ~ as.factor("bin_R"),
                        grepl("414769",microphone) ~ as.factor("146AE")
                        
                )) %>%
                pivot_longer(cols = !c(datetime,location,microphone),
                             values_to = "score",
                             names_to = "parameter")
        
        return(df.singleValue)
}

oldartemisSingleValExtract<-function(filename,location,origin_datetime,
                                  acousticParam){
        df.singleValue<-read_delim(file = filename,
                                   delim = ";",
                                   skip = 7,
                                   n_max = 4,
                                   col_names = FALSE,
                                   show_col_types = FALSE) %>%
                
                #separate single value data in first col with "," delimiter
                tidyr::separate_wider_delim(X1, ",",
                                            names=c("microphone",NA,NA,
                                                    acousticParam)) %>%
                #separate numeric data
                tidyr::separate_wider_delim(acousticParam:X8,"=",
                                            names_sep = ":") %>%
                #replace tuHMS units if exists
                dplyr::mutate(!! names(.)[3] := gsub(" tuHMS","",.[[3]])) %>%
                #separate microphone data
                tidyr::separate_wider_delim(microphone,":",
                                            names=c(NA,"microphone")) %>%
                #separate SPL data
                tidyr::separate_wider_delim(c(`X6:2`,`X7:2`,`X8:2`)," ",
                                            names_sep=":") %>%
                #remove redundant cols
                dplyr::select(!c(`X6:2:2`,`X7:2:2`,`X7:2:3`,`X7:2:4`,
                                 `X8:2:2`,`X8:2:3`,`X8:2:4`)) %>%
                #tidy columns with names
                pivot_wider(names_from = paste0(acousticParam,":1"), 
                            values_from = paste0(acousticParam,":2")) %>%
                pivot_wider(names_from = `X2:1`, values_from = `X2:2`) %>%
                pivot_wider(names_from = `X3:1`, values_from = `X3:2`) %>%
                pivot_wider(names_from = `X4:1`, values_from = `X4:2`) %>%
                pivot_wider(names_from = `X5:1`, values_from = `X5:2`) %>%
                pivot_wider(names_from = `X6:1`, values_from = `X6:2:1`) %>%
                pivot_wider(names_from = `X7:1`, values_from = `X7:2:1`) %>%
                pivot_wider(names_from = `X8:1`, values_from = `X8:2:1`) %>%
                #time column
                dplyr::mutate(datetime=origin_datetime,.before=microphone,
                              location=as.factor(location)) %>%
                `colnames<-`(gsub(" ","",colnames(.))) %>%
                #clean_names() %>% #clean up col names
                #change to numeric
                dplyr::mutate(across(!c(datetime,location,microphone),
                                     as.numeric)) %>%
                #convert microphone names
                dplyr::mutate(microphone=case_when(
                        grepl("3128104_L",microphone) ~ as.factor("bin_L"),
                        grepl("3128104_R",microphone) ~ as.factor("bin_R"),
                        grepl("445966",microphone) ~ as.factor("146AEIn"),
                        grepl("445974",microphone) ~ as.factor("146AEOut")
                )) %>%
                pivot_longer(cols = !c(datetime,location,microphone),
                             values_to = "score",
                             names_to = "parameter")
        
        return(df.singleValue)
}

stats_stars <- function(p_values) {
        
        # Use sapply to apply the conditions to each element of the vector
        result <- sapply(p_values, function(p_value) {
                if (p_value < 0.0001) {
                        return("****")
                } else if (p_value < 0.001) {
                        return("***")
                } else if (p_value < 0.01) {
                        return("**")
                } else if (p_value < 0.05) {
                        return("*")
                } else {
                        return("NS")
                }
        })
        
        return(result)
}

effect_sizes <- function(effsizes) {
        
        # Use sapply to apply the conditions to each element of the vector
        result <- sapply(effsizes, function(effsize) {
                if (effsize >= 0.14) {
                        return("L")
                } else if (effsize >= 0.06) {
                        return("M")
                } else if (effsize >= 0.01) {
                        return("S")
                } else if (effsize < 0.01) {
                        return("NA")
                }
        })
        
        return(result)
}
