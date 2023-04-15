setwd()

library(terra)

### Import ant nest coordinate data (for n = 77 nests of several Formica ant species or hybrids).
nestcoords = read.csv2("input data/FIN_33_ant_nest_coords_roundedCoord.csv")

### Import and extract Precipitation data, and derive variables for modelling.

# Precipitation data source: Daily Precipitation sum (rrday) rasters were downloaded from https://paituli.csc.fi/download.html

# This dataset is part of FMI ClimGrid, which is a gridded daily climatology dataset of Finland at a spatial resolution of 10 km x 10 km. 
# i.e., 10km_daily_precipitation/geotiff.
# The Coordinate Reference System is ETRS89 / ETRS-TM35FIN (EPSG:3067)
# see metadata here: https://etsin.fairdata.fi/dataset/0dfc9c77-e6d6-47c6-af46-f2385e0f14ff

# Here rasters for the years 1985-2020 only were retained, which cover 20 years prior to the earliest ant field samples.
# These were saved in a folder in the working directory called "FMI_daily_precipitation_1985_2020_geotiff"
# within the input data folder.

# Create SpatVector of nest coordinates in the "EPSG:3067" coord. ref. system
coords = vect(nestcoords[, c("ETRS.TM35FIN_X", "ETRS.TM35FIN_Y" )], geom = c("ETRS.TM35FIN_X", "ETRS.TM35FIN_Y"), "EPSG:3067")

#Extract daily Precipitation values at each of 77 nest sites and retain those for the months of Jan-Feb, April-May and Decemember.

#Define threshold running days in normal (not leap) years:
endFeb = 59
startApr = 91
endMay = 151
startDec = 335

winterprec = data.frame(matrix(NA, nrow = nrow(nestcoords), ncol = 0))
springprec = data.frame(matrix(NA, nrow = nrow(nestcoords), ncol = 0))

files = list.files("input data", recursive = TRUE, pattern = ".tif")
include = grep("precipitation", files)
files = files[include]
exclude = grep("tif.aux.xml", files)
files = files[-exclude]

for(i in 1:length(files))
{
  clim = rast(paste("input data/",files[i], sep = "/"))
  crs(clim) = "EPSG:3067"
  dailyclim = data.frame(extract(clim, coords)[,-1])
  if(ncol(dailyclim)==365)
  {
    winterprec = cbind(winterprec, dailyclim[, c(1:endFeb, startDec:ncol(dailyclim))])
    springprec = cbind(springprec, dailyclim[, c(startApr:endMay)])
  }
  if(ncol(dailyclim)==366)
  {
    winterprec = cbind(winterprec, dailyclim[, c(1:(endFeb+1), (startDec+1):ncol(dailyclim))])
    springprec = cbind(springprec, dailyclim[, c((startApr+1):(endMay+1))])
  }
}

# Exclude all 1985 spring values as these fall outside the required 20-year observation window 
# from the earliest ant field sampling summer (2005) 
springprec = springprec[, -c(grep("1985", names(springprec)))]

# Exclude all spring values from 2020 as the most recent ant obs. dates are from summer 2019.
springprec = springprec[, -c(grep("2020", names(springprec)))]

# Exclude Jan-Feb 1985 winter values as these fall outside required 20-year observation window 
# from the earliest ant field sampling summer (2005) 
winterprec = winterprec[, -c(1:endFeb)]

# Exclude all winter values from 2020 as the most recent ant obs. dates are from summer 2019.
winterprec = winterprec[, -c(grep("2020", names(winterprec)))]

# Exclude December values from 2019 as the most recent ant obs. dates are from summer 2019.
winterprec = winterprec[, -c(grep("rrday_2019_335", names(winterprec)):ncol(winterprec))]

# Create "output data" folder in the working directory, if it does not exist already.


### Calculate Spring seasonal precipitation indices for the years 1986-2019 ###

# Specifically, based on the daily precipitation sums in spring each year, calculated values are the seasonal
# mean and 25th and 75th percentiles of daily precipitation.

years = 1986:2019
n_vars = 3

spring_data = data.frame(nestcoords, matrix(NA, nrow = nrow(nestcoords), ncol = n_vars*length(years)))
names(spring_data)[-c(1:ncol(nestcoords))] = paste0(sort(rep(years, n_vars)), c("_mean", "_perc25", "_perc75"))


for(n in 1:length(years))
{
  year = years[n]
  seldata = springprec[, grep(year, names(springprec))]
  
  mean_col = paste0(year, "_mean")
  quantile_cols = paste0(year, c("_perc25", "_perc75"))
  
  spring_data[, mean_col] = apply(seldata, 1, mean)
  spring_data[, quantile_cols] = t(apply(seldata, 1, quantile, probs = c(0.25, 0.75)))
}

write.csv2(spring_data, "output data/SpringPrecipitation_antnests_1986-2019.csv", row.names = F)



### Calculate Winter seasonal precipitation indices for the years 1986-2019 ###

# Specifically, based on the daily precipitation sums in winter (Dec-Feb), calculated values are the seasonal
# mean and 25th and 75th percentiles of daily precipitation. 

# In the script, for each year in turn, winter combines December of the previous year with January + February of the year selected.

years = 1986:2019
n_vars = 3

winter_data = data.frame(nestcoords, matrix(NA, nrow = nrow(nestcoords), ncol = n_vars*length(years)))
names(winter_data)[-c(1:ncol(nestcoords))] = paste0(sort(rep(years, n_vars)), c("_mean", "_perc25", "_perc75"))


for(n in 1:length(years))
{
  year = years[n]
  
  seldata1 = winterprec[, grep(year-1, names(winterprec))]
  selDays1 = as.numeric(sub(paste0("rrday_", year-1, "_"), "", names(seldata1)))
  seldata1 = seldata1[, selDays1 > 300]
  
  seldata2 = winterprec[, grep(year, names(winterprec))]
  selDays2 = as.numeric(sub(paste0("rrday_", year, "_"), "", names(seldata2)))
  seldata2 = seldata2[, selDays2 < 100]
  
  combdata = cbind(seldata1, seldata2)
  
  mean_col = paste0(year, "_mean")
  quantile_cols = paste0(year, c("_perc25", "_perc75"))
  
  winter_data[, mean_col] = apply(combdata, 1, mean)
  winter_data[, quantile_cols] = t(apply(combdata, 1, quantile, probs = c(0.25, 0.75)))
}

write.csv2(winter_data, "output data/WinterPrecipitation_antnests_1986-2019.csv", row.names = F)

