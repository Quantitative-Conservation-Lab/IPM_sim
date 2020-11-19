# load libraries ####

library(here)
library(dplyr)
library(readxl)
library(lubridate)
library(hms)
library(stringr)
library(tidyr)

# load data ####

# site codes
sites <- read_excel(here("scripts", "dataSims", "IPM sim spreadsheet.xlsx", range = "E2:K5"))

# create a field for file names

# source

source(here("scripts", "dataSims", ""))

# this should load the function we need

# save data objects 