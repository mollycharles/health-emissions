library(tidyverse)
library(plyr)
library(dplyr)

options(scipen = 999)

# Read in data -------------------------------------------------------------------------------------

# Global burden of disease
GBD <- readr::read_csv("./data/IHME-GBD_2017_DATA-3ca98db6-1.csv") 

# CEDS
CEDS_BC <- readr::read_csv("./data/CEDS-master_branch-Fin_Em(2.4.20)/BC_total_CEDS_emissions.csv")
CEDS_CH4 <- readr::read_csv("./data/CEDS-master_branch-Fin_Em(2.4.20)/CH4_total_CEDS_emissions.csv")
CEDS_CO <- readr::read_csv("./data/CEDS-master_branch-Fin_Em(2.4.20)/CO_total_CEDS_emissions.csv")
CEDS_CO2 <- readr::read_csv("./data/CEDS-master_branch-Fin_Em(2.4.20)/CO2_total_CEDS_emissions.csv")
CEDS_NH3 <- readr::read_csv("./data/CEDS-master_branch-Fin_Em(2.4.20)/NH3_total_CEDS_emissions.csv")
CEDS_NMVOC <- readr::read_csv("./data/CEDS-master_branch-Fin_Em(2.4.20)/NMVOC_total_CEDS_emissions.csv")
CEDS_NOx <- readr::read_csv("./data/CEDS-master_branch-Fin_Em(2.4.20)/NOx_total_CEDS_emissions.csv")
CEDS_OC <- readr::read_csv("./data/CEDS-master_branch-Fin_Em(2.4.20)/OC_total_CEDS_emissions.csv")
CEDS_SO2 <- readr::read_csv("./data/CEDS-master_branch-Fin_Em(2.4.20)/SO2_total_CEDS_emissions.csv")

CEDS_pop <- readr::read_csv("./data/A.UN_pop_master.csv")

# Open burning (GFED)
GFED_BC <- bind_rows(readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/BC_LUC_AGRI.csv"),
                     readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/BC_LUC_BORF.csv"),
                     readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/BC_LUC_DEFO.csv"),
                     readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/BC_LUC_PEAT.csv"),
                     readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/BC_LUC_PEAT.csv"),
                     readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/BC_LUC_TEMF.csv"))

GFED_CH4 <- bind_rows(readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/CH4_LUC_AGRI.csv"),
                      readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/CH4_LUC_BORF.csv"),
                      readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/CH4_LUC_DEFO.csv"),
                      readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/CH4_LUC_PEAT.csv"),
                      readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/CH4_LUC_SAVA.csv"),
                      readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/CH4_LUC_TEMF.csv"))

GFED_CO <- bind_rows(readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/CO_LUC_AGRI.csv"),
                     readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/CO_LUC_BORF.csv"),
                     readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/CO_LUC_DEFO.csv"),
                     readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/CO_LUC_PEAT.csv"),
                     readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/CO_LUC_SAVA.csv"),
                     readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/CO_LUC_TEMF.csv"))

GFED_NH3 <- bind_rows(readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/NH3_LUC_AGRI.csv"),
                      readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/NH3_LUC_BORF.csv"),
                      readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/NH3_LUC_DEFO.csv"),
                      readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/NH3_LUC_PEAT.csv"),
                      readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/NH3_LUC_SAVA.csv"),
                      readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/NH3_LUC_TEMF.csv"))

GFED_NMVOC <- bind_rows(readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/NMVOC-bulk_LUC_AGRI.csv"),
                        readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/NMVOC-bulk_LUC_BORF.csv"),
                        readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/NMVOC-bulk_LUC_DEFO.csv"),
                        readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/NMVOC-bulk_LUC_PEAT.csv"),
                        readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/NMVOC-bulk_LUC_SAVA.csv"),
                        readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/NMVOC-bulk_LUC_TEMF.csv"))

GFED_NOx <- bind_rows(readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/NOx_LUC_AGRI.csv"),
                      readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/NOx_LUC_BORF.csv"),
                      readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/NOx_LUC_DEFO.csv"),
                      readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/NOx_LUC_PEAT.csv"),
                      readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/NOx_LUC_SAVA.csv"),
                      readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/NOx_LUC_TEMF.csv"))

GFED_OC <- bind_rows(readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/OC_LUC_AGRI.csv"),
                     readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/OC_LUC_BORF.csv"),
                     readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/OC_LUC_DEFO.csv"),
                     readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/OC_LUC_PEAT.csv"),
                     readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/OC_LUC_SAVA.csv"),
                     readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/OC_LUC_TEMF.csv"))

GFED_SO2 <- bind_rows(readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/SO2_LUC_AGRI.csv"),
                      readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/SO2_LUC_BORF.csv"),
                      readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/SO2_LUC_DEFO.csv"),
                      readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/SO2_LUC_PEAT.csv"),
                      readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/SO2_LUC_SAVA.csv"),
                      readr::read_csv("./data/1750-2015_v1.2_Emissions_bulk_em/SO2_LUC_TEMF.csv"))



# GBD's country groupings
sdi_groups <- readr::read_csv("./data/IHME_GBD_2017_SDI_2017_QUINTILES_Y2018M11D08.CSV") 
cntry_groups <- readxl::read_xlsx("./data/IHME_GBD_2017_CODEBOOK/IHME_GBD_2017_ALL_LOCATIONS_HIERARCHIES_Y2018M11D18.XLSX")

WB_high_income <- cntry_groups %>% filter(parent_id == 44575) %>%
  select(c(location_id)) %>% dplyr::mutate(WB_inc_group = "high")

WB_low_income <- cntry_groups %>% filter(parent_id == 44578) %>%
  select(c(location_id)) %>% dplyr::mutate(WB_inc_group = "low")

WB_lower_middle_income <- cntry_groups %>% filter(parent_id == 44577) %>%
  select(c(location_id)) %>% dplyr::mutate(WB_inc_group = "lower middle")

WB_upper_middle_income <- cntry_groups %>% filter(parent_id == 44576) %>%
  select(c(location_id)) %>% dplyr::mutate(WB_inc_group = "upper middle")

wb_inc_groups <- bind_rows(WB_high_income, WB_low_income, WB_lower_middle_income, WB_upper_middle_income)


G20 <- cntry_groups %>% filter(parent_id == 44586) %>%
  select(c(location_id, location_name)) 

OECD <- cntry_groups %>% filter(parent_id == 44584) %>%
  select(c(location_id, location_name)) 


# match country name, iso, ID, categories
cntry_mapping <- readr::read_csv("./data/cntry_mapping.csv") %>%
  left_join(sdi_groups, by = c("location_id" = "Location ID")) %>%
  select(-c(`2017 SDI Index Value`, `Location Name`)) %>%
  left_join(wb_inc_groups, by = "location_id") 


# CO2 emissions by country (from Global Carbon atlas; using 2016 data as 2017 and 2018 are estimates)
# get top co2 emitting countries. Add top 10 from each SDI group if not included, so all are represented
n_countries <- 50 
n_countries_SDI_group <- 10 

natl_co2_em_top <- readxl::read_xlsx("./data/natl_co2_em_2016_GC.xlsx", sheet="Data", skip = 1) %>%
  gather(-`...1`, key="country", value = "MTCO2") %>%
  left_join(cntry_mapping, by = c("country" = "Country_GCA")) %>%
  select(c(country, MTCO2, iso, location_id, `SDI Quintile`)) %>%
  na.omit() %>%
  arrange(desc(MTCO2)) %>% 
  top_n(n_countries, wt = MTCO2)

natl_co2_em_SDI_groups <- readxl::read_xlsx("./data/natl_co2_em_2016_GC.xlsx", sheet="Data", skip = 1) %>%
  gather(-`...1`, key="country", value = "MTCO2") %>%
  left_join(cntry_mapping, by = c("country" = "Country_GCA")) %>%
  select(c(country, MTCO2, iso, location_id, `SDI Quintile`)) %>%
  na.omit() %>%
  group_by(`SDI Quintile`) %>%
  arrange(desc(MTCO2)) %>% 
  top_n(n_countries_SDI_group, wt = MTCO2)

natl_co2_em <- bind_rows(natl_co2_em_top, natl_co2_em_SDI_groups)

countries_iso <- unique(natl_co2_em$iso)
countries_id <- unique(natl_co2_em$location_id)

# --------------------------------------------------------------------------------------------------

CEDS_years <- paste0("X", 1750:2018)
CEDS_years_noX <- gsub(".*X", "", CEDS_years)

GBD_years <- paste0("X", 1990:2017)
GBD_years_noX <- gsub(".*X", "", GBD_years)

GFED_years <- paste0("X", 1750:2015)
GFED_years_noX <- gsub(".*X", "", GFED_years)


# Shares of GBD causes of death/disease  ----------------------------------------------------------------

GBD_shares <- GBD %>%
  # remove unecessary variables (same across all observations)
  select(-c(sex_id, sex_name, age_id, age_name, cause_id, cause_name, metric_id, metric_name)) %>%
  # remove upper and lower for now 
  select(-c(upper, lower, rei_id)) %>%
  spread(rei_name, val) %>%
  dplyr::mutate(`Air pollution / Environmental/occupational risks` = `Air pollution` / `Environmental/occupational risks`,
                `Air pollution / All risk factors` = `Air pollution` / `All risk factors`,
                `Environ-Occup / All risk factors` = `Environmental/occupational risks` / `All risk factors`,
                `Ambient PM pollution / All risk factors` = `Ambient particulate matter pollution` / `All risk factors`) %>%
  dplyr::mutate(measure_name = if_else(measure_id == 2, "DALYs", measure_name)) %>%
  gather(-c(measure_id, measure_name, location_id, location_name, year), key = "rei_name", value = "val")

GBD_shares_grouped <- GBD %>%
  # remove unecessary variables (same across all observations)
  select(-c(sex_id, sex_name, age_id, age_name, cause_id, cause_name, metric_id, metric_name)) %>%
  # remove upper and lower for now 
  select(-c(upper, lower, rei_id)) %>%
  left_join(cntry_mapping, by = c("location_id")) %>% 
  na.omit() %>%
  group_by_at(vars(measure_id, measure_name, year, rei_name, `SDI Quintile`)) %>%
  dplyr::summarise(val = sum(val)) %>%
  ungroup() %>%
  spread(rei_name, val) %>%
  dplyr::mutate(`Air pollution / Environmental/occupational risks` = `Air pollution` / `Environmental/occupational risks`,
                `Air pollution / All risk factors` = `Air pollution` / `All risk factors`,
                `Environ-Occup / All risk factors` = `Environmental/occupational risks` / `All risk factors`,
                `Ambient PM pollution / All risk factors` = `Ambient particulate matter pollution` / `All risk factors`) %>%
  dplyr::mutate(measure_name = if_else(measure_id == 2, "DALYs", measure_name)) %>%
  gather(-c(measure_id, measure_name, `SDI Quintile`, year), key = "rei_name", value = "val")


# Process CEDS and open burning data - FFI emissions and share of total ----------------------------------------------

total_em <- function(CEDS_in, GFED_in) {
  df_out <- bind_rows(CEDS_in, GFED_in) %>%
    group_by(iso) %>%
    summarise_at(c(CEDS_years), sum, na.rm = TRUE) %>%
    select(iso, CEDS_years)
}

FFI_em <- function(CEDS_in) {
  df_out <- CEDS_in %>%
    
    # Filter out non-fossil/industrial
    filter(!fuel == "biomass", !sector == "5C_Waste-incineration") %>%
    
    # sum by country, year
    group_by(iso) %>%
    summarise_at(c(CEDS_years), sum, na.rm = TRUE) %>%
    
    # filter for years of interest
    select(iso, CEDS_years)
}

nonFFI_em <- function(CEDS_in, GFED_in) {
  df_out <- CEDS_in %>%
    
    # Filter for non-fossil/industrial
    filter(fuel == "biomass" | sector == "5C_Waste-incineration") %>%
    
    # need to add CEDS non FFI and open burning to get non-fossil
    bind_rows(GFED_in) %>%
    
    # sum by country, year
    group_by(iso) %>%
    summarise_at(c(CEDS_years), sum, na.rm = TRUE) %>%   
    
    # filter for years of interest
    select(iso, CEDS_years)
}


FFI_fraction <- function(CEDS_in, GFED_in) {
  
  FFI_em <- FFI_em(CEDS_in) %>% gather(CEDS_years, key = "year", value = "FFI_em")
  
  nonFFI_em <- nonFFI_em(CEDS_in, GFED_in)%>%  gather(CEDS_years, key = "year", value = "nonFFI_em")
  
  total_em <- total_em(CEDS_in, GFED_in) %>%
    gather(CEDS_years, key = "year", value = "total_em")
  
  nonFFI_fraction <- total_em %>%
    left_join(FFI_em, by = c("iso", "year")) %>%
    left_join(nonFFI_em, by = c("iso", "year")) %>%
    dplyr::mutate(FFI_fraction = FFI_em / total_em) %>%
    
    dplyr::mutate(year = as.integer(gsub(".*X","", year))) %>% 
    left_join(cntry_mapping, by = c("iso")) %>% 
    select(-c(Country_GBD, Country_GCA, Country_CEDS))
  
}


# Grouped by SDI quintile


total_em_grouped <- function(CEDS_in, GFED_in) {
    df <- bind_rows(CEDS_in, GFED_in) %>%
      left_join(cntry_mapping, by = "iso") %>%
      group_by(`SDI Quintile`) %>%
      summarise_at(c(CEDS_years), sum, na.rm = TRUE) %>%
      select(`SDI Quintile`, GBD_years) %>%
      na.omit()
}


FFI_em_grouped <- function(CEDS_in) {
  df_out <- CEDS_in %>%
    
    # Filter out non-fossil/industrial
    filter(!fuel == "biomass", !sector == "5C_Waste-incineration") %>%
    
    left_join(cntry_mapping, by = "iso") %>%
    
    # sum by country, year
    group_by(`SDI Quintile`) %>%
    summarise_at(c(CEDS_years), sum, na.rm = TRUE) %>%
    
    # filter for years of interest
    select(`SDI Quintile`, GBD_years) %>%
    na.omit()
}

nonFFI_em_grouped <- function(CEDS_in, GFED_in) {
  df_out <- CEDS_in %>%
    
    # Filter for non-fossil/industrial
    filter(fuel == "biomass" | sector == "5C_Waste-incineration") %>%
    
    # Add CEDS non FFI and open burning to get non-fossil
    bind_rows(GFED_in) %>%
    
    left_join(cntry_mapping, by = "iso") %>%
    
    # sum by country, year
    group_by(`SDI Quintile`) %>%
    summarise_at(c(CEDS_years), sum, na.rm = TRUE) %>%   
    
    # filter for years of interest
    select(`SDI Quintile`, GBD_years) %>%
    na.omit()
}


FFI_fraction_grouped <- function(CEDS_in, GFED_in) {
  
  FFI_em_grouped <- FFI_em_grouped(CEDS_in) 
  
  nonFFI_em_grouped <- nonFFI_em_grouped(CEDS_in, GFED_in) 
  
  total_em_grouped <- total_em_grouped(CEDS_in, GFED_in) %>%
    gather(GBD_years, key = "year", value = "total_em")
  
  FFI_em_grouped <- FFI_em_grouped %>% gather(GBD_years, key = "year", value = "FFI_em")
  nonFFI_em_grouped <- nonFFI_em_grouped %>%  gather(GBD_years, key = "year", value = "nonFFI_em")
  
  FFI_fraction <- total_em_grouped %>%
    left_join(FFI_em_grouped, by = c("SDI Quintile", "year")) %>%
    left_join(nonFFI_em_grouped, by = c("SDI Quintile", "year")) %>%
    dplyr::mutate(FFI_fraction = FFI_em / total_em) %>%
    
    dplyr::mutate(year = as.integer(gsub(".*X","", year))) 
}


# Get composite PM emissions -----------------------------------------------------------------------------------------------------------

# Conversion factors
PM_CF_BC <-	1
PM_CF_NH3	<- 0.5
PM_CF_NO2 <-	0.61
PM_CF_SO2 <-	0.83
PM_CF_OC_biomass <-	1.8
PM_CF_OC_fossil	<- 1.3


PM_emissions <- function(emissions_df, PM_conversion_factor, em_type) {
  df_out <- emissions_df %>% 
    gather(CEDS_years, key = "year", value = "value_em") %>%
    dplyr::mutate(year = as.integer(gsub(".*X","", year))) %>%
    dplyr::mutate(value_pm = value_em * PM_conversion_factor) %>%
    dplyr::mutate(em_type = em_type)
}

# Multiply by conversion factor by species and emissions type

PM_BC <- bind_rows(PM_emissions(FFI_em(CEDS_BC), PM_CF_BC, "FFI"),
                   PM_emissions(nonFFI_em(CEDS_BC, GFED_BC), PM_CF_BC, "non_FFI")) %>%
  dplyr::mutate(species = "BC") %>%
  select(-value_em) %>%
  spread(em_type, value_pm)

PM_NH3 <- bind_rows(PM_emissions(FFI_em(CEDS_NH3), PM_CF_NH3, "FFI"),
                    PM_emissions(nonFFI_em(CEDS_NH3, GFED_NH3), PM_CF_NH3, "non_FFI")) %>%
  dplyr::mutate(species = "NH3") %>%
  select(-value_em) %>%
  spread(em_type, value_pm)

PM_NOx <- bind_rows(PM_emissions(FFI_em(CEDS_NOx), PM_CF_NO2, "FFI"),
                    PM_emissions(nonFFI_em(CEDS_NOx, GFED_NOx), PM_CF_NO2, "non_FFI")) %>%
  dplyr::mutate(species = "NOx") %>%
  select(-value_em) %>%
  spread(em_type, value_pm)

PM_SO2 <- bind_rows(PM_emissions(FFI_em(CEDS_SO2), PM_CF_SO2, "FFI"),
                    PM_emissions(nonFFI_em(CEDS_SO2, GFED_SO2), PM_CF_SO2, "non_FFI")) %>%
  dplyr::mutate(species = "SO2") %>%
  select(-value_em) %>%
  spread(em_type, value_pm)

PM_OC <- bind_rows(PM_emissions(FFI_em(CEDS_OC), PM_CF_OC_fossil, "FFI"),
                   PM_emissions(nonFFI_em(CEDS_OC, GFED_OC), PM_CF_OC_biomass, "non_FFI")) %>%
  dplyr::mutate(species = "OC") %>%
  select(-value_em) %>%
  spread(em_type, value_pm)

PM_composition <- bind_rows(PM_BC, PM_NH3, PM_NOx, PM_SO2, PM_OC)  %>%
  dplyr::mutate(total_pm = FFI + non_FFI) %>% 
  left_join(cntry_mapping, by = "iso") %>%
  select(c(iso, location_id, year, species, FFI, non_FFI, total_pm, `SDI Quintile`)) %>%
  na.omit()
  
# Add to get composite PM
composite_PM <- bind_rows(PM_BC, PM_NH3, PM_NOx, PM_SO2, PM_OC) %>%
  group_by_at(vars(iso, year)) %>%
  dplyr::summarise(FFI = sum(FFI, na.rm=TRUE),
                   non_FFI = sum(non_FFI, na.rm = TRUE)) %>%
  dplyr::mutate(total_pm = FFI + non_FFI) %>%
  left_join(cntry_mapping, by = "iso") %>%
  select(c(iso, location_id, year, FFI, non_FFI, total_pm, `SDI Quintile`))

composite_PM_FFI_fraction <- composite_PM %>%
  dplyr::mutate(FFI_fraction = FFI / total_pm,
                non_FFI_fraction = non_FFI / total_pm) %>%
  select(c(iso,  location_id, year, FFI_fraction, non_FFI_fraction, `SDI Quintile`))


composite_PM_FFI_fraction_grouped <- composite_PM %>%
  group_by_at(vars(year, `SDI Quintile`)) %>%
  dplyr::summarise(FFI = sum(FFI, na.rm=TRUE),
                   non_FFI = sum(non_FFI, na.rm = TRUE),
                   total_pm = sum(total_pm, na.rm = TRUE)) %>%
  dplyr::mutate(FFI_fraction = FFI / total_pm,
                non_FFI_fraction = non_FFI / total_pm) %>%
  na.omit()

# Check country shares of PM emissions by SDI category
cntry_PM_shares <- composite_PM %>%
  dplyr::rename(cntry_pm = total_pm) %>%
  left_join(composite_PM_FFI_fraction_grouped, by = c("year", "SDI Quintile")) %>%
  select(iso, location_id, year, cntry_pm, `SDI Quintile`, total_pm) %>%
  na.omit() %>%
  dplyr::mutate(cntry_share_PM = cntry_pm / total_pm) 

# Combine with global burden of disease data
GBD_composite_PM <- GBD_shares %>% 
  left_join(composite_PM_FFI_fraction, by = c("location_id", "year")) %>% 
  na.omit() %>%
  left_join(select(cntry_PM_shares, -c(total_pm)), by = c("iso", "location_id", "year", "SDI Quintile"))

GBD_composite_PM_grouped <- GBD_shares_grouped %>% left_join(composite_PM_FFI_fraction_grouped, by = c("SDI Quintile", "year")) 



          # Get composite PM using open burning time average ---------------------------------------

                GFED_timeAvg <- function(GFED_in) { 
                  df_out <- GFED_in %>%
                    select(iso, sector, GFED_years) %>%
                    gather(GFED_years, key="year", value="value") %>%
                    group_by(iso) %>%
                    dplyr::summarise(value_timeAvg = mean(value, na.rm=TRUE))
                }

                
                nonFFI_em_BioBAvg <- function(CEDS_in, GFED_in_timeAvg) {
                  df_out <- CEDS_in %>%
                    filter(fuel == "biomass" | sector == "5C_Waste-incineration") %>%
                    group_by(iso) %>%
                    summarise_at(c(CEDS_years), sum, na.rm = TRUE) %>% 
                    select(iso, CEDS_years) %>%
                    gather(GBD_years, key="year", value = "value") %>%
                    left_join(GFED_in_timeAvg, by="iso") %>%
                    dplyr::mutate(value = value + value_timeAvg) %>%
                    select(-value_timeAvg) %>%
                    spread(year, value)
                }

          # PM emissions by species
          PM_BC_BioBAvg <- bind_rows(PM_emissions(FFI_em(CEDS_BC), PM_CF_BC, "FFI"),
                                    PM_emissions(nonFFI_em_BioBAvg(CEDS_BC, GFED_timeAvg(GFED_BC)), PM_CF_BC, "non_FFI")) %>%
            dplyr::mutate(species = "BC") %>%
            select(-value_em) %>%
            spread(em_type, value_pm)
          
          PM_NH3_BioBAvg <- bind_rows(PM_emissions(FFI_em(CEDS_NH3), PM_CF_NH3, "FFI"),
                                      PM_emissions(nonFFI_em_BioBAvg(CEDS_NH3, GFED_timeAvg(GFED_NH3)), PM_CF_NH3, "non_FFI")) %>%
            dplyr::mutate(species = "NH3") %>%
            select(-value_em) %>%
            spread(em_type, value_pm)
          
          PM_NOx_BioBAvg <- bind_rows(PM_emissions(FFI_em(CEDS_NOx), PM_CF_NO2, "FFI"),
                                      PM_emissions(nonFFI_em_BioBAvg(CEDS_NOx, GFED_timeAvg(GFED_NOx)), PM_CF_NO2, "non_FFI")) %>%
            dplyr::mutate(species = "NOx") %>%
            select(-value_em) %>%
            spread(em_type, value_pm)
          
          PM_SO2_BioBAvg <- bind_rows(PM_emissions(FFI_em(CEDS_SO2), PM_CF_SO2, "FFI"),
                                      PM_emissions(nonFFI_em_BioBAvg(CEDS_SO2, GFED_timeAvg(GFED_SO2)), PM_CF_SO2, "non_FFI")) %>%
            dplyr::mutate(species = "SO2") %>%
            select(-value_em) %>%
            spread(em_type, value_pm)
          
          PM_OC_BioBAvg <- bind_rows(PM_emissions(FFI_em(CEDS_OC), PM_CF_OC_fossil, "FFI"),
                                     PM_emissions(nonFFI_em_BioBAvg(CEDS_OC, GFED_timeAvg(GFED_OC)), PM_CF_OC_biomass, "non_FFI")) %>%
            dplyr::mutate(species = "OC") %>%
            select(-value_em) %>%
            spread(em_type, value_pm)
          
          PM_composition_BioBAvg <- bind_rows(PM_BC_BioBAvg, PM_NH3_BioBAvg, PM_NOx_BioBAvg, PM_SO2_BioBAvg, PM_OC_BioBAvg)  %>%
            dplyr::mutate(total_pm = FFI + non_FFI) %>% 
            left_join(cntry_mapping, by = "iso") %>%
            select(c(iso, location_id, year, species, FFI, non_FFI, total_pm, `SDI Quintile`)) %>%
            na.omit()
          
          # Add to get composite PM
          composite_PM_BioBAvg <- bind_rows(PM_BC_BioBAvg, PM_NH3_BioBAvg, PM_NOx_BioBAvg, PM_SO2_BioBAvg, PM_OC_BioBAvg) %>%
            group_by_at(vars(iso, year)) %>%
            dplyr::summarise(FFI = sum(FFI, na.rm=TRUE),
                             non_FFI = sum(non_FFI, na.rm = TRUE)) %>%
            dplyr::mutate(total_pm = FFI + non_FFI) %>%
            left_join(cntry_mapping, by = "iso") %>%
            select(c(iso, location_id, year, FFI, non_FFI, total_pm, `SDI Quintile`))
          
          composite_PM_FFI_fraction_BioBAvg <- composite_PM_BioBAvg %>%
            dplyr::mutate(FFI_fraction = FFI / total_pm,
                          non_FFI_fraction = non_FFI / total_pm) %>%
            select(c(iso,  location_id, year, FFI_fraction, non_FFI_fraction, total_pm, `SDI Quintile`))
          
          
          composite_PM_FFI_fraction_grouped_BioBAvg <- composite_PM_BioBAvg %>%
            group_by_at(vars(year, `SDI Quintile`)) %>%
            dplyr::summarise(FFI = sum(FFI, na.rm=TRUE),
                             non_FFI = sum(non_FFI, na.rm = TRUE),
                             total_pm = sum(total_pm, na.rm = TRUE)) %>%
            dplyr::mutate(FFI_fraction = FFI / total_pm,
                          non_FFI_fraction = non_FFI / total_pm) %>%
            na.omit()
          
          
          # Combine with global burden of disease data
          GBD_composite_PM_BioBAvg <- GBD_shares %>% 
            left_join(composite_PM_FFI_fraction_BioBAvg, by = c("location_id", "year")) %>% 
            na.omit() 
          
          GBD_composite_PM_grouped_BioBAvg <- GBD_shares_grouped %>% left_join(composite_PM_FFI_fraction_grouped_BioBAvg, by = c("SDI Quintile", "year")) 
          
# get composite PM by emissions species, with open burning (average) emissions separated
          
PM_species <- c("BC", "NH3", "NO2", "SO2", "OC")
          
          
PM_CEDS_BC <- CEDS_BC %>%
  group_by(iso) %>%
  summarise_at(c(CEDS_years), sum, na.rm = TRUE) %>%
  gather(CEDS_years, key = "year", value = "value_em") %>%
  dplyr::mutate(year = as.integer(gsub(".*X","", year))) %>%
  dplyr::mutate(value_pm = value_em * PM_CF_BC) %>%
  dplyr::mutate(em_type = "BC")
          
PM_CEDS_NH3 <- CEDS_NH3 %>%
  group_by(iso) %>%
  summarise_at(c(CEDS_years), sum, na.rm = TRUE) %>%
  gather(CEDS_years, key = "year", value = "value_em") %>%
  dplyr::mutate(year = as.integer(gsub(".*X","", year))) %>%
  dplyr::mutate(value_pm = value_em * PM_CF_NH3) %>%
  dplyr::mutate(em_type = "NH3")  

PM_CEDS_NO2 <- CEDS_NOx %>%
  group_by(iso) %>%
  summarise_at(c(CEDS_years), sum, na.rm = TRUE) %>%
  gather(CEDS_years, key = "year", value = "value_em") %>%
  dplyr::mutate(year = as.integer(gsub(".*X","", year))) %>%
  dplyr::mutate(value_pm = value_em * PM_CF_NO2) %>%
  dplyr::mutate(em_type = "NO2") 

PM_CEDS_SO2 <- CEDS_SO2 %>%
  group_by(iso) %>%
  summarise_at(c(CEDS_years), sum, na.rm = TRUE) %>%
  gather(CEDS_years, key = "year", value = "value_em") %>%
  dplyr::mutate(year = as.integer(gsub(".*X","", year))) %>%
  dplyr::mutate(value_pm = value_em * PM_CF_SO2) %>%
  dplyr::mutate(em_type = "SO2") 

PM_CEDS_OC_fossil <- CEDS_OC %>%
  filter(!fuel == "biomass") %>%
  group_by(iso) %>%
  summarise_at(c(CEDS_years), sum, na.rm = TRUE) %>%
  gather(CEDS_years, key = "year", value = "value_em") %>%
  dplyr::mutate(year = as.integer(gsub(".*X","", year))) %>%
  dplyr::mutate(value_pm = value_em * PM_CF_OC_fossil) %>%
  dplyr::mutate(em_type = "OC") 

PM_CEDS_OC_bio <- CEDS_OC %>%
  filter(fuel == "biomass") %>%
  group_by(iso) %>%
  summarise_at(c(CEDS_years), sum, na.rm = TRUE) %>%
  gather(CEDS_years, key = "year", value = "value_em") %>%
  dplyr::mutate(year = as.integer(gsub(".*X","", year))) %>%
  dplyr::mutate(value_pm = value_em * PM_CF_OC_biomass) %>%
  dplyr::mutate(em_type = "OC") 

PM_CEDS_total <- bind_rows(PM_CEDS_BC, PM_CEDS_NH3, PM_CEDS_NO2, PM_CEDS_SO2, PM_CEDS_OC_fossil, PM_CEDS_OC_bio) %>%
  group_by_at(vars(iso, year, em_type)) %>%
  dplyr::summarise(value_pm = sum(value_pm)) %>%
  spread(em_type, value_pm)
          
open_burn_avg_PM_em <- bind_rows(GFED_timeAvg(GFED_BC) %>% dplyr::mutate(species = "BC", EF = PM_CF_BC),
                                 GFED_timeAvg(GFED_NH3) %>% dplyr::mutate(species = "NH3", EF = PM_CF_NH3),
                                 GFED_timeAvg(GFED_NOx) %>% dplyr::mutate(species = "NO2", EF = PM_CF_NO2),
                                 GFED_timeAvg(GFED_SO2) %>% dplyr::mutate(species = "SO2", EF = PM_CF_SO2),
                                 GFED_timeAvg(GFED_OC) %>% dplyr::mutate(species = "OC", EF = PM_CF_OC_biomass)) %>%
  dplyr::mutate(PM_em = value_timeAvg * EF) %>%
  group_by(iso) %>%
  dplyr::summarise(open_burn = sum(PM_em, na.rm = TRUE)) 

PM_composition_open_burn_BioBAvg <- PM_CEDS_total %>%
  left_join(open_burn_avg_PM_em, by = c("iso")) %>%
  gather(c(BC, NH3, NO2, OC, SO2, open_burn), key = "emissions_type", value = "value")
          
# Get fossil fraction of emissions by species ------------------------------------------------------------------------------------------

BC_FFI_fraction <- FFI_fraction(CEDS_BC, GFED_BC)
CH4_FFI_fraction <- FFI_fraction(CEDS_CH4, GFED_CH4) 
CO_FFI_fraction <- FFI_fraction(CEDS_CO, GFED_CO) 
NH3_FFI_fraction <- FFI_fraction(CEDS_NH3, GFED_NH3) 
NMVOC_FFI_fraction <- FFI_fraction(CEDS_NMVOC, GFED_NMVOC) 
NOx_FFI_fraction <- FFI_fraction(CEDS_NOx, GFED_NOx) 
OC_FFI_fraction <- FFI_fraction(CEDS_OC, GFED_OC) 
SO2_FFI_fraction <- FFI_fraction(CEDS_SO2, GFED_SO2) 

# Grouped by SDI category

BC_FFI_fraction_grouped <- FFI_fraction_grouped(CEDS_BC, GFED_BC)
CH4_FFI_fraction_grouped <- FFI_fraction_grouped(CEDS_CH4, GFED_CH4)
CO_FFI_fraction_grouped <- FFI_fraction_grouped(CEDS_CO, GFED_CO)
NH3_FFI_fraction_grouped <- FFI_fraction_grouped(CEDS_NH3, GFED_NH3)
NMVOC_FFI_fraction_grouped <- FFI_fraction_grouped(CEDS_NMVOC, GFED_NMVOC)
NOx_FFI_fraction_grouped <- FFI_fraction_grouped(CEDS_NOx, GFED_NOx)
OC_FFI_fraction_grouped <- FFI_fraction_grouped(CEDS_OC, GFED_OC)
SO2_FFI_fraction_grouped <- FFI_fraction_grouped(CEDS_SO2, GFED_SO2)

# Merge GBD data and fossil fraction of emissions, by emissions species ----------------------------------------------------------------

GBD_BC <- GBD_shares %>% left_join(BC_FFI_fraction, by=c("location_id", "year")) %>% na.omit()
GBD_CH4 <- GBD_shares %>% left_join(CH4_FFI_fraction, by=c("location_id", "year")) %>% na.omit()
GBD_CO <- GBD_shares %>% left_join(CO_FFI_fraction, by=c("location_id", "year")) %>% na.omit()
GBD_NH3 <- GBD_shares %>% left_join(NH3_FFI_fraction, by=c("location_id", "year")) %>% na.omit()
GBD_NMVOC <- GBD_shares %>% left_join(NMVOC_FFI_fraction, by=c("location_id", "year")) %>% na.omit()
GBD_NOx <- GBD_shares %>% left_join(NOx_FFI_fraction, by=c("location_id", "year")) %>% na.omit()
GBD_OC <- GBD_shares %>% left_join(OC_FFI_fraction, by=c("location_id", "year")) %>% na.omit()
GBD_SO2 <- GBD_shares %>% left_join(SO2_FFI_fraction, by=c("location_id", "year")) %>% na.omit()

# Grouped by SDI category
GBD_BC_SDI_groups <- GBD_shares_grouped %>% left_join(BC_FFI_fraction_grouped, by = c("SDI Quintile", "year")) 
GBD_CH4_SDI_groups <- GBD_shares_grouped %>% left_join(CH4_FFI_fraction_grouped, by = c("SDI Quintile", "year")) 
GBD_CO_SDI_groups <- GBD_shares_grouped %>% left_join(CO_FFI_fraction_grouped, by = c("SDI Quintile", "year")) 
GBD_NH3_SDI_groups <- GBD_shares_grouped %>% left_join(NH3_FFI_fraction_grouped, by = c("SDI Quintile", "year")) 
GBD_NMVOC_SDI_groups <- GBD_shares_grouped %>% left_join(NMVOC_FFI_fraction_grouped, by = c("SDI Quintile", "year")) 
GBD_NOx_SDI_groups <- GBD_shares_grouped %>% left_join(NOx_FFI_fraction_grouped, by = c("SDI Quintile", "year")) 
GBD_OC_SDI_groups <- GBD_shares_grouped %>% left_join(OC_FFI_fraction_grouped, by = c("SDI Quintile", "year")) 
GBD_SO2_SDI_groups <- GBD_shares_grouped %>% left_join(SO2_FFI_fraction_grouped, by = c("SDI Quintile", "year")) 

  
# Relationship between composite PM emissions and deaths ------------------------------------------------------------

# Merge composite PM emissions, GBD data, and population 
GBD_population_PM <- GBD_composite_PM %>% left_join(CEDS_pop, by = c("iso", "year")) %>%
  dplyr::mutate(pm_pop = cntry_pm * pop)
GBD_population_PM_BioBAvg <- GBD_composite_PM_BioBAvg %>% left_join(CEDS_pop, by = c("iso", "year")) %>%
  dplyr::mutate(pm_pop = total_pm * pop)


# get slopes by country to extrapolate deaths in previous years
# use slope in early years for countries where slope changes

extrap_years <- c(1850:1989)
extrap_isos <- c("chn", "ind", "fra")

actual_years <- c(1990:2017)

fra_PM_deaths_lm <- lm(formula = val ~ pm_pop, 
                       data=filter(GBD_population_PM, iso == "fra", year <= 2002, 
                                   measure_name == "Deaths", rei_name == "Ambient PM pollution / All risk factors"))
ind_PM_deaths_lm <- lm(formula = val ~ pm_pop, 
                       data=filter(GBD_population_PM, iso == "ind", 
                                   measure_name == "Deaths", rei_name == "Ambient PM pollution / All risk factors"))
chn_PM_deaths_lm <- lm(formula = val ~ pm_pop, 
                       data=filter(GBD_population_PM, iso == "chn", year <= 1996, 
                                   measure_name == "Deaths", rei_name == "Ambient PM pollution / All risk factors"))

extrap_PM_deaths_frac <- composite_PM %>%
  left_join(CEDS_pop, by = c("iso", "year")) %>%
  dplyr::mutate(pm_pop = total_pm * pop) %>%
  na.omit() %>%
  left_join(filter(GBD_shares, measure_name == "Deaths", rei_name == "Ambient PM pollution / All risk factors"), by = c("location_id", "year")) %>%
  select(c(iso, location_id, country, year, total_pm, FFI, non_FFI, `SDI Quintile`, rei_name, pop, pm_pop, val)) %>%
  filter(year %in% extrap_years,
         iso %in% extrap_isos)

actual_PM_deaths_frac <- composite_PM %>%
  left_join(CEDS_pop, by = c("iso", "year")) %>%
  dplyr::mutate(pm_pop = total_pm * pop) %>%
  na.omit() %>%
  left_join(filter(GBD_shares, measure_name == "Deaths", rei_name == "Ambient PM pollution / All risk factors"), by = c("location_id", "year")) %>%
  select(c(iso, location_id, country, year, total_pm, FFI, non_FFI, `SDI Quintile`, rei_name, pop, pm_pop, val)) %>%
  filter(year %in% actual_years,
         iso %in% extrap_isos)

# extrapolate PM deaths fraction  for each country

extrap_PM_deaths_chn <- extrap_PM_deaths_frac %>%
  filter(iso == "chn") 
extrap_PM_deaths_chn$val <- predict(chn_PM_deaths_lm, new = extrap_PM_deaths_chn)


extrap_PM_deaths_ind <- extrap_PM_deaths_frac %>%
  filter(iso == "ind") 
extrap_PM_deaths_ind$val <- predict(ind_PM_deaths_lm, new = extrap_PM_deaths_ind)


extrap_PM_deaths_fra <- extrap_PM_deaths_frac %>%
  filter(iso == "fra") 
extrap_PM_deaths_fra$val <- predict(fra_PM_deaths_lm, new = extrap_PM_deaths_fra)


PM_deaths_frac_cntry_allyears_extrap <- bind_rows(actual_PM_deaths_frac,
                                                  extrap_PM_deaths_chn,
                                                  extrap_PM_deaths_ind,
                                                  extrap_PM_deaths_fra)
