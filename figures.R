source("./data_processing.R")

library(gganimate)
library(RColorBrewer)

em_species <- c("BC", "CH4", "CO", "NH3", "NMVOC", "NOx", "OC", "SO2")

health_indicators <- c("Air pollution / All risk factors", "Environ-Occup / All risk factors")
measures <- c("Deaths", "DALYs")

SDI_groups <- c("Low SDI", "Low-middle SDI", "Middle SDI", "High-middle SDI", "High SDI")

plot_theme <- theme_bw() + theme(text = element_text(size = 20),
                                 plot.title = element_text(size = 24, margin = margin(0, 0, 10, 0)),
                                 strip.text.x = element_text(size = 16),
                                 axis.text = element_text(size = 18))

# Composite PM (open burning average), animated ------------------------------------------------------------------------

GBD_composite_PM_BioBAvg$`SDI Quintile` <- factor(GBD_composite_PM_BioBAvg$`SDI Quintile`, levels = SDI_groups)

for (i in seq_along(health_indicators)) {
  
  health_indicator <- health_indicators[i]
  
  for (j in seq_along(measures)) {

    measure <- measures[j]

    
    plot <-
      ggplot() +
      geom_point(data=filter(GBD_composite_PM_BioBAvg, location_id %in% c(countries_id), measure_name == measure, rei_name == health_indicator),
                 aes(x=(FFI_fraction), y=val, group = location_name, color = `SDI Quintile`, fill = `SDI Quintile`,  shape = `SDI Quintile`), size = 3) +
      geom_point(data=filter(GBD_composite_PM_grouped_BioBAvg, measure_name == measure, rei_name == health_indicator),
                 aes(x=(FFI_fraction), y=val, color = `SDI Quintile`, fill = `SDI Quintile`, shape = `SDI Quintile`), size = 8) +
      transition_time(as.integer(year)) +
      labs(title = "Year: {frame_time}", x="Fossil_Industrial Fraction of Composite PM Emissions", 
           y = paste(measure, "from", health_indicator)) +
      scale_y_continuous(limits = c(0, NA),
                         labels = scales::percent_format(accuracy = 1L)) +
      scale_x_continuous(labels = scales::percent_format(accuracy = 1L)) +
      scale_color_manual(name = "Socio-Demographic Index",
                         values = brewer.pal(5, "Set1")) +   
      scale_fill_manual(name = "Socio-Demographic Index",
                         values = brewer.pal(5, "Set1")) +  
      scale_shape_manual(name = "Socio-Demographic Index",
                         values = c(21, 22, 25, 23, 24)) +
      plot_theme
    animate(plot, nframes = 100, fps = 5, width = 600, height = 400, start_pause = 5, end_pause = 20)
    anim_save(filename=paste0("./figures/", "PM", "_", measure, "_", gsub(" / ", "_", health_indicator), "_anim", "_BioBAvg", ".gif"))
    
  }
}

# Composite PM (open burning average), still with trajectory lines -----------------------------------------------------------------
# this is the main plot

for (i in seq_along(health_indicators)) {
  
  health_indicator <- health_indicators[i]
  health_indicator_label <- if_else(health_indicator == "Air pollution / All risk factors", "Air pollution", "Environ-Occup")
  
  for (j in seq_along(measures)) {

    measure <- measures[j]
    
    plot_year <- 2017
    
    
    plot <-
      ggplot() +
      # individual country points
      geom_point(data=filter(GBD_composite_PM_BioBAvg, location_id %in% c(countries_id), measure_name == measure, rei_name == health_indicator, year == plot_year),
                 aes(x=(FFI_fraction), y=val, group = location_name, color = `SDI Quintile`, fill = `SDI Quintile`, shape = `SDI Quintile`), size = 3) +
      # SDI group averages
      geom_point(data=filter(GBD_composite_PM_grouped_BioBAvg, measure_name == measure, rei_name == health_indicator, year == plot_year),
                 aes(x=(FFI_fraction), y=val, color = `SDI Quintile`, fill = `SDI Quintile`, shape = `SDI Quintile`), size = 8) +
      #small points for trail
      geom_point(data=filter(GBD_composite_PM_grouped_BioBAvg, measure_name == measure, rei_name == health_indicator, year %in% c(1991:plot_year)),
                 aes(x=FFI_fraction, y=val, color = `SDI Quintile`), shape = 13, size = 1) +
      # geom_line(data=filter(GBD_composite_PM_grouped_BioBAvg, measure_name == measure, rei_name == health_indicator, year %in% c(1991:plot_year)),
      #            aes(x=FFI_fraction, y=val, color = `SDI Quintile`), linetype = "dashed", size = 1) +
      # open point for 1990
      geom_point(data=filter(GBD_composite_PM_grouped_BioBAvg, measure_name == measure, rei_name == health_indicator, year==1990),
                 aes(x=(FFI_fraction), y=val, color = `SDI Quintile`, shape = `SDI Quintile`), size = 5) +
      labs(x="Fossil Fraction of Composite PM Emissions", 
           y = paste(measure, "from", health_indicator_label)) +
      scale_y_continuous(limits = c(0, NA),
                         labels = scales::percent_format(accuracy = 1L)) +
      scale_x_reverse(labels = scales::percent_format(accuracy = 1L)) +
      scale_color_manual(name = "Socio-Demographic Index",
                         values = brewer.pal(5, "Set1")) +   
      scale_fill_manual(name = "Socio-Demographic Index",
                         values = brewer.pal(5, "Set1")) +  
      scale_shape_manual(name = "Socio-Demographic Index",
                         values = c(21, 22, 25, 23, 24)) +
      plot_theme+
      theme(axis.text = element_text(size=20))
    ggsave(filename=paste0("./figures/", "PM", "_", measure, "_", gsub(" / ", "_", health_indicator), "_still", "_BioBAvg", ".png"), 
           height = 8.5, width = 11, units = "in")
    
  }
}

# Major countries - composite PM (open burning average), animated ------------------------------------------------------------------------

for (i in seq_along(health_indicators)) {
  
  health_indicator <- health_indicators[i]
  
  for (j in seq_along(measures)) {
    
    measure <- measures[j]
    
    plot_countries <- c("China", "India", "United States", "Indonesia", "Democratic Republic of the Congo")
    
    data_1990 <- GBD_composite_PM_BioBAvg %>% filter(year == 1990) %>% select(-year)
    
    
    plot <-
      ggplot() +
      geom_point(data=filter(GBD_composite_PM_BioBAvg, location_name %in% plot_countries, measure_name == measure, rei_name == health_indicator, year > 1990),
                 aes(x=(FFI_fraction), y=val, color = location_name, fill = location_name,  shape = location_name), size = 8) +
      geom_point(data=filter(data_1990, location_name %in% plot_countries, measure_name == measure, rei_name == health_indicator),
                 aes(x=(FFI_fraction), y=val, color = location_name, shape = location_name), size = 3.5) +
      transition_time(as.integer(year)) +
      #shadow_trail(distance = 1, size = 4) +
      shadow_mark(past = TRUE, future = FALSE, size = 1) +
      labs(title = "Year: {frame_time}", x="Fossil Fraction of Composite PM Emissions", 
           y = paste(measure, "from", health_indicator)) +
      scale_y_continuous(limits = c(0, NA),
                         labels = scales::percent_format(accuracy = 1L)) +
      scale_x_continuous(limits = c(0, NA),
                         labels = scales::percent_format(accuracy = 1L)) +
      scale_color_manual(name = "Country",
                         values = brewer.pal(length(plot_countries), "Dark2")) +
      scale_fill_manual(name = "Country",
                        values = brewer.pal(length(plot_countries), "Dark2")) +
      scale_shape_manual(name = "Country",
                         values = c(21, 22, 25, 23, 24)) +
      plot_theme
    animate(plot, nframes = 100, fps = 5, width = 600, height = 400, start_pause = 5, end_pause = 20)
    anim_save(filename=paste0("./figures/", "selected_countries_", "PM", "_", measure, "_", gsub(" / ", "_", health_indicator), "_anim", "_BioBAvg", ".gif"))
    
  }
}

# Composite PM (yearly open burning data), animated ---------------------------------------------------------

GBD_composite_PM$`SDI Quintile` <- factor(GBD_composite_PM$`SDI Quintile`, levels = SDI_groups)

for (i in seq_along(health_indicators)) {
  
  health_indicator <- health_indicators[i]
  
  for (j in seq_along(measures)) {
    
    measure <- measures[j]
    
    plot_years <- c(1990:2015)
    
    plot <-
      ggplot() +
      geom_point(data=filter(GBD_composite_PM, location_id %in% c(countries_id), year %in% plot_years, measure_name == measure, rei_name == health_indicator),
                 aes(x=FFI_fraction, y=val, group = location_name, color = `SDI Quintile`, fill = `SDI Quintile`, shape = `SDI Quintile`), size = 3) +
      geom_point(data=filter(GBD_composite_PM_grouped, year %in% plot_years, measure_name == measure, rei_name == health_indicator),
                 aes(x=FFI_fraction, y=val, color = `SDI Quintile`, fill = `SDI Quintile`, shape = `SDI Quintile`), size = 8) +
      transition_time(as.integer(year)) +
      labs(title = "Year: {frame_time}", x="Fossil Fraction of Composite PM Emissions", 
           y = paste(measure, "from", health_indicator)) +
      scale_y_continuous(limits = c(0, NA),
                         labels = scales::percent_format(accuracy = 1L)) +
      scale_x_continuous(labels = scales::percent_format(accuracy = 1L)) +
      scale_color_manual(name = "Socio-Demographic Index",
                         values = brewer.pal(5, "Set1")) +   
      scale_fill_manual(name = "Socio-Demographic Index",
                        values = brewer.pal(5, "Set1")) +  
      scale_shape_manual(name = "Socio-Demographic Index",
                         values = c(21, 22, 25, 23, 24)) +
      plot_theme
    animate(plot, nframes = 100, fps = 5, width = 600, height = 400, start_pause = 5, end_pause = 20)
    anim_save(filename=paste0("./figures/", "PM", "_", measure, "_", gsub(" / ", "_", health_indicator), "_anim", ".gif"))
    
  }
}


# Individual emissions species (yearly open burning data), animated --------------------------------------------------------------

for (i in seq_along(em_species)) {
  
# Parameters for figures:
      plot_years <- c(1990:2015)
      
      # select measure (deaths or DALYs)
      measure <- "DALYs"
      
      # select risk variable
      health_indicator <- "Air pollution / All risk factors"

  
  df <- get(paste0("GBD_", em_species[i]))
  df_grouped <- get(paste0("GBD_", em_species[i], "_SDI_groups"))
  
  df$`SDI Quintile` <- factor(df$`SDI Quintile`, levels = c("Low SDI", "Low-middle SDI", "Middle SDI", "High-middle SDI", "High SDI"))
  df$WB_inc_group <- factor(df$WB_inc_group, levels = c("low", "lower middle", "upper middle", "high"))
      
  plot <-
    ggplot() +
    geom_point(data=filter(df, location_id %in% c(countries_id), year %in% plot_years, measure_name == measure, rei_name == health_indicator),
               aes(x=(FFI_fraction), y=val, group = location_name, color = `SDI Quintile`, fill = `SDI Quintile`, shape = `SDI Quintile`), size = 3) +
    geom_point(data=filter(df_grouped, year %in% plot_years, measure_name == measure, rei_name == health_indicator),
               aes(x=(FFI_fraction), y=val, color = `SDI Quintile`, fill = `SDI Quintile`, shape = `SDI Quintile`), size = 8) +
    transition_time(as.integer(year)) +
    labs(title = paste(em_species[i]), subtitle = "Year: {frame_time}", x="Fossil Fraction of Emissions", 
         y = paste(measure, "from", health_indicator)) +
    scale_y_continuous(limits = c(0, NA),
                       labels = scales::percent_format(accuracy = 1L)) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1L)) +
    scale_color_manual(name = "Socio-Demographic Index",
                       values = brewer.pal(5, "Set1")) +   
    scale_fill_manual(name = "Socio-Demographic Index",
                       values = brewer.pal(5, "Set1")) + 
    scale_shape_manual(name = "Socio-Demographic Index",
                       values = c(21, 22, 25, 23, 24)) +
    plot_theme
  animate(plot, nframes = 100, fps = 5, width = 600, height = 400, start_pause = 10, end_pause = 30)
  anim_save(filename=paste0("./figures/", em_species[i], "_", measure, "_", gsub(" / ", "_", health_indicator), "_anim", ".gif"))
  
}


# Individual emissions species (yearly open burning data), still with trajectory lines -------------------------------------------------------

for (i in seq_along(em_species)) {
  
  # Parameters for figures:
      # select measure (deaths or DALYs)
      measure <- "DALYs"
      
      # select risk variable
      health_indicator <- "Air pollution / All risk factors"
      
      # select year
      plot_year <- 2015

  
  df <- get(paste0("GBD_", em_species[i]))
  df_grouped <- get(paste0("GBD_", em_species[i], "_SDI_groups"))
  
  
  df$`SDI Quintile` <- factor(df$`SDI Quintile`, levels = SDI_groups)

  plot <-
    ggplot() +
    geom_point(data=filter(df, location_id %in% c(countries_id), measure_name == measure, rei_name == health_indicator, year == plot_year),
               aes(x=(FFI_fraction), y=val, group = location_name, color = `SDI Quintile`, fill = `SDI Quintile`, shape = `SDI Quintile`), size = 3) +
    geom_point(data=filter(df_grouped, measure_name == measure, rei_name == health_indicator, year == plot_year),
               aes(x=(FFI_fraction), y=val, color = `SDI Quintile`, fill = `SDI Quintile`, shape = `SDI Quintile`), size = 8) +
    geom_point(data=filter(df_grouped, measure_name == measure, rei_name == health_indicator, year %in% c(1990:plot_year)),
               aes(x=(FFI_fraction), y=val, color = `SDI Quintile`, fill = `SDI Quintile`), size = 4, shape = "-") +
    # geom_text(data=filter(df_grouped, measure_name == measure, rei_name == health_indicator, year %in% c(1990:plot_year)), 
    #            mapping = aes(x=(FFI_fraction), y=val, label = year), stat = "identity") +
    labs(title = paste(em_species[i]), subtitle = plot_year, x="Fossil Fraction of Emissions", 
         y = paste(measure, "from", health_indicator)) +
    scale_y_continuous(limits = c(0, NA),
                       labels = scales::percent_format(accuracy = 1L)) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1L)) +
    scale_color_manual(name = "Socio-Demographic Index",
                        values = brewer.pal(5, "Set1")) +   
    scale_fill_manual(name = "Socio-Demographic Index",
                       values = brewer.pal(5, "Set1")) +  
    scale_shape_manual(name = "Socio-Demographic Index",
                       values = c(21, 22, 25, 23, 24)) +
    plot_theme
  ggsave(plot, filename=paste0("./figures/", em_species[i], "_", measure, "_", gsub(" / ", "_", health_indicator), "_still", ".png"), 
         height = 8.5, width = 11, units = "in")
  
}


# Graph relationship between air quality and deaths for selected countries ------------------------------------

GBD_population_PM <- GBD_population_PM %>%
  dplyr::mutate(inflection_point = if_else(iso == "chn" & year %in% c(2002, 2007), TRUE,
                                           if_else(iso == "gbr" & year == 1999, TRUE,
                                                   if_else(iso == "usa" & year %in% c(2002, 2011), TRUE,
                                                           if_else(iso == "deu" & year %in% c(1995), TRUE,
                                                                   if_else(iso == "fra" & year %in% c(2002, 2006, 2011), TRUE, FALSE))))))

GBD_population_PM_BioBAvg <- GBD_population_PM_BioBAvg %>%
  dplyr::mutate(inflection_point = if_else(iso == "chn" & year %in% c(2002, 2007), TRUE,
                                           if_else(iso == "gbr" & year == 1999, TRUE,
                                                   if_else(iso == "usa" & year %in% c(2002, 2011), TRUE,
                                                           if_else(iso == "deu" & year %in% c(1995, 2010), TRUE,
                                                                   if_else(iso == "fra" & year %in% c(2002, 2006, 2011), TRUE, FALSE))))))

#number of deaths from ambient PM pollution
ggplot(data=filter(GBD_population_PM, 
                   iso %in% c("chn", "usa", "ind", "idn", "cod", "fra", "deu", "pol", "gbr"),
                   rei_name == "Ambient particulate matter pollution",
                   measure_name == "Deaths"),
       aes(x = pm_pop/1000000, y = val, color = year)) +
  geom_point(size = 5) +
  geom_text(data=filter(GBD_population_PM, 
                        iso %in% c("chn", "usa", "ind", "idn", "cod", "fra", "deu", "pol", "gbr"),
                        rei_name == "Ambient particulate matter pollution",
                        measure_name == "Deaths",
                        inflection_point == TRUE),
            aes(label = year), size = 3, hjust=1.25, vjust=-0.5, color = "black") +
  facet_wrap(~location_name, scales = "free") +
  labs(x="Composite PM * population (millions)", y = "Deaths from ambient PM pollution") +
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  ggsave("./figures/PM_em_deaths_by_cntry.png", height = 8.5, width = 11, units = "in")

#fraction of all deaths from ambient PM pollution
ggplot(data=filter(GBD_population_PM, iso %in% c("chn", "usa", "ind", "idn", "cod", "fra", "deu", "pol", "gbr"),
                   rei_name == "Ambient PM pollution / All risk factors",
                   measure_name == "Deaths"),
       aes(x = pm_pop/1000000, y = val, color = year)) +
  geom_point(size = 5) +
  geom_text(data=filter(GBD_population_PM, iso %in% c("chn", "usa", "ind", "idn", "cod", "fra", "deu", "pol", "gbr"),
                        rei_name == "Ambient PM pollution / All risk factors",
                        measure_name == "Deaths",
                        inflection_point == TRUE),
            aes(label = year), size = 3, hjust=1.25, vjust=-0.5, color = "black") +
  facet_wrap(~location_name, scales = "free") +
  labs(x="Composite PM * population (millions)", y = "Fraction of deaths from ambient PM pollution / All risk factors") +
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  ggsave("./figures/PM_em_fraction_deaths_by_cntry.png", height = 8.5, width = 11, units = "in")

#number of deaths from ambient PM pollution, open burning average
ggplot(data=filter(GBD_population_PM_BioBAvg, 
                   iso %in% c("chn", "usa", "ind", "idn", "cod", "fra", "deu", "pol", "gbr"),
                   rei_name == "Ambient particulate matter pollution",
                   measure_name == "Deaths"),
       aes(x = pm_pop/1000000, y = val, color = year)) +
  geom_point(size = 5) +
  geom_text(data=filter(GBD_population_PM_BioBAvg, 
                        iso %in% c("chn", "usa", "ind", "idn", "cod", "fra", "deu", "pol", "gbr"),
                        rei_name == "Ambient particulate matter pollution",
                        measure_name == "Deaths",
                        inflection_point == TRUE),
            aes(label = year), size = 3, hjust=1.25, vjust=-0.25, color = "black") +
  facet_wrap(~location_name, scales = "free") +
  labs(x="Composite PM * population (millions)", y = "Deaths from ambient PM pollution") +
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  ggsave("./figures/PM_em_deaths_by_cntry_BioBAvg.png", height = 8.5, width = 11, units = "in")

#fraction of all deaths from ambient PM pollution, open burning average
ggplot(data=filter(GBD_population_PM_BioBAvg, iso %in% c("chn", "usa", "ind", "idn", "cod", "fra", "deu", "pol", "gbr"),
                   rei_name == "Ambient PM pollution / All risk factors",
                   measure_name == "Deaths"),
       aes(x = pm_pop/1000000, y = val, color = year)) +
  geom_point(size = 6) +
  geom_text(data=filter(GBD_population_PM_BioBAvg, iso %in% c("chn", "usa", "ind", "idn", "cod", "fra", "deu", "pol", "gbr"),
                        rei_name == "Ambient PM pollution / All risk factors",
                        measure_name == "Deaths",
                        inflection_point == TRUE),
            aes(label = year), size = 3, hjust=1, vjust=-0.5, color = "black") +
  facet_wrap(~location_name, scales = "free") +
  labs(x="Composite PM * population (millions)", y = "Fraction of deaths from ambient PM pollution / All risk factors") +
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  ggsave("./figures/PM_em_fraction_deaths_by_cntry_BioBAvg.png", height = 8.5, width = 11, units = "in")

# Graph composition of PM over time for countries above ----------------------------------------------------------------

#yearly open burning data
ggplot(data=filter(PM_composition,
                   iso %in% c("chn", "usa", "ind", "idn", "cod", "fra", "deu", "pol", "gbr"),
                   year %in% c(1990:2015)),
       aes(x=year, y=total_pm, fill=species)) +
  geom_area()+
  facet_wrap(~iso, scales = "free", labeller = labeller(iso = c("chn" = "China",
                                                                "usa" = "USA",
                                                                "ind" = "India",
                                                                "idn" = "Indonesia",
                                                                "cod" = "Democratic Republic of the Congo",
                                                                "fra" = "France",
                                                                "deu" = "Germany",
                                                                "pol" = "Poland", 
                                                                "gbr" = "United Kingdom"))) +
  plot_theme +
  labs(title = "Composition of PM emissions", x="Year", y="kt")+
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
        strip.text.x = element_text(size = 12)) +
  ggsave("./figures/PM_composition_by_cntry.png", height = 8.5, width = 11, units = "in")

#average open burning data
ggplot(data=filter(PM_composition_BioBAvg,
                   iso %in% c("chn", "usa", "ind", "idn", "cod", "fra", "deu", "pol", "gbr"),
                   year %in% c(1990:2017)),
       aes(x=year, y=total_pm, fill=species)) +
  geom_area()+
  facet_wrap(~iso, scales = "free", labeller = labeller(iso = c("chn" = "China",
                                                                "usa" = "USA",
                                                                "ind" = "India",
                                                                "idn" = "Indonesia",
                                                                "cod" = "Democratic Republic of the Congo",
                                                                "fra" = "France",
                                                                "deu" = "Germany",
                                                                "pol" = "Poland", 
                                                                "gbr" = "United Kingdom"))) +
  labs(title = "Composition of PM emissions", x="Year", y="kt")+
  plot_theme +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
        strip.text.x = element_text(size = 12)) +
  ggsave("./figures/PM_composition_by_cntry_BioBAvg.png", height = 8.5, width = 11, units = "in")

# separating out composite PM from open burning
ggplot(data=filter(PM_composition_open_burn_BioBAvg,
                   iso %in% c("chn", "usa", "ind", "idn", "cod", "fra", "deu", "pol", "gbr"),
                   year %in% c(1990:2017)),
       aes(x=year, y=value, fill=emissions_type)) +
  geom_area()+
  facet_wrap(~iso, scales = "free", labeller = labeller(iso = c("chn" = "China",
                                                                "usa" = "USA",
                                                                "ind" = "India",
                                                                "idn" = "Indonesia",
                                                                "cod" = "Democratic Republic of the Congo",
                                                                "fra" = "France",
                                                                "deu" = "Germany",
                                                                "pol" = "Poland", 
                                                                "gbr" = "United Kingdom"))) +
  labs(title = "Composition of PM emissions", x="Year", y="kt")+
  plot_theme +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
        strip.text.x = element_text(size = 10)) +
  ggsave("./figures/PM_composition_by_cntry_open_burn_BioBAvg.png", height = 8.5, width = 11, units = "in")

# Selected countries - Ambient PM fraction of deaths extrapolated (yearly open burning data)

p <-
  ggplot() +
  geom_point(data=filter(PM_deaths_frac_cntry_allyears_extrap, year %in% c(1850:2015)),
             aes(x=(FFI / total_pm), y=val, color = iso, fill = iso,  shape = iso), size = 8) +
  shadow_mark(past = TRUE, future = FALSE, size = 1) +
  transition_time(as.integer(year)) +
  labs(title = "Year: {frame_time}", x="Fossil_Industrial Fraction of Composite PM Emissions", 
       y = "Fraction of Deaths from Ambient PM Pollution") +
  scale_y_continuous(limits = c(0, NA),
                     labels = scales::percent_format(accuracy = 1L)) +
  scale_x_continuous(limits = c(0, NA),
                     labels = scales::percent_format(accuracy = 1L)) +
  scale_color_manual(name = "Country",
                     values = brewer.pal(3, "Dark2"),
                     labels = c("fra" = "France", "chn" = "China", "ind" = "India")) +
  scale_fill_manual(name = "Country",
                    values = brewer.pal(3, "Dark2"),
                    labels = c("fra" = "France", "chn" = "China", "ind" = "India")) +
  scale_shape_manual(name = "Country",
                     values = c(21, 22, 23),
                     labels = c("fra" = "France", "chn" = "China", "ind" = "India")) +
  plot_theme 
  animate(p, nframes = 200, fps = 5, start_pause = 5, end_pause = 20)
  anim_save(filename="./figures/selected_countries_PM_share_deaths_extrap.gif")
