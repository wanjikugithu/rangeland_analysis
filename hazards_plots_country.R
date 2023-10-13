library(tmap)
library(tmaptools)
library(raster)
library(tidyverse)
library(terra)

# population files
datadir <- "" # dir path
rp <- list.files(file.path(datadir, "rural_population"),  
                 pattern = "_rural_pop_10km.tif", full.names = TRUE)

# vv <- shapefile(file.path(datadir, "input/boundary/country_farming_system_cg_regions.shp"))
vv <- shapefile(file.path(datadir, "vector/country_farming_system_cg_regions.shp"))

vv <- vv[,c("cgregin", "NAME_EN", "ISO_A3" ,"frmng_s")]
names(vv) <- c("CG_region", "Country", "ISO_A3", "Farming_System" )
vv$Farming_System <- gsub("[[:digit:]]+","", vv$Farming_System)
vv$Farming_System <- gsub("\\.","", vv$Farming_System)
vv$Farming_System <- trimws(vv$Farming_System)

# Hazard labels with the 10 categories
haz_lab <- data.frame(HAZARD = seq(0, 10, by = 1),
                      HAZ_LAB = c("No hazard", "Extremely high climate variability (ExC)", "High climate variability (HCV)", "Thermal stress (THI)",
                                  "Climate variability (ExC+HCV) + THI", "CV + Drought (D)", "CV + Flood (F)", "CV + THI + D", "CV + THI + F",
                                  "CV + D + F", "CV + THI + D + F"),
                      colors = c('#d9d9d9','#fccde5','#bebada','#fb8072','#b3de69','#fdb462','#80b1d3','#ffffb3','#8dd3c7','#bc80bd','#ccebc5'))

# Hazards map & stacked rural population plots (by farming system)

getSummaryHazardFS <- function(iso, vv, rp, haz_lab){
  cat("Processing ", iso, "\n.....")
  
  # boundary
  v <- vv[vv$ISO_A3 == iso,] 
  v <- v[-grep("fishing", v$Farming_System),]
  v <- aggregate(v, by = "Farming_System", dissolve=TRUE)
  v$FarmSyst <- 1:nrow(v)
  fs_tab <- v@data
  
  # hazard layer
  rh <- rast("data/hazard_classes.tif")
  
  # rural pop
  pop <- grep(iso, rp, value = TRUE, ignore.case = TRUE)
  pop <- rast(pop)
  pop[pop < 100] <- NA
  
  haz <- crop(rh, pop)
  pop <- resample(pop, haz)
  haz <- mask(haz, pop)
  
  haz <- stack(haz)
  lbls <- haz_lab$HAZ_LAB
  cc <- sampleRandom(haz, length(lbls), cells=TRUE)
  haz[cc[,1]] <- 1:length(lbls)-1
  
  hmap <- tm_shape(haz) +
    tm_raster(style = "cat", labels = lbls, title = "Hazard categories",
              palette = haz_lab$colors,
              legend.is.portrait = TRUE) + 
    tm_shape(v) +
    tm_borders()+
    tm_layout(frame = FALSE,
              legend.title.size = 3,
              legend.outside = TRUE,
              legend.outside.position = "right",
              legend.outside.size = 0.4,
              legend.text.size = 1.25)
  
  tmap_save(hmap, paste0("output/",iso,"_hazards_map.png"), width=12, height=10, units="in", dpi=300)
  
  vr <- rasterize(v, haz, "FarmSyst")
  
  rr <- stack(vr, haz, stack(pop))
  names(rr) <- c( "FarmSyst", "HAZARD", "rural_pop")
  
  dd <- as.data.frame(rr)
  dd <- dd[complete.cases(dd),]
  
  tab <- dd %>%
    group_by(FarmSyst, HAZARD) %>%
    summarise(RURALPOP = sum(rural_pop, na.rm=TRUE)) %>%
    left_join(., haz_lab, by="HAZARD") %>%
    left_join(., fs_tab, by="FarmSyst") %>%
    mutate(RURALPOP=round(RURALPOP/10^6,2)) %>%
    filter(HAZARD != 0) %>%
    as_tibble()
  
  haz_lab1 <- haz_lab[haz_lab$HAZARD %in% unique(tab$HAZARD),]
  
  # rural population
  gplot2 <- ggplot(tab, aes(x=factor(Farming_System), y=RURALPOP, fill=factor(HAZ_LAB, levels=haz_lab1$HAZ_LAB))) + 
    geom_bar(stat="identity", position="stack") +
    scale_fill_manual(values = haz_lab1$colors) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(fill="Hazard category") +
    xlab("") + 
    ylab("Total rural population affected [million people]")
  print(gplot2)
  
  ggsave(gplot2, 
         filename=paste0("output/",iso,"_hazards_fs_stacked_ruralpop_hazard.png"), 
         dpi=300, height=13, width=18, units="cm")
  
  tab %>% 
    select(Farming_System, HAZ_LAB, RURALPOP) %>%
    write_csv(paste0("output/",iso,"_hazards_fs_stacked_ruralpop_hazard.csv"))
}

iso3 <- c("KEN", "ETH", "NGA")

lapply(iso3, getSummaryHazardFS, vv, rp, haz_lab)

