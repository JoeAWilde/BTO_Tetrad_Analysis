#-----------Intro-------------#

# Author: Luke Ozsanlav-Harris
# Creation: 16/05/2023

#----------------------------#

## Packages
pacman::p_load(tidyverse, data.table, sf, ggnewscale, viridis)
options(scipen = 999) # stop scientific notation



##---------------------------------##
#### 1. Load in tetrad data data ####
##---------------------------------##

## read in tetrad data and clean a bit
Abund <- fread("Data/2km_dispersal_allcounties.csv")
Abund <- Abund[duplicated(Abund$`5figgrid`)==F,] # remove some duplicates that are in Devon and Cornwall
Abund <- filter(Abund, !`5figgrid` == "5figref")
Abund <- filter(Abund, !`5figgrid` == "")

## streamline data set a bit
colnames(Abund)
Abund <- Abund %>%  dplyr::select(-c("V5", "V6", "V8", "x", "y", "6figgrid", "dum_x", "dum_y", 
                              "APHA_rel_ph&rlp", "APHA_rel_present", "logAPHA_rel_ph&rlp",           
                              "Advert_presence","Advert_details", "gamebird_breed", "gamebird_winter"))


## Load in Tetrad data and filter to the terads that we have abundance data for
Tets <- st_read("Data/SouthEnglandTetrads/SouthEngland_Tetrads.shp")
Tets <- filter(Tets, GridRef %in% c(Abund$`5figgrid`)) 

## quick plot to visualize
ggplot() +
  geom_sf(data = Tets, aes(geometry = geometry), fill = "#BDC3C7", alpha = 0.7) +
  theme_light()




##----------------------------------##
#### 2. Load in PA data (All Hab) ####
##----------------------------------##

## Read in protected area data
SPA <- st_read("Outputs/Script_0/SPA_SEng.shp")
SAC <- st_read("Outputs/Script_0/SAC_SEng.shp")
RAMSAR <- st_read("Outputs/Script_0/RAMSAR_SEng.shp")
SSSI <- st_read("Outputs/Script_0/SSSI_SEng.shp")

## quick plot to visualise
ggplot() +
  geom_sf(data = SAC, aes(geometry = geometry), fill = "#6f9969", alpha = 0.7) +
  theme_light()

## transform the crs of the tetrads so it matches the PA overlaps
Tets <- st_transform(Tets, crs = st_crs(SPA))




##---------------------------------------##
#### 3. Calculate PA overlap (All Hab) ####
##---------------------------------------##

## Firstly, just calcualte the size of each grid square
Tets <- Tets %>% mutate(GridArea = st_area(.))


##-- SPA Overlap --##

## Union the SPA areas
## This solves problems in the next step if one grid square is overlapped by multiple different SPAs
SPA_Un <- st_union(SPA)

## Calculate overlaps with SPA
SPA_lap <- st_intersection(Tets, SPA_Un) %>%  # Find the intersects between the two sets of shapes
           mutate(intersect_area = st_area(.)) %>% ## Calculate the area of the overlapping areas
           select(GridRef, intersect_area) %>% # keep only the columns needed for the join below
           st_drop_geometry() # drop geometry for the join

## Join the overlap back onto the data
Tets_SPA <- full_join(Tets, SPA_lap, by = "GridRef") %>% 
            mutate(intersect_area = ifelse(is.na(intersect_area)==T, 0, intersect_area),
                   SPA = as.numeric(intersect_area/GridArea)) %>% 
            select(-intersect_area)

## stop if join gone weird
stopifnot(nrow(Tets_SPA)==nrow(Tets))   
rm(SPA_lap)



##-- SAC Overlap --##

## Union the SAC areas
## This solves problems in the next step if one grid square is overlapped by multiple different SACs
SAC_Un <- st_union(SAC)

## Calculate overlaps with RAMSAR
SAC_lap <- st_intersection(Tets, SAC_Un) %>%  # Find the intersects between the two sets of shapes
            mutate(intersect_area = st_area(.)) %>% ## Calculate the area of the overlapping areas
            select(GridRef, intersect_area) %>% # keep only the columns needed for the join below
            st_drop_geometry() # drop geometry for the join

## Join the overlap back onto the data
Tets_SAC <- full_join(Tets_SPA, SAC_lap, by = "GridRef") %>% 
            mutate(intersect_area = ifelse(is.na(intersect_area)==T, 0, intersect_area),
                   SAC = as.numeric(intersect_area/GridArea)) %>% 
            select(-intersect_area)

## stop if join gone weird
stopifnot(nrow(Tets)==nrow(Tets_SAC)) 
rm(SAC_lap)


##-- SSSI Overlap --##

## Union the SSSI areas
## This solves problems in the next step if one grid square is overlapped by multiple different SSSIs
SSSI_Un <- st_union(SSSI)

## Calculate overlaps with RAMSAR
SSSI_lap <- st_intersection(Tets, SSSI_Un) %>%  # Find the intersects between the two sets of shapes
            mutate(intersect_area = st_area(.)) %>% ## Calculate the area of the overlapping areas
            select(GridRef, intersect_area) %>% # keep only the columns needed for the join below
            st_drop_geometry() # drop geometry for the join

## Join the overlap back onto the data
Tets_SSSI <- full_join(Tets_SAC, SSSI_lap, by = "GridRef") %>% 
            mutate(intersect_area = ifelse(is.na(intersect_area)==T, 0, intersect_area),
                   SSSI = as.numeric(intersect_area/GridArea)) %>% 
            select(-intersect_area)

## stop if join gone weird
stopifnot(nrow(Tets)==nrow(Tets_SSSI)) 
rm(SSSI_lap)


##-- RAMSAR Overlap --##

## Union the RAMSAR areas
## This solves problems in the next step if one grid square is overlapped by multiple different RAMSARs
RAMSAR_Un <- st_union(RAMSAR)

## Calculate overlaps with RAMSAR
RAMSAR_lap <- st_intersection(Tets, RAMSAR_Un) %>%  # Find the intersects between the two sets of shapes
              mutate(intersect_area = st_area(.)) %>% ## Calculate the area of the overlapping areas
              select(GridRef, intersect_area) %>% # keep only the columns needed for the join below
              st_drop_geometry() # drop geometry for the join

## Join the overlap back onto the data
Tets_RAMSAR <- full_join(Tets_SSSI, RAMSAR_lap, by = "GridRef") %>% 
                mutate(intersect_area = ifelse(is.na(intersect_area)==T, 0, intersect_area),
                       RAMSAR = as.numeric(intersect_area/GridArea)) %>% 
                select(-intersect_area)

## stop if join gone weird
stopifnot(nrow(Tets)==nrow(Tets_RAMSAR)) 
rm(RAMSAR_lap)


##-- All Protected area overlap --##

## Calculate overlaps with SPA
SPA_lap <- st_intersection(Tets, SPA_Un) %>%  # Find the intersects between the two sets of shapes
  select(GridRef) 
 
## Calculate overlaps with RAMSAR
SAC_lap <- st_intersection(Tets, SAC_Un) %>%  # Find the intersects between the two sets of shapes
  select(GridRef) 
 
## Calculate overlaps with RAMSAR
SSSI_lap <- st_intersection(Tets, SSSI_Un) %>%  # Find the intersects between the two sets of shapes
  select(GridRef) 
 
## Calculate overlaps with RAMSAR
RAMSAR_lap <- st_intersection(Tets, RAMSAR_Un) %>%  # Find the intersects between the two sets of shapes
  select(GridRef) 

## Join them all together
AllInters <- rbind(SPA_lap, SAC_lap) %>% rbind(SSSI_lap) %>% rbind(RAMSAR_lap)

## Now union the intersections but group by Grid ref
AllIntersJoin <- AllInters %>% group_by(GridRef) %>%
                summarise(geometry= st_union(geometry))

## Now calcualte the areas of these joined bits
AllIntersects <- AllIntersJoin %>% 
                 mutate(intersect_area = st_area(.)) %>% 
                  select(GridRef, intersect_area) %>% # keep only the columns needed for the join below
                  st_drop_geometry() # drop geometry for the join

## Join the overlap back onto the data
Tets_PA <- full_join(Tets_SSSI, AllIntersects, by = "GridRef") %>% 
              mutate(intersect_area = ifelse(is.na(intersect_area)==T, 0, intersect_area),
                     PA = as.numeric(intersect_area/GridArea)) %>% 
              select(-intersect_area)

## stop if join gone weird
stopifnot(nrow(Tets)==nrow(Tets_PA)) 




##-------------------------------------##
#### 5. Plot of PA overlap (All Hab) ####
##-------------------------------------##

## Get all the intersection bits again
AllIntersects <- AllIntersJoin %>% 
                  mutate(intersect_area = st_area(.)) %>% 
                  select(GridRef) 

## Add Uk outline
UK <- st_read("Data/CoastOutline/UK_Coastline.shp")
UK <- st_transform(UK, crs = st_crs(Tets)) # change crs to th same as tetrads

## Create bounding box for the map
box <- sf::st_bbox(Tets)



## plot to visualize Tetrads and shape of PAs
PAMap <- ggplot() +
  geom_sf(data = UK, aes(geometry = geometry), fill = "#ECF0F1") +
  geom_sf(data = Tets, aes(geometry = geometry), fill = "#BDC3C7", colour = "#BDC3C7", alpha = 0.6) +
  geom_sf(data = AllIntersects, aes(geometry = geometry,  fill = "#6f9969"), colour = NA, alpha = 0.9) +
  scale_fill_manual(name = expression(""),
                    values =c("#6f9969"="#6f9969"),
                    labels = c("Protected Areas")) +
  coord_sf(xlim = box[c("xmin", "xmax")], 
           ylim = box[c("ymin", "ymax")], crs = 27700, expand = T) +
  labs(x = "Longitude", y = "Latitude") +
  # set theme
  theme_light() +
  theme(legend.position = c(0.89,0.08),
        legend.box.background=element_rect(colour = "#BDC3C7"),
        legend.text = element_text(size = 13),
        axis.text=element_text(colour="black"),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(hjust=0.7),
        axis.title.y = element_text(angle=90, vjust = 0.4, size = 15),
        axis.text.y = element_text(hjust=0.7, angle=90, vjust=0.3))

## save plot
ggsave(plot = PAMap, 
       filename = "Outputs/script_1/PA_Tetrad_Map.png",
       units = "mm", width = 250, height = 160, dpi = 300,   
)



## plot to visualize heat map of PA coverage
PAMap2 <- ggplot() +
  geom_sf(data = UK, aes(geometry = geometry), fill = "#ECF0F1") +
  geom_sf(data = Tets_PA, aes(geometry = geometry, fill = PA), colour = NA, alpha = 0.8) +
  scale_fill_viridis( begin = 0, option = "C",
                       end = 1, name = "PA Coverage") +
  coord_sf(xlim = box[c("xmin", "xmax")], 
           ylim = box[c("ymin", "ymax")], crs = 27700, expand = T) +
  labs(x = "Longitude", y = "Latitude") +
  theme(axis.text=element_text(colour="black"),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(hjust=0.7),
        axis.title.y = element_text(angle=90, vjust = 0.4, size = 15),
        axis.text.y = element_text(hjust=0.7, angle=90, vjust=0.3)) +
  # set theme
  theme_light() +
  theme(legend.position = c(0.89,0.18),
        legend.box.background=element_rect(colour = "#BDC3C7"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13),
        axis.text=element_text(colour="black"),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(hjust=0.7),
        axis.title.y = element_text(angle=90, vjust = 0.4, size = 15),
        axis.text.y = element_text(hjust=0.7, angle=90, vjust=0.3))

## save plot
ggsave(plot = PAMap2, 
       filename = "Outputs/script_1/PA_Tetrad_HeatMap.png",
       units = "mm", width = 250, height = 160, dpi = 300,   
)


##-----------------------------------##
#### 6. Load in PA data (Suit Hab) ####
##-----------------------------------##

## remove some of old object, might help memory a bit
rm(SAC, SPA, SSSI, RAMSAR, SAC_Un, SPA_Un, SSSI_Un, RAMSAR_Un)

## Load in Tetrad data and filter to the terads that we have abundance data for
## load them in again just so start from a blank slate
Tets <- st_read("Data/SouthEnglandTetrads/SouthEngland_Tetrads.shp")
Tets <- filter(Tets, GridRef %in% c(Abund$`5figgrid`))

## Read in protected area data with ecologically irrelevant habitats removed
SPAr <- st_read("Outputs/Script_0/SPA_SEng_SuitHab.shp")
SACr <- st_read("Outputs/Script_0/SAC_SEng_SuitHab.shp")
SSSIr <- st_read("Outputs/Script_0/SSSI_SEng_SuitHab.shp")

## quick plot to visualise
ggplot() +
  geom_sf(data = SACr, aes(geometry = geometry), fill = "#6f9969", alpha = 0.7) +
  theme_light()

## transform the crs of the tetrads so it matches the PA overlaps
Tets <- st_transform(Tets, crs = st_crs(SPAr))




##----------------------------------------##
#### 7. Calculate PA overlap (Suit Hab) ####
##----------------------------------------##

## Firstly, just calcualte the size of each grid square
Tets <- Tets %>% mutate(GridArea = st_area(.))

##-- All Suitable protected area overlap --##

## Union the SPA areas
## This solves problems in the next step if one grid square is overlapped by multiple different SPAs
SPAr_Un <- st_union(SPAr)

## Union the SAC areas
## This solves problems in the next step if one grid square is overlapped by multiple different SACs
SACr_Un <- st_union(SACr)

## Union the SSSI areas
## This solves problems in the next step if one grid square is overlapped by multiple different SSSIs
SSSIr_Un <- st_union(SSSIr)

## Calculate overlaps with SPA
SPAr_lap <- st_intersection(Tets, SPAr_Un) %>%  # Find the intersects between the two sets of shapes
  select(GridRef) 

## Calculate overlaps with RAMSAR
SACr_lap <- st_intersection(Tets, SACr_Un) %>%  # Find the intersects between the two sets of shapes
  select(GridRef) 

## Calculate overlaps with RAMSAR
SSSIr_lap <- st_intersection(Tets, SSSIr_Un) %>%  # Find the intersects between the two sets of shapes
  select(GridRef) 


## Join them all together
AllInters2 <- rbind(SPAr_lap, SACr_lap) %>% rbind(SSSIr_lap)

## Now union the intersections but group by Grid ref
AllIntersJoin2 <- AllInters2 %>% group_by(GridRef) %>%
  summarise(geometry= st_union(geometry))

## Now calcualte the areas of these joined bits
AllIntersects2 <- AllIntersJoin2 %>% 
  mutate(intersect_area = st_area(.)) %>% 
  select(GridRef, intersect_area) %>% # keep only the columns needed for the join below
  st_drop_geometry() # drop geometry for the join


## Join the overlap back onto the data
# Tets_SSSI<- Tets
Tets_PA2 <- full_join(Tets_PA, AllIntersects2, by = "GridRef") %>% 
  mutate(intersect_area = ifelse(is.na(intersect_area)==T, 0, intersect_area),
         PA_suit = as.numeric(intersect_area/GridArea)) %>% 
  select(-intersect_area)

## stop if join gone weird
stopifnot(nrow(Tets)==nrow(Tets_PA2)) 




##--------------------------------------##
#### 8. Plot of PA overlap (Suit Hab) ####
##--------------------------------------##

## Get all the intersection bits again
AllIntersects2 <- AllIntersJoin2 %>% 
  mutate(intersect_area = st_area(.)) %>% 
  select(GridRef) 

## Add Uk outline
UK <- st_read("Data/CoastOutline/UK_Coastline.shp")
UK <- st_transform(UK, crs = st_crs(Tets)) # change crs to th same as tetrads

## Create bounding box for the map
box <- sf::st_bbox(Tets)


## plot to visualize Tetrads and shape of PAs
PAMap3 <- ggplot() +
  geom_sf(data = UK, aes(geometry = geometry), fill = "#ECF0F1") +
  geom_sf(data = Tets, aes(geometry = geometry), fill = "#BDC3C7", colour = "#BDC3C7", alpha = 0.6) +
  geom_sf(data = AllIntersects2, aes(geometry = geometry,  fill = "#6f9969"), colour = NA, alpha = 0.9) +
  scale_fill_manual(name = expression(""),
                    values =c("#6f9969"="#6f9969"),
                    labels = c("Protected Areas \n(Suitable Habitat)")) +
  coord_sf(xlim = box[c("xmin", "xmax")], 
           ylim = box[c("ymin", "ymax")], crs = 27700, expand = T) +
  labs(x = "Longitude", y = "Latitude") +
  # set theme
  theme_light() +
  theme(legend.position = c(0.89,0.08),
        legend.box.background=element_rect(colour = "#BDC3C7"),
        legend.text = element_text(size = 13),
        axis.text=element_text(colour="black"),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(hjust=0.7),
        axis.title.y = element_text(angle=90, vjust = 0.4, size = 15),
        axis.text.y = element_text(hjust=0.7, angle=90, vjust=0.3))

## save plot
ggsave(plot = PAMap3, 
       filename = "Outputs/script_1/PA_Tetrad_Map_SuitHab.png",
       units = "mm", width = 250, height = 160, dpi = 300,   
)


## plot to visualize heat map of PA coverage
PAMap4 <- ggplot() +
  geom_sf(data = UK, aes(geometry = geometry), fill = "#ECF0F1") +
  geom_sf(data = Tets_PA2, aes(geometry = geometry, fill = PA_suit), colour = NA, alpha = 0.8) +
  scale_fill_viridis( begin = 0, option = "C",
                      end = 1, name = "PA Coverage \n(Suitable Habitat)") +
  coord_sf(xlim = box[c("xmin", "xmax")], 
           ylim = box[c("ymin", "ymax")], crs = 27700, expand = T) +
  labs(x = "Longitude", y = "Latitude") +
  theme(axis.text=element_text(colour="black"),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(hjust=0.7),
        axis.title.y = element_text(angle=90, vjust = 0.4, size = 15),
        axis.text.y = element_text(hjust=0.7, angle=90, vjust=0.3)) +
  # set theme
  theme_light() +
  theme(legend.position = c(0.89,0.18),
        legend.box.background=element_rect(colour = "#BDC3C7"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13),
        axis.text=element_text(colour="black"),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(hjust=0.7),
        axis.title.y = element_text(angle=90, vjust = 0.4, size = 15),
        axis.text.y = element_text(hjust=0.7, angle=90, vjust=0.3))

## save plot
ggsave(plot = PAMap4, 
       filename = "Outputs/script_1/PA_Tetrad_HeatMap_SuitHab.png",
       units = "mm", width = 250, height = 160, dpi = 300,   
)

##--------------------------------------------------##
#### 9. Join overlaps to main data set (Suit Hab) ####
##--------------------------------------------------##

## Change column type for the join 
Tets_PA2$GridRef <- as.character(Tets_PA2$GridRef)

## join the data together
TetAb <- full_join(Abund, Tets_PA2, by = c("5figgrid" = "GridRef"))
stopifnot(nrow(TetAb)==nrow(Abund)) # check for lost data



##----------------------------------------------##
#### 10. Format data set for model (Suit Hab) ####
##----------------------------------------------##

## Add the coordinates onto the data set
## These may be useful for the modelling
coords <- as.data.frame(st_coordinates(st_centroid(TetAb$geometry)))

## bind these onto the data sets
TetAb <- cbind(TetAb, coords)

## Write out data as a shapefile
write_sf(TetAb, "Outputs/script_1/ModelData_Spatial.shp")

## write out data with the spatial part removed
TetFinal <- TetAb %>% st_drop_geometry() %>% select(-geometry)
write_csv(TetFinal, file = "Outputs/script_1/ModelData.csv")
