#-----------Intro-------------#

# Author: Luke Ozsanlav-Harris
# Creation: 16/05/2023
# Goal: Extract PAs overlapping selected BTO tetrads 
#       Remove PAs that are ecologically irrelevant to pheasants

#----------------------------#

## Packages
pacman::p_load(tidyverse, data.table, sf, terra)


##-----------------------------##
#### 1. Load in spatial data ####
##-----------------------------##

## read in tetrad data and clean a bit
Abund <- fread("Data/2km_dispersal_allcounties.csv")
Abund <- Abund[duplicated(Abund$`5figgrid`)==F,] # remove some duplicates that are in Devon and Cornwall
Abund <- filter(Abund, !`5figgrid` == "5figref")
Abund <- filter(Abund, !`5figgrid` == "")

## streamline data set, really just want to keep grid refs
Abund <- Abund %>%  select(c("County", "5figgrid"))

## Load in rough Devon outline for cropping large UK-wide data sets
Tets <- st_read("SpatialData/SouthEnglandTetrads/SouthEngland_Tetrads.shp")
Tets <- filter(Tets, GridRef %in% c(Abund$`5figgrid`)) 
rm(Abund)

## Read in shapefile of protected area data
SPA <- st_read("SpatialData/SPAs/GB_SPA_OSGB36_20220930.shp")
SAC <- st_read("SpatialData/SACs/GB_SAC_OSGB36_20191031.shp")
RAMSAR <- st_read("SpatialData/RAMSAR/UK_RAMSAR_BNG_20210308.shp")
SSSI <- st_read("SpatialData/SSSIs/Sites_of_Special_Scientific_Interest_England.shp")

## Read in UK CEH habitat data for 2007 for just the south of England
## Use this to work out dominant habitat for each PA
Hab <- rast("SpatialData/Habitat/LCM2015.tif")

## Read on conversion table, converts raster numbers in Hab into meaningful codes
ConvTab <- read_csv("Data/UKCEH_HabConversions.csv")





##-------------------------------##
#### 2. Spatial Transformation ####
##-------------------------------##

## Transform the crs of the tetrads outline to British national grid
Tets <- st_transform(Tets, crs = st_crs(SPA))

## Now crop all of the protected area data
## find the object that intersect with my outline and then just filter them by row 

## SPAs
SPA_Ind <- st_intersects(SPA, Tets, sparse = TRUE)
SPA_SEng <- SPA[(lengths(SPA_Ind) > 0) == TRUE,]
write_sf(SPA_SEng, "Outputs/Script_0/SPA_SEng.shp") # read out the data set
rm(SPA, SPA_Ind)

## SACs
SAC_Ind <- st_intersects(SAC, Tets, sparse = TRUE)
SAC_SEng <- SAC[(lengths(SAC_Ind) > 0) ==TRUE,]
write_sf(SAC_SEng, "Outputs/Script_0/SAC_SEng.shp")
rm(SAC, SAC_Ind)


## SSSIs
SSSI_Ind <- st_intersects(SSSI, Tets, sparse = TRUE)
SSSI_SEng <- SSSI[(lengths(SSSI_Ind) > 0) ==TRUE,]
write_sf(SSSI_SEng, "Outputs/Script_0/SSSI_SEng.shp")
rm(SSSI, SSSI_Ind)


## RAMSAR
RAMSAR_Ind <- st_intersects(RAMSAR, Tets, sparse = TRUE)
RAMSAR_SEng <- RAMSAR[(lengths(RAMSAR_Ind) > 0) ==TRUE,]
write_sf(RAMSAR_SEng, "Outputs/Script_0/RAMSAR_SEng.shp")
rm(RAMSAR, RAMSAR_Ind)



## quick plot to check that it has worked
ggplot() +
  geom_sf(data = Tets, aes(geometry = geometry), fill = "#BDC3C7", alpha = 0.7) +
  geom_sf(data = SPA_SEng, aes(geometry = geometry), fill = "#FA8128") +
  geom_sf(data = SAC_SEng, aes(geometry = geometry), fill = "#6f9969") +
  geom_sf(data = RAMSAR_SEng, aes(geometry = geometry), fill = "#808fe1") +
  geom_sf(data = SSSI_SEng, aes(geometry = geometry), fill = "#efc86e") +
  theme_light()




##-------------------------------##
#### 3. Subset protected areas ####
##-------------------------------##

## Calculate the dom habitat type in each PA retained so far
## Then filter out habitat types that are relevant to pheasants

## QUikc plot of Habitat raster
## see that ocean is grey (0 value)
plot(Hab)

## Convert each of the PA outlines from sf to terra
## Need to do this to wokr with the habitat raster in terra
## NOte not going to do this for RAMSAR as they are all wetlands and we want to remove them anyway
SPA_SEng_ter <- vect(SPA_SEng)
SAC_SEng_ter <- vect(SAC_SEng)
SSSI_SEng_ter <- vect(SSSI_SEng)


## create a function to take a mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

## Now extract the mode habitat type for each polygon 
SPA_HabType <- extract(x = Hab, y = SPA_SEng_ter, fun = getmode)
SAC_HabType <- extract(x = Hab, y = SAC_SEng_ter, fun = getmode)
SSSI_HabType <- extract(x = Hab, y = SSSI_SEng_ter, fun = getmode)

## add that mode habitat type back to the sf object
SPA_SEng$HabCode <- SPA_HabType$LCM2015
SAC_SEng$HabCode <- SAC_HabType$LCM2015
SSSI_SEng$HabCode <- SSSI_HabType$LCM2015

## Convert the habitat types to something meaningful
SPA_SEng <- left_join(SPA_SEng, ConvTab, by = "HabCode")
SAC_SEng <- left_join(SAC_SEng, ConvTab, by = "HabCode")
SSSI_SEng <- left_join(SSSI_SEng, ConvTab, by = "HabCode")

## Now filter out the habitat that are irrelevant to pheasants
## Mainly all the wetland and coastal ones
## These ones were kept: "Heath" "RoughGrass" "ImpGrass" "BroadWood" "AcidGrass" "ConifWood" "Arable" "HeathGrass" "NeutGrass" "CalcGrass" "Marsh"
RemHab <- c("Ocean", "LitRock", "Saltwater", "Freshwater", "InlandRock",  
            "Suburb", "LitSed", "Saltmarsh", "Bog", "Urban", "SupSed")

SPA_SEngF <- filter(SPA_SEng, !HabType %in% RemHab)
SAC_SEngF <- filter(SAC_SEng, !HabType %in% RemHab)
SSSI_SEngF <- filter(SSSI_SEng, !HabType %in% RemHab)


## Now write all these specific habitats to shape files
write_sf(SPA_SEngF, "Outputs/Script_0/SPA_SEng_SuitHab.shp")
write_sf(SAC_SEngF, "Outputs/Script_0/SAC_SEng_SuitHab.shp")
write_sf(SSSI_SEngF, "Outputs/Script_0/SSSI_SEng_SuitHab.shp")
