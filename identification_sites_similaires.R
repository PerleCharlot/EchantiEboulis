### Titre -------------------------------------
# Nom : Identification d'éboulis similaires
# Auteure : Perle Charlot
# Date de création : 07-07-2022
# Dates de modification : -2022

### Librairies -------------------------------------

library(data.table)
library(sf)
library(raster)
library(exactextractr)
library(ggplot2)
library(dplyr)
library(patchwork)

### Fonctions -------------------------------------

### Constantes -------------------------------------

# Espace de travail
wd <- getwd()
# Dossier des inputs (dans le git)
input_path <- paste0(wd,"/input/")
# Dossier des outputs (dans le git)
output_path <- paste0(wd,"/output/")
# Projection Lambert 93 (EPSG : 2154)
EPSG_2154 =  "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +units=m +no_defs "

#### Données spatiales ####

# Dossier des variables spatiales & chemins des fichiers
path_mnt <- paste0(input_path,"/mnt_25m_belledonne_cale.tif")
path_pente <- paste0(input_path,"/pente_25m.tif")
path_easting <- paste0(input_path,"/easting_25m.tif")
path_northing <- paste0(input_path,"/northing_25m.tif")
path_vect_zone_echant <- paste0(input_path,"/zones_echanti.gpkg")
path_habitat_eboulis <- paste0(input_path,"/polygones_eboulis_N2000.gpkg")

#### Tables ####

### Programme -------------------------------------

# Import couches vectorielles
zones_echanti <- st_read(path_vect_zone_echant)
masque_eboulis <- st_read(path_habitat_eboulis)

# Rasters de critères pour l'échantillonnage
stack_crit <- stack(path_mnt,path_pente,path_easting,path_northing)
# png(file=paste0(output_path,"/cartes_criteres.png"),
#     width=800, height=800)
# plot(stack_crit)
# dev.off()

# Extraction des valeurs de la stack aux zones d'échantillonnage
df_zones_echanti <- do.call(rbind,
                            exact_extract(stack_crit,zones_echanti, 
                                          include_cols=c("nom_zone","pres_senti"),
                                          include_xy = T) ) # 81 pixels par zone carrée

# Visualisation de la distribution des données
p1 = df_zones_echanti %>% ggplot(aes(x=nom_zone, y=pente_25m))+
  geom_boxplot()+
  theme(axis.title.x = element_blank())+labs(y="Pente (en °)")
p2 = df_zones_echanti %>% ggplot(aes(x=nom_zone, y=mnt_25m_belledonne_cale))+
  geom_boxplot()+
  theme(axis.title.x = element_blank())+labs(y="Altitude (en m)")
p3 = df_zones_echanti %>% ggplot(aes(x=nom_zone, y=easting_25m))+
  geom_boxplot()+
  theme(axis.title.x = element_blank())+labs(y="Easting")
p4 = df_zones_echanti %>% ggplot(aes(x=nom_zone, y=northing_25m))+
  geom_boxplot()+
  theme(axis.title.x = element_blank())+labs(y="Northing")
png(file=paste0(output_path,"/distribution_zones_echantillonnage.png"),
    width=800, height=800)
(p1 + p2)/(p3+p4)
dev.off()

# Identification des valeurs médiannes, par site
df_intervalles = df_zones_echanti %>% group_by(nom_zone) %>% summarize(borne_inf_pente = quantile(pente_25m, 0.2),
                                                      borne_sup_pente = quantile(pente_25m, 0.8),
                                                      borne_inf_alti = quantile(mnt_25m_belledonne_cale, 0.2),
                                                      borne_sup_alti = quantile(mnt_25m_belledonne_cale, 0.8),
                                                      borne_inf_east = quantile(easting_25m, 0.2),
                                                      borne_sup_east = quantile(easting_25m, 0.8),
                                                      borne_inf_north = quantile(northing_25m, 0.2),
                                                      borne_sup_north = quantile(northing_25m, 0.8))

# Extraction des valeurs de la stack en zone d'éboulis
vals <- do.call(rbind,
                exact_extract(stack_crit,masque_eboulis, 
                                          include_cols=c("lbcb"),
                                          include_xy = T) )

writeRaster(rasterFromXYZ(data.frame(vals$x, vals$y, vals$pente_25m),crs=EPSG_2154),
            paste0(output_path,"/pente_eboulis.tif"), overwrite=T)
writeRaster(rasterFromXYZ(data.frame(vals$x, vals$y, vals$mnt_25m_belledonne_cale),crs=EPSG_2154),
            paste0(output_path,"/altitude_eboulis.tif"), overwrite=T)
writeRaster(rasterFromXYZ(data.frame(vals$x, vals$y, vals$easting_25m),crs=EPSG_2154),
            paste0(output_path,"/easting_eboulis.tif"), overwrite=T)
writeRaster(rasterFromXYZ(data.frame(vals$x, vals$y, vals$northing_25m),crs=EPSG_2154),
            paste0(output_path,"/northing_eboulis.tif"), overwrite=T)


# Identification zones similaires

# Pour le secteur du lac Merlat
alti1 = df_intervalles$borne_inf_alti[df_intervalles$nom_zone == "lac merlat"]
alti2 = df_intervalles$borne_sup_alti[df_intervalles$nom_zone == "lac merlat"]
pente1 = df_intervalles$borne_inf_pente[df_intervalles$nom_zone == "lac merlat"]
pente2 = df_intervalles$borne_sup_pente[df_intervalles$nom_zone == "lac merlat"]
east1 = df_intervalles$borne_inf_east[df_intervalles$nom_zone == "lac merlat"]
east2 = df_intervalles$borne_sup_east[df_intervalles$nom_zone == "lac merlat"]
north1 = df_intervalles$borne_inf_north[df_intervalles$nom_zone == "lac merlat"]
north2 = df_intervalles$borne_sup_north[df_intervalles$nom_zone == "lac merlat"]

simil_Merlat = vals %>%                              
  dplyr::filter(dplyr::between(mnt_25m_belledonne_cale, alti1, alti2)) %>% 
  dplyr::filter(dplyr::between(pente_25m, pente1, pente2)) %>% 
  dplyr::filter(dplyr::between(easting_25m, east1, east2)) %>%
  dplyr::filter(dplyr::between(northing_25m, north1, north2))
write.csv(simil_Merlat, paste0(output_path,"/table_pixels_similaires_Merlat.csv"))

rast_simil_Merlat = rasterFromXYZ(data.frame(simil_Merlat$x, simil_Merlat$y, rep(1,dim(simil_Merlat)[1])),
                                  crs=EPSG_2154)
names(rast_simil_Merlat) = "pixel_similaire"
writeRaster(rast_simil_Merlat, paste0(output_path,"/raster_pixels_similaires_Merlat.tif"),
            overwrite=T)

# Pour le secteur des lacs Roberts
alti1 = df_intervalles$borne_inf_alti[df_intervalles$nom_zone == "lacs robert"]
alti2 = df_intervalles$borne_sup_alti[df_intervalles$nom_zone == "lacs robert"]
pente1 = df_intervalles$borne_inf_pente[df_intervalles$nom_zone == "lacs robert"]
pente2 = df_intervalles$borne_sup_pente[df_intervalles$nom_zone == "lacs robert"]
east1 = df_intervalles$borne_inf_east[df_intervalles$nom_zone == "lacs robert"]
east2 = df_intervalles$borne_sup_east[df_intervalles$nom_zone == "lacs robert"]
north1 = df_intervalles$borne_inf_north[df_intervalles$nom_zone == "lacs robert"]
north2 = df_intervalles$borne_sup_north[df_intervalles$nom_zone == "lacs robert"]

simil_Robert = vals %>%                              
  dplyr::filter(dplyr::between(mnt_25m_belledonne_cale, alti1, alti2)) %>% 
  dplyr::filter(dplyr::between(pente_25m, pente1, pente2)) %>% 
  dplyr::filter(dplyr::between(easting_25m, east1, east2)) %>%
  dplyr::filter(dplyr::between(northing_25m, north1, north2))
write.csv(simil_Robert, paste0(output_path,"/table_pixels_similaires_Robert.csv"))

rast_simil_Robert = rasterFromXYZ(data.frame(simil_Robert$x, simil_Robert$y, rep(1,dim(simil_Robert)[1])),
                                  crs=EPSG_2154)
names(rast_simil_Robert) = "pixel_similaire"
writeRaster(rast_simil_Robert, paste0(output_path,"/raster_pixels_similaires_Robert.tif"),
            overwrite=T)
