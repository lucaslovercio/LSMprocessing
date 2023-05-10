#### 0. Load packages ####
library(rgl)
library(geomorph)
library(shapes)
library(devtools)
library(ggplot2)
library(Morpho)
library(corrplot)
library(plyr)
library(dplyr)
library(rrcov)
library(Evomorph)
library(Stack)
library(gridExtra)
library(png)
library(raster)
library(sp)
library(rgdal)
library(magick)
library(stringr)
library(devtools)
install_github("marta-vidalgarcia/morpho.tools.GM")
library(morpho.tools.GM)
install_github("marta-vidalgarcia/symmetry")
library(symmetry)
library(vegan)
library(viridis)
library(ggbiplot)
library(factoextra)
library(Rvcg)
library(ggExtra)
library(Rvcg)
library(wesanderson)
library(dplyr)
library(R.matlab)
source("./R/R_functions/proliferation_import.R")
library(freesurferformats)

# 1. DATA PROCESSING ####
#### 1.1. Load data ####

land3D <- readRDS("./cluster/LM_array_correct_Dec2022.rds") # CORRECT LANDMARK ARRAY

# load("./data/ProliferationBlock_mask03.RData") # PROLIFERATION MATRIX, too big
setwd("./data/proliferation/") # PROLIFERATION MATRIX, re-sampled by Lucas
proliferation_raw <- proliferation.import()
str(proliferation_raw)
setwd("../../")


somites <- read.csv("./data/somites.csv") # SOMITES MATRIX
row.names(somites) <- somites$Name
head(somites)

atlas_e11_mesh_dec <- geomorph::read.ply("./data/E11.5_Tissues_Affine_majority_toE10.5_closed_sphere57_ascii_dec_clean.ply")
atlas_e11_mesh_dec <- Rvcg::vcgQEdecim(atlas_e11_mesh_dec, percent = 0.5)
open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700))
rgl::shade3d(atlas_e11_mesh_dec, color = "gray", alpha = 1)
writePLY("./data/E11.5_Tissues_Affine_majority_toE10.5_closed_sphere57_ascii_dec_smooth_clean.ply")


atlas_mesh <- geomorph::read.ply("./data/E10.5_atlas_simpl_2_dec_ascii.ply")
# atlas_lm <- morpho.tools.GM::pp2array()

n_land <- length(count.fields("./data/E10.5_atlas_simpl_2_dec_ascii_picked_points.pp")) - 9
raw_LM_file <- readLines("./data/E10.5_atlas_simpl_2_dec_ascii_picked_points.pp")[9:(8 + n_land)]
atlas_lm <- matrix(data = NA, nrow = n_land, ncol = 3)
for (j in 1:length(raw_LM_file)) {
  atlas_lm[j, 1] <- as.numeric(strsplit(strsplit(raw_LM_file[j], 
                                                 "x=\"")[[1]][2], "\"")[[1]][1])
  atlas_lm[j, 2] <- as.numeric(strsplit(strsplit(raw_LM_file[j], 
                                                 "y=\"")[[1]][2], "\"")[[1]][1])
  atlas_lm[j, 3] <- as.numeric(strsplit(strsplit(raw_LM_file[j], 
                                                 "z=\"")[[1]][2], "\"")[[1]][1])
}

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = lateral)
rgl::shade3d(atlas_e11_mesh_dec, color = "gray", alpha = 0.3)
rgl::plot3d(e11_atlas_LM[-22,], aspect = "iso", type = "s", size=1.2, col = "darkblue", add = T)
rgl::text3d(x = e11_atlas_LM[, 1],
            y = e11_atlas_LM[, 2],
            z=  e11_atlas_LM[, 3],
            texts = c(1:n_land),
            cex = 1.5, offset = 0.5, pos = 3)

load("./data/RGL_head_pos.rdata")
open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = lateral)
rgl::shade3d(atlas_mesh, color = "gray", alpha = 1)
rgl::plot3d(atlas_lm[-22,], aspect = "iso", type = "s", size=1.2, col = "darkblue", add = T)
# rgl::text3d(x = atlas_lm[, 1],
#             y = atlas_lm[, 2],
#             z=  atlas_lm[, 3],
#             texts = c(1:n_land),
#             cex = 1.5, offset = 0.5, pos = 3)

# So landmarks 30 & 31 are switched in the E11, use the CORRECTED FILE!!!

rgl::rgl.snapshot("./figs/E10_atlas_LMs_SIDE.png", top = TRUE)

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = frontal)
rgl::shade3d(atlas_mesh, color = "gray", alpha = 1)
rgl::plot3d(atlas_lm[-22,], aspect = "iso", type = "s", size=1.2, col = "darkblue", add = T)
# rgl::text3d(x = atlas_lm[, 1],
#             y = atlas_lm[, 2],
#             z=  atlas_lm[, 3],
#             texts = c(1:n_land),
#             cex = 1.5, offset = 0.5, pos = 3)
rgl::rgl.snapshot("./figs/E10_atlas_LMs_FRONTAL.png", top = TRUE)

frontal <- par3d()$userMatrix
lateral <- par3d()$userMatrix
rgl.close()

save(frontal, lateral, file = "./data/RGL_head_pos.rdata")

# non.sym <- c(1:5, 22)
# side.1 <- c(6:21)
# side.2 <- c(23:38)
# pairedLM <- cbind(side.1, side.2)

#### 1.2. Checking data ####
GPA <- geomorph::gpagen(land3D)
morpho.tools.GM::plotOutliers_percentile(GPA$coords, percentile = 0.95) 

# No outliers, distribution looks good
# Check we have the exact same specimen names in the different datasets

setdiff(dimnames(GPA$coords)[[3]], row.names(somites))
setdiff(row.names(somites), dimnames(GPA$coords)[[3]])
setdiff(dimnames(GPA$coords)[[3]], row.names(proliferation_raw))
setdiff(row.names(proliferation_raw), dimnames(GPA$coords)[[3]])
setdiff(row.names(proliferation_raw), row.names(somites))
setdiff(row.names(somites), row.names(proliferation_raw))

# Great, now make sure everything is in order
somites_order <- somites[match(dimnames(GPA$coords)[[3]], row.names(somites)),]
proliferation <- proliferation_raw[match(dimnames(GPA$coords)[[3]], row.names(proliferation_raw)),]

proliferation[1:2, 1:3]
# Data looks good.

saveRDS(proliferation, "./data/proliferation.rds")
# get rid of LM 22 (back of head)
land3D_no_back <- land3D[-22,,]
non.sym <- c(1:5)
side.1 <- c(6:21)
side.2 <- c(22:37)
pairedLM <- cbind(side.1, side.2)

pairedLM

GPA_Morpho <- Morpho::procSym(land3D_no_back, pairedLM = pairedLM)

saveRDS(GPA_Morpho, "./data/GPA_Morpho_corrected.rds")


# str(prolifT_mask03)


prolifT_mask03[1:40, 1:10]
dim(prolifT_mask03)

land3D
dimnames(land3D)[[3]]


#### 1.3. GPA & REGRESSIONS - Locally ####
somites$Name
dimnames(land3D)[[3]]

Somites <- somites$Tail.somite
names(Somites) <- row.names(somites)

gm_df <- geomorph.data.frame(coords_sym = GPA_Morpho$Sym, 
                             coords_asym = GPA_Morpho$Asym, 
                             Csize = GPA_Morpho$size,
                             Somites = Somites, Somites2 = Somites^2)

reg_shape_log.Csize_log.somite <- procD.lm(coords_sym ~ log(Csize) + log(Somites), data = gm_df, RRPP = TRUE)

plot(reg_shape_log.Csize_log.somite, type = "diagnostics", pch = 19, col = "blue") 
plot(reg_shape_log.Csize_log.somite, type = "PC", pch = 19, col = "blue") 
plot(reg_shape_log.Csize_log.somite, type = "regression", 
     predictor = gm_df$Csize, reg.type = "RegScore", 
     pch = 19, col = "green")
plot(reg_shape_log.Csize_log.somite, type = "regression",
     predictor = gm_df$Csize, reg.type = "RegScore",
     pch = 21, bg = "yellow")


saveRDS(gm_df, "./data/gm_df_corrected.rds")
saveRDS(reg_shape_log.Csize_log.somite, "./data/reg_shape_log.Csize_log.somite_corrected.rds")

# GRandmean
mean_shape <- reg_shape_log.Csize_log.somite$coefficients[1,] # Intercept is first coefficient in regression (grand mean)
residuals_coords <- sweep(reg_shape_log.Csize_log.somite$residuals, MARGIN = 2, STATS = mean_shape, FUN = "+")

array_reg_shape_log.Csize_log.somite <- arrayspecs(A = residuals_coords, p = dim(residuals_coords)[2]/3, k = 3)
saveRDS(array_reg_shape_log.Csize_log.somite, "./data/array_reg_shape_log.Csize_log.somite_corrected.rds")


GPA_reg_shape_log.Csize_log.somite <- procSym(array_reg_shape_log.Csize_log.somite)
saveRDS(GPA_reg_shape_log.Csize_log.somite, "./data/GPA_reg_shape_log.Csize_log.somite_corrected.rds")

GPA_reg_shape_log.Csize_log.somite$rotated

# Projected data
mean_shape_p <- reg_shape_log.Csize_log.somite$coefficients[1,] + # GrandMean (intercept)
  median(log(Somites))*reg_shape_log.Csize_log.somite$coefficients[3,] # because # coefficient 2 is Csize
residuals_coords_p <- sweep(reg_shape_log.Csize_log.somite$residuals, MARGIN = 2, STATS = mean_shape_p, FUN = "+")

array_reg_shape_log.Csize_log.somite_p <- arrayspecs(A = residuals_coords_p, p = dim(residuals_coords_p)[2]/3, k = 3)
saveRDS(array_reg_shape_log.Csize_log.somite_p, "./data/array_reg_shape_log.Csize_log.somite_projected_corrected.rds")

GPA_reg_shape_log.Csize_log.somite_p <- procSym(array_reg_shape_log.Csize_log.somite_p)
saveRDS(GPA_reg_shape_log.Csize_log.somite_p, "./data/GPA_reg_shape_log.Csize_log.somite_projected_corrected.rds")


# Submit jobs to cluster and wait. This takes 3 days. Split into 3 jobs if possbile.
# Actually doing it locally now

# 1.3.1. Raw data ####
GPA_Morpho <- readRDS("./data/GPA_Morpho_corrected.rds")
proliferation<- readRDS("./data/proliferation.rds")

# GPA_Morpho$Sym
paste0("Starting PLS between Proliferation and raw (SYM comp.) shape blocks at ", Sys.time())
PLS_raw_sym_m <- Morpho::pls2B(x = GPA_Morpho$Sym, y = proliferation, same.config = FALSE, rounds = 999) # rounds = permutations
PLS_raw_sym_g <- geomorph::two.b.pls(GPA_Morpho$Sym, proliferation)

saveRDS(PLS_raw_sym_m, "./data/PLS_raw_sym_m_corrected.rds")
saveRDS(PLS_raw_sym_g, "./data/PLS_raw_sym_g_corrected.rds")

paste0("All jobs finished at ", Sys.time())


# 1.3.2. Somite-regressed - only GRAND MEAN ####
GPA_reg_shape_log.Csize_log.somite <- readRDS("./data/GPA_reg_shape_log.Csize_log.somite_corrected.rds")

# GPA_reg_shape_log.Csize_log.somite$rotated
paste0("Starting PLS between Proliferation and regressed shape blocks at ", Sys.time())

PLS_reg_sym_m <- Morpho::pls2B(x = GPA_reg_shape_log.Csize_log.somite$rotated, y = proliferation,
                               same.config = FALSE, rounds = 999) # rounds = permutations
PLS_reg_sym_g <- geomorph::two.b.pls(GPA_reg_shape_log.Csize_log.somite$rotated, proliferation)

saveRDS(PLS_reg_sym_m, "./data/PLS_reg_sym_m_corrected.rds")
saveRDS(PLS_reg_sym_g, "./data/PLS_reg_sym_g_corrected.rds")

paste0("All jobs finished at ", Sys.time())

# While we are here, also check the correlation between Somite count and proliferation
GPA_raw <- readRDS("./data/GPA_Morpho.rds")
Somites <- somites$Tail.somite
names(Somites) <- row.names(somites)

prolif_shape_somite_DF <- geomorph.data.frame(GPA_raw, proliferation = prolifT_mask03, Somites = Somites)

prolif_DF <- cbind(prolifT_mask03, Somites)
vars_prolif <- paste0("V", 1:dim(prolif_DF)[2]-1)


reg_prolif_somite <- manova(lapply(vars_prolif[1:10], get) ~ Somites, data = prolif_DF)

summary(reg_prolif_somite)



#### 1.4. Collect RDS from cluster & load data ####
rm(list=ls())
atlas_mesh <- geomorph::read.ply("./data/E10.5_atlas_simpl_2_dec_ascii.ply")
# atlas_lm <- morpho.tools.GM::pp2array()

somites <- read.csv("./data/somites.csv")
row.names(somites) <- somites$Name
head(somites)
Somites <- somites$Tail.somite
names(Somites) <- row.names(somites)

n_land <- length(count.fields("./data/E10.5_atlas_simpl_2_dec_ascii_picked_points.pp")) - 9
raw_LM_file <- readLines("./data/E10.5_atlas_simpl_2_dec_ascii_picked_points.pp")[9:(8 + n_land)]
atlas_lm_raw <- matrix(data = NA, nrow = n_land, ncol = 3)
for (j in 1:length(raw_LM_file)) {
  atlas_lm_raw[j, 1] <- as.numeric(strsplit(strsplit(raw_LM_file[j], 
                                                 "x=\"")[[1]][2], "\"")[[1]][1])
  atlas_lm_raw[j, 2] <- as.numeric(strsplit(strsplit(raw_LM_file[j], 
                                                 "y=\"")[[1]][2], "\"")[[1]][1])
  atlas_lm_raw[j, 3] <- as.numeric(strsplit(strsplit(raw_LM_file[j], 
                                                 "z=\"")[[1]][2], "\"")[[1]][1])
}

atlas_lm <- atlas_lm_raw[-22,] # get rid of LM 22 (back of head)
load("./data/RGL_head_pos.rdata") # RGL positions

# Remember here we are using only symmetric component - GPA Morpho (that was all done in the cluster)
# We got rid of the landmark at the back of the head, maybe we need to exclude some more LMs too




# 2. Which LMs explain the most shape variation? ####

#### 2.1. Shape variance in RAW data ####
GPA_raw <- readRDS("./data/GPA_Morpho_corrected.rds")
source("./R/per_lm_variance.R")

my.variances_red <- per_lm_variance(GPA_raw$rotated)

consensus_lm <- GPA_raw$mshape
colnames(consensus_lm) <- colnames(atlas_lm)
row.names(consensus_lm) <- row.names(atlas_lm)
atlas_lm_m <- as.matrix(atlas_lm)
consensus_mesh <- tps3d(x = atlas_mesh, refmat = atlas_lm_m, tarmat = consensus_lm)

open3d(windowRect = c(20, 30, 800, 800))
shade3d(consensus_mesh, color="gray", alpha=0.95)
spheres3d(consensus_lm, col = my.variances_red$Variance_Colors, radius = 0.015, add = TRUE)

rgl::rgl.snapshot("./figs/Shape_variance_raw_LMs_FRONTAL.png", top = TRUE)
rgl::rgl.snapshot("./figs/Shape_variance_raw_LMs_SIDE.png", top = TRUE)
rgl::rgl.snapshot("./figs/Shape_variance_raw_LMs_TOP.png", top = TRUE)

# Combine the three views
frontal_lm <- image_read("./figs/Shape_variance_raw_LMs_FRONTAL.png")
frontal_lm <- image_annotate(frontal_lm, "Frontal", font = "times", location = "+80+120", size = 50)
frontal_lm <- image_crop(frontal_lm, "1000x500+0+90")
frontal_lm

side_lm <- image_read("./figs/Shape_variance_raw_LMs_SIDE.png")
side_lm <- image_annotate(side_lm, "Lateral", font = "times", location = "+80+120", size = 50)
side_lm <- image_crop(side_lm, "1000x500+0+90")
side_lm

top_lm <- image_read("./figs/Shape_variance_raw_LMs_TOP.png")
top_lm <- image_annotate(top_lm, "Dorsal", font = "times", location = "+80+120", size = 50)
top_lm <- image_crop(top_lm, "1000x500+0+90")
top_lm

stack_views <- c(frontal_lm, side_lm, top_lm)
stacked_images <- image_append(image_scale(stack_views), stack = TRUE)

image_browse(stacked_images)
image_write(stacked_images, path = "./figs/Shape_variance_raw_LMs_all_views.png", format = "png")


#### 2.2. Shape variance in RAW data (with LM 22) ####
load("./data/land3D.RData")

non.sym <- c(1:5)
side.1 <- c(6:21)
side.2 <- c(22:37)
pairedLM <- cbind(side.1, side.2)

pairedLM

GPA_Morpho <- Morpho::procSym(land3D, pairedLM = pairedLM)

my.variances_red <- per_lm_variance(GPA_Morpho$rotated)

which(my.variances_red$Per_Lm_Variance == max(my.variances_red$Per_Lm_Variance))

sort(as.numeric(my.variances_red$Per_Lm_Variance ))

consensus_lm <- GPA_Morpho$mshape
colnames(consensus_lm) <- colnames(atlas_lm_raw)
row.names(consensus_lm) <- row.names(atlas_lm_raw)
atlas_lm_m <- as.matrix(atlas_lm_raw)
consensus_mesh <- tps3d(x = atlas_mesh, refmat = atlas_lm_m, tarmat = consensus_lm)

open3d(windowRect = c(20, 30, 800, 800))
shade3d(atlas_mesh, color="gray", alpha=0.95)
spheres3d(atlas_lm_raw, col = my.variances_red$Variance_Colors, radius = 15, add = TRUE)
rgl::text3d(x = atlas_lm_raw[, 1], y = atlas_lm_raw[, 2], z=  atlas_lm_raw[, 3],
                       texts = c(1:length(atlas_lm_raw)),
                       cex = 1.5, offset = 0.5, pos = 3)
           
rgl::rgl.snapshot("./figs/Shape_variance_raw_LMs_FRONTAL_withLM22.png", top = TRUE)
rgl::rgl.snapshot("./figs/Shape_variance_raw_LMs_SIDE_withLM22.png", top = TRUE)
rgl::rgl.snapshot("./figs/Shape_variance_raw_LMs_TOP_withLM22.png", top = TRUE)

# Combine the three views
frontal_lm <- image_read("./figs/Shape_variance_raw_LMs_FRONTAL_withLM22.png")
frontal_lm <- image_annotate(frontal_lm, "Frontal", font = "times", location = "+80+120", size = 50)
frontal_lm

side_lm <- image_read("./figs/Shape_variance_raw_LMs_SIDE_withLM22.png")
side_lm <- image_annotate(side_lm, "Lateral", font = "times", location = "+80+120", size = 50)
side_lm

top_lm <- image_read("./figs/Shape_variance_raw_LMs_TOP_withLM22.png")
top_lm <- image_annotate(top_lm, "Dorsal", font = "times", location = "+80+120", size = 50)
top_lm

stack_views <- c(frontal_lm, side_lm, top_lm)
stacked_images <- image_append(image_scale(stack_views), stack = TRUE)

image_browse(stacked_images)
image_write(stacked_images, path = "./figs/Shape_variance_raw_LMs_all_views_withLM22.png", format = "png")



# #### 2.3. Shape variance in PROJECTED DATA ####
# 
# GPA_proj <- readRDS("./data/GPA_reg_shape_log.Csize_log.somite.rds")
# source("./R/per_lm_variance.R")
# 
# my.variances_red <- per_lm_variance(GPA_proj$rotated)
# 
# consensus_lm_proj <- GPA_proj$mshape
# colnames(consensus_lm_proj) <- colnames(consensus_lm)
# row.names(consensus_lm_proj) <- row.names(consensus_lm)
# consensus_mesh_proj <- tps3d(x = consensus_mesh, refmat = consensus_lm, tarmat = consensus_lm_proj)
# 
# open3d(windowRect = c(20, 30, 800, 800))
# shade3d(consensus_mesh_proj, color="gray", alpha=0.95)
# spheres3d(consensus_lm_proj, col = my.variances_red$Variance_Colors, radius = 0.015, add = TRUE)
# 
# # It looks off, I really think we need to split them into two groups. E10 and E11.5 are just too different
# rgl::rgl.snapshot("./figs/Shape_variance_proj_LMs_FRONTAL.png", top = TRUE)
# rgl::rgl.snapshot("./figs/Shape_variance_proj_LMs_SIDE.png", top = TRUE)
# rgl::rgl.snapshot("./figs/Shape_variance_proj_LMs_TOP.png", top = TRUE)
# 
# # Combine the three views
# frontal_lm <- image_read("./figs/Shape_variance_proj_LMs_FRONTAL.png")
# frontal_lm <- image_annotate(frontal_lm, "Frontal", font = "times", location = "+80+120", size = 50)
# frontal_lm <- image_crop(frontal_lm, "1000x500+0+90")
# frontal_lm
# 
# side_lm <- image_read("./figs/Shape_variance_proj_LMs_SIDE.png")
# side_lm <- image_annotate(side_lm, "Lateral", font = "times", location = "+80+120", size = 50)
# side_lm <- image_crop(side_lm, "1000x500+0+90")
# side_lm
# 
# top_lm <- image_read("./figs/Shape_variance_proj_LMs_TOP.png")
# top_lm <- image_annotate(top_lm, "Dorsal", font = "times", location = "+80+120", size = 50)
# top_lm <- image_crop(top_lm, "1000x500+0+90")
# top_lm
# 
# stack_views <- c(frontal_lm, side_lm, top_lm)
# stacked_images <- image_append(image_scale(stack_views), stack = TRUE)
# 
# image_browse(stacked_images)
# image_write(stacked_images, path = "./figs/Shape_variance_proj_LMs_all_views.png", format = "png")


# 2.4 PCA & other scatterplots:
GPA_raw <- readRDS("./data/GPA_Morpho_corrected.rds")

PCA_raw <- gm.prcomp(GPA_raw$rotated)
summary(PCA_raw)
str(PCA_raw)

# Delete file if it exists
if (file.exists("./output/PCA_shape_raw.txt")) {
  file.remove("./output/PCA_shape_raw.txt")
}
cat("PCA shape variables raw", capture.output(summary(PCA_raw)), 
    file="./output/PCA_shape_raw.txt", sep="\n", append=TRUE)


somites <- read.csv("./data/somites.csv")
row.names(somites) <- somites$Name
head(somites)

Somites <- somites$Tail.somite
names(Somites) <- row.names(somites)


col_somites <- viridis(n = n_distinct(Somites))

pdf("./figs/PCA_head_shape_raw_PC1_PC2_corrected.pdf", width = 8.25, height = 6)
palette(col_somites)
palette()
Somites <- as.factor(Somites)
plot(PCA_raw, pch = 19, col = Somites, cex = 1.75)
title("PCA of shape coordinates - RAW")
picknplot.shape(plot(PCA_raw))

dev.off()

PCA_comp <- PCA_raw
class(PCA_comp) <- "princomp"

png("./figs/PCA_head_shape_raw_scree_plot_corrected.png", width = 300, height = 300)
# pdf("./figs/PCA_head_shape_raw_scree_plot.pdf", height = 5, width = 5)
pca_scree <- fviz_eig(PCA_comp, addlabels=TRUE, hjust = -0.3,
                      barfill="darkgrey", barcolor ="black",
                      linecolor ="blue") + ylim(0, 85) + 
  theme_classic()

print(pca_scree)
dev.off()

lm_csize_somites <- lm(log(GPA_raw$size) ~ as.numeric(as.character(Somites)))
summary(lm_csize_somites)


pdf("./figs/PC1_shape_raw_Csize_corrected.pdf", width = 8.25, height = 6)
plot(PCA_raw$x[,1], log(GPA_raw$size), pch = 19, col = Somites, cex = 1.75)
# identify(x = PCA_raw$x[,1], y = log(GPA_raw$size), labels = row.names(PCA_raw$x))
title("PC1 shape ~ logCsize")
dev.off()

pdf("./figs/PC2_shape_raw_Csize_corrected.pdf", width = 8.25, height = 6)
plot(PCA_raw$x[,2], log(GPA_raw$size), pch = 19, col = Somites, cex = 1.75)
# identify(x = PCA_raw$x[,2], y = log(GPA_raw$size), labels = row.names(PCA_raw$x))
title("PC2 shape ~ logCsize")
dev.off()

Pdist <- ShapeDist(GPA_raw$rotated, GPA_raw$mshape)
names(Pdist) <- names(GPA_raw$size)

pdf("./figs/Pdist_raw_Csize.pdf", width = 8.25, height = 6)
plot(Pdist, log(GPA_raw$size), pch = 19, col = Somites, cex = 1.75)
# identify(x = Pdist, y = log(GPA_raw$size), labels = row.names(PCA_raw$x))
title("Proc dist ~ logCsize")
dev.off()

# 3. HOW MUCH SHAPE VARIATION IS EXPLAINED BY PROLIFERATION? ####
# Probably only a small amount is explained by proliferation
# Does this mean there is temporal and spatial variation in proliferation?

#### 3.1. PLS shape-proliferation - RAW SHAPE ####
PLS_raw_g <- readRDS("./cluster/PLS_raw_sym_g_corrected.rds")

summary(PLS_raw_g)
# r-PLS: 0.97 
# 
# Effect Size (Z): 2.1611
# 
# P-value: 0.001
# 
# Based on 1000 random permutations

PLS_raw_g$r.pls # 0.9700261 correlation coefficient, report this. VERY HIGH.

n_distinct(Somites)
?rainbow

# col_somites <- wes_palette("Zissou1", n_distinct(Somites), type = "continuous")
col_somites <- viridis(n = n_distinct(Somites))

row.names(PLS_raw_g$XScores)
names(Somites)

palette(col_somites)
palette()


pdf("./figs/RAW_sym_PLS/PLS_raw_sym_Dec2022_corrected.pdf", width = 8, height = 8)
# layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
Somites <- as.factor(Somites)
plot(PLS_raw_g, pch = 19, col = Somites, cex = 1.75)
# identify(x = PLS_raw_g$XScores[,1], y = PLS_raw_g$YScores[,1], labels = row.names(PLS_raw_g$XScores))
# legend("topleft", pch = 19, col = palette(), legend = levels(Skull_gdf_P0_2$genotype))
# mtext(side=3, line=3, at=0.07, adj=0, cex=1.5, "P0-P2")
# Add a simple legend
# Somites <- as.numeric(as.character(Somites))
# legend("topleft", col=palette()[c(1,length(palette()))], pch=19,
#        legend=c(range(Somites)), title = "Somites")
# mtext("r-PLS = 0.969; z = 2.155; p = 0.001", side = 3)
dev.off()

picknplot.shape(plot(PLS_raw_g))

legend_image <- as.raster(matrix(palette(), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F, xlab = '', ylab = '', main = 'Somites')
text(x=1.5, y = seq(8,20,5), labels = format(seq(8,20,5)), cex = 1.5)
rasterImage(legend_image, 0, 0, 1,1)
dev.off()

rm(PLS_raw_g)

# Morpho
PLS_raw_m <- readRDS("./cluster/PLS_raw_sym_m_corrected.rds")
# Ask Benedikt, Rebecca & Lucas if they want to have morphs (and maybe also the heatmap) showing differences between - YES

# Get the shape changes from pls2B associated with each latent variable
# -3*STD and +3*STD for both x & y, or even 5x. I personally think it would be good.
PLS_effects_raw_m <- plsCoVar(PLS_raw_m, i = 1, sdx = 3, sdy = 3) # i is the latent variable, so using the first component
# Doing x3 here, we can also do x5 no worries

palette(c("black", "red"))
deformGrid3d(PLS_effects_raw_m$x[,,1], PLS_effects_raw_m$x[,,2])##show on x - shape
# Remember
# PLS_raw_sym_m <- Morpho::pls2B(x = GPA_Morpho$Sym, y = proliferation, same.config = FALSE, rounds = 999) # rounds = permutations
# PLS_raw_sym_g <- geomorph::two.b.pls(GPA_Morpho$Sym, proliferation) 


GPA_raw <- readRDS("./data/GPA_Morpho_corrected.rds")
consensus_lm <- GPA_raw$mshape
colnames(consensus_lm) <- colnames(atlas_lm)
row.names(consensus_lm) <- row.names(atlas_lm)
atlas_lm_m <- as.matrix(atlas_lm)
consensus_mesh <- tps3d(x = atlas_mesh, refmat = atlas_lm_m, tarmat = consensus_lm)


# Actually let's do 1 for the heatmaps, it is already pretty strong effects
PLS_effects_raw_m <- plsCoVar(PLS_raw_m, i = 1, sdx = 1, sdy = 1) # i is the latent variable, so using the first component
palette(c("goldenrod", "darkblue"))
deformGrid3d(PLS_effects_raw_m$x[,,1], PLS_effects_raw_m$x[,,2])##show on x
shade3d(consensus_mesh, color = "gray", alpha = 0.4, add = TRUE)
PLS_effects_raw_m_x1_mesh <- tps3d(x = consensus_mesh, refmat = consensus_lm, tarmat = PLS_effects_raw_m$x[,,1])
PLS_effects_raw_m_x2_mesh <- tps3d(x = consensus_mesh, refmat = consensus_lm, tarmat = PLS_effects_raw_m$x[,,2])


frontal <- par3d()$userMatrix
lateral <- par3d()$userMatrix
dorsal <- par3d()$userMatrix
rgl.close()

save(frontal, lateral, dorsal, file = "./data/RGL_head_pos_PLS_shape_effects.rdata")
load("./data/RGL_head_pos_PLS_shape_effects.rdata")

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = lateral)
meshDist(PLS_effects_raw_m_x2_mesh, PLS_effects_raw_m_x1_mesh, rampcolors = c("darkblue", "blue", "white", "red", "darkred"), sign = TRUE)
rgl::rgl.snapshot("./figs/PLS_shape_effects_raw_sym_m_corrected_SIDE.png", top = TRUE)

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = dorsal)
meshDist(PLS_effects_raw_m_x2_mesh, PLS_effects_raw_m_x1_mesh, rampcolors = c("darkblue", "blue", "white", "red", "darkred"), sign = TRUE)
rgl::rgl.snapshot("./figs/PLS_shape_effects_raw_sym_m_corrected_TOP.png", top = TRUE)

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = frontal)
meshDist(PLS_effects_raw_m_x2_mesh, PLS_effects_raw_m_x1_mesh, rampcolors = c("darkblue", "blue", "white", "red", "darkred"), sign = TRUE)
rgl::rgl.snapshot("./figs/PLS_shape_effects_raw_sym_m_corrected_FRONTAL.png", top = TRUE)

# Compute the shape changes along the common axis of deformations
# this give the same results as plsCoVar. However, using common shape vectors as suggested by Mitteroecker and Bookstein (2007)
 # TO DO
rm(PLS_raw_m)

#### 3.2. PLS shape-proliferation - REGRESSED DATA ####
# Remember:
# reg_shape_log.Csize_log.somite <- procD.lm(coords_sym ~ log(Csize) + log(Somites), data = gm_df, RRPP = TRUE)
# GrandMean only

PLS_reg_g <- readRDS("./cluster/PLS_reg_sym_g_corrected.rds")

summary(PLS_reg_g)

# r-PLS: 0.606
# 
# Effect Size (Z): 1.4929
# 
# P-value: 0.075
# 
# Based on 1000 random permutations
PLS_reg_g$r.pls # 0.606285

n_distinct(Somites)
?rainbow

# col_somites <- wes_palette("Zissou1", n_distinct(Somites), type = "continuous")
col_somites <- viridis(n = n_distinct(Somites))

row.names(PLS_reg_g$XScores)
names(Somites)

palette(col_somites)
palette()



# FIX COLOURS
pdf("./figs/PLS_reg_sym_Dec2022_corrected.pdf", width = 8, height = 8)
# layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
Somites <- as.factor(Somites)
plot(PLS_reg_g, pch = 19, col = Somites, cex = 1.75)
# identify(x = PLS_reg_g$XScores[,1], y = PLS_reg_g$YScores[,1], labels = row.names(PLS_reg_g$XScores))
# legend("topleft", pch = 19, col = palette(), legend = levels(Skull_gdf_P0_2$genotype))
# mtext(side=3, line=3, at=0.07, adj=0, cex=1.5, "P0-P2")
# Add a simple legend
# Somites <- as.numeric(as.character(Somites))
# legend("topleft", col=palette()[c(1,length(palette()))], pch=19,
#        legend=c(range(Somites)), title = "Somites")
# mtext("r-PLS = 0.969; z = 2.155; p = 0.001", side = 3)
dev.off()

picknplot.shape(plot(PLS_reg_g))

legend_image <- as.raster(matrix(palette(), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F, xlab = '', ylab = '', main = 'Somites')
text(x=1.5, y = seq(8,20,5), labels = format(seq(8,20,5)), cex = 1.5)
rasterImage(legend_image, 0, 0, 1,1)
dev.off()

rm(PLS_reg_g)

# Morpho
PLS_reg_m <- readRDS("./cluster/PLS_reg_sym_m_corrected.rds")
# Ask Benedikt, Rebecca & Lucas if they want to have morphs (and maybe also the heatmap) showing differences between - YES

# Get the shape changes from pls2B associated with each latent variable
# -3*STD and +3*STD for both x & y, or even 5x. I personally think it would be good.
PLS_effects_reg_m <- plsCoVar(PLS_reg_m, i = 1, sdx = 3, sdy = 3) # i is the latent variable, so using the first component
# Doing x3 here, we can also do x5 no worries

palette(c("black", "red"))
deformGrid3d(PLS_effects_reg_m$x[,,1], PLS_effects_reg_m$x[,,2])##show on x - shape
# Remember
# PLS_reg_sym_m <- Morpho::pls2B(x = GPA_Morpho$Sym, y = proliferation, same.config = FALSE, rounds = 999) # rounds = permutations
# PLS_reg_sym_g <- geomorph::two.b.pls(GPA_Morpho$Sym, proliferation) 


GPA_reg <- readRDS("./data/GPA_Morpho_corrected.rds")
consensus_lm <- GPA_reg$mshape
colnames(consensus_lm) <- colnames(atlas_lm)
row.names(consensus_lm) <- row.names(atlas_lm)
atlas_lm_m <- as.matrix(atlas_lm)
consensus_mesh <- tps3d(x = atlas_mesh, refmat = atlas_lm_m, tarmat = consensus_lm)


# Actually let's do 1 for the heatmaps, it is already pretty strong effects
PLS_effects_reg_m <- plsCoVar(PLS_reg_m, i = 1, sdx = 1, sdy = 1) # i is the latent variable, so using the first component
palette(c("goldenrod", "darkblue"))
deformGrid3d(PLS_effects_reg_m$x[,,1], PLS_effects_reg_m$x[,,2])##show on x
shade3d(consensus_mesh, color = "gray", alpha = 0.4, add = TRUE)
PLS_effects_reg_m_x1_mesh <- tps3d(x = consensus_mesh, refmat = consensus_lm, tarmat = PLS_effects_reg_m$x[,,1])
PLS_effects_reg_m_x2_mesh <- tps3d(x = consensus_mesh, refmat = consensus_lm, tarmat = PLS_effects_reg_m$x[,,2])


load("./data/RGL_head_pos_PLS_shape_effects.rdata")

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = lateral)
meshDist(PLS_effects_reg_m_x2_mesh, PLS_effects_reg_m_x1_mesh, rampcolors = c("darkblue", "blue", "white", "red", "darkred"), sign = TRUE)
rgl::rgl.snapshot("./figs/PLS_shape_effects_reg_sym_m_corrected_SIDE.png", top = TRUE)

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = dorsal)
meshDist(PLS_effects_reg_m_x2_mesh, PLS_effects_reg_m_x1_mesh, rampcolors = c("darkblue", "blue", "white", "red", "darkred"), sign = TRUE)
rgl::rgl.snapshot("./figs/PLS_shape_effects_reg_sym_m_corrected_TOP.png", top = TRUE)

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = frontal)
meshDist(PLS_effects_reg_m_x2_mesh, PLS_effects_reg_m_x1_mesh, rampcolors = c("darkblue", "blue", "white", "red", "darkred"), sign = TRUE)
rgl::rgl.snapshot("./figs/PLS_shape_effects_reg_sym_m_corrected_FRONTAL.png", top = TRUE)

# Compute the shape changes along the common axis of deformations
# this give the same results as plsCoVar. However, using common shape vectors as suggested by Mitteroecker and Bookstein (2007)
# TO DO
rm(PLS_reg_m)


#### 3.3. PLS shape-proliferation - PROJECTED DATA ####
# Remember:
# Projected data
# reg_shape_log.Csize_log.somite <- procD.lm(coords_sym ~ log(Csize) + log(Somites), data = gm_df, RRPP = TRUE)
# mean_shape_p <- reg_shape_log.Csize_log.somite$coefficients[1,] + # GrandMean (intercept)
#   median(log(Somites))*reg_shape_log.Csize_log.somite$coefficients[3,] # because # coefficient 2 is Csize
# residuals_coords_p <- sweep(reg_shape_log.Csize_log.somite$residuals, MARGIN = 2, STATS = mean_shape_p, FUN = "+")
# array_reg_shape_log.Csize_log.somite_p <- arrayspecs(A = residuals_coords_p, p = dim(residuals_coords_p)[2]/3, k = 3)
# GPA_reg_shape_log.Csize_log.somite_p <- procSym(array_reg_shape_log.Csize_log.somite_p)

PLS_proj_g <- readRDS("./cluster/PLS_reg_sym_g_projected_corrected.rds")

summary(PLS_proj_g)
# 
# r-PLS: 0.607
# 
# Effect Size (Z): 1.5074
# 
# P-value: 0.072
# 
# Based on 1000 random permutations

PLS_proj_g$r.pls # 0.6069344

n_distinct(Somites)

# col_somites <- wes_palette("Zissou1", n_distinct(Somites), type = "continuous")
col_somites <- viridis(n = n_distinct(Somites))

row.names(PLS_proj_g$XScores)
names(Somites)

palette(col_somites)
palette()


pdf("./figs/PROJ_sym_PLS/PLS_proj_sym_Dec2022_corrected.pdf", width = 8, height = 8)
# layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
Somites <- as.factor(Somites)
plot(PLS_proj_g, pch = 19, col = Somites, cex = 1.75)
# identify(x = PLS_proj_g$XScores[,1], y = PLS_proj_g$YScores[,1], labels = row.names(PLS_proj_g$XScores))
# legend("topleft", pch = 19, col = palette(), legend = levels(Skull_gdf_P0_2$genotype))
# mtext(side=3, line=3, at=0.07, adj=0, cex=1.5, "P0-P2")
# Add a simple legend
# Somites <- as.numeric(as.character(Somites))
# legend("topleft", col=palette()[c(1,length(palette()))], pch=19,
#        legend=c(range(Somites)), title = "Somites")
mtext("r-PLS = 0.607; z = 1.507; p = 0.072", side = 3)
dev.off()

picknplot.shape(plot(PLS_proj_g))

# legend_image <- as.raster(matrix(palette(), ncol=1))
# plot(c(0,2),c(0,1),type = 'n', axes = F, xlab = '', ylab = '', main = 'Somites')
# text(x=1.5, y = seq(8,20,5), labels = format(seq(8,20,5)), cex = 1.5)
# rasterImage(legend_image, 0, 0, 1,1)
# dev.off()

rm(PLS_proj_g)

# Morpho
PLS_proj_m <- readRDS("./cluster/PLS_reg_sym_m_projected_corrected.rds")
# Ask Benedikt, Rebecca & Lucas if they want to have morphs (and maybe also the heatmap) showing differences between - YES

# Get the shape changes from pls2B associated with each latent variable
# -3*STD and +3*STD for both x & y, or even 5x. I personally think it would be good.
PLS_effects_proj_m <- plsCoVar(PLS_proj_m, i = 1, sdx = 1, sdy = 1) # i is the latent variable, so using the first component
# Doing x3 here, we can also do x5 no worries

palette(c("goldenrod", "darkblue"))
deformGrid3d(PLS_effects_proj_m$x[,,1], PLS_effects_proj_m$x[,,2])##show on x - shape

GPA_proj <- readRDS("./data/GPA_reg_shape_log.Csize_log.somite_projected_corrected.rds")
consensus_lm <- GPA_proj$mshape
colnames(consensus_lm) <- colnames(atlas_lm)
row.names(consensus_lm) <- row.names(atlas_lm)
atlas_lm_m <- as.matrix(atlas_lm)
consensus_mesh <- tps3d(x = atlas_mesh, refmat = atlas_lm_m, tarmat = consensus_lm)


# Remember
# PLS_proj_sym_m <- Morpho::pls2B(x = GPA_Morpho$Sym, y = proliferation, same.config = FALSE, rounds = 999) # rounds = permutations
# PLS_proj_sym_g <- geomorph::two.b.pls(GPA_Morpho$Sym, proliferation) 
deformGrid3d(PLS_effects_proj_m$x[,,1], PLS_effects_proj_m$x[,,2])##show on x
shade3d(consensus_mesh, color = "gray", alpha = 0.4, add = TRUE)
PLS_effects_proj_m_x1_mesh <- tps3d(x = consensus_mesh, refmat = consensus_lm, tarmat = PLS_effects_proj_m$x[,,1])
PLS_effects_proj_m_x2_mesh <- tps3d(x = consensus_mesh, refmat = consensus_lm, tarmat = PLS_effects_proj_m$x[,,2])
meshDist(PLS_effects_proj_m_x1_mesh, PLS_effects_proj_m_x2_mesh, rampcolors = c("darkblue", "blue", "white", "red", "darkred"), sign = TRUE)



# Actually let's do 1 for the heatmaps, it is already pretty strong effects
PLS_effects_proj_m <- plsCoVar(PLS_proj_m, i = 1, sdx = 1, sdy = 1) # i is the latent variable, so using the first component
palette(c("goldenrod", "darkblue"))
deformGrid3d(PLS_effects_proj_m$x[,,1], PLS_effects_proj_m$x[,,2])##show on x
shade3d(consensus_mesh, color = "gray", alpha = 0.4, add = TRUE)
PLS_effects_proj_m_x1_mesh <- tps3d(x = atlas_mesh, refmat = atlas_lm_m, tarmat = PLS_effects_proj_m$x[,,1])
PLS_effects_proj_m_x2_mesh <- tps3d(x = atlas_mesh, refmat = atlas_lm_m, tarmat = PLS_effects_proj_m$x[,,2])
meshDist(PLS_effects_proj_m_x2_mesh, PLS_effects_proj_m_x1_mesh, rampcolors = c("darkblue", "blue", "white", "red", "darkred"), sign = TRUE)

# Compute the shape changes along the common axis of deformations
# this give the same results as plsCoVar. However, using common shape vectors as suggested by Mitteroecker and Bookstein (2007)
# TO DO
rm(PLS_proj_m)



#### 3.4 PLS shape proliferation RAW DATA SUBSET 8-14 ####
somites <- read.csv("./data/somites.csv")
row.names(somites) <- somites$Name
head(somites)


length(somites$Tail.somite[which(somites$Tail.somite <= 14)]) # 20 specs
length(somites$Tail.somite[which(somites$Tail.somite > 14)]) # 20 specs

somites_8_14 <- somites[which(somites$Tail.somite <= 14),]
somites_15_20 <- somites[which(somites$Tail.somite > 14),]

prolif_8_14 <- proliferation[which(somites$Tail.somite <= 14), ]
prolif_15_20 <- proliferation[which(somites$Tail.somite > 14), ]

# get rid of LM 22 (back of head)
land3D_8_14 <- land3D[-22,,which(somites$Tail.somite <= 14)]
land3D_15_20 <- land3D[-22,,which(somites$Tail.somite > 14)]

non.sym <- c(1:5)
side.1 <- c(6:21)
side.2 <- c(22:37)
pairedLM <- cbind(side.1, side.2)

pairedLM

GPA_subset_8_14 <- Morpho::procSym(land3D_8_14, pairedLM = pairedLM)

paste0("Starting PLS (Morpho) between Proliferation and raw SUBSET 8-14 (SYM comp.) shape blocks at ", Sys.time())
PLS_raw_m_subset_8_14 <- Morpho::pls2B(x = GPA_subset_8_14$rotated, y = prolif_8_14, same.config = FALSE, rounds = 999) # rounds = permutations
saveRDS(PLS_raw_m_subset_8_14, "./PLS_raw_sym_m_corrected_subset_8_14.rds")

paste0("Starting PLS (geomorph) between Proliferation and raw SUBSET 8-14 (SYM comp.) shape blocks at ", Sys.time())
PLS_raw_g_subset_8_14 <- geomorph::two.b.pls(GPA_subset_8_14$rotated, prolif_8_14)
saveRDS(PLS_raw_g_subset_8_14, "./PLS_raw_sym_g_corrected_subset_8_14.rds")

paste0("All jobs finished at ", Sys.time())

# Plotting

PLS_raw_g_subset_8_14 <- readRDS("./PLS_raw_sym_g_corrected_subset_8_14.rds")

summary(PLS_raw_g_subset_8_14)

# r-PLS: 0.973
# 
# Effect Size (Z): 3.3413
# 
# P-value: 0.001
# 
# Based on 1000 random permutations

PLS_raw_g_subset_8_14$r.pls # 0.9726031

Somites <- somites$Tail.somite
names(Somites) <- row.names(somites)

n_distinct(Somites)

# col_somites <- wes_palette("Zissou1", n_distinct(Somites), type = "continuous")
col_somites <- viridis(n = n_distinct(Somites))

row.names(PLS_raw_g_subset_8_14$XScores)
names(Somites)


palette(col_somites[1:n_distinct(somites_8_14$Tail.somite)])
palette()


# FIX COLOURS
pdf("./figs/PLS_raw_g_subset_8_14_Feb2023.pdf", width = 8, height = 8)
# layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))

par(cex.lab = 1.25, cex.axis = 1, cex.main = 1.5)
plot(PLS_raw_g_subset_8_14, pch = 19, col = as.factor(somites_8_14$Tail.somite), 
     cex = 1.75)


# plot(PLS_raw_g_subset_8_14$XScores[,1], PLS_raw_g_subset_8_14$YScores[,1], pch = 19, col = as.factor(somites_8_14$Tail.somite), 
#      cex = 1.75, 
#      xlab  = "PLS1 - SHAPE", ylab  = "PLS2 - PROLIFERATION", 
#      main = "PLS shape-proliferation (8-14 somites)")
mtext("r-PLS = 0.973; z = 3.3413; p = 0.001", side = 3)

identify(x = PLS_raw_g_subset_8_14$XScores[,1], y = PLS_raw_g_subset_8_14$YScores[,1], 
         labels = row.names(PLS_raw_g_subset_8_14$XScores))
# legend("topleft", pch = 19, col = palette(), legend = levels(Skull_gdf_P0_2$genotype))
# mtext(side=3, line=3, at=0.07, adj=0, cex=1.5, "P0-P2")
# Add a simple legend
# Somites <- as.numeric(as.character(Somites))
legend("topleft", col = palette(), pch=19,
       legend = levels(as.factor(somites_8_14$Tail.somite)), title = "Somites")

dev.off()

# Morpho
PLS_raw_m_subset_8_14 <- readRDS("./PLS_raw_sym_m_corrected_subset_8_14.rds")
# Ask Benedikt, Rebecca & Lucas if they want to have morphs (and maybe also the heatmap) showing differences between - YES

# Get the shape changes from pls2B associated with each latent variable
# -3*STD and +3*STD for both x & y, or even 5x. I personally think it would be good.
PLS_effects_subset_8_14 <- plsCoVar(PLS_raw_m_subset_8_14, i = 1, sdx = 1, sdy = 1) # i is the latent variable, so using the first component
# Doing x3 here, we can also do x5 no worries

palette(c("goldenrod", "darkblue"))
deformGrid3d(PLS_effects_subset_8_14$x[,,1], PLS_effects_subset_8_14$x[,,2])##show on x - shape

GPA_subset_8_14 <- readRDS("./data/GPA_sym_8_14.rds")
consensus_lm <- GPA_subset_8_14$mshape


atlas_mesh <- geomorph::read.ply("./data/E10.5_atlas_simpl_2_dec_ascii.ply")
# atlas_lm <- morpho.tools.GM::pp2array()

n_land <- length(count.fields("./data/E10.5_atlas_simpl_2_dec_ascii_picked_points.pp")) - 9
raw_LM_file <- readLines("./data/E10.5_atlas_simpl_2_dec_ascii_picked_points.pp")[9:(8 + n_land)]
atlas_lm <- matrix(data = NA, nrow = n_land, ncol = 3)
for (j in 1:length(raw_LM_file)) {
  atlas_lm[j, 1] <- as.numeric(strsplit(strsplit(raw_LM_file[j], 
                                                 "x=\"")[[1]][2], "\"")[[1]][1])
  atlas_lm[j, 2] <- as.numeric(strsplit(strsplit(raw_LM_file[j], 
                                                 "y=\"")[[1]][2], "\"")[[1]][1])
  atlas_lm[j, 3] <- as.numeric(strsplit(strsplit(raw_LM_file[j], 
                                                 "z=\"")[[1]][2], "\"")[[1]][1])
}

atlas_lm <- atlas_lm[-22,]
dim(atlas_lm)
dim(consensus_lm)

colnames(consensus_lm) <- colnames(atlas_lm)
row.names(consensus_lm) <- row.names(atlas_lm)
atlas_lm_m <- as.matrix(atlas_lm)
consensus_mesh <- tps3d(x = atlas_mesh, refmat = atlas_lm_m, tarmat = consensus_lm)


# Remember
# PLS_proj_sym_m <- Morpho::pls2B(x = GPA_Morpho$Sym, y = proliferation, same.config = FALSE, rounds = 999) # rounds = permutations
# PLS_proj_sym_g <- geomorph::two.b.pls(GPA_Morpho$Sym, proliferation) 
shade3d(consensus_mesh, color = "gray", alpha = 0.4, add = TRUE)
PLS_effects_proj_m_x1_mesh <- tps3d(x = consensus_mesh, refmat = consensus_lm, tarmat = PLS_effects_subset_8_14$x[,,1])
PLS_effects_proj_m_x2_mesh <- tps3d(x = consensus_mesh, refmat = consensus_lm, tarmat = PLS_effects_subset_8_14$x[,,2])

load("./data/RGL_head_pos_PLS_shape_effects.rdata")
frontal <- par3d()$userMatrix
lateral <- par3d()$userMatrix
dorsal <- par3d()$userMatrix
rgl.close()

save(frontal, lateral, dorsal, file = "./data/RGL_head_pos_PLS_shape_effects_8_14.rdata")
load("./data/RGL_head_pos_PLS_shape_effects_8_14.rdata")


open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = lateral)
meshDist(PLS_effects_proj_m_x1_mesh, PLS_effects_proj_m_x2_mesh, rampcolors = c("darkred", "red", 
                                                                                "white", "white",
                                                                                "blue", "darkblue"), sign = TRUE)
rgl::rgl.snapshot("./figs/PLS_shape_effects_raw_subset_8_14_SIDE_ammended.png", top = TRUE)

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = dorsal)
meshDist(PLS_effects_proj_m_x1_mesh, PLS_effects_proj_m_x2_mesh, rampcolors = c("darkred", "red", 
                                                                                "white", "white",
                                                                                "blue", "darkblue"), sign = TRUE)
rgl::rgl.snapshot("./figs/PLS_shape_effects_raw_subset_8_14_DORSAL.png", top = TRUE)

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = frontal)
meshDist(PLS_effects_proj_m_x1_mesh, PLS_effects_proj_m_x2_mesh, rampcolors = c("darkred", "red", 
                                                                                "white", "white",
                                                                                "blue", "darkblue"), sign = TRUE)
rgl::rgl.snapshot("./figs/PLS_shape_effects_raw_subset_8_14_FRONTAL.png", top = TRUE)


# Add morphs of selected specimens
plot(PLS_raw_g_subset_8_14, pch = 19, col = as.factor(somites_8_14$Tail.somite), 
     cex = 1.75)
identify(x = PLS_raw_g_subset_8_14$XScores[,1], y = PLS_raw_g_subset_8_14$YScores[,1], 
         labels = row.names(PLS_raw_g_subset_8_14$XScores))

# "Dec2_E10_12_flipped"
row.names(PLS_raw_g_subset_8_14$XScores)[6]
Dec2_E10_12_flipped_LMs <- GPA_subset_8_14$rotated[,, 6]
Dec2_E10_12_flipped_mesh <- tps3d(x = consensus_mesh, refmat = consensus_lm, 
                                  tarmat = Dec2_E10_12_flipped_LMs)

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = frontal)
shade3d(Dec2_E10_12_flipped_mesh, color = "gray", alpha = 0.9, add = TRUE)
rgl::rgl.snapshot("./figs/PLS_shape_effects_Dec2_E10_12_flipped_subset_8_14_FRONTAL.png", top = TRUE)

# "Dec2_E10_9_flipped"
row.names(PLS_raw_g_subset_8_14$XScores)[14]
Dec2_E10_9_flipped_LMs <- GPA_subset_8_14$rotated[,, 14]
Dec2_E10_9_flipped_mesh <- tps3d(x = consensus_mesh, refmat = consensus_lm, 
                                 tarmat = Dec2_E10_9_flipped_LMs)

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = frontal)
shade3d(Dec2_E10_9_flipped_mesh, color = "gray", alpha = 0.9, add = TRUE)
rgl::rgl.snapshot("./figs/PLS_shape_effects_Dec2_E10_9_flipped_subset_8_14_FRONTAL.png", top = TRUE)

# "Oct27_E105_4_flipped"
row.names(PLS_raw_g_subset_8_14$XScores)[20]
Oct27_E105_4_flipped_LMs <- GPA_subset_8_14$rotated[,, 20]
Oct27_E105_4_flipped_mesh <- tps3d(x = consensus_mesh, refmat = consensus_lm, 
                                   tarmat = Oct27_E105_4_flipped_LMs)

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = frontal)
shade3d(Oct27_E105_4_flipped_mesh, color = "gray", alpha = 0.9, add = TRUE)
rgl::rgl.snapshot("./figs/PLS_shape_effects_Oct27_E105_4_flipped_subset_8_14_FRONTAL.png", top = TRUE)

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = frontal)
meshDist(Dec2_E10_12_flipped_mesh, Oct27_E105_4_flipped_mesh, rampcolors = c("darkred", "red", 
                                                                                "white", "white",
                                                                                "blue", "darkblue"), sign = TRUE)
# Checking everyone now, because there is n ewrror in 20
for (i in 1:dim(GPA_subset_8_14$rotated)[3]){
  name <- dimnames(GPA_subset_8_14$rotated)[[3]][i]
  LMs <- GPA_subset_8_14$rotated[,, i]
  mesh <- tps3d(x = consensus_mesh, refmat = consensus_lm, 
                tarmat = LMs)
  open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = frontal)
  shade3d(mesh, color = "gray", alpha = 1, add = TRUE)
  rgl.snapshot(paste0("./figs/LMs_checks/PLS_shape_effects_", name, "_subset_8_14_CHECKING_LMs.png"), top = TRUE)
  rgl.close()
  rm(name, LMs, mesh)
}


#### 3.5 PLS shape proliferation RAW DATA SUBSET 15-20 ####
somites <- read.csv("./data/somites.csv")
row.names(somites) <- somites$Name
head(somites)


length(somites$Tail.somite[which(somites$Tail.somite <= 14)]) # 20 specs
length(somites$Tail.somite[which(somites$Tail.somite > 14)]) # 20 specs

somites_8_14 <- somites[which(somites$Tail.somite <= 14),]
somites_15_20 <- somites[which(somites$Tail.somite > 14),]

load("./data/RGL_head_pos_E11.rdata")
file_LM <- "./data/E11.5_Tissues_Affine_majority_toE10.5_closed_sphere57_ascii_dec_smooth_clean_picked_points_CORRECTED.pp"
e11_atlas_LM <- pp2lm(file = file_LM)
atlas_e11_mesh <- geomorph::read.ply("./data/E11.5_Tissues_Affine_majority_toE10.5_closed_sphere57_ascii_dec_smooth_clean.ply")
atlas_e11_mesh_dec <- Rvcg::vcgQEdecim(atlas_e11_mesh, percent = 0.5)

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = lateral)
rgl::shade3d(atlas_e11_mesh_dec, color = "gray", alpha = 1)
rgl::plot3d(e11_atlas_LM, aspect = "iso", type = "s", size=1.2, col = "darkblue", add = T)
rgl.close()

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = frontal)
rgl::shade3d(atlas_e11_mesh_dec, color = "gray", alpha = 1)
rgl::plot3d(e11_atlas_LM, aspect = "iso", type = "s", size=1.2, col = "darkblue", add = T)
rgl.close()

# frontal <- par3d()$userMatrix
# lateral <- par3d()$userMatrix
# 
# save(frontal, lateral, file = "./data/RGL_head_pos_E11.rdata")

GPA_subset_15_20 <- Morpho::procSym(land3D_15_20 , pairedLM = pairedLM)

paste0("Starting PLS (Morpho) between Proliferation and raw SUBSET 15-20  (SYM comp.) shape blocks at ", Sys.time())
PLS_raw_m_subset_15_20 <- Morpho::pls2B(x = GPA_subset_15_20$rotated, y = prolif_15_20, same.config = FALSE, rounds = 999) # rounds = permutations
saveRDS(PLS_raw_m_subset_15_20 , "./PLS_raw_sym_m_corrected_subset_15_20.rds")

paste0("Starting PLS (geomorph) between Proliferation and raw SUBSET 15-20 (SYM comp.) shape blocks at ", Sys.time())
PLS_raw_g_subset_15_20 <- geomorph::two.b.pls(GPA_subset_15_20$rotated, prolif_15_20)
saveRDS(PLS_raw_g_subset_15_20, "./PLS_raw_sym_g_corrected_subset_15_20.rds")

paste0("All jobs finished at ", Sys.time())

PLS_raw_g_subset_15_20 <- readRDS("./PLS_raw_sym_g_corrected_subset_15_20.rds")

summary(PLS_raw_g_subset_15_20)

# r-PLS: 0.91
# 
# Effect Size (Z): 1.9567
# 
# P-value: 0.024
# 
# Based on 1000 random permutations

PLS_raw_g_subset_15_20$r.pls # 0.9095758

Somites <- somites$Tail.somite
names(Somites) <- row.names(somites)

n_distinct(Somites)

# col_somites <- wes_palette("Zissou1", n_distinct(Somites), type = "continuous")
col_somites <- viridis(n = n_distinct(Somites))

row.names(PLS_raw_g_subset_15_20$XScores)
names(Somites)


palette(col_somites[(1+n_distinct(somites_8_14$Tail.somite)):n_distinct(Somites)])
palette()


# FIX COLOURS
pdf("./figs/PLS_raw_g_subset_15_20_Feb2023.pdf", width = 8, height = 8)
# layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))

par(cex.lab = 1.25, cex.axis = 1, cex.main = 1.5)
plot(PLS_raw_g_subset_15_20, pch = 19, col = as.factor(somites_15_20$Tail.somite), 
     cex = 1.75)


# plot(PLS_raw_g_subset_15_20$XScores[,1], PLS_raw_g_subset_15_20$YScores[,1], pch = 19, col = as.factor(somites_15_20$Tail.somite), 
#      cex = 1.75, 
#      xlab  = "PLS1 - SHAPE", ylab  = "PLS2 - PROLIFERATION", 
#      main = "PLS shape-proliferation (8-14 somites)")
mtext("r-PLS = 0.91; z = 1.9567; p = 0.024", side = 3)

identify(x = PLS_raw_g_subset_15_20$XScores[,1], y = PLS_raw_g_subset_15_20$YScores[,1], 
         labels = row.names(PLS_raw_g_subset_15_20$XScores), pos = 2)
# legend("topleft", pch = 19, col = palette(), legend = levels(Skull_gdf_P0_2$genotype))
# mtext(side=3, line=3, at=0.07, adj=0, cex=1.5, "P0-P2")
# Add a simple legend
# Somites <- as.numeric(as.character(Somites))
legend("topleft", col = palette(), pch=19,
       legend = levels(as.factor(somites_15_20$Tail.somite)), title = "Somites")

dev.off()

# Morpho
PLS_raw_m_subset_15_20 <- readRDS("./PLS_raw_sym_m_corrected_subset_15_20.rds")
# Ask Benedikt, Rebecca & Lucas if they want to have morphs (and maybe also the heatmap) showing differences between - YES

# Get the shape changes from pls2B associated with each latent variable
# -3*STD and +3*STD for both x & y, or even 5x. I personally think it would be good.
PLS_effects_subset_15_20 <- plsCoVar(PLS_raw_m_subset_15_20, i = 1, sdx = 1, sdy = 1) # i is the latent variable, so using the first component
# Doing x3 here, we can also do x5 no worries

palette(c("goldenrod", "darkblue"))
deformGrid3d(PLS_effects_subset_15_20$x[,,1], PLS_effects_subset_15_20$x[,,2])##show on x - shape

GPA_subset_15_20 <- readRDS("./data/GPA_sym_15_20.rds")
consensus_lm <- GPA_subset_15_20$mshape


file_LM <- "./data/E11.5_Tissues_Affine_majority_toE10.5_closed_sphere57_ascii_dec_smooth_clean_picked_points_CORRECTED.pp"
e11_atlas_LM <- pp2lm(file = file_LM)
atlas_e11_mesh <- geomorph::read.ply("./data/E11.5_Tissues_Affine_majority_toE10.5_closed_sphere57_ascii_dec_smooth_clean_noneck.ply")
atlas_e11_mesh_dec <- Rvcg::vcgQEdecim(atlas_e11_mesh, percent = 0.95)

atlas_lm <- e11_atlas_LM[-22,]
dim(atlas_lm)
dim(consensus_lm)

colnames(consensus_lm) <- colnames(atlas_lm)
row.names(consensus_lm) <- row.names(atlas_lm)
atlas_lm_m <- as.matrix(atlas_lm)
consensus_mesh_big <- tps3d(x = atlas_e11_mesh_dec, refmat = atlas_lm_m, tarmat = consensus_lm*mean(GPA_subset_15_20$size))

consensus_mesh_small <- tps3d(x = consensus_mesh_big, refmat = consensus_lm*mean(GPA_subset_15_20$size), tarmat = consensus_lm)



# Remember
# PLS_proj_sym_m <- Morpho::pls2B(x = GPA_Morpho$Sym, y = proliferation, same.config = FALSE, rounds = 999) # rounds = permutations
# PLS_proj_sym_g <- geomorph::two.b.pls(GPA_Morpho$Sym, proliferation) 
deformGrid3d(PLS_effects_subset_15_20$x[,,1], PLS_effects_subset_15_20$x[,,2])##show on x
shade3d(consensus_mesh_small, color = "gray", alpha = 0.9, add = TRUE)
PLS_effects_proj_m_x1_mesh <- tps3d(x = consensus_mesh_small, refmat = consensus_lm, tarmat = PLS_effects_subset_15_20$x[,,1])
PLS_effects_proj_m_x2_mesh <- tps3d(x = consensus_mesh_small, refmat = consensus_lm, tarmat = PLS_effects_subset_15_20$x[,,2])

# load("./data/RGL_head_pos_PLS_shape_effects.rdata")
# frontal <- par3d()$userMatrix
# lateral <- par3d()$userMatrix
# dorsal <- par3d()$userMatrix
# rgl.close()
# 
# save(frontal, lateral, dorsal, file = "./data/RGL_head_pos_PLS_shape_effects_15_20.rdata")
load("./data/RGL_head_pos_PLS_shape_effects_15_20.rdata")


open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = lateral)
meshDist(PLS_effects_proj_m_x1_mesh, PLS_effects_proj_m_x2_mesh, rampcolors = c("darkred", "darkred", "red", 
                                                                                "white", "white",
                                                                                "blue", "blue", "darkblue"), sign = TRUE)
rgl::rgl.snapshot("./figs/PLS_shape_effects_raw_subset_15_20_right_SIDE.png", top = TRUE)

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = dorsal)
meshDist(PLS_effects_proj_m_x1_mesh, PLS_effects_proj_m_x2_mesh, rampcolors = c("darkred", "red", 
                                                                                "white", "white",
                                                                                "blue", "darkblue"), sign = TRUE)
rgl::rgl.snapshot("./figs/PLS_shape_effects_raw_subset_15_20_DORSAL.png", top = TRUE)

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = frontal)
meshDist(PLS_effects_proj_m_x1_mesh, PLS_effects_proj_m_x2_mesh, rampcolors = c("darkred", "darkred", "red", 
                                                                                "white", "white",
                                                                                "blue", "blue", "darkblue"), sign = TRUE)
rgl::rgl.snapshot("./figs/PLS_shape_effects_raw_subset_15_20_FRONTAL.png", top = TRUE)


# Add morphs of selected specimens
# "Feb12_E115_2_flipped" YOUNG
row.names(PLS_raw_g_subset_15_20$XScores)[4]
Feb12_E115_2_flipped_LMs <- GPA_subset_15_20$rotated[,, 4]
Feb12_E115_2_flipped_mesh <- tps3d(x = consensus_mesh_small, refmat = consensus_lm, 
                                  tarmat = Feb12_E115_2_flipped_LMs)

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = frontal)
shade3d(Feb12_E115_2_flipped_mesh, color = "gray", alpha = 1, add = TRUE)
rgl::rgl.snapshot("./figs/PLS_shape_effects_Feb12_E115_2_flipped_subset_15_20_FRONTAL.png", top = TRUE)


# "Feb12_E115_3" MID
row.names(PLS_raw_g_subset_15_20$XScores)[5]
Feb12_E115_3_LMs <- GPA_subset_15_20$rotated[,, 5]
Feb12_E115_3_mesh <- tps3d(x = consensus_mesh_small, refmat = consensus_lm, 
                                   tarmat = Feb12_E115_3_LMs)

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = frontal)
shade3d(Feb12_E115_3_mesh, color = "gray", alpha = 1, add = TRUE)
rgl::rgl.snapshot("./figs/PLS_shape_effects_Feb12_E115_3_subset_15_20_FRONTAL.png", top = TRUE)

# "May10_E11_8_flipped" OLD
row.names(PLS_raw_g_subset_15_20$XScores)[20]
May10_E11_8_flipped_LMs <- GPA_subset_15_20$rotated[,, 20]
May10_E11_8_flipped_mesh <- tps3d(x = consensus_mesh_small, refmat = consensus_lm, 
                           tarmat = May10_E11_8_flipped_LMs)

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = frontal)
shade3d(May10_E11_8_flipped_mesh, color = "gray", alpha = 1, add = TRUE)
rgl::rgl.snapshot("./figs/PLS_shape_effects_May10_E11_8_flipped_subset_15_20_FRONTAL.png", top = TRUE)


# Checking everyone now, because there is n ewrror in 20
for (i in 1:dim(GPA_subset_15_20$rotated)[3]){
  name <- row.names(PLS_raw_g_subset_15_20$XScores)[i]
  LMs <- GPA_subset_15_20$rotated[,, i]
  mesh <- tps3d(x = consensus_mesh_small, refmat = consensus_lm, 
                                    tarmat = LMs)
  open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = frontal)
  shade3d(mesh, color = "gray", alpha = 1, add = TRUE)
  rgl.snapshot(paste0("./figs/LMs_checks/PLS_shape_effects_", name, "_subset_15_20_CHECKING_LMs.png"), top = TRUE)
  rgl.close()
  rm(name, LMs, mesh)
}

#### 4. EXTRA ANALYSES ####

# 4.1. Count how many voxels are positive (they have proliferation) ####

proliferation
prolif_voxels <- vector(mode = "numeric", length = dim(proliferation)[1])

for (i in 1:length(prolif_voxels)){
  prolif_voxels[i] <- length(which(proliferation[i,] > 0))
}

names(prolif_voxels) <- row.names(proliferation)


Somites <- somites$Tail.somite
names(Somites) <- row.names(somites)

# get rid of LM 22 (back of head)

non.sym <- c(1:5,22)
side.1 <- c(6:21)
side.2 <- c(23:38)
pairedLM <- cbind(side.1, side.2)

pairedLM

GPA_all <-gpagen(land3D)

GPA_all$Csize
PCA_Morpho <- gm.prcomp(GPA_Morpho$rotated)
ggplot_prolif_somites <- as.data.frame(cbind(Somites, prolif_voxels, 
                                             GPA_Morpho$size, 
                                             PCA_Morpho$x[,1], PCA_Morpho$x[,2]))
colnames(ggplot_prolif_somites)[3:5] <- c("Csize", "PC1_shape", "PC2_shape")
sapply(ggplot_prolif_somites, class)

col_somites <- viridis(n = n_distinct(Somites))
palette(col_somites)


# Somites vs Proliferation 

plot_prolif_somites <- ggplot(ggplot_prolif_somites, aes(Somites, prolif_voxels, colour = factor(Somites))) + 
  geom_point(alpha=1, size = 3.5) + 
  scale_color_manual(values = col_somites) + 
  theme_classic() +
  ylab("PROLIFERATION (number of positive voxels)") +
  xlab("TAIL SOMITES") +
  ggtitle ("Proliferation through ontogeny") + theme(plot.title = element_text(hjust = -0.15, vjust=2.12))

plot_prolif_somites <- plot_prolif_somites + 
  geom_smooth(method = lm, formula = y ~ poly(x, 3), 
                                se = TRUE, col = "black", lwd = 0.5, lty = 2, alpha = 0.15)

svg("./figs/Proliferation_voxels_somites_Feb2023.svg", width = 7.85, height = 7)
plot_prolif_somites
dev.off()

pdf("./figs/Proliferation_voxels_somites_Feb2023.pdf", width = 7.85, height = 7)
plot_prolif_somites
dev.off()

png("./figs/Proliferation_voxels_somites_Feb2023.png", width = 785, height = 700)
plot_prolif_somites
dev.off()

# Proliferation vs PC1 shape

plot_prolif_PC1_shape <- ggplot(ggplot_prolif_somites, 
                              aes(PC1_shape, prolif_voxels, colour = factor(Somites))) + 
  geom_point(alpha=1, size = 3.5) + 
  scale_color_manual(values = col_somites) + 
  theme_classic() +
  ylab("PROLIFERATION (number of voxels)") +
  xlab("PC1 shape") +
  ggtitle ("Proliferation (# voxels) vs PC1 shape") + theme(plot.title = element_text(hjust = -0.15, vjust=2.12))

plot_prolif_PC1_shape <- plot_prolif_PC1_shape + 
  geom_smooth(method = lm, formula = y ~ poly(x, 3), 
              se = TRUE, col = "black", lwd = 0.5, lty = 2, alpha = 0.15)

svg("./figs/Proliferation_voxels_PC1_shape_Feb2023.svg", width = 7.85, height = 7)
plot_prolif_PC1_shape
dev.off()

pdf("./figs/Proliferation_voxels_PC1_shape_Feb2023.pdf", width = 7.85, height = 7)
plot_prolif_PC1_shape
dev.off()

png("./figs/Proliferation_voxels_PC1_shape_Feb2023.png", width = 785, height = 700)
plot_prolif_PC1_shape
dev.off()

# Proliferation vs PC2 shape

plot_prolif_PC2_shape <- ggplot(ggplot_prolif_somites, 
                              aes(PC2_shape, prolif_voxels, colour = factor(Somites))) + 
  geom_point(alpha=1, size = 3.5) + 
  scale_color_manual(values = col_somites) + 
  theme_classic() +
  ylab("PROLIFERATION (number of voxels)") +
  xlab("PC2 shape") +
  ggtitle ("Proliferation vs PC2 shape") + theme(plot.title = element_text(hjust = -0.15, vjust=2.12))

plot_prolif_PC2_shape 


# Somites vs PC1 shape

plot_somites_PC1_shape <- ggplot(ggplot_prolif_somites, 
                              aes(PC1_shape, Somites, colour = factor(Somites))) + 
  geom_point(alpha=1, size = 3.5) + 
  scale_color_manual(values = col_somites) + 
  theme_classic() +
  ylab("Somites") +
  xlab("PC1 shape") +
  ggtitle ("PC1 shape vs somites") + theme(plot.title = element_text(hjust = -0.15, vjust=2.12))

plot_somites_PC1_shape + 
  geom_smooth(method = lm, formula = y ~ poly(x, 4), 
              se = TRUE, col = "black", lwd = 0.5, lty = 2, alpha = 0.15)


# Somites vs PC2 shape

plot_somites_PC2_shape <- ggplot(ggplot_prolif_somites, 
                              aes(PC2_shape, Somites, colour = factor(Somites))) + 
  geom_point(alpha=1, size = 3.5) + 
  scale_color_manual(values = col_somites) + 
  theme_classic() +
  ylab("Somites") +
  xlab("PC2 shape") +
  ggtitle ("PC2 shape vs somites") + theme(plot.title = element_text(hjust = -0.15, vjust=2.12))

plot_somites_PC2_shape


# Split data into two groups

length(somites$Tail.somite[which(somites$Tail.somite <= 14)]) # 20 specs
length(somites$Tail.somite[which(somites$Tail.somite > 14)]) # 20 specs

somites_8_14 <- somites[which(somites$Tail.somite <= 14),]
somites_15_20 <- somites[which(somites$Tail.somite > 14),]

prolif_8_14 <- proliferation[which(somites$Tail.somite <= 14), ]
prolif_15_20 <- proliferation[which(somites$Tail.somite > 14), ]

# get rid of LM 22 (back of head)
land3D_8_14 <- land3D[-22,,which(somites$Tail.somite <= 14)]
land3D_15_20 <- land3D[-22,,which(somites$Tail.somite > 14)]

somites$Name[which(somites$Tail.somite <= 14)]
row.names(prolif_8_14)
dimnames(land3D_8_14)[[3]]

#### 5. Analyses on two subsets: 8-14 & 15-20 somites separately ####
# 5.1. SHAPE ####
# 5.1.1. Run new GPAs on the two subsets ####
non.sym <- c(1:5)
side.1 <- c(6:21)
side.2 <- c(22:37)
pairedLM <- cbind(side.1, side.2)

pairedLM

GPA_sym_8_14 <- Morpho::procSym(land3D_8_14, pairedLM = pairedLM)
GPA_sym_15_20 <- Morpho::procSym(land3D_15_20, pairedLM = pairedLM)

saveRDS(GPA_sym_8_14, "./data/GPA_sym_8_14.rds")
saveRDS(GPA_sym_15_20, "./data/GPA_sym_15_20.rds")


#### 5.1.2. PCA plots + density - 8-14 somites ####
PCA_8_14 <- gm.prcomp(GPA_sym_8_14$rotated)

PC1_8_14 <- PCA_8_14$x[,1]
PC2_8_14 <- PCA_8_14$x[,2]
PC3_8_14 <- PCA_8_14$x[,3]
PC4_8_14 <- PCA_8_14$x[,4]
PC5_8_14 <- PCA_8_14$x[,5]

summary(PCA_8_14)



scores_8_14 <- as.data.frame(cbind(somites_8_14$Tail.somite, 
                                   PC1_8_14, PC2_8_14, PC3_8_14, PC4_8_14, PC5_8_14))
head(scores_8_14)
sapply(scores_8_14, class)
scores_8_14[2:6] <- sapply(scores_8_14[2:6], as.numeric)
colnames(scores_8_14)[1] <- "Somites"
scores_8_14$Somites <- as.factor(scores_8_14$Somites)
levels(scores_8_14$Somites)

Somites <- somites$Tail.somite
names(Somites) <- row.names(somites)

col_somites <- viridis(n = n_distinct(Somites))



n_distinct(somites_8_14$Tail.somite)
n_distinct(Somites)
palette(col_somites[1:n_distinct(somites_8_14$Tail.somite)])
palette()


# PC1 vs PC2 Shape

plot1_2 <- ggplot(scores_8_14, aes(PC1_8_14, PC2_8_14, colour = factor(somites_8_14$Tail.somite))) + 
  geom_point(alpha=1, size = 3.5) + 
  scale_color_manual(values = col_somites[1:n_distinct(somites_8_14$Tail.somite)]) + 
  # geom_density2d(aes(colour = factor(Somites)), linewidth=0) + 
  # stat_density2d(aes(fill = ..level..), geom="polygon", alpha=0.1)+ylim(-0.08,0.08) + 
  guides(colour=guide_legend(title="Somites")) + 
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme_classic() +
  xlab(paste0("PC1 (", round(PCA_8_14$d[1]/sum(PCA_8_14$d)*100, digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(PCA_8_14$d[2]/sum(PCA_8_14$d)*100, digits = 2), "%)")) +
  ggtitle ("PC1 vs PC2 - Shape 8-14 somites") + theme(plot.title = element_text(hjust = -0.15, vjust=2.12))
plot1_2 <- ggMarginal(plot1_2, scores_8_14, PC1_8_14, PC2_8_14, type = c("density"), 
                      margins = c("both", "x", "y"), size = 3, groupColour = TRUE, groupFill = TRUE)

plot1_2



svg("./figs/Raw_shape_PC1_PC2_8-14_somites.svg", width = 7.85, height = 7)
plot1_2
dev.off()

pdf("./figs/Raw_shape_PC1_PC2_8-14_somites.pdf", width = 7.85, height = 7)
plot1_2
dev.off()

png("./figs/Raw_shape_PC1_PC2_8-14_somites.png", width = 785, height = 700)
plot1_2
dev.off()

#### 5.1.3. PCA plots + density - 15-20 somites ####
PCA_15_20 <- gm.prcomp(GPA_sym_15_20$rotated)

PC1_15_20 <- PCA_15_20$x[,1]
PC2_15_20 <- PCA_15_20$x[,2]
PC3_15_20 <- PCA_15_20$x[,3]
PC4_15_20 <- PCA_15_20$x[,4]
PC5_15_20 <- PCA_15_20$x[,5]

summary(PCA_15_20)



scores_15_20 <- as.data.frame(cbind(somites_15_20$Tail.somite, 
                                    PC1_15_20, PC2_15_20, PC3_15_20, PC4_15_20, PC5_15_20))
head(scores_15_20)
sapply(scores_15_20, class)
scores_15_20[2:6] <- lapply(scores_15_20[2:6], as.numeric)
colnames(scores_15_20)[1] <- "Somites"
scores_15_20$Somites <- as.factor(scores_15_20$Somites)
levels(scores_15_20$Somites)

Somites <- somites$Tail.somite
names(Somites) <- row.names(somites)

col_somites <- viridis(n = n_distinct(Somites))


tail(col_somites, n = n_distinct(somites_15_20$Tail.somite))

n_distinct(somites_15_20$Tail.somite)
n_distinct(Somites)
palette(tail(col_somites, n = n_distinct(somites_15_20$Tail.somite)))
palette()


# PC1 vs PC2 Shape

plot1_2 <- ggplot(scores_15_20, aes(PC1_15_20, PC2_15_20, colour = factor(somites_15_20$Tail.somite))) + 
  geom_point(alpha=1, size = 3.5) + 
  scale_color_manual(values = tail(col_somites, n = n_distinct(somites_15_20$Tail.somite))) + 
  # geom_density2d(aes(colour = factor(Somites)), linewidth=0) + 
  # stat_density2d(aes(fill = ..level..), geom="polygon", alpha=0.1)+ylim(-0.08,0.08) + 
  guides(colour=guide_legend(title="Somites")) + 
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme_classic() +
  xlab(paste0("PC1 (", round(PCA_15_20$d[1]/sum(PCA_15_20$d)*100, digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(PCA_15_20$d[2]/sum(PCA_15_20$d)*100, digits = 2), "%)")) +
  ggtitle ("PC1 vs PC2 - Shape 8-14 somites") + theme(plot.title = element_text(hjust = -0.15, vjust=2.12))
plot1_2 <- ggMarginal(plot1_2, scores_15_20, PC1_15_20, PC2_15_20, type = c("density"), 
                      margins = c("both", "x", "y"), size = 3, groupColour = TRUE, groupFill = TRUE)

plot1_2



svg("./figs/Raw_shape_PC1_PC2_15-20_somites.svg", width = 7.85, height = 7)
plot1_2
dev.off()

pdf("./figs/Raw_shape_PC1_PC2_15-20_somites.pdf", width = 7.85, height = 7)
plot1_2
dev.off()

png("./figs/Raw_shape_PC1_PC2_15-20_somites.png", width = 785, height = 700)
plot1_2
dev.off()




# 5.2. PROLIFERATION - Somites added as a variable ####

# 5.2.1. Proliferation 8-14 ####
prolif_8_14_age <- cbind(somites_8_14$Tail.somite, prolif_8_14)
PCA_prolif_8_14_age <- prcomp(prolif_8_14_age) # cannot use princomp(), as we have more variables than specimens

summary(PCA_prolif_8_14_age) # pretty good actually! 22.19% PC1, 11.59% PC2

PC1_prolif_8_14_age <- PCA_prolif_8_14_age$x[,1]
PC2_prolif_8_14_age <- PCA_prolif_8_14_age$x[,2]
PC3_prolif_8_14_age <- PCA_prolif_8_14_age$x[,3]
PC4_prolif_8_14_age <- PCA_prolif_8_14_age$x[,4]
PC5_prolif_8_14_age <- PCA_prolif_8_14_age$x[,5]


scores <- as.data.frame(cbind(somites_8_14$Tail.somite, 
                              PC1_8_14, PC2_8_14, PC3_8_14, PC4_8_14, PC5_8_14, 
                              PC1_prolif_8_14_age, PC2_prolif_8_14_age, PC3_prolif_8_14_age, PC4_prolif_8_14_age, PC5_prolif_8_14_age))
head(scores)
sapply(scores, class)
scores[2:11] <- sapply(scores[2:11], as.numeric)
colnames(scores)[1] <- "Somites"
scores$Somites <- as.factor(scores$Somites)
levels(scores$Somites)


Somites <- somites$Tail.somite
names(Somites) <- row.names(somites)

col_somites <- viridis(n = n_distinct(Somites))



n_distinct(somites_8_14$Tail.somite)
n_distinct(Somites)
palette(col_somites[1:n_distinct(somites_8_14$Tail.somite)])
palette()

# PC1 vs PC2 Proliferation (+ somites variable)

plot1_2_prolif_8_14_age <- ggplot(scores, aes(PC1_prolif_8_14_age, PC2_prolif_8_14_age, colour = factor(somites_8_14$Tail.somite))) + 
  geom_point(alpha=1, size = 3.5) + 
  scale_color_manual(values = col_somites[1:n_distinct(somites_8_14$Tail.somite)]) + 
  # geom_density2d(aes(colour = factor(Somites)), linewidth=0) + 
  # stat_density2d(aes(fill = ..level..), geom="polygon", alpha=0.1)+ylim(-0.08,0.08) + 
  guides(colour=guide_legend(title="Somites")) + 
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme_classic() +
  xlab(paste0("PC1 (", round(PCA_prolif_8_14_age$sdev[1]^2/sum(PCA_prolif_8_14_age$sdev^2)*100, digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(PCA_prolif_8_14_age$sdev[2]^2/sum(PCA_prolif_8_14_age$sdev^2)*100, digits = 2), "%)")) +
  ggtitle ("PC1 vs PC2 - Proliferation (+ somites variable)") + theme(plot.title = element_text(hjust = -0.15, vjust=2.12))
plot1_2Prolif <- ggMarginal(plot1_2_prolif_8_14_age, scores, 
                            PC1_prolif_8_14_age, PC2_prolif_8_14_age, type = c("density"), margins = c("both", "x", "y"), 
                            size = 3, groupColour = TRUE, groupFill = TRUE)

plot1_2Prolif



svg("./figs/Proliferation_PC1_PC2_8_14.svg", width = 7.85, height = 7)
plot1_2Prolif
dev.off()

pdf("./figs/Proliferation_PC1_PC2_8_14.pdf", width = 7.85, height = 7)
plot1_2Prolif
dev.off()

png("./figs/Proliferation_PC1_PC2_8_14.png", width = 785, height = 700)
plot1_2Prolif
dev.off()


# PC1 shape vs PC1 Proliferation (+ somites variable)

plot1_2_prolif_8_14_age <- ggplot(scores, aes(PC1_8_14, PC1_prolif_8_14_age, colour = factor(somites_8_14$Tail.somite))) + 
  geom_point(alpha=1, size = 3.5) + 
  scale_color_manual(values = col_somites[1:n_distinct(somites_8_14$Tail.somite)]) + 
  # geom_density2d(aes(colour = factor(Somites)), linewidth=0) + 
  # stat_density2d(aes(fill = ..level..), geom="polygon", alpha=0.1)+ylim(-0.08,0.08) + 
  guides(colour=guide_legend(title="Somites")) + 
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme_classic() +
  xlab(paste0("PC1 shape (", round(PCA_8_14$d[1]/sum(PCA_8_14$d)*100, digits = 2), "%)")) +
  ylab(paste0("PC1 proliferation (", round(PCA_prolif_8_14_age$sdev[1]^2/sum(PCA_prolif_8_14_age$sdev^2)*100, digits = 2), "%)")) +
  ggtitle ("PC1 shape vs PC1 proliferation") + theme(plot.title = element_text(hjust = -0.15, vjust=2.12))
plot1_1Prolif_shape <- ggMarginal(plot1_2_prolif_8_14_age, scores, 
                                  PC1_prolif_8_14_age, PC2_prolif_8_14_age, type = c("density"), margins = c("both", "x", "y"), 
                                  size = 3, groupColour = TRUE, groupFill = TRUE)

plot1_1Prolif_shape



svg("./figs/PC1_shape_PC1_prolif_8_14.svg", width = 7.85, height = 7)
plot1_1Prolif_shape
dev.off()

pdf("./figs/PC1_shape_PC1_prolif_8_14.pdf", width = 7.85, height = 7)
plot1_1Prolif_shape
dev.off()

png("./figs/PC1_shape_PC1_prolif_8_14.png", width = 785, height = 700)
plot1_1Prolif_shape
dev.off()

# PC2 shape vs PC2 Proliferation (+ somites variable)

plot1_2_prolif_8_14_age <- ggplot(scores, aes(PC2_8_14, PC2_prolif_8_14_age, colour = factor(somites_8_14$Tail.somite))) + 
  geom_point(alpha=1, size = 3.5) + 
  scale_color_manual(values = col_somites[1:n_distinct(somites_8_14$Tail.somite)]) + 
  # geom_density2d(aes(colour = factor(Somites)), linewidth=0) + 
  # stat_density2d(aes(fill = ..level..), geom="polygon", alpha=0.1)+ylim(-0.08,0.08) + 
  guides(colour=guide_legend(title="Somites")) + 
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme_classic() +
  xlab(paste0("PC2 shape (", round(PCA_8_14$d[2]/sum(PCA_8_14$d)*100, digits = 2), "%)")) +
  ylab(paste0("PC2 proliferation (", round(PCA_prolif_8_14_age$sdev[2]^2/sum(PCA_prolif_8_14_age$sdev^2)*100, digits = 2), "%)")) +
  ggtitle ("PC2 shape vs PC2 proliferation") + theme(plot.title = element_text(hjust = -0.15, vjust=2.12))
plot2_2Prolif_shape <- ggMarginal(plot1_2_prolif_8_14_age, scores, 
                                  PC1_prolif_8_14_age, PC2_prolif_8_14_age, type = c("density"), margins = c("both", "x", "y"), 
                                  size = 3, groupColour = TRUE, groupFill = TRUE)

plot2_2Prolif_shape



svg("./figs/PC2_shape_PC2_prolif_8_14.svg", width = 7.85, height = 7)
plot2_2Prolif_shape
dev.off()

pdf("./figs/PC2_shape_PC2_prolif_8_14.pdf", width = 7.85, height = 7)
plot2_2Prolif_shape
dev.off()

png("./figs/PC2_shape_PC2_prolif_8_14.png", width = 785, height = 700)
plot2_2Prolif_shape
dev.off()



# 5.2.2. Proliferation 15-20 ####
prolif_15_20_age <- cbind(somites_15_20$Tail.somite, prolif_15_20)
PCA_prolif_15_20_age <- prcomp(prolif_15_20_age) # cannot use princomp(), as we have more variables than specimens

summary(PCA_prolif_15_20_age) # pretty good actually! 22.19% PC1, 11.59% PC2

PC1_prolif_15_20_age <- PCA_prolif_15_20_age$x[,1]
PC2_prolif_15_20_age <- PCA_prolif_15_20_age$x[,2]
PC3_prolif_15_20_age <- PCA_prolif_15_20_age$x[,3]
PC4_prolif_15_20_age <- PCA_prolif_15_20_age$x[,4]
PC5_prolif_15_20_age <- PCA_prolif_15_20_age$x[,5]


scores <- as.data.frame(cbind(somites_15_20$Tail.somite, 
                              PC1_15_20, PC2_15_20, PC3_15_20, PC4_15_20, PC5_15_20, 
                              PC1_prolif_15_20_age, PC2_prolif_15_20_age, PC3_prolif_15_20_age, PC4_prolif_15_20_age, PC5_prolif_15_20_age))
head(scores)
sapply(scores, class)
scores[2:11] <- sapply(scores[2:11], as.numeric)
colnames(scores)[1] <- "Somites"
scores$Somites <- as.factor(scores$Somites)
levels(scores$Somites)


Somites <- somites$Tail.somite
names(Somites) <- row.names(somites)

col_somites <- viridis(n = n_distinct(Somites))



n_distinct(somites_15_20$Tail.somite)
n_distinct(Somites)
palette(col_somites[1:n_distinct(somites_15_20$Tail.somite)])
palette()

# PC1 vs PC2 Proliferation (+ somites variable)

plot1_2_prolif_15_20_age <- ggplot(scores, aes(PC1_prolif_15_20_age, PC2_prolif_15_20_age, colour = factor(somites_15_20$Tail.somite))) + 
  geom_point(alpha=1, size = 3.5) + 
  scale_color_manual(values = tail(col_somites, n = n_distinct(somites_15_20$Tail.somite))) + 
  # geom_density2d(aes(colour = factor(Somites)), linewidth=0) + 
  # stat_density2d(aes(fill = ..level..), geom="polygon", alpha=0.1)+ylim(-0.08,0.08) + 
  guides(colour=guide_legend(title="Somites")) + 
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme_classic() +
  xlab(paste0("PC1 (", round(PCA_prolif_15_20_age$sdev[1]^2/sum(PCA_prolif_15_20_age$sdev^2)*100, digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(PCA_prolif_15_20_age$sdev[2]^2/sum(PCA_prolif_15_20_age$sdev^2)*100, digits = 2), "%)")) +
  ggtitle ("PC1 vs PC2 - Proliferation (+ somites variable)") + theme(plot.title = element_text(hjust = -0.15, vjust=2.12))
plot1_2Prolif <- ggMarginal(plot1_2_prolif_15_20_age, scores, 
                            PC1_prolif_15_20_age, PC2_prolif_15_20_age, type = c("density"), margins = c("both", "x", "y"), 
                            size = 3, groupColour = TRUE, groupFill = TRUE)

plot1_2Prolif



svg("./figs/Proliferation_PC1_PC2_15_20.svg", width = 7.85, height = 7)
plot1_2Prolif
dev.off()

pdf("./figs/Proliferation_PC1_PC2_15_20.pdf", width = 7.85, height = 7)
plot1_2Prolif
dev.off()

png("./figs/Proliferation_PC1_PC2_15_20.png", width = 785, height = 700)
plot1_2Prolif
dev.off()


# PC1 shape vs PC1 Proliferation (+ somites variable)

plot1_2_prolif_15_20_age <- ggplot(scores, aes(PC1_15_20, PC1_prolif_15_20_age, colour = factor(somites_15_20$Tail.somite))) + 
  geom_point(alpha=1, size = 3.5) + 
  scale_color_manual(values = tail(col_somites, n = n_distinct(somites_15_20$Tail.somite))) + 
  # geom_density2d(aes(colour = factor(Somites)), linewidth=0) + 
  # stat_density2d(aes(fill = ..level..), geom="polygon", alpha=0.1)+ylim(-0.08,0.08) + 
  guides(colour=guide_legend(title="Somites")) + 
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme_classic() +
  xlab(paste0("PC1 shape (", round(PCA_15_20$d[1]/sum(PCA_15_20$d)*100, digits = 2), "%)")) +
  ylab(paste0("PC1 proliferation (", round(PCA_prolif_15_20_age$sdev[1]^2/sum(PCA_prolif_15_20_age$sdev^2)*100, digits = 2), "%)")) +
  ggtitle ("PC1 shape vs PC1 proliferation") + theme(plot.title = element_text(hjust = -0.15, vjust=2.12))
plot1_1Prolif_shape <- ggMarginal(plot1_2_prolif_15_20_age, scores, 
                                  PC1_prolif_15_20_age, PC2_prolif_15_20_age, type = c("density"), margins = c("both", "x", "y"), 
                                  size = 3, groupColour = TRUE, groupFill = TRUE)

plot1_1Prolif_shape



svg("./figs/PC1_shape_PC1_prolif_15_20.svg", width = 7.85, height = 7)
plot1_1Prolif_shape
dev.off()

pdf("./figs/PC1_shape_PC1_prolif_15_20.pdf", width = 7.85, height = 7)
plot1_1Prolif_shape
dev.off()

png("./figs/PC1_shape_PC1_prolif_15_20.png", width = 785, height = 700)
plot1_1Prolif_shape
dev.off()

# PC2 shape vs PC2 Proliferation (+ somites variable)

plot1_2_prolif_15_20_age <- ggplot(scores, aes(PC2_15_20, PC2_prolif_15_20_age, colour = factor(somites_15_20$Tail.somite))) + 
  geom_point(alpha=1, size = 3.5) + 
  scale_color_manual(values = tail(col_somites, n = n_distinct(somites_15_20$Tail.somite))) + 
  # geom_density2d(aes(colour = factor(Somites)), linewidth=0) + 
  # stat_density2d(aes(fill = ..level..), geom="polygon", alpha=0.1)+ylim(-0.08,0.08) + 
  guides(colour=guide_legend(title="Somites")) + 
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme_classic() +
  xlab(paste0("PC2 shape (", round(PCA_15_20$d[2]/sum(PCA_15_20$d)*100, digits = 2), "%)")) +
  ylab(paste0("PC2 proliferation (", round(PCA_prolif_15_20_age$sdev[2]^2/sum(PCA_prolif_15_20_age$sdev^2)*100, digits = 2), "%)")) +
  ggtitle ("PC2 shape vs PC2 proliferation") + theme(plot.title = element_text(hjust = -0.15, vjust=2.12))
plot2_2Prolif_shape <- ggMarginal(plot1_2_prolif_15_20_age, scores, 
                                  PC1_prolif_15_20_age, PC2_prolif_15_20_age, type = c("density"), margins = c("both", "x", "y"), 
                                  size = 3, groupColour = TRUE, groupFill = TRUE)

plot2_2Prolif_shape



svg("./figs/PC2_shape_PC2_prolif_15_20.svg", width = 7.85, height = 7)
plot2_2Prolif_shape
dev.off()

pdf("./figs/PC2_shape_PC2_prolif_15_20.pdf", width = 7.85, height = 7)
plot2_2Prolif_shape
dev.off()

png("./figs/PC2_shape_PC2_prolif_15_20.png", width = 785, height = 700)
plot2_2Prolif_shape
dev.off()

# compare the variance-covariance matrix between (1) 8-14 somites and 
# (2) 15-20 somites.
# We would expect the covariance to be positive.
# We can look at the negative covariation within the matrix & see when & where it arises


cov_8_14 <- cov(prolif_8_14)
# Error: vector memory exhausted (limit reached?)
prolif_15_20 <- proliferation[which(somites$Tail.somite > 14), ]


# 5.3. SHAPE COVARIANCE  ####

# 5.3.1. variance covariance comparisons between the two ages

GPA_sym_8_14_coords <- two.d.array(GPA_sym_8_14$rotated)
dimnames(GPA_sym_8_14$rotated)[[3]]
row.names(GPA_sym_8_14_coords)

GPA_sym_15_20_coords <- two.d.array(GPA_sym_15_20$rotated)
dimnames(GPA_sym_15_20$rotated)[[3]]
row.names(GPA_sym_15_20_coords)

matrix_all_shape <- rbind(GPA_sym_8_14_coords, GPA_sym_15_20_coords)
somites_all_shape <- rbind(somites_8_14, somites_15_20)

somites_all_shape$Group[which(somites_all_shape$Tail.somite <= 14)] <- "8_14"
somites_all_shape$Group[which(somites_all_shape$Tail.somite > 14)] <- "15_20"
somites_all_shape$Group <- as.factor(somites_all_shape$Group)


# Covariance matrix of each age group
CVC_shape <- cov.group(matrix_all_shape, groups = somites_all_shape$Group)




prop.vcv.test(n = c(20, 20), CVC_shape[,,"8_14"], CVC_shape[,,"15_20"]) # not good sample size
# Warning message:
#   In prop.vcv.test(n = c(20, 20), CVC_shape[, , "8_14"], CVC_shape[,  :
#   The sample size is not very large compared to the number of relative eigenvalues.

# So perhaps what we could focus on here is to see which components differ between the two groups
covW(matrix_all_shape, groups = somites_all_shape$Group)

# Try with PCA
summary(PCA_8_14) # first 5 PCs explain 68.33% variance
summary(PCA_15_20) # first 5 PCs explain 65.74% variance

matrix_PCA_shape <- rbind(PCA_8_14$x[,1:5], PCA_15_20$x[,1:5])
CVC_shape <- cov.group(matrix_PCA_shape, groups = somites_all_shape$Group)
prop.vcv.test(n = c(50, 50), CVC_shape[,,"8_14"], CVC_shape[,,"15_20"]) # not good sample size
# I had to fake the sample size from 20 to 50

group <- somites$Tail.somite
group[which(somites$Tail.somite <= 14)] <- "8_14"
group[which(somites$Tail.somite > 14)] <- "15_20"


covW(as.matrix(PCA$x[,1:5]), groups = as.factor(group))

CVC_shape <- cov.group(matrix_PCA_shape, groups = somites_all_shape$Group)

cov(PCA_8_14$x[,1:5])
cov(PCA_15_20$x[,1:5])



# 5.4. PROLIFERATION COVARIANCE  ####
# 5.4.1. CORRELATION GROWTH WITH PROLIFERATION ####
PCA_prolif <- prcomp(proliferation) 
cor_prolif_age <- cor(as.numeric(somites$Tail.somite),as.matrix(PCA_prolif$x))

summary(cor_prolif_age)
# PC1              PC2               PC3               PC4                PC5              PC6               PC7               PC8          
# Min.   :0.9382   Min.   :0.05814   Min.   :0.06865   Min.   :-0.08238   Min.   :0.2076   Min.   :-0.1373   Min.   :0.05949   Min.   :-0.02347  
# 1st Qu.:0.9382   1st Qu.:0.05814   1st Qu.:0.06865   1st Qu.:-0.08238   1st Qu.:0.2076   1st Qu.:-0.1373   1st Qu.:0.05949   1st Qu.:-0.02347  
# Median :0.9382   Median :0.05814   Median :0.06865   Median :-0.08238   Median :0.2076   Median :-0.1373   Median :0.05949   Median :-0.02347  
# Mean   :0.9382   Mean   :0.05814   Mean   :0.06865   Mean   :-0.08238   Mean   :0.2076   Mean   :-0.1373   Mean   :0.05949   Mean   :-0.02347  
# 3rd Qu.:0.9382   3rd Qu.:0.05814   3rd Qu.:0.06865   3rd Qu.:-0.08238   3rd Qu.:0.2076   3rd Qu.:-0.1373   3rd Qu.:0.05949   3rd Qu.:-0.02347  
# Max.   :0.9382   Max.   :0.05814   Max.   :0.06865   Max.   :-0.08238   Max.   :0.2076   Max.   :-0.1373   Max.   :0.05949   Max.   :-0.02347

# 4.2. COVARIANCE PCs ####
# So now I am looking at which PCs have a negative covariance with PC1
covariance_prolif_PCs <- cov(PCA_prolif$x)

neg_PCs_cov <- which(covariance_prolif_PCs[,1] < 0)

# Display PCs from more negative
covariance_prolif_PCs[neg_PCs_cov[order(covariance_prolif_PCs[neg_PCs_cov,1], decreasing = FALSE)],1]

# PC5          PC37          PC15          PC10          PC14          PC36          PC13          PC39          PC34          PC19 
# -8.947302e-10 -6.688379e-10 -6.564299e-10 -5.431958e-10 -5.153979e-10 -4.344479e-10 -4.284663e-10 -4.271620e-10 -3.522512e-10 -3.489692e-10 
# PC29          PC20          PC27          PC25          PC18          PC40          PC33           PC9 
# -2.337671e-10 -2.264332e-10 -1.382544e-10 -7.741168e-11 -7.224664e-11 -4.248393e-11 -8.405715e-12 -4.374831e-12 

# Plot PC5

# Well that is great, but not very biologically meaningful. Need to check directly the proliferation matrix
# And of course the first variable is number of somites
prolif_clean_cov <- cbind(as.numeric(somites$Tail.somite), proliferation)

prolif_clean_cov[1:10, 1:10]


covariance_prolif_raw <- cov(prolif_clean_cov)
# Error: vector memory exhausted (limit reached?)
# Ugh. Too big as I expected.
dim(prolif_clean_cov)[2] # 1889479

# Check a small subset just for fun (n = 10000)
covariance_prolif_raw_1_10000 <- cov(prolif_clean_cov[,1:10000])

min(covariance_prolif_raw_1_10000[,1]) # -0.6791252

neg_raw_cov <- which(covariance_prolif_raw_1_10000[,1] < 0)

covariance_prolif_raw_1_10000[neg_raw_cov[order(covariance_prolif_raw_1_10000[neg_raw_cov,1], decreasing = FALSE)],1]


covariance_prolif_raw_100001_20000 <- cov(prolif_clean_cov[,c(1, 10001:20000)])

min(covariance_prolif_raw_100001_20000[,1]) # -0.7906083

# So ideally we would do it for the whole proliferation matrix, then check which 
# are the variables that have the smallest numbers (most negative)

# That should indicate that there is no proliferation in those areas overtime.

