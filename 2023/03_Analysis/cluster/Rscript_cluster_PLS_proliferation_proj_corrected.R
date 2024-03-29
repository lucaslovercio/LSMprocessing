# I added this in case you have multiple packages and you are not sure if they have been installed in the cluster
packages_needed <- c("geomorph", "Morpho")
packages_to_install <- packages_needed[!(packages_needed %in% installed.packages()[,"Package"])]
if(length(packages_to_install)) install.packages(packages_to_install)

# Beware X11 cannot be loaded in the cluster
suppressPackageStartupMessages(suppressWarnings(library(geomorph)))
suppressPackageStartupMessages(suppressWarnings(library(Morpho)))

paste0("Starting to load Proliferation file at ", Sys.time())

# 3. Somite-regressed - PROJECTED TO MEDIAN SOMITE ####
GPA_reg_shape_log.Csize_log.somite_p <- readRDS("./GPA_reg_shape_log.Csize_log.somite_projected_corrected.rds")
proliferation_raw <- readRDS("./proliferation.rds")
proliferation <- proliferation_raw[match(dimnames(GPA_reg_shape_log.Csize_log.somite_p$rotated)[[3]], row.names(proliferation_raw)),]


# GPA_reg_shape_log.Csize_log.somite$rotated

paste0("Starting PLS (Morpho) between Proliferation and regressed-projected shape blocks at ", Sys.time())

PLS_reg_sym_p_m <- Morpho::pls2B(x = GPA_reg_shape_log.Csize_log.somite_p$rotated, y = proliferation, 
                                 same.config = FALSE, rounds = 999) # rounds = permutations
saveRDS(PLS_reg_sym_p_m, "./PLS_reg_sym_m_projected_corrected.rds")

paste0("Starting PLS (geomorph) between Proliferation and regressed-projected shape blocks at ", Sys.time())
PLS_reg_sym_p_g <- geomorph::two.b.pls(GPA_reg_shape_log.Csize_log.somite_p$rotated, proliferation)

saveRDS(PLS_reg_sym_p_g, "./PLS_reg_sym_g_projected_corrected.rds")

paste0("All jobs finished at ", Sys.time())
