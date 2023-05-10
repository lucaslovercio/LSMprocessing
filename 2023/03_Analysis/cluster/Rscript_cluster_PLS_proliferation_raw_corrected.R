# I added this in case you have multiple packages and you are not sure if they have been installed in the cluster
packages_needed <- c("geomorph", "Morpho")
packages_to_install <- packages_needed[!(packages_needed %in% installed.packages()[,"Package"])]
if(length(packages_to_install)) install.packages(packages_to_install)

# Beware X11 cannot be loaded in the cluster
suppressPackageStartupMessages(suppressWarnings(library(geomorph)))
suppressPackageStartupMessages(suppressWarnings(library(Morpho)))

paste0("Starting to load Proliferation file at ", Sys.time())

# 1. Raw data ####
GPA_Morpho <- readRDS("./GPA_Morpho_corrected.rds")
proliferation_raw <- readRDS("./proliferation.rds")
proliferation <- proliferation_raw[match(dimnames(GPA_Morpho$rotated)[[3]], row.names(proliferation_raw)),]

# GPA_Morpho$rotated
paste0("Starting PLS (Morpho) between Proliferation and raw (SYM comp.) shape blocks at ", Sys.time())
PLS_raw_sym_m <- Morpho::pls2B(x = GPA_Morpho$rotated, y = proliferation, same.config = FALSE, rounds = 999) # rounds = permutations
saveRDS(PLS_raw_sym_m, "./PLS_raw_sym_m_corrected.rds")

paste0("Starting PLS (geomorph) between Proliferation and raw (SYM comp.) shape blocks at ", Sys.time())
PLS_raw_sym_g <- geomorph::two.b.pls(GPA_Morpho$rotated, proliferation)
saveRDS(PLS_raw_sym_g, "./PLS_raw_sym_g_corrected.rds")

paste0("All jobs finished at ", Sys.time())
