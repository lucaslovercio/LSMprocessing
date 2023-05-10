
### PROLIFERATION IMPORT ####
proliferation.import <- function(dir = NULL, ID = NULL){
  if (is.null(dir) == TRUE) {
    path <- getwd()
  }
  else {
    path <- dir
  }
  proliferation_list <- dir(pattern = "*.mat")
  n <- length(R.matlab::readMat(proliferation_list[1])$valuesProliferation)
  proliferation_matrix <- matrix(data = NA, nrow = length(proliferation_list), ncol = n)
  row.names(proliferation_matrix) <- gsub(pattern = "_proliferationValuesFront_mask03_downsample05.mat", replacement = "", x = proliferation_list)
  n_fields <- vector("numeric", length = length(proliferation_list))
  for (i in 1:length(proliferation_list)){
    proliferation_matrix[i,] <- R.matlab::readMat(proliferation_list[i])$valuesProliferation
  }
  if (length(unique(n_fields)) != 1) {
    stop("Specimens have different number of proliferation rows. Re-check the data.")
  }
  return(proliferation_matrix)
}

