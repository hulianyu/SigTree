#### Load ODS files ####
# Get the path of the current working directory
my_dir <- getwd()
# Get the paths of the subdirectories
sub_dirs <- c("ODS")
# Get the paths of TXT files in the subdirectories
txt_paths <- list.files(sub_dirs, pattern = "\\.txt$", full.names = TRUE)
# Read all TXT files
for (txt_path in txt_paths) {
  file_name <- gsub(".txt$", "", basename(txt_path))
  assign(file_name,  read.csv(txt_path, header=FALSE, colClasses="numeric"))
}

# [1]lenses [2]lung_cancer [3]soybean_small [4]zoo [5]dna_promoter [6]hayes_roth [7]lymphography [8]heart_disease [9]solar_flare
# [10]primary_tumor [11]dermatology [12]house_votes [13]balance_scale [14]credit_approval 
# [15]breast_cancer_wisconsin [16]mammographic_mass [17]tic_tac_toe [18]car
data_obj_list <-c("lenses", "lung_cancer", "soybean_small", "zoo", "dna_promoter", "hayes_roth",
                  "lymphography", "heart_disease", "solar_flare", "primary_tumor", "dermatology", "house_votes", "balance_scale", "credit_approval",
                  "breast_cancer_wisconsin", "mammographic_mass", "tic_tac_toe", "car")
ODSs <- paste("ODS_", data_obj_list, sep = "")

library(R.matlab)
source('CUBT_MI.R')
CUBT_MI_pi <- list()
CUBT_MI_times <- numeric(length(data_obj_list))
CUBT_MI_Depth <- matrix(0, nrow = 18, ncol = 3)

for (i in 1:18){
data <- get(ODSs[i])
print(i)
data[] <- lapply(data, as.character)
data <- as.matrix(data)

start_time <- proc.time()
# construct the maximal tree
my_cubt_maximal = CUBT_Cat_growth(data)
# prune it 
my_cubt_prune = CUBT_Cat_prune(my_cubt_maximal,data)
# join leaves
my_cubt_join = CUBT_Cat_join(my_cubt_prune,data)
# assign cluster labels to leaf nodes
pi_Node <- rep(NA, nrow(data))
leaf_node = rownames(my_cubt_join[["frame"]])[my_cubt_join[["frame"]][,1] == "<leaf>"]
K = length(leaf_node)
for (k in 1:K){
  pi_Node[my_cubt_join[["who"]][[leaf_node[k]]]] = k
}
end_time <- proc.time()

## plot the final tree
# Plot.CUBT(my_cubt_join, type="uniform")
# Text.CUBT(my_cubt_join)

# save partitions, running times
CUBT_MI_pi[[data_obj_list[i]]] <- pi_Node
CUBT_MI_times[i] <- end_time[3] - start_time[3]
# save Leaf number, Tree Depth, Avg Leaf Depth
CUBT_MI_Depth[i,1] <- K
leaf_growth_depths = tree.depth(as.numeric(row.names(my_cubt_join$frame)))[my_cubt_join[["frame"]][,1] == "<leaf>"]
CUBT_MI_Depth[i,2] <- max(leaf_growth_depths)
num_growth_leaf = length(leaf_growth_depths)
CUBT_MI_Depth[i,3] <- sum(leaf_growth_depths)/num_growth_leaf
}
writeMat(con = "CUBT_MI_pi.mat", CUBT_MI_pi = CUBT_MI_pi)
writeMat(con = "CUBT_MI_times.mat", CUBT_MI_times = CUBT_MI_times)
writeMat(con = "CUBT_MI_Depth.mat", CUBT_MI_Depth = CUBT_MI_Depth)