load("intval_n500_train5_rawdist_fpg2024.RData")
vr2 <- valid_res
load("intval_n500_train5_rawdist_2024.RData")
for (i in seq_along(valid_res)){
valid_res[[i]][,16:20] <- vr2[[i]]  
}
save(valid_res,
     file="intval_n500_train5_rawdist_2024.RData")

load("intval_n500_train5_rawdist_fpg2025.RData")
vr2 <- valid_res
load("intval_n500_train5_rawdist_2025.RData")
for (i in seq_along(valid_res)){
  valid_res[[i]][,19:24] <- vr2[[i]]  
}
save(valid_res,
     file="intval_n500_train5_rawdist_2025.RData")
