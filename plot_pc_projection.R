# adjust study PCs to same scale as the ref.eigenvec
# to adjust, computedPCs <- projectedPC/(-sqrt(eigenval)/2)
# make new set of PCs
basque_cPCs <- data.frame(basque_pcs$superpop)
names(basque_cPCs)[1]<- "superpop"
for (i in 1:10){
  basque_cPCs[1+i] <- basque_pcs[,paste0("PC", i)]/(-sqrt(eigen_vals[i,])/2)
  names(basque_cPCs)[1+i] <- paste0("PC", i)
}
