library(sommer)
########################
#####-----input ----####
########################
### read the map file
map_v4 <- data.table::fread("~/work/AGPv4_Annotations/map_v4.txt")
### read the marker data
markers <- readRDS("~/save/PhD_work/Genotypes/GenotypeMat.RDS")
markers <- markers-1 #changing to -1, 0, and 1 as required by A.mat in the sommer package
### read evolutionary constraints data
evo_cons <- read.csv("bloc_synthese2.csv",sep = ";")
evo_cons$block_id_new <- paste0("Chr_",evo_cons$ChrName,"_block_",evo_cons$Bloc)
########################
####--setting loop params #####
########################
loop.start <- 1
loop.end <- nrow(evo_cons)
selected.kinships <- list()
remain.kinships <- list()
no.snp.found <- vector()
snp.count <- data.frame(matrix(data=NA, nrow=nrow(evo_cons),ncol=2))
colnames(snp.count) <- c("Block","Count")
#########################
########################
########################


####start of different loops####
## i will guide GO terms selection
for (i in loop.start:loop.end){
  selected.block.id <- evo_cons$block_id_new[i]
  sel_evo_cons <- evo_cons[evo_cons$block_id_new==selected.block.id,]
  
  snps_in_data <- map_v4[1,]
  snps_in_data[1,] <- NA
  for(j in 1:nrow(sel_evo_cons)){
    step1 <- map_v4[map_v4$chr==sel_evo_cons$ChrName[j],]
    snps_present <- step1[(step1$pos>=sel_evo_cons$BlocStart[j])& (step1$pos<=sel_evo_cons$BlocStop[j]),]
    snps_in_data <- rbind(snps_in_data,snps_present)
  }
  snps_in_data <- na.omit(snps_in_data)
  
  ##setting kinship matrices
  
  selected_snps <- unique(snps_in_data$snp)
  snp.count$Block[i] <- selected.block.id
  snp.count$Count[i] <- length(selected_snps)
  
  if (length(selected_snps)==0) {no.snp.found[i]=selected.block.id} else{
    remain_snps <- setdiff(as.character(colnames(markers)),selected_snps)
    
    selected.kinships[[selected.block.id]] <- A.mat(X = markers[,colnames(markers)%in%selected_snps],min.MAF = 0.01)
    #selected.kinships[[selected.block.id]] <- statgenGWAS::kinship(markers[,match(selected_snps,colnames(markers))],method = "vanRaden")
    remain.kinships[[selected.block.id]] <- A.mat(X = markers[,colnames(markers)%in% remain_snps],min.MAF = 0.01)

    print(i)
  }
}

no.snp.found <- na.omit(no.snp.found)
save(no.snp.found,selected.kinships,remain.kinships,snp.count,file = "kinships.Rdata")
