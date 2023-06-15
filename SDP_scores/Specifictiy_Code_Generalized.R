args <- commandArgs(trailingOnly = TRUE)

compute_specificity <- function(files, distances, output_name) {

  MAIN <- read.csv(files[1])
  
  rem_subfams <- list()
  for (i in 2:length(files)){
    rem_subfams[[names[i]]] <-  read.csv(files[i])
  }
  
  weight <- 1/distances
  
  total_pos <- length(MAIN[,1])
  aas <- colnames(MAIN[,2:21])
  pos <- 1:total_pos
  specificity <- matrix(0,total_pos,(1+3*length(files)))
  specificity[,1] <- pos
    
  for (i in 1:total_pos) { 
      main <- as.numeric(MAIN[i,2:21])
      main <- main/sum(main)
      
      aa_main <- aas[which(main==max(main))]
      obs_aas <- c()
      keep_rem_subfams <- c()
      for (j in 2:length(names)){
        subfam <- matrix(as.numeric(unlist(rem_subfams[[names[j]]])),total_pos,21)
        subfam <- subfam[i,2:21]
        subfam <- subfam/sum(subfam)
        
        aa_subfam <- aas[which(subfam==max(subfam))]
        keep_rem_subfams <- rbind(keep_rem_subfams, subfam)
        
        if (length(aa_subfam)>1){
          if (length(intersect(aa_subfam, aa_main))>=1){
            aa_subfam <- intersect(aa_subfam, aa_main)[1]
          } else {
            aa_subfam <- aa_subfam[1]
          }
        }
        obs_aas <- c(obs_aas, aa_subfam)
      }
      
      if (length(aa_main)>1){
        x1 <- c()
        for (f in 1:length(aa_main)){
          x1 <- c(x1, length(which(obs_aas==aa_main[f])))
        }
        ind_c <- which(x1==max(x1))[1]
      } else {
        ind_c <- 1
      }
      
      aa_main <- aa_main[ind_c]
      aa_main_ind <- which(main==max(main))
      score_main <- max(main)
      others_main <- exp(sum(main))-max(main)
      

      specificity[i,3:4] <- c(aa_main,  score_main)
      k <- 5
      c <- matrix(0, 1, (length(names)-1))
      weight_j <- c()
      for (j in 2:length(names)){
        subfam_j <- keep_rem_subfams[(j-1),]
        aa_subfam_j <- obs_aas[(j-1)]
        
        score_subfam_j <- max(subfam_j)
        main_subfam_j <- max(subfam_j[aa_main_ind])
        others_subfam_j <- exp(sum(subfam_j))-exp(main_subfam_j)
        
        specificity[i,k:(k+2)] <- c(aa_subfam_j, score_subfam_j, main_subfam_j) 
        k <- k+3
        
        val <- (score_subfam_j<=0.5)*(-1) + (score_subfam_j>0.5)*(1)

        c[(j-1)] <-  ((aa_main==aa_subfam_j)*exp(score_subfam_j)*(val==1) 
                      - (aa_main==aa_subfam_j)*others_subfam_j*(val==-1)
                      - (aa_main!=aa_subfam_j)*others_subfam_j) 
        weight_j <- c(weight_j, main_subfam_j)
      }
    
      weight_j <- max(weight_j)
      specificity[i,2] <- (exp(score_main) - (1-weight_j)*sum(weight*c))^score_main
      
            
  }
   
  nms <-  c("Pos", "SDP Score", "Main AA", "Score")
  for (i in 2:length(names)){
    nms <- c(nms, paste("AA_", names[i], sep = ""), "Score", paste("Main AA_", names[i], sep = ""))
  }
  colnames(specificity) <- nms
  write.csv(specificity, quote = F, file = paste(output_name, ".csv", sep = ""))   
  
}


