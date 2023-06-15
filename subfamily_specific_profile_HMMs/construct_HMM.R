construct_HMM <- function(folder_name, main_fam, close_fam, rest_fam, cfunvals_file, bound1_conslevel, bound2_conslevel, bound_blosum){ 
  
    cfunvals <- read.csv(sprintf("%s/%s", folder_name, cfunvals_file))
    
    conservation_main <- cfunvals[, which(colnames(cfunvals)==main_fam)]
    conservation_main_aa <- cfunvals[, match(sprintf("%s_aa", main_fam), colnames(cfunvals))]
    
    conservation_close <- cfunvals[, match(close_fam, colnames(cfunvals))]
    conservation_close_aa <- cfunvals[, match(sprintf("%s_aa", close_fam), colnames(cfunvals))]
    
    conservation_rest <- cfunvals[, match(rest_fam, colnames(cfunvals))]
    conservation_rest_aa <- cfunvals[, match(sprintf("%s_aa", rest_fam), colnames(cfunvals))]
  
    logodd2prob <- function(A) {
      A <- exp(A)
      A <- A / (1 + A) 
      return(A)
    }
    logodd2dist <- function(A) {
      A <- exp(-A)
      A <- A / (1 + A) 
      return(A)
    }
    amino_acids <- c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")
    data("BLOSUM62")
    blosum_mat <- BLOSUM62[amino_acids, amino_acids]
    Z_dist <- logodd2dist(blosum_mat)
    Z_prob <- logodd2prob(blosum_mat)
    
    for (ii in 1:20){
      Z_prob[ii,] <- Z_prob[ii,]/max(Z_prob[ii,])
    }
    aa_ord <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
    
    main_hmm <- read.table(sprintf("%s/%s.hmm", folder_name, main_fam), fill = TRUE, skip = 16)
    hmm_original <- readLines(sprintf("%s/%s.hmm", folder_name, main_fam))
    
    inds <- seq(4,length(main_hmm[,1]),4)
    inds2 <- seq(5,length(main_hmm[,1]),4)
    
    profile_upd<-main_hmm
    profile_upd[2,] <- c("",main_hmm[2,1:20])
    for (i in c(inds, inds2)){
      profile_upd[i, ] <- c("",main_hmm[i,1:20])
    }
    
    for (i in seq(length(profile_upd[,1]),1,-1)){
      if (!is.na(profile_upd[i,1]) & profile_upd[i,1] == "//"){
        length_mtx <- (i-1)
        break
      }
    }
    inds3 <- seq(7,length_mtx,4)
    add_later <- profile_upd[inds3,]
    profile_upd <- profile_upd[-inds3,]
    
    match_states <- as.numeric(add_later[1:length(add_later[,1]),1])
    max_pos <- max(match_states)
    insertion_states <- setdiff(1:max_pos, match_states)
    length_mtx <- length_mtx-length(add_later[,1])
    match_states_place <- seq(6,length_mtx,3)
    
    match_vec <- matrix(0,length(match_states),2)
    match_vec[,1] <- match_states
    match_vec[,2] <- match_states_place
    
    numpos <- length(cfunvals[,1])
    weight <- matrix(0,numpos,1)
    score <- matrix(0,numpos,1)
    spec_pos <- c()
    
    
    ##### REPRESENTATIVE FROM "MAIN"
    main_sub <- matrix(0,1,numpos)
    main_sub[1,] <- conservation_main
    main_sub_aa <- matrix(0,1,numpos)
    main_sub_aa[1,] <- conservation_main_aa
    
    ##### REPRESENTATIVE FROM "CLOSE"
    
    if (length(close_fam)==1){
      close <- matrix(0,1,numpos)
      close[1,] <- conservation_close
      close[which(close==0)] <- 1
      close_aa <- matrix(0,1,numpos)
      close_aa[1,] <- conservation_close_aa
    } else {
      close <- matrix(0,1,numpos)
      close_aa <- c()
      change1 <- matrix(1,1,numpos)
      change <- matrix(1,1,numpos)
      
      for (i in 1:numpos){
        aas <- conservation_close_aa[i,]
        aas_num <- aa_to_num(aas)
        vals <- conservation_close[i,]
        vals[which(vals==0)]<-1
        unq <- length(unique(aas_num))
        
        if (unq == 1){
          aa_sub <- unique(aas_num)
          close[i] <- sum(vals)/length(aas_num)
          ind_mx <- which(vals==max(vals))
          if (aa_sub==21){
            tr <- "-"
          } else {
            tr <- num_to_aa(aa_sub)
          }
          close_aa <- c(close_aa, tr)
          
        } else if (unq >= 2){
          aa_sub <- unique(aas_num)
          val_sub <- matrix(0,1,unq)
          for (j in 1:unq){
            if (aa_sub[j]==21){
              aa_c <- "-"
              val_sub[j] <- 0
            } else {
              aa_c <- num_to_aa(aa_sub[j])
              num_sub <- which(aas==aa_c)
              val_sub[j] <- max(vals[num_sub])
            }
          }
          
          aa_mains<-aa_to_num(main_sub_aa[i])
          if (is.element(aa_mains,aa_sub) & aa_mains!=21){
            plc <- which(aa_sub==aa_mains)
            close[i] <- val_sub[plc]
            ind_mx <- plc
            if (num_to_aa(aa_sub[ind_mx[1]])==21){
              aa_c <- "-"
            } else {
              aa_c <- num_to_aa(aa_sub[ind_mx[1]])
            }
            close_aa <- c(close_aa, aa_c)
            
          } else {
            close[i] <- sum(val_sub)/length(val_sub)
            ind_mx <- which(val_sub==max(val_sub))
            if (num_to_aa(aa_sub[ind_mx[1]])==21){
              aa_c <- "-"
            } else {
              aa_c <- num_to_aa(aa_sub[ind_mx[1]])
            }
            close_aa <- c(close_aa, aa_c)
          }
        } 
        
      }
      close_aa <- unlist(close_aa)
    }
    
    
    ##### REPRESENTATIVE FROM "REST"
    
    if (length(rest_fam)==1){
      rest <- matrix(0,1,numpos)
      rest[1,] <- conservation_rest
      rest[which(rest==0)] <- 1
      rest_aa <- matrix(0,1,numpos)
      rest_aa[1,] <- conservation_rest_aa
    } else {
      rest <- matrix(0,1,numpos)
      rest_aa <- c()
      change <- matrix(1,1,numpos)
      for (i in 1:numpos){
        aas <- conservation_rest_aa[i,]
        aas_num <- aa_to_num(aas)
        vals <- conservation_rest[i,]
        vals[which(vals==0)]<-1
        unq <- length(unique(aas_num))
        if (unq == 1){
          aa_sub <- unique(aas_num)
          rest[i] <- sum(vals)/length(aas_num)
          ind_mx <- which(vals==max(vals))
          if (aa_sub==21){
            tr <- "-"
          } else {
            tr <- num_to_aa(aa_sub)
          }
          rest_aa <- c(rest_aa, tr)
          
        } else if (unq >= 2){
          aa_sub <- unique(aas_num)
          val_sub <- matrix(0,1,unq)
          for (j in 1:unq){
            if (aa_sub[j]==21){
              tr <- "-"
              val_sub[j] <- 0
            } else {
              tr <- num_to_aa(aa_sub[j])
              num_sub <- which(aas==tr)
              val_sub[j] <- max(vals[num_sub])
            }
          }
          
          aa_mains<-aa_to_num(main_sub_aa[i])
          if (is.element(aa_mains,aa_sub) & aa_mains!=21){
            
            plc <- which(aa_sub==aa_mains)
            rest[i] <- val_sub[plc]
            ind_mx <- plc
            if (num_to_aa(aa_sub[ind_mx[1]])==21){
              tr <- "-"
            } else {
              tr <- num_to_aa(aa_sub[ind_mx[1]])
            }
            rest_aa <- c(rest_aa, tr)
            change[i] <- length(aa_sub)-as.numeric(is.element(21,aa_sub))
          } else {
            rest[i] <- sum(val_sub)/length(val_sub)
            ind_mx <- which(val_sub==max(val_sub))
            if (num_to_aa(aa_sub[ind_mx[1]])==21){
              tr <- "-"
            } else {
              tr <- num_to_aa(aa_sub[ind_mx[1]])
            }
            rest_aa <- c(rest_aa, tr)
            change[i] <-  length(aa_sub)-as.numeric(is.element(21,aa_sub))
          }
        } 
      }
      rest_aa <- unlist(rest_aa)
    }
    
    
    #################### DETERMINE 4 TYPES & WEIGHTS ##############################
    
    overall <- cbind(t(main_sub),t(close),t(rest))
    type <- matrix(0,1,numpos)
    for (k in 1:numpos){
      main_sub_obs <- main_sub_aa[k]
      close_obs <- close_aa[k]
      rest_obs <- rest_aa[k]
      
      if (close_obs=="-" && rest_obs=="-"){
        if (main_sub_obs =="-"){
          weight[k] <- 1/sum(overall[k,])
          type[k] <- 2
        } else {
          weight[k] <- (sum(overall[k,]))
          type[k] <- 1
        }
      } else if (main_sub_obs =="-"){
        if (close_obs=="-") {
          weight[k] <- 1/(sum(overall[k,]))
          type[k] <- 2
        } else if (rest_obs=="-") {
          weight[k] <- 1/(sum(overall[k,]))
          type[k] <- 2
        } else {
          weight[k] <- (sum(overall[k,]))
          type[k] <- 1
        }
        
      } else if (main_sub_obs!=close_obs && main_sub_obs!=rest_obs && rest_obs!=close_obs) {
        if (aa_to_num(main_sub_obs)!=21 && aa_to_num(close_obs)!=21 && aa_to_num(rest_obs)!=21) {
          if (sum(overall[k,1]>=bound1_conslevel)==1 && sum(overall[k,2:3]>=bound2_conslevel)==2){
            weight[k] <- (sum(overall[k,]))
            type[k] <- 1
          } else if (Z_prob[aa_to_num(main_sub_obs),aa_to_num(close_obs)]<=bound_blosum && Z_prob[aa_to_num(main_sub_obs),aa_to_num(rest_obs)]<=bound_blosum){
            weight[k] <- (sum(overall[k,]))
            type[k] <- 1
          } else {
            weight[k] <- (sum(overall[k,]))
            type[k] <- 4
          }
        } else if (aa_to_num(close_obs) ==21) {
          if (sum(overall[k,1]>=bound1_conslevel)==1 && sum(overall[k,2:3]>=bound2_conslevel)==2){
            weight[k] <- (sum(overall[k,]))
            type[k] <- 1
          } else if (Z_prob[aa_to_num(main_sub_obs),aa_to_num(rest_obs)]<=bound_blosum){
            weight[k] <- (sum(overall[k,]))
            type[k] <- 1
          } else {
            weight[k] <- (sum(overall[k,]))
            type[k] <- 4
          }
        } else if (aa_to_num(rest_obs) ==21) {
          if (sum(overall[k,1]>=bound1_conslevel)==1 && sum(overall[k,2:3]>=bound2_conslevel)==2){
            weight[k] <- (sum(overall[k,]))
            type[k] <- 1
          } else if (Z_prob[aa_to_num(main_sub_obs),aa_to_num(close_obs)]<=bound_blosum){
            weight[k] <- (sum(overall[k,]))
            type[k] <- 1
          } else {
            weight[k] <- (sum(overall[k,]))
            type[k] <- 4
          }
        }
        
      } else if (main_sub_obs==close_obs) {
        weight[k] <- 1/sum(overall[k,])
        type[k] <- 2
        
      } else if (main_sub_obs != close_obs && main_sub_obs==rest_obs) {
        if (close_obs=="-"){
          weight[k] <- (sum(overall[k,]))
          type[k] <- 3
        } else {
          if (Z_prob[aa_to_num(main_sub_obs),aa_to_num(close_obs)]>bound_blosum && overall[k,2]<bound2_conslevel){
            weight[k] <- 1/(sum(overall[k,]))
            type[k] <- 2
          } else {
            weight[k] <- (sum(overall[k,]))
            type[k] <- 3
          }
        }
      } else if (main_sub_obs!=close_obs && main_sub_obs!=rest_obs && rest_obs==close_obs) {
        if (sum(overall[k,1]>=bound1_conslevel)==1 && sum(overall[k,2:3]>=bound2_conslevel)==2){
          weight[k] <- (sum(overall[k,]))
          type[k] <- 1
        } else if (Z_prob[aa_to_num(main_sub_obs),aa_to_num(close_obs)]<=bound_blosum){
          weight[k] <- (sum(overall[k,]))
          type[k] <- 1
        } else {
          weight[k] <- 1/(sum(overall[k,]))
          type[k] <- 2
        }
      }
    }
    
    weg <- cbind(weight, t(type), 1:numpos)
    type1 <- weg[which(weg[,2]==1),]
    min_1 <- min(type1[,1])
    tyy <- type1[which(type1[,1]<3),1]
    max_1 <- max(tyy)
    type1 <- cbind(type1, type1[,1]/min_1)
    
    min_w1 <- mean(type1[,4])
    type4 <- weg[which(weg[,2]==4),]
    max_4 <- max(type4[,1])
    min_4 <- min(type4[,1])
    type4 <- cbind(type4, type4[,1]/(max_4)*min_w1)
    
    min_w1 <- min(type4[,4])
    type3 <- weg[which(weg[,2]==3),]
    max_3 <- max(type3[,1])
    min_3 <- min(type3[,1])
    type3 <- cbind(type3, type3[,1]/(max_3)*min_w1)
    
    min_w3 <- 0.2
    
    type2 <- weg[which(weg[,2]==2),]
    max_2 <- max(type2[,1])
    min_2 <- min(type2[,3])
    type2 <- cbind(type2, (type2[,1]/(max_2))*min_w3)
    
    type1[,4]<-type1[,4]*2
    type2[,4]<-type2[,4]
    type3[,4]<-type3[,4]
    type4[,4]<-type4[,4]
    
    weight_upd <- rbind(type1,type2,type3, type4)
    ord <- sort(weight_upd[,3], index.return=T)
    weight_upd <- weight_upd[ord$ix,]
    
    overall <- cbind(overall, weight_upd[,4])
    main_fam_aa <- main_sub_aa
    weights <- cbind(1:numpos, weight_upd[,4])
    
    
    
    #################### UPDATE PROFILE PROBABILITIES ##############################
    
    profile_weighted <- profile_upd
    residue_place <-  seq(6,length_mtx,3)
    profile_weighted[, 22]<-""
    profile_weighted[residue_place, 23:27]<-add_later[,1:5]
    
    
    for (i in match_states){
      if (i<max(weights[,1])){
        s <- match_vec[match_vec[,1]==i,2]
        spec_weight <- weights[which(weights[,1]==i),2]
        probs<-matrix(0,1,20)
        for (j in 2:21){
          if (profile_weighted[s,j] != "*"){
            if (aa_ord[(j-1)]==main_fam_aa[i]){
              profile_weighted[s,j] <- -log(exp(-as.numeric(unlist(profile_weighted[s,j])))*spec_weight)
            } else {
              profile_weighted[s,j] <- -log(exp(-as.numeric(unlist(profile_weighted[s,j])))*spec_weight)
            }
          }
        }
      }
    }
    profile_weighted_prob <- profile_weighted[6:length_mtx, 2:21]
    
    
    ########## WRITE THE HMMER MODEL TO THE FILE #####################################################
    
    position_txt <- matrix(0, 20, 2)
    position_txt[1,] <- c(11, 17)
    for (i in 2:20){
      position_txt[i, ] <- c(position_txt[(i-1),2]+3, position_txt[(i-1),2]+9)
    }
    
    txt<-""
    for (i in 1:(length(hmm_original)-22)){
      txt <- hmm_original[i+21]
      for (j in 1:20){
        substr(txt, position_txt[j,1], position_txt[j,2]) <- profile_weighted_prob[i,j]
      }
      hmm_original[(i+21)] <- txt
    }
    ln <- length(hmm_original)-1
    for (j in 1:length(hmm_original)){
      for (i in 1:(nchar(hmm_original[j])-6)){
        if (substr(hmm_original[j], i, i+6)== "*     *"){
          substr(hmm_original[j], i, i+6) <- "*      "
        }  
      }
    }
    hmm_original[ln+1] <- "//"
    hmm_original[2] <- sprintf("NAME Updated Profile %s", main_fam)
    writeLines(hmm_original, sprintf("%s/%s_UpdatedProfiles.hmm", folder_name, main_fam))


}

