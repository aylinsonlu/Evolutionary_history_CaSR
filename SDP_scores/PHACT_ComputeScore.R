library(ape)
library(tidytree)
library(stringr)
library(dplyr)
library(bio3d)
args = commandArgs(trailingOnly=TRUE)

aa_to_num <- function(aa) {
  amino_acids <- c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")
  num <- sapply(aa, function(a){ifelse(sum(amino_acids %in% a) == 1, as.numeric(which(amino_acids %in% a)), 21)})
  # num <- ifelse(sum(amino_acids %in% aa) == 1, as.numeric(which(amino_acids %in% aa)), 21)
  return(num)
}

num_to_aa <- function(num) {
  amino_acids <- c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")
  aa <- ifelse(num == 21, 21, amino_acids[num])
  return(aa)
}


compute_score <- function(file_nwk, file_nwk_main, file_rst, file_fasta, output_name, parameters, main_root) {
  
  nodes <- unlist(str_split(nodes, pattern = ","))
 
  # Read tree file
  tr_org_main <- read.tree(file_nwk_main)
  tr_org <- read.tree(file_nwk)
  x <- read.table(file = file_rst, sep = '\t', header = TRUE, fill = TRUE)
  colnames(x)[4:ncol(x)] <- gsub("p_", replacement = "", x = colnames(x)[4:ncol(x)], fixed = TRUE )
  x[,1] <- str_remove(x[,1], "Node")  
  
  # Tree_info: node-node, node-leaf connections
  tree_info_main <- as.data.frame(as_tibble(tr_org_main))
  tree_info <- as.data.frame(as_tibble(tr_org))
  
  # Read fasta file, MSA
  fasta <- read.fasta(file = file_fasta)
  msa <- fasta$ali
  
  # connections_1: Parent node, connections_2: connected node/leaf
  connections_1 <- tree_info$parent
  connections_2 <- tree_info$node
  
  # Names of leaves
  names_all <- tr_org[["tip.label"]]
  msa <- msa[names_all, ]
  # Number of total leaves&nodes
  num_leaves_main <- length(tr_org_main[["tip.label"]])
  num_nodes_main <- tr_org_main[["Nnode"]]
  num_leaves <- length(tr_org[["tip.label"]])
  num_nodes <- tr_org[["Nnode"]]
  
  # Distance between leaves & nodes
  dd_node <- dist.nodes(tr_org)
  dist_leaf <- dd_node[1:num_leaves, 1:num_leaves]
  dist_node <- dd_node[(num_leaves+1):(num_leaves + num_nodes), (num_leaves+1):(num_leaves + num_nodes)]
  
  # Human position (leaf & node)
  nodes_raxml_main <- as.numeric(gsub(pattern = "Node", replacement = "", x = tree_info_main[num_leaves_main+1:num_nodes_main, "label"])) #Node or Branch
  names(nodes_raxml_main) <- tree_info_main[num_leaves_main+1:num_nodes_main, "node"]
  
  nodes_raxml <- as.numeric(gsub(pattern = "Node", replacement = "", x = tree_info[num_leaves+1:num_nodes, "label"])) #Node or Branch
  names(nodes_raxml) <- tree_info[num_leaves+1:num_nodes, "node"]
  nodes_raxml[1] <- as.numeric(substr(main_root,5,nchar(main_root)))
  
  root <- 1
  # Total number of positions from ancestralProbs file
  total_pos <- max(x$Site)
  
  # Chosen positions (all or some)
  positions <- 1:total_pos
  score_all <- matrix(0, total_pos, 21)
  
  # Connections between leaves & nodes
  chosen_leaves <- tree_info[1:num_leaves,c("parent", "node")]
  # Connections between nodes & nodes
  chosen_nodes <- tree_info[(num_leaves+2):(num_leaves +num_nodes),c("parent", "node")]
  leaf_names <- tree_info$label
  
  
  #########################################3
  
  parameters <- unlist(str_split(parameters, pattern = ","))
  
  score_norm <- t(mapply(function(ps, parameter){position_score(ps, x, msa, main_root, num_nodes_main, num_leaves_main, num_nodes, num_leaves, total_pos, nodes_raxml, parameter, chosen_leaves, chosen_nodes, d_n, d_l, root, tree_info)}, rep(positions, length(parameters)), rep(parameters, each = length(positions))))
  score_norm_with_leaf <- matrix(unlist(score_norm[ ,1]), nrow = length(positions) * length(parameters), ncol = 20, byrow = TRUE)
  score_norm_with_leaf <- cbind(rep(positions, length(parameters)), score_norm_with_leaf)
  colnames(score_norm_with_leaf) <- c("Pos/AA", num_to_aa(1:20))
  
  print_wl <- lapply(1:length(parameters), function(p){
    score_to_print_wl <- score_norm_with_leaf[(positions + length(positions)*(p - 1)), ]
    write.csv(score_to_print_wl, sprintf("%s.csv", output_name), row.names = FALSE, quote = FALSE)
  })

}

position_score <- function(ps, x, msa, main_root, num_nodes_main, num_leaves_main, num_nodes, num_leaves, total_pos, nodes_raxml, parameter, chosen_leaves, chosen_nodes, d_n, d_l, root, tree_info) {
  position <- ps
  b1 <- position + total_pos*(0:(num_nodes_main-1))
  TT <- x[b1,]
  
  node_info <- as.numeric(TT[,1])
  sort_node_info <- sort(node_info, decreasing = F, index.return=T)
  TT <- TT[sort_node_info$ix,]
  
  label_tr <- tree_info$label[(num_leaves+1):(num_leaves+num_nodes)]
  label_tr[1] <- main_root
  
  nms_TT <- TT$Node
  elim <- c()
  for (i in 1:num_nodes_main){
    if (!is.element(paste("Node",TT$Node[i], sep = ""),label_tr)){
      elim <- c(elim, i)
    }
  }
  
  TT <- TT[-elim,]
  matrix_prob <- matrix(0, num_nodes, 20)
  
  names_TT <- TT$Node
  ord <- matrix(0,num_nodes,1)
  for (i in 1:num_nodes){
    if (paste("Node",TT$Node[i], sep = "")==main_root){
      if (length(which(is.na(tree_info$branch.length)==1))>=1){
        ord[i] <- tree_info$node[which(is.na(tree_info$branch.length)==1)]
      } else {
        ord[i] <- tree_info$node[which(tree_info$branch.length==0)]
      }
    } else {
      ord[i] <- tree_info$node[which(tree_info$label==paste("Node",TT$Node[i], sep = ""))]
    }
  }
  
  aa <- sort(ord, decreasing = F, index.return = T)
  TT <- TT[aa$ix,]
  
  probs <- data.matrix((TT[, (4:ncol(TT))]))
  rownames(probs) <- NULL
  rr <- aa_to_num(colnames(x)[4:ncol(TT)])
  matrix_prob[,rr] <- probs
  
  position_vec <- msa[, ps]
  
  position_num <- aa_to_num(position_vec)
  prob_leaves <- matrix(0, num_leaves, 20)
  prob_leaves[cbind(which(position_num <= 20), position_num[which(position_num <= 20)])] <- 1
  
  gaps <- which(position_num == 21)
  
  diff_leaves <- matrix(0, num_leaves, 20)
  diff_leaves <- prob_leaves - matrix_prob[(chosen_leaves$parent - num_leaves), ]
  
  diff_nodes <- matrix(0, num_nodes-1, 20)
  diff_nodes <- matrix_prob[chosen_nodes[,2] - num_leaves, ] - matrix_prob[chosen_nodes[,1] - num_leaves, ]
  
  ################## weights
  weight_leaf <- matrix(1,1,num_leaves)
  weight_node <- matrix(1,1,num_nodes)
  ####################
  
  score <- matrix(0,1,20)
  
  s1 <- sapply(1:20, function(ii){
    dif_pr <- diff_nodes[1:(num_nodes-1),ii]
    dif_pr[dif_pr<0] <- 0
    # CHANGED 29.03
    score[ii] <<- score[ii] + sum(1 * dif_pr)
  })
  
  vect_root <- matrix_prob[root,]
  
  score_without_leaf <- score
  score_without_leaf <- score_without_leaf + vect_root
  score <- score + vect_root
  
  s2 <- sapply(1:20, function(ii){
    diff_lf <- diff_leaves[1:num_leaves,ii]
    diff_lf[gaps] <-  0
    diff_lf[diff_lf<0] <- 0
    
    s1 <- sum(weight_leaf[(1:length(diff_lf))] * diff_lf[(1:length(diff_lf))])
    score[ii] <<- score[ii] + s1
  })
  
  score_norm <- score
  score_wol_norm <- score_without_leaf
  
  kk1 <- sum(score_norm)
  kk2 <- sum(score_wol_norm)
  
  sums1 <- kk1-max(score_norm)
  sums2 <- kk2-max(score_wol_norm)
  
  sums1_upd <- (-(length(which(score_norm<0.0001))*0.1)/20+0.1)*(sums1)
  sums2_upd <- (-(length(which(score_wol_norm<0.0001))*0.1)/20+0.1)*(sums2)
  scores <- list()
  scores$score_with_leaf <- (score_norm*0.9 + sums2_upd)/(num_nodes+num_leaves)
  scores$score_without_leaf <- (score_wol_norm*0.9 + sums2_upd)/num_nodes
  return(scores)
}

