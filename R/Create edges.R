create edges<-
function(d,relationship,beta){
       reln <- list(relationship)
       beta <- list(beta)
        ##finding EntryId from pathway
       rel <- pathway@edges
       x <- getEntryID(rel)
       ##finding genesnames of nodes
       ###finding numbers of nodes
       drn <- matrix(rownames(d))
       nods <- pathway@nodes
       nod <- nods
       lnods <- length(nods)
       j <- 0
       for (i in 1:lnods) {            
         if (getName(nods[[i]]) == "undefined") {
           nod[[i - j]] <- NULL
           j <- j + 1
    
         }
       }
       nods <- nod
       dd <- NULL
       # D<-matrix(nrow = 150,ncol= 2)
       for (i in 1:nrow(x)) {
         if (length(nods[[x[i,1]]]) == 0 || length(nods[[x[i,2]]]) == 0) next
         NN1 <- matrix(getName(nods[[x[i,1]]]))    
         NN2 <- matrix(getName(nods[[x[i,2]]]))
         num1 <- nrow(NN1)
         num2 <- nrow(NN2)
         if (length(rel[i]$relation@subtype) == 0) {
           beta_I <- 0.5
         }
         else {
           SN <- rel[i]$relation@subtype[[1]]@name 
           a <- which(reln[[1]] == SN,arr.ind = TRUE)
           beta_I <- beta[[1]][a]
         }

         for (j in 1:num1) {
           a <- which(drn == NN1[[j,1]],arr.ind = TRUE)
           if (length(a) == 0) next
           da <- t(matrix(d[a[[1]],]))
    
           for (h in 1:num2) {
             b <- which(drn == NN2[[h,1]],arr.ind = TRUE)
             if (length(b) == 0) {next}
             db <- t(matrix(d[b[[1]],]))
             dx <- as.numeric(da)
             dy <- as.numeric(db)
             corxy <- cor(dx,dy,method = "pearson")
             corxy <- abs(corxy)
             dc <- matrix(ncol = ncol(d))
             for (k in 1:ncol(d)) {
               dc[1,k][[1]] <- beta_I * corxy *(da[1,k][[1]]/2 + db[1,k][[1]]/2)
             }
             genenames <- paste(NN1[j,1],NN2[h,1],sep = " & ")
             rownames(dc) <- genenames
             dd <- rbind(dd,dc)
           }
        }
       }
       colnames(dd) <- colnames(d)
       dd <- dd[rowSums(dd == 0) == 0,]
       ddx <- dd
       return(dd)
}       
