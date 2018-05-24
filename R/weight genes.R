weight genes<-
function(d_ori,labels){
        d_new <- scale(d_ori,center = T,scale = T) #standardization of data set
        labels <- as.matrix(labels,nrow = 1,ncol = ncol(d_new))
        labels <- matrix(labels,nrow = 1,ncol = ncol(d_new))
        d_new <- rbind(d_new,labels)
        d_new <- d_new[,order(d_new[nrow(d_new),])] #sort the samples by labels
        stn <- which(d_new[nrow(d_new),] == "or",arr.ind = T)
        strm <- which(d_new[nrow(d_new),] == "rm",arr.ind = T)
        ln <- append(stn,strm)
        nor_d <- d_new[,ln]
        lnor <- ncol(nor_d)
        ntu_d <- d_new[,-ln]
        ltu <- ncol(ntu_d)
        d_new <- cbind(nor_d,ntu_d)
        d_new <- d_new[1:nrow(d_new) - 1,]
        d_new1 <- matrix(as.numeric(d_new),nrow = nrow(d_new),ncol = ncol(d_new))
        rownames(d_new1) <- rownames(d_new)
        colnames(d_new1) <- colnames(d_new1)
        d_new <- d_new1
        lsam <- ncol(d_new)
        tv1 <- matrix(nrow = nrow(d_new),ncol = 1)
        lab <- c(rep(1,lnor),rep(0,ltu))
        colnames(tv1) <- c("t-values")
        for (i in 1:nrow(d_new)) {        #using the "t value" as weights
          d1 <- d_new[i,1:lnor];d2 <- d_new[i,(lnor + 1):lsam]
          d1 <- as.numeric(d1);d2 <- as.numeric(d2)
          t <- t.test(d1,d2,alternative = "two.sided")
          corxy <- cor(as.numeric(d_new[i,]) ,lab,method = "pearson")
          t <- abs(t$statistic[[1]])
          d_new[i,] <- d_new[i,] * t ** 2 * corxy ** 2  
          tv1[i,1] <- t
        }
        d<-d_new
        return(d)
}
