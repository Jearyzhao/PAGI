Getpathwayactivity<-
function(d,dd){
        ddx <- rbind(d,dd)
        datax <- t(ddx)
        px <- prcomp(datax,scale. = F)
        datax1 <- predict(px)
        pa <- data.frame(features = datax1[,1],labels = lab)
        return(pa)
}        
