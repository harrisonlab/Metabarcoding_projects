options(stringsAsFactors = FALSE)
snps <- read.table("snps",header=T,sep="\t")
 
m27_trio <- snps[,c(1:6,13:14,7:8)]

m27_trio <- remove_impossible(m27_trio)
m27_trio <- remove_hom_child(m27_trio)
m27_trio <- remove_all_het(m27_trio)

phased_m27 <- t(as.data.frame(phaser(m27_trio[,5:10])))
colnames(phased_m27) <- c("m9.m27","m13.m27")
 
remove_impossible <- function(X) {
	X <- X[
		(
			(X[,9]==X[,5]| X[,9]==X[,6])|
			(X[,10]==X[,5]| X[,10]==X[,6])
		) &
		( 	
			(X[,9]==X[,7]| X[,9]==X[,8]) |
			(X[,10]==X[,7]| X[,10]==X[,8])
		)
	,]
	return(X)
}

remove_hom_child <- function(X) {
	X <- X[
		X[,9]!=X[,10]
	,]
	return(X)
}

remove_all_het <- function(X) {
	X <- X[
		(X[,5]==X[,6])|
		(X[,7]==X[,8])|
		(X[,9]==X[,10])
	,]
	return(X)
}

phaser <- function(X) {
	apply(X,1,function(x) 
	{
		count=0
		if(x[5]==x[4]){count=count+1}
		if(x[5]==x[3]){count=count+2}
		if(x[5]==x[2]){count=count+4}
		if(x[5]==x[1]){count=count+8}
		if(x[5]==1) {
			if(count==1|count==2|count==3|count==7|count==11){
				return(c(0,1))
			} else {
				return(c(1,0))
			}
		} else {
			if(count==1|count==2|count==3|count==7|count==11){
				return(c(1,0))
			} else {
				return(c(0,1))
			}
		}
	})
}

