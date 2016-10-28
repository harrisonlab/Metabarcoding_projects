# produces a correlation matrix 

correr <- function(x) {
	l <- (length(x)-2)
	count <- (0.5*length(x))*((length(x)-3))+1
	
	x1 <- matrix(rep(seq((l+1),2),l),l,l)
	x2 <- t(matrix(rev(x1),l,l))
	x3 <- x2 + matrix(rep(seq(l,1),l),l,l)
	
	vx0 <- x[rep(1,count)]
	vx1 <- x[x1[lower.tri(x1,T)]]
	vx2 <- x[x2[lower.tri(x2,T)]]
	vx3 <- x[x3[lower.tri(x3,T)]]
	df <- cbind(vx0,vx1,vx2,vx3)
	cor.tot <- apply(df,1,function(y) cor(y[1:2],y[3:4]))
	sum(cor.tot)/count
}


x1 <- matrix(seq(1,5),5,5)
x1[lower.tri(x1,T)]


#x1 <- t(matrix(seq(1,5),5,5)) +  matrix(seq(4,0),5,5)

#mydist <- x1[lower.tri(x1,T)]

correr2 <- function(x) {
	y <- x
	mycorr <- numeric(0)
	count<-1
	while (length(x) >2) {
		mycorr[count] <- cor(y[1:(length(x)-1)],x[2:length(x)])
		count <- count +1
		x<-x[-1]
	}
	return(mycorr)
}
