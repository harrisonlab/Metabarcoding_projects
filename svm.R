# SUPPORT VECTOR MACHINE
library(e1071)
# Train SVM

Xm <- t(Sub.Set.New$A)
resp <- Sub.Set.New$targets$Condition
svm1 <- svm(Xm, resp, type="C-classification", kernel="linear")

# Estimate training error
trpred <- predict(svm1, Xm)
sum( trpred != resp)
table(trpred, resp)

# 10 fold cross-validation
trcv <- svm(Xm, resp, type="C-classification", kernel="linear",cross=10)
summary(trcv)

# Load and Normalise test set
targets.test <- readTargets("targets.mvx.txt")
RG.test <- read.maimages(
	targets.test, 
	path="../Data", 
	columns = list(
		G = "gMedianSignal", 
		Gb = "gBGMedianSignal", 
		R = "gProcessedSignal",
		Rb = "gIsPosAndSignif"
	),
	other.columns=list(
		"gIsPosAndSignif",
		"gIsWellAboveBG",
		"gIsFeatNonUnifOL",
		"gIsBGNonUnifOL",
		"gIsFeatPopnOL",
		"gIsBGPopnOL",
		"IsManualFlag",
		"gIsSaturated"
	),
	annotation = c("Row", "Col","FeatureNum", "ControlType","ProbeName","SystematicName")
)

RG.test <- backgroundCorrect(RG.test, method="normexp", offset=50)
RG.test$G <- normalizeBetweenArrays(RG.test$G, method="quantile")
RG.test$G <- log2(RG.test$G)
RG.test.filt <- RG.test[filt,]
RG.test.filt <- RG.test.filt[RG.test.filt$genes$ControlType==0,]
RG.test.ess <- new("MAList", list(targets=RG.test.filt$targets, genes=RG.test.filt$genes, source=RG.test.filt$source,  A=RG.test.filt$G))
E.avg.test <- avereps(RG.test.ess, ID=RG.test.ess$genes$ProbeName)
E.avg.test <- avereps(E.avg.test, ID=E.avg.test$genes$SystematicName)

# Predict classes of test set
Sub.Set.Test <- E.avg.test[E.avg.test$genes$ProbeName %in% results$ProbeName,]
Xmtr <- t(Sub.Set.Test$A)
tepred <- predict(svm1, Xmtr)
sum(tepred != Sub.Set.Test$targets$Condition)
table(tepred, Sub.Set.Test$targets$Condition)


func <- function(cluster=10,m=1) {
	subset1 <- E.avg[E.avg$genes$ProbeName %in% cas.comp.sub$ProbeName,1:4] 
	subset1$genes$Group <- f.ward(subset1$A,cluster,m)
	clustered <- merge(cas.comp.sub, subset1$genes[,c("ProbeName","Group")],by.x="ProbeName", by.y="ProbeName")
	results <- do.call(rbind, lapply(split(clustered <- clustered[order(clustered$adj.P.Val,decreasing=TRUE),], clustered$Group), tail, 1)) 
	Sub.Set.New <- E.avg[E.avg$genes$ProbeName %in% results$ProbeName,]
	Xm <- t(Sub.Set.New$A)
	resp <- Sub.Set.New$targets$Condition
	svm1 <- svm(Xm, resp, type="C-classification", kernel="linear")
	trpred <- predict(svm1, Xm)
	print("Estimated training error:")
	print (sum( trpred != resp))
	trcv <- svm(Xm, resp, type="C-classification", kernel="linear",cross=10)
	print(summary(trcv))
	Sub.Set.Test <- E.avg.test[E.avg.test$genes$ProbeName %in% results$ProbeName,]
	Xmtr <- t(Sub.Set.Test$A)
	tepred <- predict(svm1, Xmtr)
	sum(tepred != Sub.Set.Test$targets$Condition)
	vtab <- table(tepred, Sub.Set.Test$targets$Condition)
	
	return(vtab)
}


#leave one out cross validation
Sub.Set.New <- E.avg[E.avg$genes$ProbeName %in% results$ProbeName,]
Xm <- t(Sub.Set.New$A[,c(1:3,5:7,9:11,13:15)])
resp <- Sub.Set.New$targets$Condition[,c(1:3,5:7,9:11,13:15)]
svm1 <- svm(Xm, resp, type="C-classification", kernel="linear")
trpred <- predict(svm1, Xm)
sum( trpred != resp)
table(trpred, resp)
trcv <- svm(Xm, resp, type="C-classification", kernel="linear",cross=10)
summary(trcv)
#test data
Xmtr <- t(Sub.Set.New$A[,c(4,8,12,16)])
tepred <- predict(svm1, Xmtr)
sum(tepred != Sub.Set.New$targets$Condition[c(4,8,12,16)])
table(tepred, Sub.Set.New$targets$Condition[c(4,8,12,16)])

