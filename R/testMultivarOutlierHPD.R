# multivariate_multimodal_HDF_06-19-18





#################################################
# dealing with correlated params

x<-rnorm(100,1,1)
y<-(x*1.5)+rnorm(100)
pIntX<-highestPostDensity(x,prob=0.8)
pIntY<-highestPostDensity(y,prob=0.8)

# plot it
xLims<-c(min(c(x,pIntX)),max(c(x,pIntX)))
yLims<-c(min(c(y,pIntY)),max(c(y,pIntY)))
plot(x,y,xlim=xLims,ylim=yLims)
rect(pIntX[1],pIntY[1],pIntX[2],pIntY[2])

# that doesn't look good

#####
# with PCA

x<-rnorm(100,1,1)
y<-(x*3)+rnorm(100)
pcaRes<-princomp(cbind(x,y))
PC<-pcaRes$scores

pIntPC1<-highestPostDensity(PC[,1],prob=0.8)
pIntPC2<-highestPostDensity(PC[,2],prob=0.8)

# plot it
xLims<-c(min(c(PC[,1],pIntPC1)),max(c(PC[,1],pIntPC1)))
yLims<-c(min(c(PC[,2],pIntPC2)),max(c(PC[,2],pIntPC2)))
plot(PC[,1],PC[,2],xlim=xLims,ylim=yLims)
rect(pIntPC1[1],pIntPC2[1],pIntPC1[2],pIntPC2[2])


#####################
x<-rnorm(100,1,1)
y<-(x*0.8)+rnorm(100)
z<-sample(c(rnorm(50,3,2),rnorm(50,30,3)))
a<-sample(c(rnorm(50,3,2),rnorm(50,10,3)))+x^2

#plot(x,y)
#plot(x,z)
#plot(x,a)
data<-cbind(x,y,z,a)

# actual outlier, but maybe not obvious if PCA isn't applied
outlier<-c(2,0.6,10,8)

# this one should be fine
outlier2<-c(1,0,30,5)

# this one should be an obvious outlier no matter what
outlier3<-c(3,-2,20,18)



testMultivarOutlierHPD<-function(dataMatrix,outlier,prob,pca=TRUE,...){
	dataAll<-rbind(data,outlier)
	# use PCA or not
	if(pca){
		pcaRes<-princomp(dataAll)
		varia<-pcaRes$scores
	}else{
		# use raw
		varia<-dataAll
		}
	# seperate out the outlier to be tested
	variaOut<-varia[nrow(varia),]
	# now remove outlier from variable sample to get HPD from
	varia<-varia[-nrow(varia),]
	#
	# get HPDs
	HPD<-lapply(1:ncol(varia),function(i)
		 highestPostDensity(varia[,i],prob=prob,...))
	#
	# test if outlier is within intervals listed for each variable
	withinHPD<-sapply(1:ncol(varia),function(x) 
		any(sapply(1:nrow(HPD[[x]]),function(y) 
			HPD[[x]][y,1]<=variaOut[x] & HPD[[x]][y,2]>=variaOut[x]
			))
		)
	res<-all(withinHPD)
	return(res)
	}

# this one should appear to be an outlier (but only if PCA is applied)
testMultivarOutlierHPD(dataMatrix=data,outlier=outlier,prob=0.8)
testMultivarOutlierHPD(dataMatrix=data,outlier=outlier,prob=0.8,pca=FALSE)

# this one should be within the 80% area
testMultivarOutlierHPD(dataMatrix=data,outlier=outlier2,prob=0.8)
testMultivarOutlierHPD(dataMatrix=data,outlier=outlier2,prob=0.8,pca=FALSE)

# this one should be outside the 80% area
testMultivarOutlierHPD(dataMatrix=data,outlier=outlier3,prob=0.8)
testMultivarOutlierHPD(dataMatrix=data,outlier=outlier3,prob=0.8,pca=FALSE)

