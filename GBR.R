# Activate our dataset

# Activate our libraries
library(pracma)
library(kernlab)
library(pca3d)
library(plot3D)
library(plotly)
library(ggplot2)
library(diffusionMap)
library(destiny)
library(Matrix)
library(matrixStats)
library(rdist)
library(rgl)
library(gridExtra)
library(ggbiplot)

GBR <- read.csv("~/Desktop/Graduate Biology Research.csv")
v<-c(NA,NA,0.87,0.98,1.09,1.33,0.87,0.98,1.09,1.33,0.87,0.98,1.09,1.33)
GBR <- rbind(v,GBR)
GBR$Gene_ID <- NULL
GBR$length <- NULL
mrow<-rowSums(GBR)/12
Abb<-GBR-mrow
vAbb<-rowSums(Abb^2)
Aba<-Abb/sqrt(vAbb)
rowSums(Aba^2)
head(GBR)
GBR=scale(GBR,center=TRUE,scale=TRUE)
w1<-c(0.87,0.98,1.09,1.33,0.87,0.98,1.09,1.33,0.87,0.98,1.09,1.33)
w2<-c(2311,2727,2857,3541,2311,2727,2857,3541,2311,2727,2857,3541)
w3<-c(342.9,348.8,362.8,396.8,342.9,348.8,362.8,396.8,342.9,348.8,362.8,396.8)
w1=scale(w1,center=TRUE,scale=TRUE)
w2=scale(w2,center=TRUE,scale=TRUE)
w3=scale(w3,center=TRUE,scale=TRUE)
g=GBR%*%w1
g=GBR%*%w2
g=GBR%*%w3
print(g)

g1_row <- vector() #Stores row index of g matrix with values < -.3
g1_val <- vector() #Stores row value of g matrix when row value < -.3
g1 <- data.frame(g1_row, g1_val)

for(r in 1:nrow(g)){
    g_val = g[r, 1]
    if(g_val < -0.3){
      g1_row <- r
      g1_val <- g_val
      new_row <- c(g1_row, g1_val)
      g1 <- rbind(g1, new_row)
    }
}
g1_names <- c("G_INDEX", "G_VALUE")
names(g1) <- g1_names
print(g1)


g2_row <- vector() #Stores row index of g matrix with values > .3
g2_val <- vector() #Stores row value of g matrix when row value > .3
g2 <- data.frame(g2_row, g2_val)

for(r in 1:nrow(g)){
  g_val = g[r, 1]
  if(g_val > 0.3){
    g2_row <- r
    g2_val <- g_val
    new_row <- c(g2_row, g2_val)
    g2 <- rbind(g2, new_row)
  }
}
g2_names <- c("G_INDEX", "G_VALUE")
names(g2) <- g2_names
print(g2)


threshold_one <- data.frame() # subset of GBR with dot product values < -.3 
for (row in 1:nrow(g1)){
  gbr_index = g1[row,1]
  threshold_one = rbind(threshold_one, GBR[gbr_index, 1:12])
}
threshold_names <- colnames(GBR)
names(threshold_one) <- threshold_names
print(threshold_one)

threshold_two <- data.frame() # subset of GBR with dot product values > .3 
for (row in 1:nrow(g2)){
  gbr_index = g2[row,1]
  threshold_two = rbind(threshold_two, GBR[gbr_index, 1:12])
}

names(threshold_two) <- threshold_names
print(threshold_two)

# Merging g1 and g2 to get combined subset of GBR
g_indices = rbind(g1["G_INDEX"], g2["G_INDEX"])
g_indices <- g_indices[order(g_indices$G_INDEX), ]

threshold_three <- data.frame() # subset of GBR with dot product values > .3 
for (index in g_indices){
  gbr_index = index
  threshold_three = rbind(threshold_three, GBR[gbr_index, 1:12])
}
names(threshold_three) <- threshold_names
print(threshold_three)

a=prcomp(threshold_one)
b=prcomp(threshold_two)
c=prcomp(threshold_three)
summary(a)
summary(b)
summary(c)

print(a$rotation)
class(a$rotation)
print(a$x)

library(plotly)
plot_ly(x=a$x[,"PC4"], y=a$x[,"PC7"],z=a$x[,"PC8"], type="scatter3d", mode="markers", color=a$x[,"PC4"])
a1=ggbiplot(a,choices = c(1,2))
a2=ggbiplot(a,choices = c(3,4))
a3=ggbiplot(a,choices = c(5,6))
a4=ggbiplot(a,choices = c(7,8))
a5=ggbiplot(a,choices = c(9,10))
a6=ggbiplot(a,choices = c(11,12))

b1=ggbiplot(b,choices = c(1,2))
b2=ggbiplot(b,choices = c(3,4))
b3=ggbiplot(b,choices = c(5,6))
b4=ggbiplot(b,choices = c(7,8))
b5=ggbiplot(b,choices = c(9,10))
b6=ggbiplot(b,choices = c(11,12))

c1=ggbiplot(c,choices = c(1,2))
c2=ggbiplot(c,choices = c(3,4))
c3=ggbiplot(c,choices = c(5,6))
c4=ggbiplot(c,choices = c(7,8))
c5=ggbiplot(c,choices = c(9,10))
c6=ggbiplot(c,choices = c(11,12))

grid.arrange(a1,a2,a3,a4,a5,a6,nrow=3)
grid.arrange(b1,b2,b3,b4,b5,b6,nrow=3) 
grid.arrange(c1,c2,c3,c4,c5,c6,nrow=3)


# PCA

par(mar=c(1,1,1,1))
GBR<-GBR[1:200,2:14]
Mpca<-prcomp(GBR, center = TRUE, scale. = TRUE)
pca3d(Mpca, group = GBR$length, legend = "topright", biplot = TRUE, biplot.vars = 2)
pca3d(Mpca, group = GBR$length, legend = 'bottomright', show.shadows = TRUE)

# KPCA

train<-as.matrix(GBR[,3:6])
gau=kpca(train,kernel="rbfdot",sigma=1,features=2)
lap=kpca(train,kernel="laplacedot",sigma=1,features=2)
bessel=kpca(train,kernel="besseldot",sigma=1,features=2)
anova=kpca(train,kernel="anovadot",sigma=1,features=2)
plot3d(rotated(gau),col=as.integer(train),xlab = "Gaussian",ylab = "",zlab = "")
plot3d(rotated(lap),col=as.integer(train),xlab = "Laplacian",ylab = "",zlab = "")
plot3d(rotated(bessel),col=as.integer(train),xlab = "Bessel",ylab = "",zlab = "")
plot3d(rotated(anova),col=as.integer(train),xlab = "Anova",ylab = "",zlab = "")

# Diffusion Maps

dis <- rdist(GBR, metric = "hamming")
dif <- diffuse(dis, neigen = 15)
plot3d(dif$X[,2:9],xlab = "Diffusion Biology",ylab="",zlab = "")
