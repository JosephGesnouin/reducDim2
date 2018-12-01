library(dimRed)
library(reshape2)
library(ggplot2)
library(FactoMineR)
library(tm)
library(stringr)
pen.tra = read.table("/Users/jzk/Documents/M2/reducDim/pendigits.tra", sep = ",")
pen.tes = read.table("/Users/jzk/Documents/M2/reducDim/pendigits.tes", sep = ",")
pen = rbind(pen.tra, pen.tes)

names(pen) = c(paste0(c("X", "Y"), rep(1:8, each = 2)), "digit")
class(pen$digit)

pen$digit = factor(pen$digit)
class(pen$digit)
table(pen$digit)
dat=pen
ggplot(pen, aes(digit, fill = digit)) + 
  geom_bar() +
  guides(fill = FALSE)

tab = apply(pen[-17], 2, function(x) {
  return (c(Mean = mean(x),
            Median = median(x),
            StdDev = sd(x)))
})
round(t(tab), 2)

pen.X = pen[seq(1, 15, 2)]
names(pen.X)
pen.X.melt = melt(pen.X)
knitr::kable(head(pen.X.melt))
ggplot(pen.X.melt, aes(variable, value)) +
  geom_boxplot()

ggplot(melt(pen[seq(2, 16, 2)]), aes(variable, value)) + geom_boxplot()

apply(pen[-17], 2, tapply, pen$digit, mean)

agg = aggregate(. ~ digit, pen, mean)
knitr::kable(agg, digits = 2)
pen.melt = melt(pen, id.vars = 17)
head(pen.melt)
ggplot(pen.melt, aes(digit, value, color = digit)) + 
  geom_boxplot() +
  facet_wrap(~ variable)

first.0 = subset(pen, digit == 0)[1,]
x = unlist(first.0[seq(1,15, by = 2)])
y = unlist(first.0[seq(2,16, by = 2)])
par(mar = c(0, 0, 0, 0) + .1)
plot(x, y, type = "b", pch = " ", axes = FALSE, col = "gray70")
text(x, y, 1:8)

############PCA
res = PCA(pen, quali.sup = 17, graph = FALSE)
eig = data.frame(comp = 1:16,
                 setNames(res$eig, c("eigenvalue", "percent", "cum.percent")))
res2 = data.frame(res$ind$coord, digit = pen$digit)
ggplot(res2, aes(Dim.1, Dim.2, color = digit)) + 
  geom_point()

ggplot(res2, aes(Dim.1, Dim.2, color = digit)) +
  geom_point() +
  facet_wrap(~ digit) +
  guides(color = FALSE)
############MDS
library(kernlab)
# Nonmetric MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name

library(MASS)
d <- dist(pen) # euclidean distances between the rows
fit <- isoMDS(d, k=2) # k is the number of dim
fit # view results

# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Nonmetric MDS", type="n")
text(x, y, labels = row.names(mydata), cex=.7)


####ISOMAP
library(dimRed)
emb2 <- embed(pen, "Isomap", .mute = NULL, knn = 10)
plot(emb2, type = "2vars")

####Lle 
library(lle)
neighbours <- find_nn_k(pen, k=15)
neighbours[1:6, 1:6]
weights <- find_weights(neighbours, pen, m=2, reg=2)
weights$wgts[1:6, 1:6]
library(scatterplot3d)
# the 3-D graph of original data
scatterplot3d(x=pen[,1], y=pen[,2], z=pen[,3], color=pen[,2])
k5 <- lle(pen, m=2, k=5, reg=2, ss=FALSE, id=TRUE, v=0.9 )
plot(k5$Y, main="K=5 data", xlab=expression(y[1]), ylab=expression(y[2]))

# K=15
k15 <- lle(pen, m=2, k=15, reg=2, ss=FALSE, id=TRUE, v=0.9 )
plot(k15$Y, main="K=15 data", xlab=expression(y[1]), ylab=expression(y[2]))
# K=40
k40 <- lle(X, m=2, k=40, reg=2, ss=FALSE, id=TRUE, v=0.9 )
plot(k40$Y, main="K=40 data", xlab=expression(y[1]), ylab=expression(y[2]))

#########LTSA 
library(ltsa)
library(Rdimtools)
## 1. use 10%-connected graph
output1 <- do.ltsa(pen,ndim=2)

## 2. use 25%-connected graph
output2 <- do.ltsa(X,ndim=2,type=c("proportion",0.25))

## 3. use 50%-connected graph
output3 <- do.ltsa(X,ndim=2,type=c("proportion",0.50))

## Visualize three different projections
par(mfrow=c(1,3))
plot(output1$Y[,1],output1$Y[,2],main="10%")
plot(output2$Y[,1],output2$Y[,2],main="25%")
plot(output3$Y[,1],output3$Y[,2],main="50%")


######t-SNE
library(caret)  
library(Rtsne)

data_tsne=pen
set.seed(9)  
tsne_model_1 = Rtsne(as.matrix(data_tsne), check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.5, dims=2)
## getting the two dimension matrix
d_tsne_1 = as.data.frame(tsne_model_1$Y)
## plotting the results without clustering
ggplot(d_tsne_1, aes(x=V1, y=V2)) +  
  geom_point(size=0.25) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  xlab("") + ylab("") +
  ggtitle("t-SNE") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_colour_brewer(palette = "Set2")

## use the S4 Class directly:
isomap <- Isomap()
emb <- isomap@fun(dat, isomap@stdpars)
## or simpler, use embed():
samp <- sample(nrow(dat), size = 200)
emb2 <- embed(dat[samp], "Isomap", .mute = NULL, knn = 10)
emb3 <- emb2@apply(dat[-samp])
plot(emb2, type = "2vars")
plot(emb3, type = "2vars")


kamada_kawai <- KamadaKawai()
kk <- kamada_kawai@fun(dat, kamada_kawai@stdpars)
plot(kk@data@data)


## use the S4 class directly:
kpca <- kPCA()
emb <- kpca@fun(dat, kpca@stdpars)
## simpler, use embed():
emb2 <- embed(dat, "kPCA")
plot(emb, type = "2vars")


leim <- LaplacianEigenmaps()
emb <- leim@fun(dat, leim@stdpars)
plot(emb@data@data)


## directy use the S4 class:
lle <- LLE()
emb <- lle@fun(dat, lle@stdpars)
## using embed():
emb2 <- embed(dat, "LLE", knn = 45)
plot(emb, type = "2vars")
plot(emb2, type = "2vars")


mds <- MDS()
emb <- mds@fun(dat, mds@stdpars)
## use embed():
emb2 <- embed(dat, "MDS", d = function(x) exp(stats::dist(x)))
plot(emb, type = "2vars")
plot(emb2, type = "2vars")

nmds <- nMDS()
emb <- nmds@fun(dat, nmds@stdpars)
## using embed()
emb2 <- embed(dat, "nMDS", d = function(x) exp(dist(x)))
plot(emb, type = "2vars")
plot(emb2, type = "2vars")


set.seed(4646)
factorization <- embed(dat, "NNMF")
proj_dat <- factorization@apply(dat)
plot(proj_dat@data[, 1], proj_dat@data[, 2])
# project new values:
nn_proj <- predict(factorization, iris[1:7, 1:4])
nn_proj


pca <- PCA()
emb <- pca@fun(dat, pca@stdpars)
## using embed()
emb2 <- embed(dat, "PCA")
plot(emb, type = "2vars")
plot(emb@inverse(emb@data), type = "3vars")

pca_l1 <- PCA_L1()
emb <- pca_l1@fun(dat, pca_l1@stdpars)
## using embed()
emb2 <- embed(dat, "PCA_L1")
plot(emb, type = "2vars")
plot(emb@inverse(emb@data), type = "3vars")
