library(igraph)
library(netdep)
library(e1071)
library(doParallel)
library(ggplot2)
library(gridExtra)
library(ggplotify)

internet <- read_graph("as-22july06.gml", format="gml") # internet snap shot
hep <- read_graph("hep-th.gml", format="gml") # co-authorships among the scientists posting on High-Energy Theory E-Print Archive


# sampling 400 sub-graph from internet data
set.seed(100)
internet.subgraphs <- list()
for (i in 1:400) {
  n.nodes <- sample(50:200, 1)
  internet.subgraphs[[i]] <- snowball.sampling(internet, n.nodes)$subG
  if(i %% 100==0) {print(as.character(i))}
}
summary(sapply(internet.subgraphs, gsize))

# sampling 400 sub-graph from hep data
set.seed(300)
hep.subgraphs <- list()
for (i in 1:400) {
  n.nodes <- sample(50:200, 1)
  repeat {
    error <- tryCatch(hep.subgraphs[[i]] <- snowball.sampling(hep, n.nodes)$subG,
                      error= function(e) return(T))
    if (is.list(error)) {break}
  }
  if(i %% 100 == 0) {print(as.character(i))}
}


# plot of sub-graphs
igraph.options(vertex.size=3, vertex.label=NA, edge.arrow.size=0.5)
set.seed(100)
par(mfrow=c(2,3))
for (i in sample(1:400, 6)) {
  plot(internet.subgraphs[[i]], layout=layout.kamada.kawai, vertex.label=NA)
}
for (i in sample(1:400, 6)) {
  plot(hep.subgraphs[[i]], layout=layout.kamada.kawai, vertex.label=NA)
}
par(mfrow=c(1,1))

degree.centrality <- function(igraph) {
  # returns the vector of degree centralities for each vertex
  N <- gorder(igraph)
  adj.matrix <- get.adjacency(igraph)
  result <- numeric(N)
  for (i in 1:N) {
    result[i] <- sum(adj.matrix[i,][-i]) / N
  }
  return(result)
}

nc <- detectCores()
registerDoParallel(nc-1)

# calculate statistics of degree centralities (Case 2)
D.stats <- foreach(i=1:400, .combine=rbind, .packages=c("igraph", "e1071")) %dopar% {
  D.INT <- degree.centrality(internet.subgraphs[[i]])
  row.INT <- c(mean(D.INT), var(D.INT), skewness(D.INT, type=3), kurtosis(D.INT, type=3))
  row.INT[is.nan(row.INT)] <- 0
  D.HEP <- degree.centrality(hep.subgraphs[[i]])
  row.HEP <- c(mean(D.HEP), var(D.HEP), skewness(D.HEP, type=3), kurtosis(D.HEP, type=3))
  row.HEP[is.nan(row.HEP)] <- 0
  df <- as.data.frame(rbind(row.INT, row.HEP))
  names(df) <- c("D.mean", "D.var", "D.skew", "D.kurt")
  df$graph <- c("INT", "HEP")
  return(df)
}
row.names(D.stats) <- 1:800
D.stats$graph <- as.factor(D.stats$graph)

# calculate statistics of clustering coefficients (Case 3)
CC.stats <- foreach(i=1:400, .combine=rbind, .packages=c("igraph", "e1071")) %dopar% {
  CC.INT <- transitivity(internet.subgraphs[[i]], type="local", isolates="zero") # clustering coefficient
  row.INT <- c(mean(CC.INT), var(CC.INT), skewness(CC.INT, type=3), kurtosis(CC.INT, type=3))
  row.INT[is.nan(row.INT)] <- 0
  CC.HEP <- transitivity(hep.subgraphs[[i]], type="local", isolates="zero")
  row.HEP <- c(mean(CC.HEP), var(CC.HEP), skewness(CC.HEP, type=3), kurtosis(CC.HEP, type=3))
  row.HEP[is.nan(row.HEP)] <- 0
  df <- as.data.frame(rbind(row.INT, row.HEP))
  names(df) <- c("CC.mean", "CC.var", "CC.skew", "CC.kurt")
  df$graph <- c("INT", "HEP")
  return(df)
}
row.names(CC.stats) <- 1:800
CC.stats$graph <- as.factor(CC.stats$graph)

# calculate statistics of betweenness centralities (Case 4)
B.stats <- foreach(i=1:400, .combine=rbind, .packages=c("igraph", "e1071")) %dopar% {
  B.INT <- betweenness(internet.subgraphs[[i]], directed=F)
  row.INT <- c(mean(B.INT), var(B.INT), skewness(B.INT, type=3), kurtosis(B.INT, type=3))
  row.INT[is.nan(row.INT)] <- 0
  B.HEP <- betweenness(hep.subgraphs[[i]], directed=F)
  row.HEP <- c(mean(B.HEP), var(B.HEP), skewness(B.HEP, type=3), kurtosis(B.HEP, type=3))
  row.HEP[is.nan(row.HEP)] <- 0
  df <- as.data.frame(rbind(row.INT, row.HEP))
  names(df) <- c("B.mean", "B.var", "B.skew", "B.kurt")
  df$graph <- c("INT", "HEP")
  return(df)
}
row.names(B.stats) <- 1:800
B.stats$graph <- as.factor(B.stats$graph)

# calculate statistics of closeness centralities (Case 5)
C.stats <- foreach(i=1:400, .combine=rbind, .packages=c("igraph", "e1071")) %dopar% {
  C.INT <- closeness(internet.subgraphs[[i]], mode="all", normalized=T)
  row.INT <- c(mean(C.INT), var(C.INT), skewness(C.INT, type=3), kurtosis(C.INT, type=3))
  row.INT[is.nan(row.INT)] <- 0
  C.HEP <- closeness(hep.subgraphs[[i]], mode="all", normalized=T)
  row.HEP <- c(mean(C.HEP), var(C.HEP), skewness(C.HEP, type=3), kurtosis(C.HEP, type=3))
  row.HEP[is.nan(row.HEP)] <- 0
  df <- as.data.frame(rbind(row.INT, row.HEP))
  names(df) <- c("C.mean", "C.var", "C.skew", "C.kurt")
  df$graph <- c("INT", "HEP")
  return(df)
}
row.names(C.stats) <- 1:800
C.stats$graph <- as.factor(C.stats$graph)

# calculate statistics of eigenvector centralities (Case 6)
E.stats <- foreach(i=1:400, .combine=rbind, .packages=c("igraph", "e1071")) %dopar% {
  E.INT <- eigen_centrality(internet.subgraphs[[i]], scale=F)$vector
  row.INT <- c(mean(E.INT), var(E.INT), skewness(E.INT, type=3), kurtosis(E.INT, type=3))
  row.INT[is.nan(row.INT)] <- 0
  E.HEP <- eigen_centrality(hep.subgraphs[[i]], scale=F)$vector
  row.HEP <- c(mean(E.HEP), var(E.HEP), skewness(E.HEP, type=3), kurtosis(E.HEP, type=3))
  row.HEP[is.nan(row.HEP)] <- 0
  df <- as.data.frame(rbind(row.INT, row.HEP))
  names(df) <- c("E.mean", "E.var", "E.skew", "E.kurt")
  df$graph <- c("INT", "HEP")
  return(df)
}
row.names(E.stats) <- 1:800
E.stats$graph <- as.factor(E.stats$graph)

# data frame consists of statistics of degree centralities and clusterinc coefficients (Case 1)
D.CC.stats <- cbind(D.stats[-5], CC.stats)

stats.list <- list(D.CC.stats, D.stats, CC.stats, B.stats, C.stats, E.stats)

get.ND <- function(stats, train.ratio=0.6, seed=100) {
  # function that calculates normalized distances for internet and hep data
  n <- nrow(stats)
  n.train <- round(n*train.ratio)
  set.seed(seed)
  train.index <- sample(1:n, n.train)
  stats.train <- stats[train.index,]
  factor.index <- ncol(stats)
  INT.index <- which(stats.train$graph == "INT")
  train.INT <- stats.train[INT.index, -factor.index]
  train.HEP <- stats.train[-INT.index, -factor.index]
  mean.vector <- list(INT=apply(train.INT, 2, mean), HEP=apply(train.HEP, 2, mean))
  
  mean.INT.mat <- matrix(mean.vector$INT, n-n.train, factor.index-1, byrow=T)
  ND.INT.i <- abs(stats[-train.index, -factor.index] - mean.INT.mat) / mean.INT.mat
  mean.HEP.mat <- matrix(mean.vector$HEP, n-n.train, factor.index-1, byrow=T)
  ND.HEP.i <- abs(stats[-train.index, -factor.index] - mean.HEP.mat) / mean.HEP.mat
  
  ND.test <- data.frame(ND.INT=apply(ND.INT.i, 1, sum), ND.HEP=apply(ND.HEP.i, 1, sum))
  ND.test$graph <- stats$graph[-train.index]
  return(ND.test)
}

ND.list <- lapply(stats.list, get.ND) # ND with 4 statistics for 6 cases

# ND without kurtosis for 6 cases
ND.list.without.kurt <- list()
for (i in 1:6) {
  if(i == 1) {
    ND.list.without.kurt[[i]] <- get.ND(stats.list[[i]][,-c(4,8)])
    }
  else {
    ND.list.without.kurt[[i]] <- get.ND(stats.list[[i]][,-4])
    }
}

# ND without skewness and kurtosis for 6 cases
ND.list.without.skew.kurt <- list()
for (i in 1:6) {
  if(i == 1) {
    ND.list.without.skew.kurt[[i]] <- get.ND(stats.list[[i]][,-c(3,4,7,8)])
  }
  else {
    ND.list.without.skew.kurt[[i]] <- get.ND(stats.list[[i]][,-c(3,4)])
  }
}

# plots of result with 4 statistics for 6 cases
plot.list1 <- list()
for (i in 1:6) {
  plot.title <- paste("case", as.character(i), "with all statistics")
  plot.list1[[i]] <- as.grob(ggplot(data=ND.list[[i]], aes(x=ND.INT, y=ND.HEP, shape=graph, color=graph)) + 
                             geom_point() + ggtitle(plot.title))
}
marrangeGrob(plot.list1, nrow=2, ncol=3)

# plots of result without kurtosis for 6 cases
plot.list2 <- list()
for (i in 1:6) {
  plot.title <- paste("case", as.character(i), "without kurtosis")
  plot.list2[[i]] <- as.grob(ggplot(data=ND.list.without.kurt[[i]], aes(x=ND.INT, y=ND.HEP, shape=graph, color=graph)) + 
                               geom_point() + ggtitle(plot.title))
}
marrangeGrob(plot.list2, nrow=2, ncol=3)

# plots of result without skewness and kurtosis for 6 cases
plot.list3 <- list()
for (i in 1:6) {
  plot.title <- paste("case", as.character(i), "without skewness & kurtosis")
  plot.list3[[i]] <- as.grob(ggplot(data=ND.list.without.skew.kurt[[i]], aes(x=ND.INT, y=ND.HEP, shape=graph, color=graph)) + 
                               geom_point() + ggtitle(plot.title))
}
marrangeGrob(plot.list3, nrow=2, ncol=3)
