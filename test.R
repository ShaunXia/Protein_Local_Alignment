library(reshape2)
library(ggplot2)
#rnorm(4)
mat = matrix( round(runif(4,0,1)),2)
m = melt(mat)
g = ggplot(m, aes(x=Var1, y=Var2, fill=value,group=factor(Var2),colour=factor(Var2)))+
  xlab('X-labels')+
  ylab("Y-labels")

g1=g+geom_tile(); print(g1)