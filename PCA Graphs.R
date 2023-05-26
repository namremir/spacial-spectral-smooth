library(tidyverse)
library(devtools)
install_github("kassambara/factoextra", force=T)
library(factoextra)
library(cowplot)
library(viridis)

data(iris)
iris_species <- iris$Species
iris <- select(iris, -c(5))
iris.pca <- prcomp(iris, scale = F)
iris_pca <- data.frame(iris.pca$x)

t = seq(0, 1, length.out = 101)
iris_pca <- mutate(iris_pca, Species = iris_species)

iris_pca <- iris_pca %>%
  mutate(Row=row_number(), .before=1) %>% 
  uncount(101) %>% 
  add_column(T=rep(t, times=nrow(iris)))

for (i in 1:5){
  for (j in 1:5){
    tempvec = iris_pca$PC1*cos(2*pi*(iris_pca$T+(i-1)/5)) + iris_pca$PC2*sin(2*pi*(iris_pca$T+(i-1)/5)) + 
              iris_pca$PC3*cos(4*pi*(iris_pca$T+(j-1)/5)) + iris_pca$PC4*sin(4*pi*(iris_pca$T+(j-1)/5))
    iris_pca <- mutate(iris_pca, tempcol = tempvec)
    names(iris_pca)[names(iris_pca) == "tempcol"] <- paste("Andrews",i-1,j-1,sep="-")
  }
}

myplots <- list()
for (i in 1:25){
    tempplot <- ggplot(iris_pca) +
      geom_line(aes_q(x=iris_pca$T, y=iris_pca[,7+i], group=iris_pca$Row, color=iris_pca$Species)) +
      scale_color_brewer(palette="Paired", name="Species") +
      labs(x="", y="") +
      theme_minimal()
    
    myplots[[i]] <- tempplot
}

myplots[[3]]
plot_grid(plotlist=myplots)
