#we are going to follow this tutorial: https://joey711.github.io/phyloseq/plot_ordination-examples.html

##1.Load packages
library(phyloseq)
library(ggplot2)

##2.load the phyloseq object from dada2 analysis 
ps <- readRDS("ps.rds")

##3.Normalize the counts (convert to proportions)
ps.prop <- transform_sample_counts(ps, function(x) x/sum(x))

##4.ordinate using bray-curtis distance and PCoA method 
ps.ord <- ordinate(ps.prop, "PCoA", "bray")

##5.visualization
plot_ordination(ps.prop, ps.ord, type="samples") +
  stat_ellipse(type="t", linetype=2, level=0.95) +
  geom_point(size=3) +
  labs(title="PCoA with Bray Curtis distance") + 
  theme_minimal()

