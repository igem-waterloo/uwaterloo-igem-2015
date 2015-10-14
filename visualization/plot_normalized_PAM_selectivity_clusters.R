# Visualize the PAM selectivity data from PyRosetta molecular dynamics and
# Kleinstiver et al. (2015) experiments by plotting the cluster centroids of
# the normalized selectivity scores
#
# This is the graphic used on the poster- it is agreed to be somewhat unclear
# and was replaced in the presentation and on the wiki
library(ggplot2)    # lib for plotting

# The CSV contains 11 columns, including 'SOURCE', either known (from
# Kleinstiver et al.) or simulated and either 'wt' (NGG-binding) or 'eqr'
# (NGAG-binding); CLUSTER_NORM, the normalized score of the cluster centroid;
# and TOTAL_PAMS, the number of 4-nucleotide PAM sequences in the cluster.
clusters <- read.csv("./normalized_PAM_selectivity_clusters.csv")

clusters.graph <-
  ggplot(clusters, aes(log(CLUSTER_NORM + exp(1)), SOURCE)) +
  geom_point(aes(size=TOTAL_PAMs, colour="fbb67a")) +
  scale_x_continuous(limits=c(0.99, 1.33)) + 
  scale_size(range = c(4, 12)) 

ggsave(file="PAM_selectivity_clusters.svg", clusters.graph, width=7, height=3, bg="transparent")