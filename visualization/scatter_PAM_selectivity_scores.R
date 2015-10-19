# Visualize the PAM selectivity data from PyRosetta molecular dynamics and
# Kleinstiver et al. (2015) experiments by plotting a scatter of the two scores
#
# This is the graphic used in the presentation and on the wiki, replacing the
# cluster-based graph on the poster
library(ggplot2)    # lib for plotting
library(plyr)       # lib to manipulate data frames

# The CSV contains 5 columns: the PAM sequence, then scores for 'wt' (NGG-
# binding) Cas9 from PyRosetta and Kleinstiver et al. and then scores for 'eqr'
# (NGAG-binding) Cas9, again from two sources.
pam.selectivity <- read.csv("./normalized_PAM_selectivity_scores.csv")

# Set up theme
scatter.theme <- theme(
  # Set up transparent background
  panel.background = element_blank(),
  plot.background = element_rect(fill = "transparent",colour = NA),
  # Set up axis
  axis.text = element_text(size = 14, colour ="#DDDDDD", family="Calibri"),
  axis.title.x =  element_text(size = 18, colour ="#DDDDDD", family="Calibri"),
  axis.title.y =  element_text(size = 18, colour ="#DDDDDD", family="Calibri"),
  axis.ticks = element_line(color="#DDDDDD"),
  # Set up text
  panel.grid.major.y = element_line(colour = "#DDDDDD"),
  panel.grid.minor.y = element_blank(),
  panel.grid.major.x = element_line(colour = "#DDDDDD"),
  panel.grid.minor.x = element_blank()
)


# Plot wild type
plot <- ggplot(pam.selectivity, aes(x=PyRosetta_wt, y=Kleinstiver_wt)) +
  geom_point(colour="#fbb67a33", size=9) +
  stat_smooth(method="lm", se=FALSE, colour="#DDDDDD", size=0.9) + 
  scale_y_continuous("Kleinstiver Empirical Selectivity", expand = c(0.005,0.005)) + 
  scale_x_continuous("PyRosetta Molecular Binding Score", expand = c(0.005,0.005)) + 
  scatter.theme
ggsave(file="scatter_PAM_selectivity.wt.svg", plot, width=6, height=4, bg="transparent")

plot <- ggplot(pam.selectivity, aes(x=PyRosetta_eqr, y=Kleinstiver_eqr)) +
  geom_point(colour="#fbb67a33", size=9) +
  stat_smooth(method="lm", se=FALSE, colour="#DDDDDD", size=0.9) + 
  scale_y_continuous("Kleinstiver Empirical Selectivity", expand = c(0.005,0.005)) + 
  scale_x_continuous("PyRosetta Molecular Binding Score", expand = c(0.005,0.005)) + 
  scatter.theme
ggsave(file="scatter_PAM_selectivity.eqr.svg", plot, width=6, height=4, bg="transparent")

# Check correlation with permutation (these stats were manually added to figures)
permute.cor <- function() {
  tmp <- data.frame( PyRosetta_wt=pam.selectivity$PyRosetta_wt,
                    Kleinstiver_wt=sample( pam.selectivity$Kleinstiver_wt ))
  cor(tmp$PyRosetta_wt, tmp$Kleinstiver_wt)
}
cor.permuted <- c( cor(pam.selectivity$PyRosetta_wt, pam.selectivity$Kleinstiver_wt), 
          replicate(9999, permute.cor()) )

cor.test(pam.selectivity$PyRosetta_wt, pam.selectivity$Kleinstiver_wt)

cor.test(pam.selectivity$PyRosetta_eqr, pam.selectivity$Kleinstiver_eqr)
cor.permuted <- c( cor(pam.selectivity$PyRosetta_eqr, pam.selectivity$Kleinstiver_eqr), 
                   replicate(9999, permute.cor()) )
