# Create time series of P6/time from 1000 simualtions of the targeting model,
# with an exponential decay fit.
library(ggplot2)      # lib for plotting
library(reshape2)     # lib for reshaping data (tall-narrow <-> short-wide)

# The CSV contains 10 columns (time and the 10 genome regions defined in
# .../models/targeting/init_genome_camv.py. In the simulation only P6 has been
# targeted, so other columns are not converted to percent intact
data_1000 <- read.csv("./gene_decay_P6_targeting.csv")
data_1000$gene_P6 <- data_1000$gene_P6 / 1000
data_1000$time_h <- data_1000$time / (60*60)

# Fit an exponential decay model
model <- nls(gene_P6 ~ I(exp(time*b)),
             data = data_1000, start = list(b=0), trace = T)
decay <- coef(model)
data_1000$fit_P6 <- predict(model)
data_1000_melted <- melt(data_1000, measure.vars = c("gene_P6", "fit_P6"), 
                    id.vars=c("time"), variable.name="Simulation",
                    value.name="gene_P6")
data_1000_melted$time_h <- data_1000$time / (60*60)

# Plot with exponential fit
ggplot(data=data_1000_melted, aes(x=time_h, y=gene_P6)) +
       geom_line(aes(colour=Simulation)) + 
       ylab("% Gene P6 Active (1000 Simulations)") + xlab("Time (hours)") +
       scale_x_continuous(breaks = seq(0, 18, 2)) +
       theme(panel.background = element_blank()) + 
       scale_colour_manual(values = c("#373b40","#92cc78"))

# Data written to a file and cleaned up in Igor Pro for final figures
write.csv(data_1000[, c("time_h", "gene_P6", "fit_P6")], 
          file="./fit_gene_decay_P6_targeting.csv", row.names=FALSE)
