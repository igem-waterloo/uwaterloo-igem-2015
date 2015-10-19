# Create histograms of sgRNA exchange flow cytometry data
# The CSV contains 3 trials (10000 measurements from 3 separate days) for each
# of: the positive RFP control, the dCas9 knockdown with modified sgRNA
# structure and the dCas9 knockdown with the original sgRNA structure
library(ggplot2)      # lib for plotting
library(scales)       # lib for log scales in the plot ticks
library(extrafont)    # lib to include Calibri in plots

data <- read.csv("./sgrna-exchange-flow-data.csv")
data$Date <- factor(data$Date)

# Set up beautiful theme 
density.theme <- theme(
    # Nice font, matching colours for text+ticks
    axis.title =  element_text(size = 30, colour ="white", family="Calibri"),
    axis.text = element_text(size = 24, colour ="#DDDDDD", family="Calibri"),
    axis.ticks = element_line(colour ="#DDDDDD"),
    # Set grid to nice colours
    panel.grid.major = element_line(colour = "#DDDDDD"),
    panel.grid.minor.y = element_line(colour = "#DDDDDD"),
    panel.grid.minor.x = element_blank(),
    # Set up transparent background
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )

rfp.control <- 
  ggplot(data, aes(RFP.Control, group=Date)) +
  geom_density(alpha=0.4, lwd=0.8, adjust=0.2, fill="#DDDDDD", linetype="blank") +
  # Natural log scale for x axis, set breaks to display as powers of e
  scale_x_continuous("RFP Intensity", limits=c(exp(3), exp(12.5)), expand = c(0,0),
                     trans = "log",
                     breaks = trans_breaks("log", function(x) exp(x), n=5),
                     labels = trans_format("log", math_format(e^.x))) + 
  # Hardcoded breaks
  scale_y_continuous("Normalized Frequency", limits=c(0.0, 0.95), expand = c(0,0),
                     breaks=c(0.0, 0.2, 0.4, 0.6, 0.8)) + 
  annotation_logticks(sides = "b", color="#DDDDDD") + 
  density.theme

sgrna.control <- 
  ggplot(data, aes(sgRNA.Control, group=Date)) +
  geom_density(alpha=0.4, lwd=0.8, adjust=0.2, fill="#BA93C8", linetype="blank") +
  # Natural log scale for x axis, set breaks to display as powers of e
  scale_x_continuous("RFP Intensity", limits=c(exp(3), exp(12.5)), expand = c(0,0),
                     trans = "log",
                     breaks = trans_breaks("log", function(x) exp(x), n=5),
                     labels = trans_format("log", math_format(e^.x))) + 
  # Hardcoded breaks
  scale_y_continuous("Normalized Frequency", limits=c(0.0, 0.95), expand = c(0,0),
                     breaks=c(0.0, 0.2, 0.4, 0.6, 0.8)) + 
  annotation_logticks(sides = "b", color="#DDDDDD")  + 
  density.theme

sgrna.exchange <-   
  ggplot(data, aes(sgRNA.Modified, group=Date)) +
  geom_density(alpha=0.4, lwd=0.8, adjust=0.2, fill="#79bcc7", linetype="blank") +
  # Natural log scale for x axis, set breaks to display as powers of e
  scale_x_continuous("RFP Intensity", limits=c(exp(3), exp(12.5)), expand = c(0,0),
                     trans = "log",
                     breaks = trans_breaks("log", function(x) exp(x), n=5),
                     labels = trans_format("log", math_format(e^.x))) + 
  # Hardcoded breaks
  scale_y_continuous("Normalized Frequency", limits=c(0.0, 0.9), expand = c(0,0),
                     breaks=c(0.0, 0.2, 0.4, 0.6, 0.8)) + 
  annotation_logticks(sides = "b", color="#DDDDDD")  + 
  density.theme


ggsave(file="rfp.control.test.svg", plot=rfp.control, width=10, height=8, bg="transparent")
ggsave(file="sgrna.control.test.svg", plot=sgrna.control, width=10, height=8, bg="transparent")
ggsave(file="sgrna.exchange.test.svg", plot=sgrna.exchange, width=10, height=8, bg="transparent")
