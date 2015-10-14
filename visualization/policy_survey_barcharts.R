# Filled bar graphs showing response to two survey questions that were asked on
# two similarly-worded surveys, in which one used language about 'GMO' and
# another language about 'CRISPR/Cas9'. The questions graphed are about
# participant knowledge of the technology and how comfortable they would feel
# eating food made using GMO or CRISPR/Cas9 technology
library(ggplot2)      # lib for plotting

# The CSV contains data from four survey questions, though only 2 questions are
# graphed by this script. There are 4 columns: Question, Response, GMO (showing
# counts of partictpants who gave that response on the GMO survet)
data <- read.csv("./policy_survey_responses.csv")

# Set up theme
barchart.theme <- theme(
  # Set up transparent background
  panel.background = element_blank(),
  plot.background = element_rect(fill = "transparent",colour = NA),
  legend.background = element_blank(),
  # Titles
  legend.title = element_text(size=40, colour ="#FFFFFF", family="Calibri"),
  legend.text = element_text(size=40, colour ="#DDDDDD", family="Calibri"),
  legend.position="bottom",
  # Axis text
  axis.text.x = element_text(size = 40, colour ="#DDDDDD", family="Calibri"),
  axis.text.x = element_blank(),
  axis.title = element_blank(),
  axis.ticks = element_blank()
)

# Data needs to be formatted with columns for questions and responses listed
# in the rows
# Question: #"Have you heard of [GMOs/CRISPR-Cas9]?"
knowledge <- data[which(data$Question=="Knowledge"),]
knowledge.toplot <- data.frame(
  Technology = c(rep("GMO", 60), rep("CRISPR Cas9",61)),
  Response = c(
      rep("Yes", knowledge[which(knowledge$Response=="Yes"), "GMO"]), 
      rep("No", knowledge[which(knowledge$Response=="No"), "GMO"]), 
      rep("Yes", knowledge[which(knowledge$Response=="Yes"), "Gene.Editing"]), 
      rep("No", knowledge[which(knowledge$Response=="No"), "Gene.Editing"]))
  )
knowledge.plot <-
  ggplot(knowledge.toplot, aes(Technology, fill=Response)) +
  geom_bar() +
  scale_fill_manual(values = c("#fbb67a", "#92cc78")) + 
  scale_y_discrete("Percent", expand = c(0,0), breaks=NULL) + 
  scale_x_discrete(expand = c(0,0)) + 
  barchart.theme
ggsave(file="policy_survey_barchart.knowledge.svg", knowledge.plot, width=3, height=7, bg="transparent")

# Question: "How do you feel about eating food with [GM/CRISPR-Cas9] ingredients?"
#"Comfortable", "Neutral", "Avoid", "No Comment"
food <- data[which(data$Question=="Food"),]
food.toplot <- data.frame(
  Technology = c(rep("GMO", 60), rep("CRISPR Cas9",61)),
  "Response" = c(
      rep("Comfortable",food[which(food$Response=="Comfortable"), "GMO"]),
      rep("Neutral",food[which(food$Response=="Neutral"), "GMO"]),
      rep("Avoid",food[which(food$Response=="Avoid"), "GMO"]),
      rep("No Comment",food[which(food$Response=="No comment"), "GMO"]),
      rep("Comfortable",food[which(food$Response=="Comfortable"), "Gene.Editing"]),
      rep("Neutral",food[which(food$Response=="Neutral"), "Gene.Editing"]),
      rep("Avoid",food[which(food$Response=="Avoid"), "Gene.Editing"]),
      rep("No Comment",food[which(food$Response=="No comment"), "Gene.Editing"]))
  )

# Reorder so not alphabetical by label
food.toplot$Response <- relevel(food.toplot$Response, "Neutral")
food.toplot$Response<- relevel(food.toplot$Response, "Avoid")
food.toplot$Response <- relevel(food.toplot$Response, "No Comment")

food.plot <-
  ggplot(food.toplot, aes(Technology, fill=Response)) +
  geom_bar() +
  scale_fill_manual(values = c("#DDDDDD", "#fbb67a", "#79bcc7", "#92cc78")) + 
  scale_y_discrete("Percent", expand = c(0,0), breaks=NULL) + 
  scale_x_discrete(expand = c(0,0)) +
  barchart.theme
  
ggsave(file="policy_survey_barchart.food.svg", plot=food.plot, width=8, height=7, bg="transparent")

