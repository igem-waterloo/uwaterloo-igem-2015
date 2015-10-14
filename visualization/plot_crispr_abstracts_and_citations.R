# Line graph showing PubMed entries mentioning CRISPR vs. iGEM abstracts
# mentioning CRISPR
library(ggplot2)      # lib for plotting
library(plyr)         # lib to manipulate data frames
library(jsonlite)     # lib to fetch JSON data from Aalto-Helsinki database
library(httr)         # lib to open connection to fetch online data
library(RISmed)       # lib to fetch data from Pub Med


## PubMed Data
# Fetch PubMed articles on CRISPR, based on the tutorial
# https://freshbiostats.wordpress.com/2013/12/03/analysis-of-pubmed-search-results-using-r/
query = "(CRISPR OR CRISPRi OR Cas9 OR dCas9)"
crispr.search <- EUtilsSummary(query, type="esearch", db = "pubmed",
                               mindate=2000, maxdate=2014, retmax=30000)
QueryCount(crispr.search)
crispr.records <- EUtilsGet(crispr.search)
years <- YearAccepted(crispr.records)
crispr.pubm.freq <- count(years)

# Fetch total number of publications/year for normalizing
total <- NULL
for (i in 2000:2014){
  peryear <- EUtilsSummary("", type="esearch", db="pubmed", mindate=i, maxdate=i)
  total[i] <- QueryCount(peryear)
}
year <- 2000:2014
total.pubs.count<- as.data.frame(cbind(year,total[year]))
crispr.pubm.freq$freq <- apply(crispr.pubm.freq, 1, function(x) {x["freq"] / total.pubs.count[total.pubs.count$year == x[1], "V2"]} )

## iGEM Data
# Fetch abstracts stored in the 2014 Aalto-Helsinki iGEM Team-Seeker Tool
igem.abstracts <- fromJSON("https://raw.githubusercontent.com/iGEM-QSF/iGEM-Team-Seeker/master/app/teams.json")

# Find abstracts that match CRISPR or Cas9 projects
crispr.igem.rows <- unique(c(
  grep("CRISPR|CRISPRi|Cas9|dCas9", igem.abstracts$Description, ignore.case=T),
  grep("CRISPR|CRISPRi|Cas9|dCas9", igem.abstracts$Abstract, ignore.case=T)))
crispr.igem.freq <- count(igem.abstracts[crispr.igem.rows, "Year"])

# Merge for plotting
crispr.freq <- rbind(crispr.pubm.freq, crispr.igem.freq)
colnames(crispr.freq) <- c("Year", "Value")
crispr.freq$Source <- c(rep("PubMed", nrow(crispr.pubm.freq)), rep("iGEM", nrow(crispr.igem.freq)))
crispr.freq <- crispr.freq[!is.na(crispr.freq$Year),  ]
crispr.freq$Value <- unlist(crispr.freq$Value)
crispr.freq[crispr.freq$Source=="PubMed","Value"] <- crispr.freq[crispr.freq$Source=="PubMed","Value"]*100000

# CRISPR per 100000 PubMed entries
plot <-
  ggplot(crispr.freq, aes(x=Year, y=Value, colour=Source)) + geom_line(size=1) +
  scale_colour_manual(values = c("#92cc78", "#79bcc7")) + 
  scale_y_continuous("abstracts | citations/100,000 PubMed Entries", expand = c(0,0), breaks=seq(0,30,5)) + 
  scale_x_continuous("Year", expand = c(0,0), breaks=seq(2003,2014)) +
  theme(
    # Set up transparent background
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA),
    legend.background = element_blank(),
    axis.text = element_text(size = 14, colour ="#DDDDDD", family="Calibri"),
    axis.title.x =  element_text(size = 16, colour ="#DDDDDD", family="Calibri"),
    axis.title.y =  element_text(size = 14, colour ="#DDDDDD", family="Calibri"),
    panel.grid.major.y = element_line(colour = "#DDDDDD"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=18, colour ="#DDDDDD", family="Calibri"),
    legend.position="bottom"
  )
ggsave(file="crispr.pubmed.svg", plot, width=7, height=4, bg="transparent")
