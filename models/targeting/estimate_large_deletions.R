# data from the large multiplexed deletions of Table 5 in
#http://www.nature.com/ncomms/2015/150218/ncomms7244/fig_tab/ncomms7244_T1.html
# which gives % large deletions at 3 days (4320 min) and at 10 days (14400 min)
percent_intact_4320 <- rep(100, 10) - c(6.6, 10.3, 11.9, 12.4, 16.1, 1.3, 13.2,
                                      6.8, 22.5, 26.4)
percent_intact_14400 <- rep(100, 9) - c(9.3, 14, 14.4, 13.3, 16.9, 11, 7.1,
                                        20.9, 24.7)

quantile_4320 <- quantile(percent_intact_4320)
quantile_14400 <- quantile(percent_intact_14400)

# Find data in lower quartile, adding 3 data points for time = 0
lower_4320 <- percent_intact_4320[percent_intact_4320 <= quantile_4320["25%"]]
lower_14400 <- percent_intact_4320[percent_intact_4320 <= quantile_4320["25%"]]

del_data <- data.frame(minutes = c(rep(0.00001,3), rep(4320,length(lower_4320)),
                                   rep(14400, length(lower_14400))),
                       percent_intact = c(rep(100, 3), lower_4320, lower_14400),
                       quartile = rep("lower quartile",
                                      3 +length(lower_4320) + length(lower_14400)));

# Add equivalent data from upper quartile
upper_4320 <- percent_intact_4320[percent_intact_4320 >= quantile_4320["75%"]]
upper_14400 <- percent_intact_4320[percent_intact_4320 >= quantile_4320["75%"]]

del_data <- rbind(del_data,
            data.frame(minutes = c(rep(0.00001,3), rep(4320,length(upper_4320)),
                                   rep(14400, length(upper_14400))),
                      percent_intact = c(rep(100, 3), upper_4320, upper_14400),
                      quartile = rep("upper quartile",
                                      3 +length(upper_4320) + length(upper_14400))))

# Plot with exponential fits
library(ggplot2)
ggplot(data=del_data, aes(x=minutes, y=percent_intact)) + 
  ylab("% large deletion targets intact") + xlab("Time (minutes)") + 
  ggtitle("Exponential Fit to Lower and Upper Quartiles of DSB Data") +
  geom_point() + geom_smooth(method = "glm", family = gaussian(link="log")) +
  facet_grid(. ~ quartile)

# Find decay constant via exponential models
model_lower <- nls(percent_intact ~ I(exp(a+minutes*b)),
                   data = del_data[which(del_data$quartile == "lower quartile"),],
                   start = list(a=0, b=0), trace = T)
model_upper <- nls(percent_intact ~ I(exp(a+minutes*b)),
                   data = del_data[which(del_data$quartile == "upper quartile"),],
                   start = list(a=0, b=0), trace = T)

decay_lower <- coef(model_lower)[2]
decay_upper <- coef(model_upper)[2]



