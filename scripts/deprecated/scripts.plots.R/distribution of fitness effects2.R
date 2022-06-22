dfe <- rexp(5000,rate=3.5)-1
dfe <- dfe[dfe<=1]
dat<- with(density(dfe), data.frame(x, y))

library(ggplot2)
ggplot(data = dat, mapping = aes(x = x, y = y)) +
  geom_line()+
  geom_area(mapping = aes(x = ifelse(x < 0 , x, 0)), fill = "dark green") +
  geom_area(mapping = aes(x = ifelse(x > 0 , x, 0)), fill = "purple") +
  xlim(-1, 1) +
  ylim(0,3) +
  ggtitle("Figure 1: Distribution of Fitness Effects") +
  xlab("Selection Coefficeint") +
  ylab("Proportion")

# examine the looks of our mutations
hist(dfe,breaks=100, xlim=c(-1,1))
abline(v=0)

