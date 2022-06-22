dfe <- rexp(5000,rate=3.5)-1
dfe <- dfe[dfe<=1]
dat <- with(density(dfe), data.frame(x, y))



ggplot(data = dat, mapping = aes(x = x, y = y)) +
  geom_line()+
  geom_area(mapping = aes(x = ifelse(x < 0 , x, 0)), fill = "red") +
  geom_area(mapping = aes(x = ifelse(x > 0 , x, 0)), fill = "blue") +
  xlim(-1, 1) +
  ylim(0,3) +
  xlab("selection coefficeint") +
  ylab("density")
