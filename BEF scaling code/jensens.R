library(dplyr)

s <- 1:40
f <- function(x, b=0.53) x^b


fun_lost <- function(s, ch, b=0.53)
  (f(s) - (f(s-ch,b)+f(s+ch,b))/2)/f(s)

####Plot Jensen's inequality
m <- 20 #mean richness
ch <- 10 #deviation from mean - loss or gain

par(mfrow=c(1,2))

#
plot(s, f(s), type="l", xlab="Number of Species", ylab="Productivity", xlim=c(-5, 50))
segments(m-ch, f(m-ch), m+ch, f(m+ch), col="red")
points(m-ch, f(m-ch), col="red", pch=19)
points(m+ch, f(m+ch), col="red", pch=19)
segments(20,0,20,f(m), lty=2, col="blue")
text(20.0, 3, "Average Richness", pos=4)
text(m-ch+0.2, 4, "Function After\nLoss", pos=2)
text(m+ch+0.5, 5.3, "Function After\nGain", pos=4)

plot(s, f(s), type="l", xlim=c(19,21), ylim=c(4.5, 5.1), xlab="Number of Species", ylab="Productivity")
segments(m-ch, f(m-ch), m+ch, f(m+ch), col="red")
points(m, f(m), col="black", pch=19)
points(m, (f(m-ch)+f(m+ch))/2, col="red", pch=19)
segments(m, (f(m-ch)+f(m+ch))/2, m, f(m), lty=2, col="blue", lwd=2)
text(20.01, f(m)-0.04, "Loss of Productivity", pos=4)

par(mfrow=c(1,1))

ggplot() +
  geom_line(mapping=aes(x=s, y=f(s))) +
  geom_segment(aes(x=m-ch, xend=m+ch, y=f(m-ch), yend=f(m+ch)), color="red")+
  geom_segment(aes(x=m, xend=m, y=(f(m-ch)+f(m+ch))/2, yend=f(m)), color="blue", lty=1)+
  geom_point(aes(x=m, y=f(m)), color="black") +
  geom_point(aes(x=m, y=(f(m-ch)+f(m+ch))/2)) +
  xlim(c(m-ch,m+ch)) +
  ylim(c(f(m-ch), f(m+ch)))

#Plot loss of function by initial diversity
library(ggplot2)

jensens_frame <- data.frame(expand.grid(change = c(3,5,9), 
                                        initial_s = 5:30, 
                                        b=c(0.25, 0.47, 0.53))) %>%
  mutate(function_lost = fun_lost(initial_s, change, b)) %>%
  mutate(b=paste0("b = ", b)) %>%
  mutate(change=paste0("± ", change, " species")) 

ggplot(data=jensens_frame, aes(x=initial_s, y=function_lost)) +
  geom_line(lwd=1.2) +
  facet_grid(change ~ b) +
  theme_bw(base_size=14) +
  ylab("Propotion of Productivity Lost") +
  xlab("Mean Diversity")

##### Look at this another way

jensens_loss_frame <- data.frame(expand.grid(change = seq(1,10,.01), 
                                        initial_s = 30, 
                                        b=c(0.25, 0.47, 0.53))) %>%
  mutate(function_lost = fun_lost(initial_s, change, b)) %>%
  mutate(b=paste0("b = ", b)) 

ggplot(data=jensens_loss_frame, aes(x=change, y=function_lost)) +
  geom_line() +
  facet_wrap(~b, scale="free_y")+
  theme_bw(base_size=14) +
  ylab("Propotion of Productivity Lost\n") +
  xlab("Gain/Loss of Species") +
  scale_x_continuous(breaks=seq(0,10,2))
