library(iNEXT)

A_asympt<-iNEXT(results,
				q = 0, endpoint = 10000, knots = 100, se = F)
A_1_asympt<-iNEXT(results,
                q = 1, endpoint = 10000, knots = 100, se = F)
A_2_asympt<-iNEXT(results,
                q = 2, endpoint = 10000, knots = 100, se = F)

png(filename = "asymp_plot.png",
    width = 4*2, height = 3, units = "in", pointsize = 10,
    bg = "white", res = 1000)
par(mar = c(4, 3, 3, 2) , mfrow = c(1, 3))

X<-c(0, 10000)
Y<-c(000, 1000)
plot(X, Y,
     type = "n",
	 xlab = "",
	 ylab = "")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted", lwd = 0.5)
plot.iNEXT(A_asympt,
           show.legend = T, add = T, SITE = host_names, show.main = T)


plot(X ,Y,
     type = "n", xlab = "", ylab = "")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted", lwd = 0.5)
plot.iNEXT(A_1_asympt,
           show.legend = F, add = T, SITE = host_names, show.main = T)
 

plot(X, Y,
     type = "n", xlab = "", ylab = "")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted", lwd = 0.5)
plot.iNEXT(A_2_asympt,
           show.legend = F, add = T, SITE = host_names, show.main = T)

dev.off()
