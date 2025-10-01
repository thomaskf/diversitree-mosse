library(ape)
library(diversitree)

lambda <- sigmoid.x
mu <- constant.x
drift <- 0
diffusion <- 0.001
sd <- 0.001
# args <- list(lambda=1:4, mu=5, drift=6, diffusion=7)
pars <- c(0.1,0.2,0,2.5,0.03, drift, diffusion)#y0,y1,x0,r,mu,drift,diffusion

phy <- read.tree(text="(A:0.4,B:0.4);")
types <- c(1,2)
names(types) <- phy$tip.label
states <- c(0.15,0.1)
names(states) <- phy$tip.label
states.sd <- rep(0.01,length(phy$tip.label))
names(states.sd) <- phy$tip.label

# not normalized
# Q <- t(matrix(c(-3/4,1/4,1/4,1/4,
#                1/4,-3/4,1/4,1/4,
#                1/4,1/4,-3/4,1/4,
#                1/4,1/4,1/4,-3/4),4,4))

Q <- t(matrix(c(-1,1/3,1/3,1/3,
                1/3,-1,1/3,1/3,
                1/3,1/3,-1,1/3,
                1/3,1/3,1/3,-1),4,4))


control.C.1 <- list(xmax=1)
lik.C.1 <- make.mosse(phy, types,ntypes=4,states,states.sd,sigmoid.x, constant.x, Q,
                      control.C.1)
lik.C.1(pars)  

