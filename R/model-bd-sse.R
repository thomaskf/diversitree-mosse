## Models should provide:
##   1. make
##   2. info
##   3. make.cache, including initial tip conditions
##   4. initial.conditions(init, pars,t, idx)
##   5. rootfunc(res, pars, ...)

## Common other functions include:
##   stationary.freq
##   starting.point
##   branches

## 1: make
make.bd.sse <- function (tree, sampling.f=NULL, unresolved=NULL, times=NULL, control=list()) {
	#this generates 1) the order of nodes to calculate, 2) the children of each node, 3) the starting time and bl of each branch, 4) what kind of calculation, 1=D*, 2=D
    control <- check.control.bd(control, times)
	cache <- make.cache.bd.ode(tree, sampling.f, unresolved)
	all.branches <- make.all.branches.bd.sse(cache)
	rootfunc <- function (res, pars, t=tree$root.depth, condition.surv) {
  		vals <- res$vals
  		lq <- res$lq
  		d.root <- vals[2] #the root must be in state R
  
  		if ( condition.surv ) { #condition the likelihood on that the survivial of the root
    		b <- pars[1]
    		mu <- pars[2]
    		e.root <- vals[1]
    		er.root <- (mu-mu*exp((b-mu)*t)) / (mu-b*exp((b-mu)*t))
    		d.root <- d.root / (b * (1 - e.root) * (1 - er.root))
  		}

  		loglik <- log(d.root) + sum(lq)

  		loglik
  	}
	ll <- function (pars, t=tree$root.depth, condition.surv=T) {
		ans <- all.branches(pars)
		rootfunc(ans, pars, t, condition.surv)
	}
	class(ll) <- c("bdsse","dtlik","function")
	ll
}

## 2: info
make.info.bd.sse <- function (phy) {
	list(name="bd.sse",
		 ## Parameters:
		 np=2L,
		 argnames=default.argnames.bd(),
		 ## Variables:
		 ny=2L,
		 k=1L,
		 idx.e=1L,
		 idx.d=2L,
		 ## Phylogeny:
		 phy=phy,
		 ## Inference:
		 ml.default="subplex",
		 mcmc.lowerzero=TRUE)
	}

######################################################################
## Extra core stuff:
#this is similar to make.all.branches.dtlik
make.all.branches.bd.sse <- function (cache) {
	branches <- function(y, len, pars, t0, idx) {
        b <- pars[1]
		mu <- pars[2]
		r <- b-mu
		z <- exp(len * r)
		e0 <- y[1]
  		A <- (z * r * r)/(z * b - mu + (1-z)*b*e0)^2
  		ret <- rbind((mu + z*(e0 - 1)*mu - b*e0) / (mu + z*(e0 - 1)*b - b*e0),A * y[2])

        lq <- log(ret[2,])
        ret[2,] <- 1
  		list(lq, ret)
		}
	function (pars) {
		cache$y <- initial.tip.bd.ode(cache)
		all.branches.matrix(pars,cache,initial.conditions=initial.conditions.bd.ode,branches)
	}
}
