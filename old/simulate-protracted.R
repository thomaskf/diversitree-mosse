make.tree.protracted <- function(pars, max.taxa=Inf, max.t=Inf, x0=1,
                            single.lineage=TRUE) {
  #pars <- c(b,mu,lambda)   
  extinct <- FALSE
  split   <- FALSE
  parent <- 0
  species <- 1
  
  r <- sum(pars)
  
  len <- 0
  t <- 0
  hist <- list()
  
  if ( single.lineage ) {
  	states <- x0
    n.taxa <- lineages <- 1
    start <- 0    
  } else {
    stop("Nope.")    
  }
  
  while ( n.taxa <= max.taxa && n.taxa > 0 ) {
    ## When does an event happen?
    r.n <- r * n.taxa
    r.tot <- sum(r.n)
    dt <- rexp(1, r.tot)
    t <- t + dt

    if ( t > max.t ) {
      dt <- dt - (t - max.t)
      len[lineages] <- len[lineages] + dt
      t <- max.t
      break
    }

    len[lineages] <- len[lineages] + dt

    ## Pick a lineage for an event to happen to:
    lineage.i <- sample(n.taxa, 1)
    lineage <- lineages[lineage.i]

    type <- sample(3, 1, FALSE, pars)

    if ( type == 1 ) {
      ## Speciating:
      if ( n.taxa == max.taxa )
        ## Don't add this one.
        break
      new.i <- length(extinct) + 1:2
      split[lineage] <- TRUE
      extinct[new.i] <- split[new.i] <- FALSE
      states[new.i] <-c(0,0)
      parent[new.i] <- lineage
      species[new.i] <- species[lineage]
      start[new.i] <- t      
      len[new.i] <- 0 
      n.taxa <- n.taxa + 1      
      lineages <- which(!split & !extinct)
    } else if ( type == 2 ) {
      ## Extinct
      extinct[lineage] <- TRUE
      lineages <- which(!split & !extinct)
      n.taxa <- n.taxa - 1
    } else {
      ## become a distinct species:
      species[lineage] <- max(species)+1
      hist[[length(hist)+1]] <- c(lineage, t, 0, 1)
      states[lineage] <- 1
    }
  }

  info <- data.frame(idx=seq_along(extinct), len=len, parent=parent, species=species,
                     start=start, state=states, extinct=extinct,
                     split=split)

  hist <- as.data.frame(do.call(rbind, hist))
  if ( nrow(hist) == 0 )
    hist <- as.data.frame(matrix(NA, 0, 4))
  names(hist) <- c("idx", "t", "from", "to")
  hist$x0 <- info$start[match(hist$idx, info$idx)]
  hist$tc <- hist$t - hist$x0
  
  attr(info, "t") <- t
  attr(info, "hist") <- hist  
  info
}

tree.protracted <- function(pars, max.taxa=Inf, max.t=Inf,
                       include.extinct=FALSE, x0=1,sampling=FALSE) {
  info <- make.tree.protracted(pars, max.taxa, max.t, x0)
  phy <- me.to.ape.protracted(info[-1,], info$state[1], root.depth=attr(info, "t"))
  if (!include.extinct) {phy <- prune.protracted(phy)}
  if (sampling) {
  	species.label <- unique(phy$species)
  	n.species <- length(species.label)
  	tip.species <- sapply(species.label,function(i) is.element(phy$species,i))
  	n.tip.species <- colSums(tip.species)
  	keep.tips <- sapply(1:n.species,function (i) ifelse(n.tip.species[i]==1,which(tip.species[,i]),sample(which(tip.species[,i]),1)))
  	to.drop <- rep(TRUE,length(phy$tip.label))
  	to.drop[keep.tips] <- FALSE
	sampling.f <- 1/n.tip.species
	names(sampling.f) <- phy$tip.label[keep.tips]
	phy <- prune.protracted(phy,to.drop)
	phy$sampling.f <- sampling.f[phy$tip.label]
  }
  phy
}

me.to.ape.protracted <- function(x, root.state, root.depth) {
  if ( nrow(x) == 0 )
    return(NULL)
  Nnode <- sum(!x$split) - 1
  n.tips <- sum(!x$split)
    
  x$idx2 <- NA
  x$idx2[!x$split] <- 1:n.tips
  x$idx2[ x$split] <- order(x$idx[x$split]) + n.tips + 1

  i <- match(x$parent, x$idx)
  x$parent2 <- x$idx2[i]
  x$parent2[is.na(x$parent2)] <- n.tips + 1

  tip.label <- ifelse(subset(x, !split)$extinct,
                      sprintf("ex%d", 1:n.tips),
                      sprintf("sp%d", 1:n.tips))
  node.label <- sprintf("nd%d", 1:Nnode)

  x$name <- NA
  x$name[!x$split] <- tip.label
  ## More useful, but I don't want to clobber anything...
  x$name2 <- c(tip.label, node.label)[x$idx2]

  tip.state <- x$state[match(1:n.tips, x$idx2)]
  names(tip.state) <- tip.label

  node.state <- x$state[match(1:Nnode + n.tips, x$idx2)]
  names(node.state) <- node.label
  node.state["nd1"] <- root.state

  species <- x$species[match(1:n.tips, x$idx2)]
  
  hist <- attr(x, "hist")
  if ( !is.null(hist) ) {
    hist$idx2 <- x$idx2[match(hist$idx, x$idx)]
    hist$name2 <- x$name2[match(hist$idx, x$idx)]
    if ( nrow(hist) > 0 )
      hist <- hist[order(hist$idx2),]
  }

  phy <- reorder(structure(list(edge=cbind(x$parent2, x$idx2),
                                Nnode=Nnode,
                                tip.label=tip.label,
                                species=species,
                                tip.state=tip.state,
                                node.label=node.label,
                                node.state=node.state,
                                edge.length=x$len,
                                root.depth=root.depth,
                                orig=x,
                                hist=hist),
                           class="phylo"))
  phy$edge.state <- x$state[match(phy$edge[,2], x$idx2)]
  phy
}

prune.protracted <- function (phy, to.drop = NULL) 
{
    if (is.null(to.drop)) 
        to.drop <- subset(phy$orig, !split)$extinct
    if (sum(!to.drop) < 2) {
        NULL
    }
    else if (any(to.drop)) {
        phy2 <- drop.tip.fixed(phy, phy$tip.label[to.drop])
        phy2$orig <- phy2$orig[!phy2$orig$extinct, ]
        phy2$tip.state <- phy2$tip.state[!to.drop]
        phy2$species <- phy2$species[!to.drop]
        phy2$node.state <- phy2$node.state[phy2$node.label]
        phy2$hist <- prune.hist(phy, phy2)
        phy2
    }
    else {
        phy
    }
}