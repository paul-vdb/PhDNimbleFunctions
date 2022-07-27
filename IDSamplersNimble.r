sampler_myX <- nimbleFunction(
    name = 'sampler_myX',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        xlim <- control$xlim
        ylim <- control$ylim
        calcNodesAll <- model$getDependencies(target)
        calcNodesNoIDs <- model$getDependencies(target)
        # if(model$getDistribution(target) != 'dX') stop('myX sampler error')
        nodeIndex <- as.numeric(gsub('^X\\[([[:digit:]]+), .*$', '\\1', target))
        zNode <- paste0('z[', nodeIndex, ']')
        dNode <- paste0('d2[', nodeIndex, ', 1:', control$J, ']')
        copyNodes <- c(target, dNode)
    },
    run = function() {
        newx <- runif(1, xlim[1], xlim[2])
        newy <- runif(1, ylim[1], ylim[2])
        if(model[[zNode]] == 0) {
            logMHR <- 0
            jump <- decide(logMHR)
            if(jump) {
                model[[target]] <<- c(newx, newy)
                model$calculate(copyNodes)
                nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodes, logProb = FALSE)
            }
            return()
        }
        anyID <- any(model[['ID']] == nodeIndex)
        if(anyID) { lpcurrent <- model$getLogProb(calcNodesAll)
        } else { lpcurrent <- model$getLogProb(calcNodesNoIDs) }
        model[[target]] <<- c(newx, newy)
        if(anyID) { lpprop <- model$calculate(calcNodesAll)
        } else { lpprop <- model$calculate(calcNodesNoIDs) }
        logMHR <- lpprop - lpcurrent
        jump <- decide(logMHR)
        if(jump) {
            if(anyID) {
                nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesAll, logProb = TRUE)
            } else {
                nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoIDs, logProb = TRUE)
            }
        } else {
            if(anyID) {
                nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesAll, logProb = TRUE)
            } else {
                nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoIDs, logProb = TRUE)
            }
        }
    },
    methods = list( reset = function() {} )
)

# This is my new sampler for multimodal data such as this "label-switching" problem.
# The proposal distribution is a mixture distribution of a uniform and normal.
# The scale is for the normal standard deviation and the
# 'temp' is the temperature, which really is just the mixing proportion of how often to jump.
# If we are running hot, that means mostly uniform, cold is mostly RW.
# Proposal is symmetric so cancels out in MH update.
sampler_myJAM <- nimbleFunction(
    name = 'sampler_myJAM',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        xlim <- extractControlElement(control, 'xlim', c(0,1))		
        ylim <- extractControlElement(control, 'ylim', c(0,1))		
		scale <- extractControlElement(control, 'scale', 1)		
		temp <- extractControlElement(control, 'temp', 0.5)		
		nv <- extractControlElement(control, 'occ', 1)

        calcNodes <- model$getDependencies(target)
        nodeIndex <- as.numeric(gsub(".*?([0-9]+).*", '\\1', target))
		if(nv == 1){ 
			zNode <- paste0('z[', nodeIndex, ']')
		}else{
			v <- as.numeric(sub(".*,\\D*(\\d+).*", "\\1", target))
			zNode <- paste0('z[', nodeIndex,',', v, ']')
		}		
    },
    run = function() {
		p <- runif(1,0,1)
		if(p < temp)
		{
			# Turn up the heat and try a big step
			newX <- c(runif(1, xlim[1], xlim[2]), runif(1, ylim[1], ylim[2]))
		}else{
			# Do a local move.
			newX <- model[[target]] + rnorm(2, 0, scale)
		}
		if(model[[zNode]] == 0) {
            logMHR <- 0
            jump <- decide(logMHR)
            if(jump) {
                model[[target]] <<- newX
                model$calculate(calcNodes)
                nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = FALSE)
            }
            return()
        }
        lpcurrent <- model$getLogProb(calcNodes)
        model[[target]] <<- newX
        lpprop <- model$calculate(calcNodes)
        logMHR <- lpprop - lpcurrent
        jump <- decide(logMHR)
        if(jump) {
			nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        } else {
			nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        }
    },
    methods = list( reset = function() {} )
)


# Small changes here.
# Tried to make it multi-session general but needs some work.
# Seems to be fine for now.
sampler_myBinary <- nimbleFunction(
    name = 'sampler_myBinary',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
        isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)
        calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
        calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
        nodeIndex <- as.numeric(sub("\\D*(\\d+).*", "\\1", target))
		nv <- extractControlElement(control, 'Noccasion', 1)		
		occ <- extractControlElement(control, 'IDoccasion', 1)
		n_obs <- length(model[['ID']])
		if(nv == 1){ 
			IDMatch <- 1:n_obs
		}else{
			v <- as.numeric(sub(".*,\\D*(\\d+).*", "\\1", target))
			IDMatch <- which(occ == v)
		}
    },
    run = function() {
        if(model[[target]] == 1){
			if(any(model[['ID']][IDMatch] == nodeIndex)) return()
        }
		currentLogProb <- model$getLogProb(calcNodes)
        model[[target]] <<- 1 - model[[target]]
        otherLogProb <- model$calculate(calcNodes)
        acceptanceProb <- 1/(exp(currentLogProb - otherLogProb) + 1)
        jump <- (!is.nan(acceptanceProb)) & (runif(1,0,1) < acceptanceProb)
        if(jump) {
            nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        } else {
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        }
    },
    methods = list( reset = function() { } )
)

######################################################
# Basic sampler for categorical distribution
# I've added the pID here of z as a prior so
# I can do a very quick check on pID[k] = 0
# via calculate(target).
# This generalizes really nicely.
#######################################################
sampler_myCategorical <- nimbleFunction(
    name = 'sampler_myCategorical',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
        isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)
        calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
        calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
        M <- control$M
        probs <- numeric(M)
        logProbs <- numeric(M)
    },
    run = function() {
        currentValue <- model[[target]]
        logProbs[currentValue] <<- model$getLogProb(calcNodes)
        for(k in 1:M) {
            if(k != currentValue) {
                model[[target]] <<- k
                logProbPrior <- model$calculate(target)
                if(logProbPrior == -Inf) {
                    logProbs[k] <<- -Inf
                } else {
                    if(is.nan(logProbPrior)) {
                        logProbs[k] <<- -Inf
                    } else {
                        logProbs[k] <<- logProbPrior + model$calculate(calcNodesNoSelf)
                        if(is.nan(logProbs[k])) logProbs[k] <<- -Inf
                    }
                }
            }
        }
        logProbs <<- logProbs - max(logProbs)
        probs <<- exp(logProbs)
        newValue <- rcat(1, probs)
        if(newValue != currentValue) {
            model[[target]] <<- newValue
            model$calculate(calcNodes)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        } else {
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        }
    },
    methods = list( reset = function() { } )
)

######################################################
# Basic sampler for categorical distribution
# sped up for when z is not specifically included.
# Also adding constraints:
# AnimalLink - Observations cannot link to animal index k,
# 	This information may be because we know that animals 1-10 were
#	collared and this observation was not collared.
# cannotlink - Pairwise constraints on the detections themselves that cannot 
# 	be assigned to the same animal.
#######################################################
sampler_mySPIM <- nimbleFunction(
    name = 'sampler_mySPIM',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
		calcNodesZ <- model$getDependencies('z')
        targetNodesAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)		
        isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)
        calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
        calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
        M <- control$M
        probs <- numeric(M)
        logProbs <- numeric(M)		
		cannotlink <- extractControlElement(control, 'cannotlink', 'identity')	# n x n matrix where a 1 indicates i cannot link to detection j.
		nodeIndex <- as.numeric(gsub('[^[:digit:]]', '', targetNodesAsScalar[1]))
		n_grp <- length(targetNodesAsScalar)
    },
    run = function() {
        currentValue <- model[["ID"]][nodeIndex]
        logProbs[currentValue] <<- model$getLogProb(calcNodes)
        for(k in 1:M) {
			model[[target]] <<- k
			logProbPrior <- model$calculate(target)
			if(logProbPrior == -Inf) {
				logProbs[k] <<- -Inf
			} else {
				# Add in the pairwise constraint. 1 means cannot link.
				no_link <- sum(cannotlink[model[['ID']] == k, nodeIndex])	# This is the check to see if there are cannot links.
				if(no_link > 0)
				{
					logProbs[k] <<- -Inf
				}else {					
					values(model, targetNodesAsScalar) <<- rep(k, n_grp)
					logProbs[k] <<- model$calculate(calcNodes)
					if(is.nan(logProbs[k])) logProbs[k] <<- -Inf
				}	
			}
		}
        logProbs <<- logProbs - max(logProbs)
        probs <<- exp(logProbs)
        newValue <- rcat(1, probs)		
        if(newValue != currentValue) {
			values(model, targetNodesAsScalar) <<- rep(newValue, n_grp)	
            model$calculate(calcNodes)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        } else {
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        }
    },
    methods = list( reset = function() { } )
)

sampler_myCategoricalBernoulli <- nimbleFunction(
    name = 'sampler_myCategoricalBernoulli',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
        isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)
        calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
        calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]

		nodeIndex <- as.numeric(gsub('[^[:digit:]]', '', target))

        M <- control$M
		omega <- control$omega
        probs <- numeric(M)
        logProbs <- numeric(M)
    },
    run = function() {
		ID_j <- model[['ID']][omega == omega[nodeIndex]]
        currentValue <- model[[target]]
        logProbs[currentValue] <<- model$getLogProb(calcNodes)
		
        for(k in 1:M) {
            if(k != currentValue) {
				if(any(ID_j == k)){
					logProbs[k] <<- -Inf
				}else{	
					model[[target]] <<- k
					logProbPrior <- model$calculate(target)
					if(logProbPrior == -Inf) {
						logProbs[k] <<- -Inf
					} else {
						if(is.nan(logProbPrior)) {
							logProbs[k] <<- -Inf
						} else {
							logProbs[k] <<- logProbPrior + model$calculate(calcNodesNoSelf)
							if(is.nan(logProbs[k])) logProbs[k] <<- -Inf
						}
					}
				}
			}	
        }
        logProbs <<- logProbs - max(logProbs)
        probs <<- exp(logProbs)
        newValue <- rcat(1, probs)
        if(newValue != currentValue) {
            model[[target]] <<- newValue
            model$calculate(calcNodes)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        } else {
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        }
    },
    methods = list( reset = function() { } )
)

######################################################
## Nimble Implementation of the CRP
## This is Neal's Algorithm 8 for SCR.
#######################################################
sampler_myCRP <- nimbleFunction(
    name = 'sampler_myCRP',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
		
		calcNodesZ <- model$getDependencies('z')
		nodesAll <- c(calcNodesZ, calcNodes)
		nodesAll <- nodesAll[!duplicated(nodesAll)]

        M <- control$M
		m <- control$m
        probs <- numeric(M)
        logProbs <- numeric(M)		
		nodeIndex <- as.numeric(gsub('[^[:digit:]]', '', target))
    },
    run = function() {
		keep <- 0
		
		currentValue <- model[[target]]
		n <- sum(model[['ID']] == currentValue)
		if(n == 1){
			keep <- -1
			logProbs[currentValue] <<- model$getLogProb(calcNodes) + log(model[['alpha']]/m)
		}else{
			logProbs[currentValue] <<- model$getLogProb(calcNodes) + log(n-1)
		}
        for(k in 1:M) {
			if(k != currentValue)
			{
				if(model[['z']][k] != 0)
				{
					n <- sum(model[['ID']] == k)
					model[[target]] <<- k
					logProbs[k] <<- model$calculate(calcNodes) + log(n)
				}else{
					logProbs[k] <<- -Inf
				}
			}
		}
		
		a <- m + keep
		for(k in 1:a)
		{
			newVal <- rcat(1, 1-model[['z']])
			model[[target]] <<- newVal
			logProbs[newVal] <<- model$calculate(calcNodesNoSelf) + log(model[['alpha']]/m)			
		}
		
        logProbs <<- logProbs - max(logProbs)
        probs <<- exp(logProbs)
        newValue <- rcat(1, probs)
		model[[target]] <<- newValue
        if(newValue != currentValue) {
			model[['z']][newValue] <<- 1
			if(keep == -1) model[['z']][currentValue] <<- 0			
            model$calculate(nodesAll)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = nodesAll, logProb = TRUE)
        } else {
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        }
    },
    methods = list( reset = function() { } )
)

sampler_myAlpha <- nimbleFunction(
    name = 'sampler_myAlpha',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
		shape <- extractControlElement(control, 'shape', 1)		
		rate <- extractControlElement(control, 'rate', 1)		
		nObs <- length(model[['ID']])
    },
    run = function() {
		# This is described in the dirichletprocess package vignette and from West 1992
		# We assume a Gamma prior with parameters alpha_0 (user input).
		# Draw random gamma eta:
		alphai <- model[['alpha']]
		eta <- rbeta(1, alphai + 1, nObs)
		nGroups <- sum(model[['z']])
		
		# Calculate mixing proportions:
		pi1 <- shape + nGroups - 1
		pi2 <- nObs*(rate - log(eta))
		piSelect <- pi1/(pi1 + pi2)

		# Sample from posterior alpha:
		if(runif(1) < piSelect){
			model[[target]] <<- rgamma(1, shape + nGroups, rate - log(eta))
		}else{
			model[[target]] <<- rgamma(1, shape + nGroups - 1, rate - log(eta))
		}
		model$calculate(calcNodes)
		
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    },
    methods = list( reset = function() { } )
)

## The point of this sampler is to only sample animals that exist.
## Quick return if the animal doesn't exist for an efficient MCMC algorithm.
##--------------------------------------------------------------------------
sampler_myzs_CRP <- nimbleFunction(
    name = 'sampler_myzs_BD',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
		scale <- extractControlElement(control, 'scale', 1)		
        calcNodes <- model$getDependencies(target)
        nodeIndex <- as.numeric(gsub(".*?([0-9]+).*", '\\1', target))		
		xlim <- control$xlim
		ylim <- control$ylim
    },
    run = function() {
		if(model[['z']][nodeIndex] == 0) 
		{
				model[[target]][1] <<- runif(1, xlim[1], xlim[2])
				model[[target]][2] <<- runif(1, ylim[1], ylim[2])
				model$calculate(calcNodes)
				jump <- TRUE
		}else{
			newX <- model[[target]] + rnorm(2, 0, scale)
			lpcurrent <- model$getLogProb(calcNodes)
			model[[target]] <<- newX
			lpprop <- model$calculate(calcNodes)
			logMHR <- lpprop - lpcurrent
			jump <- decide(logMHR)
		}
		if(jump) {
			nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        } else {
			nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        }
    },
    methods = list( 
		reset = function() {} 
	)
)

## Fast computation of ESA for semi-complete likelihood and conditional!!
## Assuming Bernoulli capture history. In particular this is for acoustics stevenson 2015
RESA_C <- function(d2mask, sigma, g0, area, detfn, nmask) {

		ESA <- 0
		for(i in 1:nmask)
		{
			if(detfn == 1) # Half-normal
			{
				p.mask <- g0*exp(-d2mask[i,]/(2*sigma^2))
			}else{	# Hazard half-normal
				p.mask <- (1-exp(-g0*exp(-d2mask[i,]/(2*sigma^2))))
			}
			p.detect <- 1-prod(1-p.mask)
			ESA <- ESA + p.detect*area	
		}
   		return(ESA)
}

## Fast computation of ESA for semi-complete likelihood and conditional!!
## This one is for a repeated caller at a fixed location from stevenson 2021
RESA_A <- function(d2mask, lambda, sigma, g0, area, detfn, nmask) {

		ESA <- 0
		for(i in 1:nmask)
		{
			if(detfn == 1) # Half-normal
			{
				p.mask <- g0*exp(-d2mask[i,]/(2*sigma^2))
			}else{	# Hazard half-normal
				p.mask <- (1-exp(-g0*exp(-d2mask[i,]/(2*sigma^2))))
			}
			p.detect <- 1-prod(1-p.mask)
			
			ESA <- ESA + (1-exp(-p.detect*lambda))*area	# Assumes lambda is total expected calls.
		}
   		return(ESA)
}

ESA_C <- nimbleRcall(function(d2mask = double(2), sigma=double(0), 
								g0=double(0), area = double(0), detfn = double(0), 
								nmask = double(0)){}, 
		Rfun = 'RESA_C',
		returnType = double(0)
		)
ESA_A <- nimbleRcall(function(d2mask = double(2), lambda = double(0), sigma=double(0), 
								g0=double(0), area = double(0), detfn = double(0), 
								nmask = double(0)){}, 
		Rfun = 'RESA_A',
		returnType = double(0)
		)
















# These samplers can be ignored for now....
#########################################################


####################################
# Generalized for Multi Session:
####################################
sampler_myIDZ <- nimbleFunction(
    name = 'sampler_myIDZ',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
		# Find those that are dependent and not dependent for ID
        isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)
        calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
        calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]

        M <- control$M
		nv <- extractControlElement(control, 'Noccasion', 1)		
		v <- extractControlElement(control, 'occasion', 1)
		if(nv == 1){ 
			zNodes <- paste0('z[1:', M, ']')
			Hkv <- 'Hk'
		}else{
			zNodes <- paste0('z[1:', M,',', v,']')
			Hkv <- paste0('Hk[1:', M, ',', v,']')
		} 
		calcNodesZ <- model$getDependencies(zNodes)
		# Find those that are dependent and not dependent for z
        calcNodesZNoSelf <- model$getDependencies(zNodes, self = FALSE)		
        isStochCalcNodesZNoSelf <- model$isStoch(calcNodesZNoSelf)
        calcNodesZNoSelfDeterm <- calcNodesZNoSelf[!isStochCalcNodesZNoSelf]
        calcNodesZNoSelfStoch <- calcNodesZNoSelf[isStochCalcNodesZNoSelf]

        nodes1 = c(target, zNodes)
        nodes2 = c(calcNodesNoSelfDeterm, calcNodesZNoSelfDeterm)
        nodes3 = c(calcNodesNoSelfStoch, calcNodesZNoSelfStoch)
		nodesAll <- c(calcNodesZ, calcNodes)
		nodesAll <- nodesAll[!duplicated(nodesAll)]

        probs <- numeric(M)
        logProbs <- numeric(M)
	    nodeIndex <- as.numeric(gsub('[^[:digit:]]', '', target))
    },
    run = function() {
		psi <- model[['psi']]
        currentValue <- model[[target]]
		logProbs[currentValue] <<- model$getLogProb(calcNodes)
		n_currentValue <- sum(model[['ID']] == currentValue)
		if(n_currentValue == 1)
		{
			logProbs[currentValue] <<- logProbs[currentValue] + log(psi) - model[[Hkv]][currentValue]
			model[[zNodes]][currentValue] <<- 0
		}
        for(k in 1:M) {
			if(k != currentValue){
				# Start with check if it can match and assign 0 to 1 prob or 1 already prob.
				if(model[[zNodes]][k] == 0){
					logProbs[k] <<- log(psi)-model[[Hkv]][k]				
				}else{
					if(sum(model[['ID']] == k) == 0){
						logProbs[k] <<- -Inf
					}else{
						logProbs[k] <<- log(1-psi)
					}
				}
				# If it is not a zero-inflated value then find the full conditional.
				if(logProbs[k] != -Inf)
				{
					model[[target]] <<- k
					logProbs[k] <<- logProbs[k] + model$calculate(calcNodes)
				}
				
				# If it's a bad number make it -Inf.
				if(is.nan(logProbs[k])){
						logProbs[k] <<- -Inf
				}
			}
        }
        logProbs <<- logProbs - max(logProbs)
        probs <<- exp(logProbs)
        newValue <- rcat(1, probs)
        if(newValue != currentValue) {
			model[[target]] <<- newValue
            model$calculate(calcNodes)						
			# If I added a new zero animal, or killed an existing one, let's make sure it's updated for z.
			if(model[[zNodes]][newValue] == 0 | model[[zNodes]][currentValue] == 0){
				model[[zNodes]][newValue] <<- 1
				model$calculate(calcNodesZ)
				nimCopy(from = model, to = mvSaved, row = 1, nodes = nodesAll, logProb = TRUE)
			}else{
				model[[zNodes]][currentValue] <<- 1
				nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
				nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
				nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)				
			}
        } else {
				# Make sure we have the current value switched on!
				model[[zNodes]][currentValue] <<- 1
				model[[target]] <<- currentValue				
				nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
				nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
				nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
		}
	},
    methods = list( reset = function() { } )
)









# Joint ID Z sampler with constraints.
sampler_myIDZSPIM <- nimbleFunction(
    name = 'sampler_myIDZSPIM',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
		calcNodesZ <- model$getDependencies('z')
        targetNodesAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)		
        isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)
        calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
        calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
        M <- control$M
        probs <- numeric(M)
        logProbs <- numeric(M)		
		cannotlink <- extractControlElement(control, 'cannotlink', 'identity')	# n x n matrix where a 1 indicates i cannot link to detection j.
		AnimalLink <- extractControlElement(control, 'AnimalLink', -1)	# A zero represents an impossible animal match.
        if(AnimalLink == -1) AnimalLink <- numeric(length = M, value = 1)
		nodeIndex <- as.numeric(gsub('[^[:digit:]]', '', targetNodesAsScalar[1]))
		n_grp <- length(targetNodesAsScalar)
    },
    run = function() {
		psi <- model[['psi']]
        currentValue <- model[["ID"]][nodeIndex]
		logProbs[currentValue] <<- model$getLogProb(calcNodes)
		n_currentValue <- sum(model[['ID']] == currentValue)
		no_link <- 0
		if(n_currentValue == n_grp)
		{
			model[['z']][currentValue] <<- 0
			logProbs[currentValue] <<- logProbs[currentValue] + log(psi)-model[['Hk']][currentValue]
		}
        for(k in 1:M) {
		if(AnimalLink[k] == 0) {
			logProbs[k] <<- -Inf
		}else {	
				if(k != currentValue){
					if(model[['z']][k] == 1) {
						nk <- sum(model[['ID']] == k)
						if(nk == 0)
						{
							logProbs[k] <<- -Inf
						}else {
							no_link <- sum(cannotlink[model[['ID']] == k, nodeIndex])	# This is the check to see if there are cannot links.
							if(no_link > 0)
							{
								logProbs[k] <<- -Inf
							}else {
								values(model, targetNodesAsScalar) <<- rep(k, n_grp) 					
								logProbs[k] <<- model$calculate(calcNodes) + log(1-psi)
								if(is.nan(logProbs[k])) logProbs[k] <<- -Inf
							}
						}
					}else {
						values(model, targetNodesAsScalar) <<- rep(k, n_grp)
						logProbs[k] <<- model$calculate(calcNodes) + log(psi)-model[['Hk']][k]
						if(is.nan(logProbs[k])) logProbs[k] <<- -Inf
					}
				}
			}
		}
		# Note that logProbs of z=1 and nk=0 is -Inf, or it had better be!y
        logProbs <<- logProbs - max(logProbs)
        probs <<- exp(logProbs)
        newValue <- rcat(1, probs)
        if(newValue != currentValue) {
			values(model, targetNodesAsScalar) <<- rep(newValue, n_grp) ##  replace with this			
			if(model[['z']][newValue] == 0){
				model[['z']][newValue] <<- 1
				model$calculate(calcNodesZ)
			}
            model$calculate(calcNodes)	# I've made ID independent of z so this shouldn't double effort.
            nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        } else {
			model[['z']][currentValue] <<- 1
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        }
    },
    methods = list( reset = function() { } )
)

# If I put a uniform prior on sigma_toa^2 this is the conjugate prior sampler...
sampler_mySigmaToa <- nimbleFunction(
    name = 'sampler_mySigmaToa',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodesAll <- model$getDependencies(target)
		mi <- control$mi
		n <- length(mi)
		J <- control$J
    },
    run = function() {
		alpha <- 1
		ssq <- 0
		mdiff <- 0
		tdiff <- numeric(J)
		for(i in 1:n){
			if(mi[i] > 1){
				tdiff <- (model[['toa']][i,1:J] - model[['expTime']][model[['ID']][i],1:J])*model[['y']][i,1:J]
				mdiff <- sum(tdiff[1:J])/mi[i]
				ssq <- ssq + sum(model[['y']][i,1:J]*(tdiff[1:J] - mdiff)^2)
				alpha = alpha + (mi[i]-1)/2
			}
		}
	model[[target]]	<<- 1/sqrt(rgamma(1, shape = alpha, rate = ssq/2))
	model$calculate(calcNodesAll)
	nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesAll, logProb = TRUE)
    },
    methods = list( reset = function() {} )
)

# Sample sigma... test:
sampler_mySigma <- nimbleFunction(
    name = 'sampler_mySigma',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodesAll 	<- model$getDependencies(target)
        scale           <- extractControlElement(control, 'scale', 1)
    },
    run = function() {
	
		logProb <- -sum(model[['Hk']]*model[['z']]) + sum(log(model[['pobs']]))
		model[[target]]	<<- model[[target]] + rnorm(1, 0, sd = scale)
		model$calculate(calcNodesAll)
		logProbProposal <- -sum(model[['Hk']]*model[['z']]) + sum(log(model[['pobs']]))
		logMHR <- logProbProposal - logProb
		jump <- decide(logMHR)

		if(jump) {
				nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesAll, logProb = TRUE)
		} else {
				nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesAll, logProb = TRUE)
		}
    },
    methods = list( reset = function() {} )
)

sampler_myCat <- nimbleFunction(
    name = 'sampler_myCat',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
        isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)
        calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
        calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
        k <- control$M
        probs <- numeric(k)
        logProbs <- numeric(k)
	    nodeIndex <- as.numeric(gsub('[^[:digit:]]', '', target))
    },
    run = function() {
		currentValue <- model[[target]]
		trap <- model[['omega']][nodeIndex]
		probs <<- model[['pkj']][1:k,trap]*model[['z']]
		probs <<- probs/max(probs)
        newValue <- rcat(1, probs)
        if(newValue != currentValue) {
            model[[target]] <<- newValue
            model$calculate(calcNodes)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        } else {
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        }
    },
    methods = list( reset = function() { } )
)

