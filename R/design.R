#' Conversion to numeric
#' 
#' Simply an abbreviation for \code{as.numeric(as.character())}
#' 
#' 
#' @param x Something to be converted to numeric
#' @author Zach Kurtz
#' @export as.num
as.num = function(x) as.numeric(as.character(x))



#' Put LLM vector into a LLM design matrix
#' 
#' Given the output of \code{pop.to.counts}, i.e., a 1-row dataframe with as
#' many columns as there are observable capture patterns, create a design
#' matrix for log-linear modeling
#' 
#' 
#' @param densi The capture pattern counts
#' @param rasch Logical: if TRUE, the square of the number of captures is
#' included in the design matrix as the column \code{ck2}
#' @return A design matrix
#' @author Zach Kurtz
#' @export string.to.array
string.to.array = function(densi, rasch = FALSE){ #densi = ydens[1,]
	## Put the data into cannonical poisson regression form
	# Get the design matrix:
	k = nchar(colnames(densi)[1])
	dat = make.design.matrix(k=k, rasch = rasch)
	# Populate the design matrix with the counts of each capture pattern
	cpatterns = as.character(apply(dat[,1:k], 1, paste, collapse = ""))
	dat$c = rep(0, nrow(dat)) #initialize counts to zero
	for(pat in colnames(densi)) dat$c[which(cpatterns == pat)] = densi[,pat] 
	return(dat)
}



#' Construct standard LLM design matrix.
#' 
#' Makes a design matrix for a (local) log-linear model
#' 
#' 
#' @param k The number of lists
#' @param order.max The maximum number of lists to include in interaction terms
#' @param rasch Logical: if TRUE, include a column for the square of the number
#' of captures
#' @return A design matrix as a data frame
#' @author Zach Kurtz
#' @export make.design.matrix
#' @import combinat
make.design.matrix = function(k = 3, order.max = k-1, rasch = FALSE){
	# k is the number of lists; include all list interactions up to order order.max
	# make the first-order design matrix
	dat = expand.grid( rep( list(c(0,1) ), k) )[-1,]
	names(dat) = apply(cbind(rep("c", k), 1:k), 1, paste, collapse = "")
	if(order.max > 1){ # Put in the interactions up to order.max
		for(ord in 2:order.max){ #ord = 2
			com.m = t(combn(k, m = ord))
			com.f = apply(com.m, 1, paste, collapse = "")
			com.f = paste("c", com.f, sep = "")
			for(t in 1:nrow(com.m)){ #t = 1
				dat[,com.f[t]] = apply(dat[, com.m[t,]], 1, prod)
			}
		} # end for(ord
	} # end if(order.max > 1)
	if(rasch) dat$ck2 = rowSums(dat)^2 
	return(dat)
}



#' Capture patterns to design matrix
#' 
#' Tally the multinomial counts for a vector of capture patterns, and put the
#' result into design matrix
#' 
#' 
#' @param y A vector of capture patterns (each capture pattern is a single
#' string of binary digits)
#' @return A design matrix, including capture pattern counts, that is ready for
#' log-linear analysis
#' @author Zach Kurtz
#' @export y.string.to.y.glm
y.string.to.y.glm = function(y){
	# Given the set of capture patterns (strings) return a glm-ready data.frame
	#   for traditional log-linear modeling (no covariates here)
	splity = function(k) strsplit(y[k], "")[[1]]
	dat = data.frame(t(sapply(1:length(y), splity)))
	names(dat) = paste("c", seq(1:nchar(y[1])), sep = "")
	dat$count = rep(1, nrow(dat))
	dat = xtabs(count ~ c1 + c2 + c3, data = dat)
	dat = as.data.frame(dat)
	dat = apply(dat, 2, as.num)
	dat = as.data.frame(dat)
	dat$c = dat$Freq; dat$Freq = NULL
	dat = dat[-1,]
	return(dat)
}



#' Generate all observable capture patterns
#' 
#' Returns a vector of binary strings representing all possible capture
#' patterns excluding complete noncapture
#' 
#' 
#' @param k The number of lists
#' @author Zach Kurtz
#' @export patterns.possible
patterns.possible = function(k) as.character(apply(expand.grid( rep( list(c(0,1) ), k) )[-1,], 1, paste, collapse = ""))



#' Put CRC data into LLM vector
#' 
#' Essentially, this does \code{table(y)}, where \code{y} is a vector of
#' capture patterns, with the exception that here every capture pattern that is
#' observable (in principle) is included, with a possibly zero count.
#' 
#' 
#' @param y A vector of capture patterns
#' @author Zach Kurtz
#' @export pop.to.counts
pop.to.counts = function(y){
	# Given the full vector of individually observed capture patterns, return 1-row df 
	#   of capture-pattern counts in standard form
	cts = table(y)
	k = nchar(y[1])
	pats = patterns.possible(k)
	out = data.frame(matrix(0, ncol = length(pats), nrow =1))
	names(out) = pats
	out[names(cts)] = cts
	return(out)
}



#' Template for capture-pattern counts
#' 
#' Make a matrix with one row and as many columns as there are observeable
#' capture patterns, with column names corresponding to the possible capture
#' patterns.
#' 
#' By default, the matrix gets populated with uniform value summing to 1 across
#' all the entries.
#' 
#' @param k The number of lists
#' @author Zach Kurtz
#' @export make.patterns.template
make.patterns.template = function(k=3){
	# Set up a row vector in the form of multinomial capture-pattern probabilities
	x = patterns.possible(k)
	n.cases = length(x)
	y = data.frame(matrix(rep(1/n.cases, n.cases), nrow = 1))
	names(y) = x
	return(y)
}



#' Generate a universe of hierarchical model specifications.
#' 
#' Creates a list of all hierarchical sets of log-linear modelling terms that
#' include all main effects and exclude the highest order interaction.
#' Currently implemented only for k=3 lists.
#' 
#' Needs generalization to k>3 very much!
#' 
#' @param k The number of lists
#' @param rasch If TRUE, include the Rasch model (simplest Rasch model in
#' Darroch et. al. 1993).
#' @return A list of character vectors
#' @author Zach Kurtz
#' @export make.hierarchical.term.sets
make.hierarchical.term.sets = function(k=3, rasch = FALSE){
	if(k!=3) stop("In make.hierarchical.term.sets, we have not yet generalized any case other than k=3 lists")
	term.sets = list()
	term.sets[[1]] = c("c1", "c2", "c3")
	term.sets[[2]] = c("c1", "c2", "c3", "c12")
	term.sets[[3]] = c("c1", "c2", "c3", "c13")
	term.sets[[4]] = c("c1", "c2", "c3", "c23")
	term.sets[[5]] = c("c1", "c2", "c3", "c12", "c13")
	term.sets[[6]] = c("c1", "c2", "c3", "c12", "c23")
	term.sets[[7]] = c("c1", "c2", "c3", "c13", "c23")
	term.sets[[8]] = c("c1", "c2", "c3", "c12", "c13", "c23")
	if(rasch) term.sets[[9]] = c("c1", "c2", "c3", "ck2")
	return(term.sets)
}


# Put data into a somewhat archaic standard form:


#' Format the CRC data
#' 
#' Put a CRC dataset into a standard form that is understood by functions like
#' \code{lllcrc}.
#' 
#' 
#' @param x The data frame of CRC data, such as the output of the simulation
#' function \code{poptop}.
#' @param y A character string indicating which column of \code{x} contains the
#' capture patterns.
#' @param cont.vars A character vector of variable names for continuous
#' variables.
#' @param categ.vars A character vector of variable names for categorical
#' variables.
#' @return A data frame that is an expanded version of the argument \code{x}.
#' Variables are renamed using an "x.con...." or "x.dis...." for continuous and
#' discrete variables, respectively.  Each discrete/categorical variable with L
#' levels is expanded into L separate binary variables.
#' @author Zach Kurtz
#' @export formatdata
formatdata = function(x, y, cont.vars = NULL, categ.vars = NULL){ 
	if(is.null(categ.vars)&is.null(cont.vars)){
		stop("You need to specify at least one of categ.vars or cont.vars")
	}
	# Set up a list of variable keys
	var.keys = list()

	# Standard the names of continuous covariates to be in the form x.con.i for the ith covariate
	if(length(cont.vars)>0){
		for(i in 1:length(cont.vars)){
			xvar = paste("x.con.", i, sep = "")
			var.keys[xvar] = cont.vars[i]
			temp = x[,cont.vars[i]]
			x[,cont.vars[i]] = NULL
			x[,xvar] = temp
		}
	}
	# Standardize the names of categorical (discrete) covariates to be in the form x.con.i.j for level j of covariate i
	if(length(categ.vars)>0){
		for(i in 1:length(categ.vars)){
			temp = x[, categ.vars[i]]
			xuniq = unique(temp)
			x[,categ.vars[i]] = NULL
			var.keys[paste("x.dis.", i, sep = "")] = categ.vars[i]
			for(j in 1:length(xuniq)){
				xvar = paste("x.dis.", i, ".", xuniq[j], sep = "")
				x[, xvar] = 0+(temp == xuniq[j])
			}				
		}
	}
	# Split y into a multinomial outcome:
	ty = x[,y]; x[,y] = NULL; y = ty; ty = NULL
	k = nchar(y[1])
	ynames = apply(expand.grid( rep( list(c(0,1) ), k) )[-1,], 1, paste, collapse = "")
	for(yn in ynames) x[, paste("y", yn, sep = "")] = 0+(y == yn)
	# Affix variable keys:
	attributes(x)$var.keys = var.keys
	return(x)
}



#' Collapse CRC data through micro post-stratification
#' 
#' Rounding continuous covariates creates "micro-post-strata" and therefore
#' tends to reduce the number of distinct covariate vectors.  After rounding,
#' the data is collapsed so that there is exactly one row for each distinct
#' covariate vector, and a column called \code{mct} (for multinomial cell
#' count) is appended with that contains the number of records corresponding to
#' each row.
#' 
#' Continuous variables will be divided by \code{rounding.scale}, then rounded
#' to the nearest whole number, and then multiplied by \code{rounding.scale}.
#' The net effect is to round to the nearest multiple of \code{rounding.scale}
#' 
#' @param dat The data in a matrix form
#' @param round.vars A vector of names of variables to be rounded for the
#' purpose of collapsing the data.
#' @param rounding.scale A vector of scalars that determines how much each
#' corresponding variable in \code{round.vars} is to be rounded.  For example,
#' the first variable \code{round.vars[1]} will be divided by
#' \code{rounding.scale[1]}, then rounded to the nearest whole number, and then
#' multiplied by \code{rounding.scale[1]}.  The net effect is to round to the
#' nearest multiple of \code{rounding.scale[1]}.
#' @return Another matrix, just like the input \code{dat} except that there are
#' fewer rows there is the new column \code{mct}
#' @author Zach Kurtz
#' @export micro.post.stratify
micro.post.stratify = function(dat, round.vars = NULL, rounding.scale = NULL){
	if(!is.null(round.vars)){
		for(i in 1:length(round.vars)){
			temp = dat[,round.vars[i]]
			temp = round(temp/rounding.scale[i])
			dat[,round.vars[i]] = temp*rounding.scale[i]
		}
	}
	# Collapse the data over the set of unique covariate vectors
	x.covs = names(dat)[substr(names(dat), 1,1) == "x"]
	y.bits = names(dat)[substr(names(dat), 1,1) == "y"]
	x.part = paste(x.covs, collapse = "+")
	y.part = paste(y.bits, collapse = ",")
	agg.form = paste("out = aggregate(cbind(", y.part, ") ~ ", x.part, ", data = dat, FUN = 'sum')", sep = "")
	out = eval(parse(text = agg.form))
	out$mct = rowSums(out[,y.bits])
	#out = data.matrix(out)
	return(out)
}
