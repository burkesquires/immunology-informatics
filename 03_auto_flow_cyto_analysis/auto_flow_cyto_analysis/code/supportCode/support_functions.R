####################################################################
### FLOW CYTOMETRY PRE-PROCESSING FUNCTIONS:
### Author: Radina Droumeva
### Date Last Modified: May 31, 2013
###
### The following functions can be called independently or through flowPrep(...).
####################################################################


####################################################################
# Remove scatter margin events
# Args:
#   frame: a flowFrame
# Value:
#   A flowFrame with its margin events removed by automatically detecting
#     scatter channels and removing cells in the high-range expression values.
cleanMargins <- function(frame){
  scatter.chans <- which(grepl("FS|SS", colnames(frame), ignore.case=TRUE))
  for (chan in scatter.chans){
    boundary.vals <- as.numeric(frame@parameters@data[chan, 
                                c("minRange", "maxRange")])
    boundary.vals[1] <- max(c(0, min(exprs(frame)[, chan])))
    largest.val <- quantile(exprs(frame)[, chan], 0.99)
    boundary.vals[2] <- min(c(largest.val, boundary.vals[2]))
    frame <- frame[which.between(exprs(frame)[, chan], boundary.vals)]
  }
  # Remove bad forced compensation margin events
  for (chan in setdiff(1:ncol(frame), scatter.chans)){
    min.val <- min(exprs(frame)[, chan]) + .Machine$double.eps * 2
    max.val <- max(exprs(frame)[, chan]) - .Machine$double.eps * 2
    frame <- frame[which(exprs(frame)[, chan] != min.val)]
    frame <- frame[which(exprs(frame)[, chan] != max.val)]
  }
  return (frame)
}


####################################################################
# Compensate, optionally specify spillover matrix
# Args:
#   frame: flowFrame
#   comp.mat: optionally a compensation matrix. If not specified,
#     the function tries to look for any keyword matching "SPILL" for it.
# Value:
#   A compensated flowFrame.
flowComp <- function(frame, comp.mat = NULL){
  if (is.null(comp.mat)){
    word <- which(grepl("SPILL", names(frame@description), ignore.case=TRUE))
    word <- word[which(unlist(lapply(frame@description[word], is.matrix)))][1]
    if (length(word) > 0){
      comp.mat <- frame@description[[word]]
    }
  }
  if (!is.null(comp.mat)){
    try(frame <- compensate(frame, comp.mat))
  } else {
    warning("No spill-over matrix found!")
  }
  return (frame)
}


####################################################################
# Remove doublets by gating high doublet channel cells
# Args:
#   frame: flowFrame
# Value:
#   A flowFrame object with doublets removed. A density-based gating is used
#     to remove high scatter-width cells, with 95th percentile used as quality check.
#     If no scatter width channel is found, an ellipsoid gate is fit over the
#     forward scatter height - area plane.
removeDoublets <- function(frame){
  channels <- colnames(frame)
  chan <- channels[which(grepl("FS.*W|SS.*W", channels, ignore.case = TRUE))][1]
  if (!is.na(chan)){
    doublet.gate <- deGate(frame, chan)
    doublet.gate <- max(doublet.gate, quantile(exprs(frame)[, chan], 
                                      0.95, na.rm=TRUE))
    frame <- frame[which(exprs(frame)[, chan] < doublet.gate)]
  } else {
    chans <- channels[which(grepl("FS.*H|FS.*A", channels, ignore.case = TRUE))]
    if (length(chans) == 2){
      frame <- getflowFrame(flowDensity(frame, channels = chans, scale = 0.9,
            position = c(TRUE, TRUE), gates = c(-Inf, -Inf), ellip.gate = TRUE))
    } else {
      warning("No scatter width channel found, only one of height and area available. Doublets not removed.")
    }
  }
  return (frame)
}


####################################################################
# Gate lymphocytes
# Args:
#   frame: flowFrame 
#   chans: 2-vector of forward and side scatter channels, such as c("FSC-A", "SSC-A")
#   ellipse.scale: 'scale' parameter in flowDensity, how tight the ellipsoid gate will be
#   fsc.gate: 2-vector used as a guideline for the lower and upper FSC gate
#   ssc.gate: simiarly, for the side scatter channel
gateLymphocytes <- function(frame, chans = NULL, ellipse.scale=0.975, 
                            fsc.gate = NULL, ssc.gate = NULL){
                            print(identifier(frame))
  channels <- colnames(frame)
  # Automatically detect one forward scatter and one side scatter channel if not provided.
  # Note that although flowPrep does this already, a user could call this function
  # outside of flowPrep and make use of the auto-channel detection.
  if (is.null(chans)){
    chans <- channels[which(grepl("FS", channels, ignore.case=TRUE))[1]]
    chans <- c(chans, channels[which(grepl("SS", channels, ignore.case = TRUE))[1]])
  } else if (is.numeric(chans)){
    chans <- channels[chans]
  }
  # Define default gates to be used as a quality check if none are supplied
  if (is.null(fsc.gate)){
    min.fsc <- min(exprs(frame)[, chans[1]])
    # If no cells in the lowest range (first 1/20th) of scatter values are present, 
    # it is likely the debris has already been removed manually.
    f.rng <- as.numeric(frame@parameters@data[which(colnames(frame) == chans[1]),
                        c('minRange', 'maxRange')])
    if (min.fsc > f.rng[1] + 1/20*(f.rng[2] - f.rng[1])){
      fsc.gate <- c(min.fsc, quantile(exprs(frame)[, chans[1]], 0.95, na.rm=TRUE))
    } else {
      # The default FSC gate should be the smaller of 15th quantile (i.e. 15%
      # debris), and 15% of the theoretical range. It is possible some data set
      # has much more than 15% debris, but this is staying on the conservative side.
      # Even if way more than 15% of debris exists, the gate shouldn't be much
      # further on the FSC value scale.      
      fsc.low.default <- min(quantile(exprs(frame)[, chans[1]], 0.15, na.rm=TRUE),
                             f.rng[1] + 0.15*(f.rng[2] - f.rng[1]))
      fsc.gate <- c(fsc.low.default, quantile(exprs(frame)[, chans[1]], 0.95,
                    na.rm=TRUE))
    }
  }
  # For side scatter, for quality checking simply use 1st and 95th percentiles.
  if (is.null(ssc.gate)){
    ssc.gate <- as.numeric(quantile(exprs(frame)[, chans[2]], c(0.01, 0.95), 
                           na.rm=TRUE))
  }
  # Calculate automatic gates for FSC and compare with default/passed gate.
  # If there is a difference of > 25% from the default/passed value (likely no debris)
  # try using standard deviation method instead and compare again.
  # If in the end the automated calculation still diverges from the default/passed
  # by > 25%, revert to default/passed.
  fsc.cuts <- deGate(frame, chans[1], all.cut=TRUE)
  fsc.low <- fsc.cuts[1]
  if (abs(fsc.low - fsc.gate[1])/fsc.gate[1] > 0.25){
    fsc.low <- deGate(frame, chans[1], sd.threshold=TRUE, n.sd=2, upper=FALSE)
    if (abs(fsc.low - fsc.gate[1])/fsc.gate[1] > 0.25){
      fsc.low <- fsc.gate[1]
    }
  }
  fsc.high <- tail(fsc.cuts, 1)
  if (abs(fsc.high - fsc.gate[2])/fsc.gate[2] > 0.25){
    fsc.high <- deGate(frame, chans[1], sd.threshold=TRUE, n.sd=2, upper=TRUE)
    if (abs(fsc.high - fsc.gate[2])/fsc.gate[2] > 0.25){
      fsc.high <- fsc.gate[2]
    }
  }
  fsc.gate <- c(fsc.low, fsc.high)
  
  # Since most debris is removed by FSC, simply use 1st percentile or passed
  # low-end SSC gate.
  ssc.low <- ssc.gate[1]
  first.gating <- flowDensity(obj=frame, channels=chans,
                 position=c(TRUE, TRUE), gates=c(fsc.gate[1], ssc.low))
  ssc.high <- deGate(getflowFrame(first.gating), chans[2], upper=TRUE)
  if (abs(ssc.high - ssc.gate[2])/ssc.gate[2] > 0.25){
    ssc.high <- ssc.gate[2]
  }
  ssc.gate <- c(ssc.low, ssc.high)
  a <- try({second.gating <- flowDensity(obj = first.gating, 
                                  channels = chans, position = c(FALSE, FALSE), 
                                  gates = c(fsc.gate[2], ssc.gate[2]), 
                                  ellip.gate=TRUE, scale = ellipse.scale)})
  # In rare cases the ellipsoid gate fails -- in such a case, simply turn it off.
  if (is( a, 'try-error')){
  second.gating <- flowDensity(obj = first.gating, 
                                  channels = chans, position=c(FALSE, FALSE), 
                                  gates = c(fsc.gate[2], ssc.gate[2]), 
                                  ellip.gate=FALSE)
  }
  return (getflowFrame(second.gating))
}


####################################################################
# Apply logicle transform to logarithmically scaled channels
# Args:
#   frame: single flowFrame
#   chans: character/numeric vector of channels to be transformed.
#     NOTE: if using flowPrep, chans can be automatically generated using
#           the P#DISPLAY keywords, but if calling outside of flowPrep,
#           chans must be specified. If it is NULL, no transformation is performed.
#   use.estimateL: if TRUE, the estimate logicle is used only if estimate.transform
#     is specified. Otherwise the generic Logicle is used.
#   estimate.transform: an object returned by 'estimateLogicle'.
#   M, A, bins, w: parameters for the generic Logicle transform.
transformChannels <- function(frame, chans=NULL, use.estimateL=FALSE, 
                              estimate.transform=NULL, 
                              M=4, A=1, bins=0, W=0.5){
                              print(identifier(frame))
  if (is.null(chans)){
    warning("No channels will be transformed.")
    return (frame)
  }
  if (is.numeric(chans)){
    chans <- colnames(frame)[chans]
  }
  if (use.estimateL && !is.null(estimate.transform)){
     frame <- transform(frame, estimate.transform)
  } else {
    for (chan in chans){
      # Separate Logicle for each channel since parameter T is channel-specific
      l <- Logicle::create(T=frame@parameters@data[which(colnames(frame)==chan), 
                           "maxRange"], M=M, A = A, bins=bins, W=W)
      exprs(frame)[, chan] <- Logicle::scale(l, exprs(frame)[, chan])
    }
  }
  return (frame)
}


####################################################################
### HELPER FUNCTIONS
####################################################################

####################################################################
# Finds the indices of the points in 1D data which lie between the two values in 2-vector v.
# Args:
#   data: numeric vector.
#   v: a numeric vector. Typically length 2, otherwise its smallest and largets
#     values are used.
# Value:
#   a vector of indices describing which values of 'data' fall in the v-range.
which.between <- function(data, v){
  data <- as.vector(data)
  if (!is.numeric(data) || !is.numeric(v)){
    return (NULL)
  }
  return (intersect(which(data >= min(v)), which(data < max(v))))
}


####################################################################
# Create a global frame which is a random sample of cells from a flowSet fs. 
# Args:
#   length: gives the number of samples used
#   r: the weight of each frame. Roughly, if this is 1, the expected number of
#     cells in the final sampled frame is about the same as that of a single frame
#     in the flow set. If r is 2, the final returned frame should have roughly
#     twice as many cells as any individual frame.
# Value:
#   a flowFrame is returned which contains randomly sampled cells from the flow set.
getGlobalFrame <- function(fs, length = 10, r = 2){
  if (is(fs, 'flowFrame')){
    return (frame)
  }
  n <- length(fs)
  sample.n <- min(length, n)
  global.frame <- fsApply(fs[sample(n, sample.n)], 
                          function(frame){
                            m <- nrow(frame)
                            frame <- frame[sample(m, min(m, m*r/sample.n))]
                            return (frame)
                          })
  global.frame <- as(global.frame, 'flowFrame')
  return (global.frame)
}


####################################################################
# Simple outlier check for gates. A use case is when a gate is calculated for
# the same stain for each of many samples, and a quality check is needed to 
# minimize the number of errors due to density peaks not being caught by flowDensity, or
# other issues. The gate values can be quality checked by setting any gate which is
# 'sd.coeff' standard deviations away from the average gate to the average.
#
# Assuming a normal distribution of the errors for the gate, an unbiased
# estimator of the standard deviation of a vector of values x is:
# s = sd(x)/c4(N),
# where sd(x) is the usual sqrt(1/(N-1)*sum(x_i - average(x), i = 1..N),
# N is the sample size,
# c4(N) is a correction term in terms of the Gamma function given by:
# c4(N) = sqrt(2/(N-1)) * Gamma(N/2)/Gamma((N-1)/2)
#
# (http://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation)
averageGates <- function(x, sd.coeff = 2){
# Args:
#   x: a numeric vector of values
# Value:
#   x with outlying values replaced by the median
  x <- as.vector(x)
  m <- median(x)
  N <- length(x)
  c4 <- sqrt(2/(N-1)) * gamma(N/2)/gamma((N-1)/2)
  sdev <- sd(x)/c4
  outliers <- which(abs(x - m)/sdev > sd.coeff)
  if (any(outliers)){
    cat("Changing outliers", outliers, ".\n")
    x[outliers] <- median(x[-outliers])
  }
  return (x)
}

# Convert colour name to hexadecimal, and optionally add transparency via the 
# alpha parameter. For example, 'col2hex('red', '55') will generate a transaprent red.
# Useful when placing a legend over a plot, can make the background colour of the
# legend col2hex("white", 70) for example.
col2hex <- function(colour, alpha = "FF"){
    colour <- col2rgb(colour)
    hex <- rgb(red = colour['red', ]/255, 
             green = colour['green', ]/255, 
              blue = colour['blue', ]/255)
    hex <- paste(hex, alpha, sep = "")
    return (hex)
}

# Remove bad characters from a string
validString <- function(str){
# This function takes in a string and returns an altered string that would be a valid file name
# Args:
#   str: a character vector that may contain bad characters
# Value:
#   a string with bad characters substituted with underscore.
  chars <- unlist(strsplit(str, ""))
  valid <- grep(pattern = "[a-zA-Z0-9_.]+", chars)
  bad <- setdiff(1:length(chars), valid)
  chars[bad] <- "_"
  return (paste(chars, collapse = ""))
}

