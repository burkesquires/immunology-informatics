####################################################################
### FLOW CYTOMETRY PRE-PROCESSING: flowPreP
### MARGIN EVENT REMOVAL, COMPENSATION, LOGICLE TRANSFORMATION, 
### DOUBLET REMOVAL, LYMPHOCYTE GATING
### Author: Radina Droumeva
### Date Last Modified: May 31, 2013
####################################################################


####################################################################
# Main preprocessing function
# Args:
#   fs: a 'flowSet' object
#   apply.comp: boolean, should the flow set be compensated. Default is FALSE.
#   comp.mat: optional spillover matrix. If unspecified and apply.comp = TRUE,
#     the embedded matrix is applied. If it does not exist, no compensation is
#     applied and a warning is printed to the screen.
#   lympho.chans: a two vector (character or numeric) identifying one forward 
#     and one side scatter channel to be used in the gating out of debris.
#     It is automatically detected if not specified.
#   ellipse.scale: to be passed onto flowDensity in the lymphocyte gate.
#   plot.for.lympho: boolean, if TRUE, a locator plot is displayed for a random
#     frame from the flow set for the scatter channels, and the user can use
#     six mouse clicks to outline the boundary/gate around the cells of interest-
#     these could be lymphocytes or some other population. The resulting gate is
#     then used as a quality check for the automatic lymphocyte gate on the rest
#     of the flow frames.
#   rem.doublets: boolean, should an attempt be made to remove doublets.
#     The approach looks for a scatter width channel and gates out the higher
#     population. If there is no scatter width channel, then an ellipsoid gate
#     is fit over the FSC-H -- FSC-A plane.
#   trans.chans: a vector of channels which require transformation. If not 
#     specified, the P1DISPLAY, P2DISPLAY, etc. keywords are used to detect
#     which channels are on a logarithmic scale and should be transformed.
#   use.estimateL: boolean, if TRUE, the estimate logicle transform is estimated
#     using a global sampling of all cells in the flow set. Otherwise, the
#     standard Logicle transform is used with parameters M, A, bins, W below.
#   M, A, bins, W: parameters for the Logicle transform. ?Logicle for info.

# Testing parameters copy paste convinience:
if (FALSE){
  apply.comp = FALSE
  comp.mat = NULL
  lympho.chans = NULL
  ellipse.scale = 0.975
  plot.for.lympho = TRUE
  plot.preproc = FALSE
  rem.doublets = TRUE
  trans.chans = NULL
  use.estimateL = TRUE
  M = 4
  A = 1
  bins = 0
  W = 0.5
}
flowPrep <- function(fs, apply.comp = FALSE, comp.mat = NULL, lympho.chans = NULL,
           ellipse.scale = 0.975, plot.for.lympho = TRUE, plot.preproc = FALSE, 
           rem.doublets = TRUE, trans.chans = NULL, 
           use.estimateL = FALSE, M = 4, A = 1, bins = 0, W = 0.5){
  # TO DO: input validity checks
  
  # This check enables this for use on a single frame
  is.flowset <- is(fs, 'flowSet')
  if (!is.flowset){
    fs <- as(fs, 'flowSet')
  }
  
  # Remove margin events
  fs <- fsApply(fs, cleanMargins)
  
  # Apply compensation if requested
  if (apply.comp){
    fs <- fsApply(fs, flowComp, comp.mat)
  }
  
  channels <- colnames(fs[[1]])
  # Locate channels which require a transformation
  if (is.null(trans.chans)){
    display.keywords <- sapply(1:ncol(fs[[1]]), function(x) 
                      fs[[1]]@description[[paste("P", x, "DISPLAY", sep = "")]])
    if (is.null(unlist(display.keywords))){
      warning("No keyword found to detect transformation channels!")
    } else {
      locate.null.entries <- which(unlist(lapply(display.keywords, is.null)))
      display.keywords[locate.null.entries] <- "NULL"
      trans.chans <- which(grepl("LOG|LG", unlist(display.keywords), 
                           ignore.case = TRUE))
    }
  } 
  if (is.numeric(trans.chans)){
      trans.chans <- channels[trans.chans]
  }
  # el stands for Estimate Logicle
  el <- NULL
  if (!is.null(trans.chans)){
    if (use.estimateL){
    # Define the transform using a global subset of cell from the whole flow set
      global.frame <- getGlobalFrame(fs)
      el <- estimateLogicle(global.frame, trans.chans)
    } else {
      el <- NULL
    }
    print("Transforming channels...")
    fs <- fsApply(fs, transformChannels, chans = trans.chans,
                                use.estimateL = use.estimateL, 
                                estimate.transform = el, 
                                M = M, A = A, bins = bins, W = W)
  }
  # Keep track of flow set before lymphocyte gating for plotting QA purposes:
  fs.start <- fs
  # Attempt to remove doublets
  if (rem.doublets){
    fs <- fsApply(fs, removeDoublets)
  }

  # Automatically detect one forward scatter and one side scatter channel if not provided
  if (is.null(lympho.chans)){
    lympho.chans <- channels[which(grepl("FS", channels, ignore.case = TRUE))[1]]
    lympho.chans <- c(lympho.chans, channels[which(grepl("SS", channels, 
                      ignore.case = TRUE))[1]])
  } else if (is.numeric(lympho.chans)){
    lympho.chans <- channels[lympho.chans]
  }
  # Initialize 'fsc.gate' and 'ssc.gate' to NULL
  fsc.gate <- ssc.gate <- NULL
  # Try to plot one frame at random and get 6 manual points as a scatter plot gate
  if (plot.for.lympho){
    try({
      dev.new()
      # Plot a global subset of the flow set for lymphocyte gating on screen
      plot.frame <- getGlobalFrame(fs)
      plotDens(plot.frame, lympho.chans, main = "Select cells of interest in 6 clicks\n(ensure lowest and highest FSC and SSC points are included)")
      coords <- locator(n = 6, type = "l")
      # Define a reference scatter gate for the cells of interest using user's clicks
      fsc.gate <- c(min(coords$x), max(coords$x))
      ssc.gate <- c(min(coords$y), max(coords$y))
    })
  }
  # Gate out debris and fit ellipsoid gate on the two scatter channels
  print("Gating lymphocytes...")
  fs <- fsApply(fs, gateLymphocytes, chans = lympho.chans, 
                ellipse.scale = ellipse.scale, 
                fsc.gate = fsc.gate, ssc.gate = ssc.gate)
  # l and w define the plotting region length and width
  l <- ceiling(sqrt(length(fs)))
  w <- ceiling(length(fs)/l)
  # Plot gated lymphocytes [during testing stage]
  if (plot.preproc == TRUE){
    par(mfrow = c(l, w), mar = c(3,3,3,2), mgp = c(2, 1, 0))
    for (i in 1:length(fs)){
      plotDens(fs.start[[i]], lympho.chans, main = identifier(fs[[i]]), 
               devn = FALSE)
      ch <- chull(exprs(fs[[i]])[, lympho.chans[1]],
                  exprs(fs[[i]])[, lympho.chans[2]])
      ch <- c(ch, ch[1])
      lines(exprs(fs[[i]])[ch, lympho.chans], lwd = 2, lty = "dashed")
    }
  }
  if (!is.flowset){
    fs <- as(fs, 'flowFrame')
  }
  print("Done.")
  return (fs)
}

