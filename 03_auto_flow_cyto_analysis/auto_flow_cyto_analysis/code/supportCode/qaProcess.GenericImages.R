qaProcess.GenericImages <- function(
  filenames,
  frameIDs,
  outdir="QAReport",
  name="Gating Images",
  flags = NULL,
  width=200,
  pdf = FALSE,
  summary.graph = NULL,
  ...)
{
# This function generates a qaProcess object using images, which can be used to generate an HTML report about a flowset
# Args:
#   filenames: a vector of paths to images which will be displayed in a column of the report
#   frameIDs: a vector of strings which correspond to each of the frames (should match sampleNames(flowset))
#   outdir: the QA report directory
#   name: your choice of string describing the QA process
#   flags: a vector of frameIDs corresponding to frames that should "fail" this QA procedure
#   width (default is 200): width for the HTML images. The html image extension defaults to the extension of the first 
#       image in filenames, or to "jpeg" if the extension is not one of png, jpeg, pdf or eps
#   summary.graph: an optional file path for this QA process's summary image. If NULL, the first image in filenames is
#       used.
# Value:
#   returns a qaProcess object
#

  if(!is.null(dim(filenames)) || !is.null(dim(frameIDs)))
    stop("filenames and frameIDs must be vectors")

  if(length(filenames) != length(frameIDs))
    stop("Number of images passed and number of frame IDs do not match")

  # Assess image file extension
  file.extension <- strsplit(filenames[1], "\\.")[[1]]
  file.extension <- file.extension[length(file.extension)]
  if(!is.element(file.extension, c("png", "jpeg", "pdf", "eps")))
    file.extension <- "jpeg"
  
  # Generate new directory for the images to be associated with the HTML report later
  unique.id <- format.hexmode(as.integer(Sys.time())/runif(1)*proc.time()["elapsed"])
  frame.names <- frameIDs
  num.frames <- length(frame.names)
  
  # Create a summary graph by copying the first image in the filenames vector
  summary.graph <- ifelse(is.null(summary.graph), filenames[1], summary.graph)

  super.dir <- unlist(strsplit(summary.graph, "/"))
  super.dir <- paste(paste(super.dir[-length(super.dir)], collapse="/"), "/", sep="")
  if(file.extension == "pdf")
  {
    system(paste('rm ', super.dir, "*.jpg", sep=""))
  }
  summary.qa.graph <- qaGraph(fileName = summary.graph, 
                              imageDir = super.dir, 
	                            width    = width, pdf=pdf)
	      
  # Generate individual qaProcessFrame objects:
  frameProcesses <- list()
  counter <- 0
  bar <- txtProgressBar(min = 0, max = length(file.names), style = 3)
  super.dir <- unlist(strsplit(filenames[1], "/"))
  super.dir <- paste(paste(super.dir[-length(super.dir)], collapse="/"), "/", sep="")
  for(i in seq_len(num.frames))
  {
    file.names <- NULL
    aggregator.list <- aggregatorList()
    frame.id <- frame.names[i]

    # In order for the HTML report to be built properly, each image has to be copied into a fresh directory that was
    # created above. Images are copied one by one and their new destination is used instead:
    pass <- ifelse(is.element(frame.id, flags), FALSE, TRUE) # check if the frame id was in the "flags" vector
    aggregator.list[[1]] <- binaryAggregator(pass) # create a binary aggregator (pass or fail) associated with the image
    summary.aggregator <- binaryAggregator(pass) # create a summary aggregator associated with the frame overall

    # Create the necessary flowQ object to store the image information:
    fGraphs <- qaGraphList(imageFiles = filenames[i], imageDir = super.dir, width = width, pdf = pdf)
    # Create individual qaProcessFrame objects:
    frameProcesses[[frame.id]] <- qaProcessFrame(frameID = frame.id,
                                       summaryAggregator = summary.aggregator,
                                       frameAggregators  = aggregator.list,
                                       frameGraphs       = fGraphs)
    counter <- counter + 1
    setTxtProgressBar(bar, counter)
  }

  # Finally compile everything in a single qaProcess object and return:
  output <- qaProcess(id             = unique.id, 
                      name           = name,
                      type           = "images",
                      summaryGraph   = summary.qa.graph,
                      frameProcesses = frameProcesses)
  return(output)
}

