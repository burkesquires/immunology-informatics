qaProcess.GenericNumber <- function(numbers, frameIDs, outdir, name="generic", cutoff=-Inf, twoclass=FALSE)
{
# This function generates a qaProcess object using numbers, which can be used to generate an HTML report about a flowset
# Args:
#   numbers: a vector of numbers which will be displayed in a column of the report
#   frameIDs: a vector of strings which correspond to each of the frames (should match sampleNames(flowset))
#   outdir: the QA report directory
#   name: your choice of string describing the QA process
#   cutoff: a numeric value, such that the numbers below this value will be displayed in red, indicating a "fail" 
#           (default -Inf)
#   twoclass (default FALSE): if TRUE, kmeans will be applied to the numbers to separate them into two arbitrary groups
#       (e.g. positive/negative classes) and different colours will be used in the summary plot. The threshold value,
#       taken as the maximum value in the lower class, will also be plotted in the summary plot.
# Value:
#   returns a qaProcess object
#

  if(!is.null(dim(numbers)) || !is.null(dim(frameIDs)))
    stop("numbers and frameIDs must be vectors")

  if(length(numbers) != length(frameIDs))
    stop("Number of numeric values passed and number of frame IDs do not match")

  if(!file.exists(outdir))
    dir.create(outdir, recursive=TRUE)
  names(numbers) <- frameIDs
  
  # Create directory for images to be used in the HTML report, as well as a single summary image for this QA process:
  unique.id <- format.hexmode(as.integer(Sys.time())/runif(1)*proc.time()["elapsed"])
  image.dir <- paste(outdir, gsub(" ", "_", name), "/", sep="")
  dir.create(image.dir, recursive=TRUE, showWarnings=FALSE)
  summary.file <- file.path(image.dir, "summary.pdf")

  # Use k-means to define two arbitrary classes of numbers (e.g. if you have really low positive counts vs. really
  # high positive counts, you may be interested in the cut off value, and you may want to visualize this)
  # A bar plot with blue and purple bars will display the two classes.
  colours <- rep("blue", length(numbers))
  if(twoclass && length(unique(numbers)) > 1)
  {
    km <- kmeans(x = numbers, centers = 2, nstart = 50)
    indices <- km$cluster
    threshold <- min(max(numbers[which(indices == 1)]), max(numbers[which(indices == 2)]))
    smaller.cluster <- which.min(km$centers[,1])
    bigger.cluster <- which.max(km$centers[,1])
    colours[which(indices == smaller.cluster)] <- "blue"
    colours[which(indices == bigger.cluster)] <- "purple"
  }
  colours[which(numbers < cutoff)] <- "red"
  
  # Generate and save the summary image:
  pdf(file = summary.file, width = max(5,length(numbers)/10), height=4)
  par(cex = 1.2, mar = c(3, 3, 1, 1))
  b <- barplot(numbers, border = "gray", col = colours, names.arg = frameIDs, 
          cex.names = 0.5, ylab=name, space = 0, axisnames=FALSE, xaxt="n")
  text(b, 0, labels = frameIDs, srt=90,  xpd=TRUE, cex=.5, pos=1, offset = 1, col=colours)
  if(twoclass && length(unique(numbers)) > 1)
  {
    abline(h = threshold, lwd = 3, lty = "dashed", col = "gray")
    text(2, threshold, sprintf("%.2f", threshold), pos=3, offset=0.5)
    title(main = paste("Two-class separating threshold: ", sprintf("%.2f", threshold), sep = ""))
  } else {
    title(main = paste("Average: ", sprintf("%.2f", mean(numbers)), sep = ""))
    abline(h = mean(numbers), lwd = 3, lty = "dashed", col = "gray")
  }
  dev.off()
  summary.graph <- qaGraph(fileName = summary.file, imageDir = image.dir)
  
  # Formally define individual qaProcessFrame objects. Note the use of a numeric aggregator with a "passed" slot
  # defined by comparing a value to the cutoff value. If using the default cutoff = -Inf, all frames will "pass."
  frame.processes <- lapply(frameIDs, 
    function(i) 
    {
      qaProcessFrame(i, summaryAggregator = numericAggregator(numbers[i] > cutoff, x = numbers[i]))
    })
  names(frame.processes) <- frameIDs
  
  # Compile all the qaFrameProcess objects into one qaProcess object and return:
  qa.process <- qaProcess(id = unique.id, name = name, type = "number",
                          summaryGraph = summary.graph, frameProcesses=frame.processes)
  return(qa.process)
}


