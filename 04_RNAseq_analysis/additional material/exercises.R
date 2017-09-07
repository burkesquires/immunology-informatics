#highlight, command-enter

2+2
1+1
2 + 4 * 5      
# Order of operations
log (10)       
# Natural logarithm with base e=2.7182
log10(5)      
# Common logarithm with base 10
5^2             
# 5 raised to the second power
5/8             
# Division
sqrt (16)      
# Square root
abs (3-7)     
# Absolute value 
pi                
# 3.14
round(pi,0)  
# Round pi to a whole number
round(pi,1)  
# Round pi to 1 decimal place
round(pi,4)  
# Round pi to 4 decimal places
floor(15.9)   
# Rounds down
ceiling(15.1)  
# Rounds up

x = c(1,2,3)
x
median(x)

x <- c(1:3, 7:8)
#[1] 1 2 3 7 8

x <- c(1:100)
x
x <- runif(100,0,1)
hist(x)
histinfo<-hist(x)
histinfo
str(histinfo)
histinfo$counts
 [1]  8 14 11 14  7  9  7 12 12  6
histinfo$counts[1]
[1] 8

#create some normally distributed data yourself. 
#In R, you can generate normal data this way using the rnorm() function:
BMI<-rnorm(n=1000, m=24.2, sd=2.2)
hist(BMI)
histinfo<-hist(BMI)
histinfo
str(histinfo)
summary(histinfo)
#function that combines its arguments to form a vector

A = c(2, 4, 3, 1, 5, 7)
sort(A)

#A = matrix( 
#+   c(2, 4, 3, 1, 5, 7), # the data elements 
#+   nrow=2,              # number of rows 
#+   ncol=3,              # number of columns 
#+   byrow = TRUE)        # fill matrix by rows 

A = matrix( c(2, 4, 3, 1, 5, 7), nrow=2, ncol=3, byrow = TRUE)
heatmap(A)

# see links: http://www.inside-r.org/packages/cran/gplots/docs/heatmap.2
#https://stat.ethz.ch/pipermail/bioconductor/2008-July/023300.html
#http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/r/heatmap/
#http://www.molecularecologist.com/2013/08/making-heatmaps-with-r-for-microbiome-analysis/

#data_matrix <- data.matrix(data)
#data_matrix <- data.matrix(A)
install.packages("RColorBrewer")
library("RColorBrewer")
#display all colour schemes
display.brewer.all()
heatmap(data_matrix,col=brewer.pal(9,"Blues"))

A[2, 3]      
# element at 2nd row, 3rd column 
[1] 7 

A[2, ]       
# the 2nd row 
[1] 1 5 7 

A[ ,3]       
# the 3rd column 
[1] 3 7 

A[ ,c(1,3)]  
# the 1st and 3rd columns 
     [,1] [,2] 
[1,]    2    3 
[2,]    1    7 

dat <- data.frame(xvar = 1:20 + rnorm(20,sd=3),
                  yvar = 1:20 + rnorm(20,sd=3),
                  zvar = 1:20 + rnorm(20,sd=3))
head(dat)

# Plot the points using the vectors xvar and yvar
plot(dat$xvar, dat$yvar)

#abline: This function adds one or more straight lines through the current plot
## Setup up coordinate system (with x == y aspect ratio):
plot(c(-2,3), c(-1,5), type = "n", xlab = "x", ylab = "y", asp = 1)
## the x- and y-axis, and an integer grid
abline(h = 0, v = 0, col = "gray60")
text(1,0, "abline( h = 0 )", col = "gray60", adj = c(0, -.1))
abline(h = -1:5, v = -2:3, col = "lightgray", lty = 3)
abline(a = 1, b = 2, col = 2)
text(1,3, "abline( 1, 2 )", col = 2, adj = c(-.1, -.1))

## Simple Regression Lines:
require(stats)
sale5 <- c(6, 4, 9, 7, 6, 12, 8, 10, 9, 13)
plot(sale5)
abline(lsfit(1:10, sale5))
abline(lsfit(1:10, sale5, intercept = FALSE), col = 4) 
# less fitting

#data frames 
myFirstDataframe = data.frame(       
# Press Enter to start a new line.
name=c("Bob","Fred","Barb","Sue","Jeff"),
age=c(21,18,18,24,20), hgt=c(70,67,64,66,72),
wgt=c(180,156,128,118,202),
race=c("Cauc","Af.Am","Af.Am","Cauc","Asian"),
year=c("Jr","Fr","Fr","Sr","So"),
SAT=c(1080,1210,840,1340,880))    
# End with double close parenthesis. Why?

myFirstDataframe
  name age hgt wgt  race year  SAT
1  Bob  21  70 180  Cauc   Jr 1080
2 Fred  18  67 156 Af.Am   Fr 1210
3 Barb  18  64 128 Af.Am   Fr  840
4  Sue  24  66 118  Cauc   Sr 1340
5 Jeff  20  72 202 Asian   So  880
str(myFirstDataframe)
summary(myFirstDataframe)

#functions
hello<-function() cat("Hello, class!\n")
hello()

hello<-function() {
	cat(paste("Hello, ",system("whoami",T),"!\n",sep="",collapse=""))
}
hello()

hello<-function() {
	cat(paste("Hello, ",system("whoami",T),"!\n",sep="",collapse=""))
	today<-as.list(unlist(strsplit(system("date",T)," ")))
	names(today)<-c("weekday","month","day","time","timezone","year")
	return(today)
}
hello()

