#
# Example of ploting R graphs with error bars
# Modelled precisely on Justin Zobel's example in jgraph
#
# Authour: andrew.turpin@rmit.edu.au
# Date:    Fri Apr 20 16:38:28 EST 2007
#

    # uncomment this for pdf output to file
#pdf("ebars.pdf", width=7, height=5, family="Times", paper="a4")
#postscript("ebars.ps", width=7, height=5, family="Times", paper="a4")

    # A little function I wrote to draw a point at (x,y)
    # and a vertical error bar from (x, y-len) to (x, y+len)
    # assumes plotting axis already set up
    # Note the use of ... to pass on unspecified arguments
plotPandB <- function(x, y, len, sfrac=0.01, ...) {
    points(x, y, ...)
    arrows(x, y, x, y+len, length=par("fin")[1] * sfrac, angle=90)
    arrows(x, y, x, pmax(0, y-len), length=par("fin")[1] * sfrac, angle=90)
}

    # Read in the 4 data sets into four list elements
    # Only read the first 12 rows as in jz's original eg
d <- list(1:4)
for (i in 1:4)
    d[[i]] <- read.table(
                paste("ebars-d",i,".data",sep=""),
                nrows = 12              
              )

    #set up the axis
plot(0,0,
    xlim=c(0, 12), 
    ylim=c(0, 0.045), 
    type="n",
    xlab = "Level",
    ylab = "Av. probability at level"
)

    # draw points (x offsets as per jz's original jgraph)
    #         square     x     circle   triangle
markType <- c(22,       120,    19,     24)
backGs   <- c("black", "white", "gray", "white")
cols     <- c("black", "gray" , "gray", "gray")
offsets  <- c(-0.24,   -0.08,   +0.08,  +0.24)

for(i in 1:4) 
    plotPandB(1:length(d[[i]][,1]) + offsets[i], d[[i]][,1], d[[i]][,2], 
                pch = markType[i], 
                bg  = backGs[i], 
                col = cols[i],
                sfrac = 0.005
    )

    # plonk on a legend
legend(8.5, 0.025, 
    c("Major", "Upper", "Lower", "Minor"),
    pch=markType,
    col = cols,
    pt.bg = backGs
)

#dev.off()  # This is only needed if you use pdf/postscript in interactive mode
