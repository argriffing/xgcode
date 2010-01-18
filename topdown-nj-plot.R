# Do it like this:
# R CMD BATCH thisfilename.R
#
# Here are some obsolete instructions from before I started using 'mar'.
# This old way is cleaner in the sense that I don't have to
# manually tweak the margin numbers, but it is worse in the sense
# that I have to post-process the resulting tikz file:
# "Then afterwards, modify the resulting .tikz file as follows.
# First, remove the white box drawn around the border.
# Second, remove the multiple lines of code that define clipping boxes.
# These modifications allow LaTeX to infer the correct margins.
# Without these modifications, there is too much whitespace."
#
# In the LaTeX caption that goes with this plot,
# it should be mentioned that 
# nj is neighbor joining and
# rs is recursive splitting

# read the table
mytable <- read.table('nj-like-1000.table')

# load the tikz library
require('tikzDevice')

# Begin writing a tikz figure.
# For the old method where I used the default margins
# and manually hacked the .tikz output to remove the bounding box,
# the best settings were width=3.5 and height=4.
# Note that only sizes 10pt, 11pt, and 12pt are valid.

# Try one way of doing this.
options(tikzDocumentDeclaration = '\\documentclass[12pt]{article}')
tikz('nj-like-1000.tikz', width=3.25, height=4)

# This alternative way did not work.
#tikz('nj-like-1000.tikz', documentDeclaration='\\documentclass[12pt]{article}', width=3.25, height=4)

# Define the margins.
# The syntax is mar=c(bottom, left, top, right).
# The default values are c(5,4,4,4) + 0.1.
# I need some extra space at the bottom and the left
# to show the ticks, the tick labels, and the axis labels.
# Also I need some space to the right so the 10000 tick label
# doesn't get cut off.
# My attempts to leave room for LaTeX (through TikZ) to
# use the 12 point font for Systematic Biology journal failed.
par(mar=c(4.5,4,0,1))

# Initialize the plot.
# Turn off axes and annotations (axis labels) so we can customize them.
plot(mytable$first.split.informative, pch=1, type='o',
	xlab='sequence length', ylab='number correct', ylim=c(0, 1000),
	axes=FALSE, ann=FALSE)

# make the x axis with a logarithmic scale
axis(1, at=c(1, 15, 29), lab=c('100', '1000', '10000'))

# make the y axis
axis(2, las=1, at=c(0, 200, 400, 600, 800, 1000))

# draw a box around the plot
#box()

# label the x and y axes
title(xlab='sequence length (base pairs)')
title(ylab='number correct (whole tree or first split)')

# draw more lines onto the plot
lines(mytable$nsuccesses.both, pch=2, type='o')
lines(mytable$nsuccesses.nj.only, pch=3, type='o')
lines(mytable$nsuccesses.topdown.only, pch=4, type='o')
lines(mytable$nsuccesses.neither, pch=5, type='o')

# draw the legend for the plot
legend('right', c('first split', 'nj and rs', 'nj only', 'rs only', 'neither'), pch=c(1, 2, 3, 4, 5))

# stop writing the tikz figure
dev.off()

