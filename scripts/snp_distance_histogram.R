# Create a histogram of the # of bp truncated from the end of the repeat element
library(Cairo)
require(lattice)
filename="alu_element_histogram.jpg"
Cairo(file=filename,type="jpeg",units="in",width=14,height=11,dpi=72)
histogram(~conMiss|alufamily,data=alu_ranges,type="density",breaks=50)
dev.off()
