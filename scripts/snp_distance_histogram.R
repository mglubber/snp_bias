require(lattice)
jpeg(filename="alu_element_histogram.jpg")
histogram(~distance|alufamily,data=distance_vs_family,type="density",breaks=50)
dev.off()
