

# A treemap R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0002376","immune system process",0.600,4.2083,0.925,0.000,"immune system process"),
c("GO:0006952","defense response",0.568,9.0329,0.453,0.000,"defense response"),
c("GO:0034097","response to cytokine",0.136,5.0177,0.356,0.314,"defense response"),
c("GO:0009607","response to biotic stimulus",0.342,4.6946,0.545,0.340,"defense response"),
c("GO:0045087","innate immune response",0.148,7.9245,0.263,0.504,"defense response"),
c("GO:0035456","response to interferon-beta",0.004,4.2757,0.366,0.557,"defense response"),
c("GO:0070887","cellular response to chemical stimulus",1.007,4.5560,0.380,0.603,"defense response"),
c("GO:0007155","cell adhesion",0.544,4.8069,0.925,0.000,"cell adhesion"),
c("GO:0021527","spinal cord association neuron differentiation",0.003,4.6108,0.908,0.000,"spinal cord association neuron differentiation"),
c("GO:0022610","biological adhesion",0.550,4.7595,0.925,0.000,"biological adhesion"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap_SMBvsMPI.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
tmPlot(
	stuff,
	index = c("representative","description"),
	vSize = "abslog10pvalue",
	type = "categorical",
	vColor = "representative",
	title = "REVIGO Gene Ontology treemap",
	inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
	lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
	bg.labels = "#CCCCCCAA",     # define background color of group labels
												       # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
	position.legend = "none"
)

dev.off()
