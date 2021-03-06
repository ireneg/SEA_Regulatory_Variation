

# A plotting R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800


# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );
# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );

# SMBvsMPI -----------------------------------------------------------------

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0002376","immune system process", 0.600,-0.159,-7.386, 4.886,-4.2083,0.925,0.000),
c("GO:0006952","defense response", 0.568, 4.989,-1.157, 4.863,-9.0329,0.453,0.000),
c("GO:0007155","cell adhesion", 0.544,-4.167,-2.976, 4.844,-4.8069,0.925,0.000),
c("GO:0021527","spinal cord association neuron differentiation", 0.003,-0.542, 7.009, 2.543,-4.6108,0.908,0.000),
c("GO:0022610","biological adhesion", 0.550,-4.294, 2.431, 4.849,-4.7595,0.925,0.000),
c("GO:0034097","response to cytokine", 0.136, 5.957, 1.240, 4.242,-5.0177,0.356,0.314),
c("GO:0009607","response to biotic stimulus", 0.342, 6.980,-0.885, 4.643,-4.6946,0.545,0.340),
c("GO:0045087","innate immune response", 0.148, 5.519,-2.057, 4.279,-7.9245,0.263,0.504),
c("GO:0035456","response to interferon-beta", 0.004, 5.498, 2.217, 2.692,-4.2757,0.366,0.557),
c("GO:0070887","cellular response to chemical stimulus", 1.007, 5.494, 0.505, 5.111,-4.5560,0.380,0.603));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0.6) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$log10_p_value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ]; 
p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 );
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);
p1 <- p1 + ggtitle("Sumba vs Mappi")


# MTWvsMPI -----------------------------------------------------------------

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0006952","defense response", 0.568,-6.153,-0.130, 4.863,-6.6498,0.546,0.000),
c("GO:0007155","cell adhesion", 0.544, 4.888, 0.136, 4.844,-10.5528,0.527,0.000),
c("GO:0022610","biological adhesion", 0.550, 1.010,-6.348, 4.849,-11.1062,0.736,0.000),
c("GO:0032693","negative regulation of interleukin-10 production", 0.003, 0.993, 6.601, 2.587,-4.8182,0.696,0.000),
c("GO:0006955","immune response", 0.337,-5.788, 1.714, 4.635,-5.8069,0.413,0.340),
c("GO:0006954","inflammatory response", 0.110,-5.726,-1.220, 4.151,-5.1746,0.555,0.491));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p2 <- ggplot( data = one.data );
p2 <- p2 + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0.6) ) + scale_size_area();
p2 <- p2 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$log10_p_value), 0) );
p2 <- p2 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
p2 <- p2 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ]; 
p2 <- p2 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 );
p2 <- p2 + labs (y = "semantic space x", x = "semantic space y");
p2 <- p2 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p2 <- p2 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p2 <- p2 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);
p2 <- p2 + ggtitle("Mentawai vs Mappi")
