library("phyloseq")
library("dplyr")
library("ggpubr")
library("Matrix")
library("reshape2")
library("vegan")
library("ggplot2")
library("dplyr")
library("microbiome")


dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")


input <- function(inputfile) {
  pfix = prefix()
  if (length(pfix) != 0) {
     prefix <- paste(pfix, "/", sep="")
  }
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
   # Need to get the three files
#paste(pfix, toString(parameters["inputfile",2]), sep="")
   otu.path <<- paste(pfix, toString(parameters["otufile", 2]), sep="")
   tree.path <<- paste(pfix, toString(parameters["tree", 2]), sep="")
   map.path <<- paste(pfix, toString(parameters["mapping", 2]), sep="")
   colors.path <<- paste(pfix, toString(parameters["colors", 2]), sep="")
physeq1 <<- read_csv2phyloseq(otu.file=otu.path, taxonomy.file=tree.path, metadata.file=map.path, sep=",")
   distmeth <<- toString(parameters["distance", 2])
   diffmeth <<- toString(parameters["differential", 2])
   column <<- toString(parameters["column", 2])
   allGroupsColors <<- readLines(file(colors.path, "r"))
}

run <- function() {
relab_genera <<- transform_sample_counts(physeq1, function(x) x / sum(x) * 100)
abrel_bray <- phyloseq::distance(relab_genera, method = distmeth)
abrel_bray <- as.matrix(abrel_bray)
sub_dist <- list()
#print(column)
#print((sample_data(relab_genera))[["Description"]])
groups_all <- sample_data(relab_genera)[[column]]
#groups_all <- attributes(sample_data(relab_genera))[column]#sample_data(relab_genera)$as.name(column)
for (group in groups_all) {
     row_group <- which(groups_all == group)
     sample_group <- sample_names(relab_genera)[row_group]
     sub_dist[[group]] <- abrel_bray[ sample_group, sample_group]
     sub_dist[[group]][!lower.tri(sub_dist[[group]])] <- NA
}
braygroups<- melt(sub_dist)
df.bray <<- braygroups[complete.cases(braygroups), ]
df.bray$L1 <<- factor(df.bray$L1, levels=names(sub_dist))
}


output <- function(outputfile) {
	pdf(paste(outputfile,"pdf",sep="."))#,  width = 10*300,        # 5 x 300 pixels
y <- ggplot(df.bray, aes(x=L1, y=value, colour=L1)) +
     geom_jitter() +
   geom_boxplot(alpha=0.6) +
   theme(legend.position="none") +
   scale_color_manual(values = allGroupsColors) +
   ylab("Beta Diversity") +
   theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12), axis.text.y=element_text(size=12))
ord = ordinate(relab_genera, method=diffmeth, distance = distmeth)
z <- plot_ordination(relab_genera, ord, color = column, label="Name") +
geom_point(size=4) +
   scale_color_manual(values = allGroupsColors) +
stat_ellipse()
#stat_ellipse(aes(group=Description))
write.csv(y$data, paste(outputfile,"csv",sep="."))
print(y)
print(z)
#samples <- data.frame(sample_data(relab_genera))
#adonis(abrel_bray ~ Test, data = samples)
dev.off()
}
