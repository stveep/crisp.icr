#' Add a column with reading frame calculated from the number of deleted/inserted bases. 0 = native reading frame.
#'
add_reading_frame <- function(x) {x$frame = (x$Deletions+x$Insertions) %% 3; x}


#' Aggregates mutation data for plotting(generated as a report with the -r flag from the mutations tool: github.com/stveep/bio-sam-mutation).
#' Also creates a column with normalised read counts.
aggregate_mutation_data <- function(y) {
	y.ag <- aggregate(y$Reads, by=list(y$File,y$frame,y$Allele),FUN=mean)
	names(y.ag) <- c("File","Frame","Allele","Reads")
	totals <- by(y$Reads,y$File,sum)
	y.ag$total = totals[y.ag$File]
	y.ag$norm = y.ag$Reads/y.ag$total
	y.ag
}


#! Plot
reading_frame_plot <- function(y) {
	library('ggplot2')
	y.ag = aggregate_mutation_data(y)
	ggplot(y.ag,aes(x= factor(Frame),y=Reads,fill=factor(Allele))) + 
    	  geom_bar(stat="identity",aes(width=0.5))                       + 
    	  facet_grid(~ File)                                             + 
    	  labs(x="Frameshift")                                           + 
    	  theme(legend.position="none", text=element_text(size = 18))    +  
    	  scale_fill_manual(values=cbPalette)
}

norm_reading_frame_plot <- function(y) {
	library('ggplot2')
	y.ag = aggregate_mutation_data(y)
	ggplot(y.ag,aes(x= factor(Frame),y=norm,fill=factor(Allele))) +
	  geom_bar(stat="identity",aes(width=0.5))                     +
	  facet_grid(~ File)                                           +
	  labs(x="Frameshift",y="Fraction of reads")                   +
	  theme(legend.position="none", text=element_text(size = 18))  +  
	  scale_fill_manual(values=cbPalette)
}

