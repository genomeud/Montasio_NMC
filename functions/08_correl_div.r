#Modified by Fabio Marroni on 2023/04/30

# Run with --help flag for help.
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-D", "--diversity"), type="character", default="/projects/populus/ep/share/marroni/collaborations/lattoinnesto/05_tables/diversity.txt",
              help="Table with diversity indices [default= %default]", metavar="character"),
  make_option(c("-I", "--indices"), type="character", default="obstaxa,simpson,shannon,pielou", 
              help="comma separated list of diversity indices to plot (names should exactly match those in the file specified by the 'D' parameter [default= %default]", metavar="character"),
  make_option(c("-O", "--outdir"), type="character", default="/projects/populus/ep/share/marroni/collaborations/lattoinnesto/06_chempar/", 
              help="output dir [default= %default]", metavar="character"),
  make_option(c("-F", "--usefdr"), type="logical", default=FALSE, 
              help="Should pvalue b corrected for multiple testing using fdr? [default= %default]", metavar="character"),
  make_option(c("-G", "--genescorr"), type="numeric", default=20, 
              help="Number of genes for which to compute the correlation [default= %default]", metavar="character"),
  make_option(c("-C", "--chempar"), type="character", default="/projects/populus/ep/share/marroni/collaborations/lattoinnesto/06_chempar/Treatment.txt", 
              help="Chemical parameter file (sample names in first column) [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$diversity)) {
  stop("WARNING: No diversity specified with '-D' flag.")
} else {  cat ("diversity is ", opt$diversity, "\n")
  diversity <- opt$diversity  
  }

if (is.null(opt$indices)) {
  stop("WARNING: No indices specified with '-I' flag.")
} else {  cat ("indices is ", opt$indices, "\n")
  indices <- opt$indices  
  }

if (is.null(opt$chempar)) {
  stop("WARNING: No chempar specified with '-C' flag.")
} else {  cat ("chempar is ", opt$chempar, "\n")
  chempar <- opt$chempar  
  }

if (is.null(opt$usefdr)) {
  stop("WARNING: No usefdr specified with '-r' flag.")
} else {  cat ("usefdr is ", opt$usefdr, "\n")
  usefdr <- opt$usefdr  
  }

if (is.null(opt$genescorr)) {
  stop("WARNING: No genescorr specified with '-V' flag.")
} else {  cat ("genescorr is ", opt$genescorr, "\n")
  genescorr <- opt$genescorr  
  }

  if (is.null(opt$outdir)) {
  stop("WARNING: No output directory specified with '-O' flag.")
} else {  cat ("Output dir is ", opt$out, "\n")
  outdir <- opt$outdir 
  }


correl<-function() 
{
	library(corrplot)
    library(data.table)
	library(openxlsx)
	library("RColorBrewer")
	library("ggplot2")
	countdata<-fread(diversity,data.table=F)
	chemdata<-fread(chempar,data.table=F)
	#We sort based on the first row, which is assumed to be the sample
	#Well, actually, I guess this liine is useless
	chemdata<-chemdata[order(chemdata[,1]),]
	#Account for situations in which the same genus or species is associated to different taxa (it is a rare event, but it happens) 
	row.names(countdata)<-countdata[,1]
	countdata[,1]<-NULL
	#Select indices specified by user
	keepind<-unlist(strsplit(indices,","))
	countdata<-countdata[,keepind]
	chemdata<-aggregate(chemdata[,2:ncol(chemdata)],by=list(chemdata[,1]),FUN="mean")
	row.names(chemdata)<-chemdata[,1]
	chemdata[,1]<-NULL
	mydiv<-names(countdata)
	mychem<-names(chemdata)
	tcount<-merge(countdata,chemdata,by.x="row.names",by.y="row.names",sort=F)
	row.names(tcount)<-tcount[,1]
	tcount[,1]<-NULL
	#Filter to remove values with "crazy" distributions (i.e. too many zeros or zero variance).
	#Count how many zeros
	#We ask that we have less than 60% of 0s
	zerocount<-apply(apply(tcount,1,">",0),1,sum,na.rm=T)
	tcount<-tcount[,zerocount/nrow(tcount)>0.6]
	#Remove values with zero variance
    myvar<-apply(tcount,2,var,na.rm=T)
    tcount<-tcount[,myvar>0]	
	mydiv<-mydiv[mydiv%in%names(tcount)]
	mychem<-mychem[mychem%in%names(tcount)]
	mychem<-mychem[!mychem%in%c("Area TOTALE")]

	if(length(mychem)==1) cdf<-data.frame(Diversity=mydiv,Parameter=mychem,stringsAsFactors=F)
	if(length(mychem)>1) 
	{
	cdf<-data.frame(expand.grid(mydiv,mychem),stringsAsFactors=F)
	names(cdf)<-c("Diversity","Parameter")
	}
	cdf$Diversity<-as.character(cdf$Diversity)
	cdf$Parameter<-as.character(cdf$Parameter)
	myfile<-paste0(outdir,basename(chempar))
	pdf(gsub(".txt","_allcor_div.pdf",myfile))
	for(aaa in 1:nrow(cdf))
	{
		pp<-cor.test(tcount[,cdf$Diversity[aaa]],tcount[,cdf$Parameter[aaa]],method="spearman")
		cdf$rho[aaa]<-round(pp$estimate,3)
		cdf$pvalue[aaa]<-signif(pp$p.value,3)
		myxlim<-c(min(tcount[,cdf$Diversity[aaa]]),max(tcount[,cdf$Diversity[aaa]]))
		myylim<-c(min(tcount[,cdf$Parameter[aaa]]),max(tcount[,cdf$Parameter[aaa]]))
		plot(tcount[,cdf$Diversity[aaa]],tcount[,cdf$Parameter[aaa]],xlab=cdf$Diversity[aaa],ylab=cdf$Parameter[aaa])
		text(myxlim[2]-0.1*(myxlim[2]-myxlim[1]),myylim[2]-0.05*(myylim[2]-myylim[1]),bquote(rho == .(cdf$rho[aaa])))
		text(myxlim[2]-0.1*(myxlim[2]-myxlim[1]),myylim[2]-0.1*(myylim[2]-myylim[1]),bquote(italic(p-value) == .(cdf$pvalue[aaa])))
	}
	cdf$fdr<-p.adjust(cdf$pvalue,method="fdr")
	write.table(cdf,gsub(".txt","_correlations_div.txt",myfile),sep="\t",quote=F,row.names=F)
	cormat<-corp<-corf<-matrix(NA,nrow=length(mychem),ncol=length(mydiv),dimnames=list(mychem,mydiv))
	#I run this second loop instead of using just one because I want to use fdr, and this is computed AFTER finisheing the first loop
	for(aaa in 1:nrow(cdf))
	{
		cormat[cdf$Parameter[aaa],cdf$Diversity[aaa]]<-cdf$rho[aaa]
		corf[cdf$Parameter[aaa],cdf$Diversity[aaa]]<-cdf$fdr[aaa]
		corp[cdf$Parameter[aaa],cdf$Diversity[aaa]]<-cdf$pvalue[aaa]
	}
	
	dev.off()
	col2 = colorRampPalette(c('red', 'white', 'blue3')) 
	#Only select classes that have at least one significant result
	if(usefdr) corp<-corf
	sig<-corp<=0.05
	if(nrow(corp)>1&ncol(corp)>1)
	{
	#Aesthetic changes to names
	colnames(cormat)[colnames(cormat)=="obstaxa"]<-"N. of species"
	colnames(cormat)[colnames(cormat)=="simpson"]<-"Simpson"
	colnames(cormat)[colnames(cormat)=="shannon"]<-"Shannon"
	colnames(cormat)[colnames(cormat)=="pielou"]<-"Pielou"
	row.names(cormat)[row.names(cormat)=="TreatTemp"]<-"Treatment\nTemperature [°C]"
	row.names(cormat)[row.names(cormat)=="TreatTime (sec)"]<-"Treatment Time [s]"
	row.names(cormat)[row.names(cormat)=="IncTemp"]<-"Incubation\nTemperature [°C]"
	row.names(cormat)[row.names(cormat)=="IncTime (min)"]<-"Incubation Time [min]"
	colnames(corp)<-colnames(cormat)
	row.names(corp)<-row.names(cormat)
	mytitle<-paste0("Correlation between ", gsub(".txt","",basename(chempar)), " and diversity indices")
	pdf(gsub(".txt","_div.pdf",myfile))
	corrplot(cormat, p.mat = corp, 
         tl.col="black", tl.cex=0.8,
         diag = TRUE, type = 'full',
         col=col2(256),
		 addgrid.col=NA,
         insig = 'label_sig', sig.level = c(0.05), cl.align.text ="l", cl.cex =0.7,
		 mar=c(0,0.5,1.2,0.5),oma=c(0,0,0),
         pch.cex = 0.9, pch.col = 'black')
	dev.off()
	}
}
correl()
