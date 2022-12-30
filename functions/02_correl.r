#Modified by Fabio Marroni on 2022/06/05

# Run with --help flag for help.
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-A", "--abundance"), type="character", default="",
              help="Abundance file in bracken output format [default= %default]", metavar="character"),
  make_option(c("-O", "--outdir"), type="character", default="", 
              help="output dir [default= %default]", metavar="character"),
  make_option(c("-r", "--removeme"), type="character", default="L471,L593", 
              help="Comma separate list of samples to be removed, e.g. because of low yield [default= %default]", metavar="character"),
  make_option(c("-F", "--usefdr"), type="logical", default=FALSE, 
              help="Should pvalue b corrected for multiple testing using fdr? [default= %default]", metavar="character"),
  make_option(c("-G", "--genescorr"), type="numeric", default=20, 
              help="Number of genes for which to compute the correlation [default= %default]", metavar="character"),
  make_option(c("-C", "--chempar"), type="character", default="", 
              help="Chemical parameter file (sample names in first column) [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$abundance)) {
  stop("WARNING: No abundance specified with '-I' flag.")
} else {  cat ("abundance is ", opt$abundance, "\n")
  abundance <- opt$abundance  
  }

if (is.null(opt$chempar)) {
  stop("WARNING: No chempar specified with '-C' flag.")
} else {  cat ("chempar is ", opt$chempar, "\n")
  chempar <- opt$chempar  
  }

if (is.null(opt$removeme)) {
  stop("WARNING: No removeme specified with '-r' flag.")
} else {  cat ("removeme is ", opt$removeme, "\n")
  removeme <- opt$removeme  
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
	countdata<-fread(abundance,data.table=F)
	countdata$taxonomy_id<-NULL
	names(countdata)[grep("perc",names(countdata))]<-paste0("L",unlist(lapply(strsplit(names(countdata)[grep("perc",names(countdata))],"-"),"[",2)))
	#Patch to sum genera with same name (but different taxon ID): this is basically used only for unculutured
	if(length(unique(countdata$name))<nrow(countdata))
	{
	countdata<-aggregate(countdata[,2:ncol(countdata)],by=list(countdata$name),FUN="sum")
	names(countdata)[1]<-"name"
	}
	rownames(countdata)<-countdata$name
	countdata$name<-NULL
	tot<-rowSums(countdata)
	genescorr<-min(genescorr,nrow(countdata))
	countdata<-countdata[order(tot,decreasing=T),][1:genescorr,]
	chemdata<-fread(chempar,data.table=F)
	mynames<-c("Sample",names(chemdata)[2:ncol(chemdata)])
	chemdata<-aggregate(chemdata[,2:ncol(chemdata)],by=list(chemdata[,1]),FUN="mean")
	names(chemdata)<-mynames
	tcount<-data.frame(t(countdata),check.names=F)
	tcount$Unmapped<-NULL
	tcount<-merge(tcount,chemdata,by.x="row.names",by.y="Sample",sort=F)
	row.names(tcount)<-tcount[,1]
	tcount[,1]<-NULL
	#Filter to remove values with "crazy" distributions.
	#Count how many zeros
	#We ask that we have less than 50% of 0s
	zerocount<-apply(apply(tcount,1,">",0),1,sum,na.rm=T)
	tcount<-tcount[,zerocount/nrow(tcount)>0.5]
	myprop<-mynames[2:length(mynames)]
	#Remove unwanted parameters (and in the future translate names in English)
	myprop<-myprop[!myprop%in%c("Area TOTALE")]
	#Only use retained species and measures to build the database for correlation analysis
	mygenes<-row.names(countdata)[row.names(countdata)%in%names(tcount)]
	myprop<-myprop[myprop%in%names(tcount)]
	if(length(myprop)==1) cdf<-data.frame(Species=mygenes,Parameter=myprop,stringsAsFactors=F)
	if(length(myprop)>1) 
	{
	cdf<-data.frame(expand.grid(mygenes,myprop),stringsAsFactors=F)
	names(cdf)<-c("Species","Parameter")
	}
	cdf$Species<-as.character(cdf$Species)
	cdf$Parameter<-as.character(cdf$Parameter)
	myfile<-paste0(outdir,basename(chempar))
	pdf(gsub(".txt","_allcor.pdf",myfile))
	for(aaa in 1:nrow(cdf))
	{
		pp<-cor.test(tcount[,cdf$Species[aaa]],tcount[,cdf$Parameter[aaa]],method="spearman")
		cdf$rho[aaa]<-round(pp$estimate,3)
		cdf$pvalue[aaa]<-signif(pp$p.value,3)
		myxlim<-c(min(tcount[,cdf$Species[aaa]]),max(tcount[,cdf$Species[aaa]]))
		myylim<-c(min(tcount[,cdf$Parameter[aaa]]),max(tcount[,cdf$Parameter[aaa]]))
		plot(tcount[,cdf$Species[aaa]],tcount[,cdf$Parameter[aaa]],xlab=cdf$Species[aaa],ylab=cdf$Parameter[aaa])
		text(myxlim[2]-0.1*(myxlim[2]-myxlim[1]),myylim[2]-0.05*(myylim[2]-myylim[1]),bquote(rho == .(cdf$rho[aaa])))
		text(myxlim[2]-0.1*(myxlim[2]-myxlim[1]),myylim[2]-0.1*(myylim[2]-myylim[1]),bquote(italic(p-value) == .(cdf$pvalue[aaa])))
	}
	cdf$fdr<-p.adjust(cdf$pvalue,method="fdr")
	write.table(cdf,gsub(".txt","_correlations.txt",myfile),sep="\t",quote=F,row.names=F)
	cormat<-corp<-corf<-matrix(NA,nrow=length(myprop),ncol=length(mygenes),dimnames=list(myprop,mygenes))
	#I run this second loop instead of using just one because I want to use fdr, and this is computed AFTER finisheing the first loop
	browser()	
	for(aaa in 1:nrow(cdf))
	{
		cormat[cdf$Parameter[aaa],cdf$Species[aaa]]<-cdf$rho[aaa]
		corf[cdf$Parameter[aaa],cdf$Species[aaa]]<-cdf$fdr[aaa]
		corp[cdf$Parameter[aaa],cdf$Species[aaa]]<-cdf$pvalue[aaa]
	}
	
	dev.off()
	col2 = colorRampPalette(c('darkred', 'white', 'blue3')) 
	#Only select classes that have at least one significant result
	if(usefdr) corp<-corf
	sig<-corp<=0.05
	# tsig<-apply(sig,1,"sum")
	# cormat<-cormat[tsig>0,]
	# corp<-corp[tsig>0,]
	#Only select large values
	# large<-abs(cormat)>0.2
	# tlarge<-apply(large,1,"sum")
	# cormat<-cormat[tlarge>0,]
	# corp<-corp[tlarge>0,]
	if(nrow(corp)>1&ncol(corp)>1)
	{
	mytitle<-paste0("Correlation between ", gsub(".txt","",basename(chempar)), " and species abundance")
	pdf(gsub(".txt",".pdf",myfile))
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
