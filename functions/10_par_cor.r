#Modified by Fabio Marroni on 2022/12/30

# Run with --help flag for help.
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-p", "--pH"), type="character", default="",
              help="pH values file [default= %default]", metavar="character"),
  make_option(c("-O", "--outfile"), type="character", default="", 
              help="output dir [default= %default]", metavar="character"),
  make_option(c("-F", "--usefdr"), type="logical", default=FALSE, 
              help="Should pvalue b corrected for multiple testing using fdr? [default= %default]", metavar="character"),
  make_option(c("-C", "--chempar"), type="character", default="", 
              help="Chemical parameter file (sample names in first column) [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$pH)) {
  stop("WARNING: No pH specified with '-I' flag.")
} else {  cat ("pH is ", opt$pH, "\n")
  pH <- opt$pH  
  }

if (is.null(opt$usefdr)) {
  stop("WARNING: No usefdr specified with '-r' flag.")
} else {  cat ("usefdr is ", opt$usefdr, "\n")
  usefdr <- opt$usefdr  
  }

  if (is.null(opt$outfile)) {
  stop("WARNING: No output directory specified with '-O' flag.")
} else {  cat ("Output dir is ", opt$out, "\n")
  outfile <- opt$outfile 
  }

if (is.null(opt$chempar)) {
  stop("WARNING: No chempar specified with '-C' flag.")
} else {  cat ("chempar is ", opt$chempar, "\n")
  chempar <- opt$chempar  
  }

correl<-function() 
{
    library(data.table)
    pHdat<-fread(pH,data.table=F)
    pHdat<-aggregate(pHdat[,2],by=list(pHdat$Campione),FUN="mean")
    names(pHdat)<-c("Name","pH")
    chemdat<-fread(chempar,data.table=F)
    fulldat<-merge(pHdat,chemdat,sort=F)
    chemparnames<-names(fulldat)[-1]
    #This would perform an n by n correlation, but I want only pH against everything else, and I also want p-value.
    #cor(fulldat[,2:ncol(fulldat)])
    
	cdf<-data.frame(expand.grid(chemparnames,chemparnames,stringsAsFactors=F))
    cdf$rho<-cdf$pvalue<-NA
    for(aaa in 1:nrow(cdf))
    {
    tcor<-cor.test(fulldat[,cdf$Var1[aaa]],fulldat[,cdf$Var2[aaa]],method="spearman",exact=F)
    cdf$rho[aaa]<-tcor$estimate
    cdf$pvalue[aaa]<-tcor$"p.value"
    }
    write.table(cdf,outfile,sep="\t",row.names=F,quote=F)
}
correl()
