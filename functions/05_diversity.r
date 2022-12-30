#Modified by Fabio Marroni on 2022/07/30

# Run with --help flag for help.
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-A", "--abundance"), type="character", default="",
              help="Abundance file in bracken output format [default= %default]", metavar="character"),
  make_option(c("-O", "--outfile"), type="character", default="", 
              help="output dir [default= %default]", metavar="character"),
  make_option(c("-R", "--removeme"), type="character", default="", 
              help="Comma delimited list of organisms to remove [default= %default]", metavar="character"),
  make_option(c("-m", "--mintot"), type="numeric", default=1, 
              help="Only keep species that have a total % abundance of mintot [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$abundance)) {
  stop("WARNING: No abundance specified with '-I' flag.")
} else {  cat ("abundance is ", opt$abundance, "\n")
  abundance <- opt$abundance  
  }

if (is.null(opt$mintot)) {
  stop("WARNING: No mintot specified with '-m' flag.")
} else {  cat ("mintot is ", opt$mintot, "\n")
  mintot <- opt$mintot  
  }

if (is.null(opt$removeme)) {
  stop("WARNING: No outfile specified with '-T' flag.")
} else {  cat ("removeme is ", opt$removeme, "\n")
  removeme <- opt$removeme  
  }

  if (is.null(opt$outfile)) {
  stop("WARNING: No output file specified with '-O' flag.")
} else {  cat ("Output file is ", opt$out, "\n")
  outfile <- opt$outfile 
  }

est_div<-function()
{
library(data.table)
library(vegan)
library(entropart)
species<-fread(abundance,data.table=F)
#Patch to sum genera with same name (but different taxon ID): this is basically used only for unculutured
if(length(unique(species$name))<nrow(species))
{
species<-aggregate(species[,2:ncol(species)],by=list(species$name),FUN="sum")
names(species)[1]<-"name"
}

toremove<-unlist(strsplit(removeme,","))
species<-species[!species$name%in%removeme,]
species<-species[order(species$Total,decreasing=TRUE),]
#This was for checking that I correctly removed Streptococcus thermophilus when launching from shell.
#print(head(mdat))

row.names(species)<-species$name
species$name<-species$taxonomy_id<-NULL
species<-species[species$Total>=mintot,]
species$Total<-NULL
names(species)<-paste0("L",unlist(lapply(strsplit(names(species),"-"),"[",2)))
shannon<-diversity(t(species))
obstaxa<-specnumber(t(species))
simpson<-diversity(t(species),index="simpson")
fisher<-fisher.alpha(t(species))
chao1<-estimateR(t(species))[2,]
ACE<-estimateR(t(species))[4,]
pielou<-diversity(t(species))/log(specnumber(t(species)))
good<-apply(t(species),1,Coverage,CheckArguments=FALSE)
myres<-data.frame(Sample=names(species),shannon=shannon,obstaxa=obstaxa,simpson=simpson,fisher=fisher,chao1=chao1,ACE=ACE,pielou=pielou,good=good)
write.table(myres,outfile,sep="\t",quote=F,row.names=F)
}
est_div()

