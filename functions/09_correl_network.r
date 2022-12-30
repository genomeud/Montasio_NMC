#Modified by Fabio Marroni on 2022/06/05

# Run with --help flag for help.
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-A", "--abundance"), type="character", default="/projects/populus/ep/share/marroni/collaborations/lattoinnesto/05_tables_silva/genus_bracken.txt",
              help="Abundance file in bracken output format [default= %default]", metavar="character"),
  make_option(c("-O", "--outdir"), type="character", default="/projects/populus/ep/share/marroni/collaborations/lattoinnesto/06_chempar_silva/", 
              help="output dir [default= %default]", metavar="character"),
  make_option(c("-g", "--outgraph"), type="character", default="/projects/populus/ep/share/marroni/collaborations/lattoinnesto/06_chempar_silva/sugar_aroma_netcorr.pdf", 
              help="output dir [default= %default]", metavar="character"),
  make_option(c("-r", "--removeme"), type="character", default="", 
              help="Comma separated list of samples to be removed, e.g. because of low yield [default= %default]", metavar="character"),
  make_option(c("-F", "--usefdr"), type="logical", default=FALSE, 
              help="Should pvalue be corrected for multiple testing using fdr? [default= %default]", metavar="character"),
  make_option(c("-G", "--genescorr"), type="numeric", default=20, 
              help="Number of genes for which to compute the correlation [default= %default]", metavar="character"),
  make_option(c("-a", "--aroma"), type="character", default="/projects/populus/ep/share/marroni/collaborations/lattoinnesto/06_chempar/aromi.txt", 
              help="File with data regarding aromatic compounds concentration (sample names in first column) [default= %default]", metavar="character"),
  make_option(c("-S", "--sugar_acid"), type="character", default="/projects/populus/ep/share/marroni/collaborations/lattoinnesto/06_chempar/Acid_Sugars_mean.txt", 
              help="File with data regarding sugar and acid concentration (sample names in first column) [default= %default]", metavar="character"),
  make_option(c("-C", "--chempar"), type="character", default="/projects/populus/ep/share/marroni/collaborations/lattoinnesto/06_chempar/Treatment.txt", 
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

if (is.null(opt$aroma)) {
  stop("WARNING: No aroma specified with '-a' flag.")
} else {  cat ("aroma is ", opt$aroma, "\n")
  aroma <- opt$aroma  
  }

if (is.null(opt$sugar_acid)) {
  stop("WARNING: No sugar_acid specified with '-S' flag.")
} else {  cat ("sugar_acid is ", opt$sugar_acid, "\n")
  sugar_acid <- opt$sugar_acid  
  }

  if (is.null(opt$outdir)) {
  stop("WARNING: No output directory specified with '-O' flag.")
} else {  cat ("Output dir is ", opt$out, "\n")
  outdir <- opt$outdir 
  }

  if (is.null(opt$outgraph)) {
  stop("WARNING: No output graph specified with '-g' flag.")
} else {  cat ("graph dir is ", opt$outgraph, "\n")
  outgraph <- opt$outgraph 
  }


#Develop a function to draw a three layered network
three_layered<-function(x,Lev1,Lev2,Lev3)
{
x[Lev2,2]<-0.5
x[Lev1,2]<-0
x[Lev3,2]<-1
x[,2]<-jitter(x[,2])
x
}

three_layered_vert<-function(x,Lev1,Lev2,Lev3,jitter=F)
{
y<-x
y[,1]<-x[,2]
y[,2]<-x[,1]
x<-y
x[Lev1,1]<-0
x[Lev2,1]<-0.5
x[Lev3,1]<-1
x[Lev1,2]<-seq(0,1,length.out=length(Lev1))
x[Lev2,2]<-seq(0,1,length.out=length(Lev2))
x[Lev3,2]<-seq(0,1,length.out=length(Lev3))
if(jitter) x[,1]<-jitter(x[,1])
x
}

correl<-function() 
{
	library(igraph)
	library(WGCNA)
	#library(corrplot)
    library(data.table)
	#library(openxlsx)
	library("RColorBrewer")
	library("ggplot2")
	countdata<-fread(abundance,data.table=F)
	countdata$taxonomy_id<-countdata$Total<-NULL
	names(countdata)[grep("perc",names(countdata))]<-paste0("L",unlist(lapply(strsplit(names(countdata)[grep("perc",names(countdata))],"-"),"[",2)))
	#Patch to sum genera with same name (but different taxon ID): this is basically used only for unculutured
	if(length(unique(countdata$name))<nrow(countdata))
	{
	countdata<-aggregate(countdata[,2:ncol(countdata)],by=list(countdata$name),FUN="sum")
	names(countdata)[1]<-"name"
	}
	rownames(countdata)<-countdata$name
	countdata$name<-NULL
	#Eventually remove unwanted species
	toremove<-unlist(strsplit(removeme,","))
	countdata<-countdata[!row.names(countdata)%in%toremove,]
	#Remove contamination from Homo
	countdata<-countdata[!row.names(countdata)%in%c("Homo","Homo sapiens"),]
	genescorr<-min(genescorr,nrow(countdata))
	tot<-rowSums(countdata)
	countdata<-countdata[order(tot,decreasing=T),][1:genescorr,]
	tcount<-data.frame(t(countdata),check.names=F)
	bnames<-row.names(countdata)
	aromadata<-fread(aroma,data.table=F)
	aronames<-c("Sample",names(aromadata)[2:ncol(aromadata)])
	aromadata<-aggregate(aromadata[,2:ncol(aromadata)],by=list(aromadata[,1]),FUN="mean")
	names(aromadata)<-aronames
	tcount<-merge(tcount,aromadata,by.x="row.names",by.y="Sample",sort=F)
	row.names(tcount)<-tcount[,1]
	tcount[,1]<-NULL
	sugardata<-fread(sugar_acid,data.table=F)
	sugarnames<-c("Sample",names(sugardata)[2:ncol(sugardata)])
	sugardata<-aggregate(sugardata[,2:ncol(sugardata)],by=list(sugardata[,1]),FUN="mean")
	names(sugardata)<-sugarnames
	chemdata<-fread(chempar,data.table=F)
	chemdata<-chemdata[,names(chemdata)!="TreatTime (sec)"]
	sugardata<-merge(sugardata,chemdata,by="Sample",sort=F)
	sugardata<-sugardata[,grep("Acid",names(sugardata),invert=T)]
	sugarnames<-names(sugardata)
	tcount<-merge(tcount,sugardata,by.x="row.names",by.y="Sample",sort=F)
	row.names(tcount)<-tcount[,1]
	tcount[,1]<-NULL
	tcount$TOTALE<-NULL

	#Remove parameters/bacteria with more than 50% empty data (too much missing data do not allow reliable correlation estimate)
	zerocount<-apply(apply(tcount,1,">",0),1,sum,na.rm=T)
	tcount<-tcount[,zerocount/nrow(tcount)>0.5]
	allcor<-cor(tcount,method="spearman",use="p")
	allp<-corPvalueStudent(allcor,nSamples=nrow(tcount))
	#Only keep significant correlations
	#allcor[allcor<0.5]<-0
	allcor[allp>=0.05]<-0
	#Only keep correlation between bacteria and parameters (i.e. set to zero all correlation between parameters)
	allcor[row.names(allcor)%in%c(aronames,sugarnames),row.names(allcor)%in%c(aronames,sugarnames)]<-0
	#Also remove correl between bacteria. We only want bacteria vs param
	allcor[row.names(allcor)%in%bnames,row.names(allcor)%in%bnames]<-0
	mynet<-graph_from_adjacency_matrix(allcor,weighted=T,mode="undirected",diag=F)
	#Remove isolated nodes
	Isolated = which(degree(mynet)==0)	
	mynet<-delete.vertices(mynet, Isolated)
	#Only keep bacterial names that are retained in the network
	bnames<-bnames[bnames%in%V(mynet)$name]
	#Attribute size proportional to abundance
	V(mynet)$size<-rep(14,length(V(mynet)$name))
	#Add a fixed amount to size
	V(mynet)$size[V(mynet)$name%in%bnames]<-V(mynet)$size[V(mynet)$name%in%bnames]+sqrt(colSums(tcount[,bnames]))
	#We attribute different colors to bacteria and parameters
	V(mynet)$color<-"lightblue"
	V(mynet)$frame.color<-"lightblue"
	V(mynet)$color[V(mynet)$name%in%aronames]<-"orange"
	V(mynet)$frame.color[V(mynet)$name%in%aronames]<-"orange"
	V(mynet)$color[V(mynet)$name%in%sugarnames]<-"red"
	V(mynet)$frame.color[V(mynet)$name%in%sugarnames]<-"red"
	V(mynet)$label.color<-"black"
	#We define names of colors for positive and negativ correlation
	poscol<-"purple"
	negcol<-"black"
	#We plot the colors with some transparency
	E(mynet)$color[E(mynet)$weight > 0 ] <- rgb(col2rgb(poscol)[1],col2rgb(poscol)[2],col2rgb(poscol)[3],max=255,alpha=125)
	E(mynet)$color[E(mynet)$weight <= 0 ] <- rgb(col2rgb(negcol)[1],col2rgb(negcol)[2],col2rgb(negcol)[3],max=255,alpha=125)
	pp<-rank(abs(E(mynet)$weight))
	E(mynet)$width<-seq(2,3,length.out=length(pp))[pp]
	L1<-layout_with_fr(mynet)
	L2<-layout_on_grid(mynet)
	L3<- layout_as_star(mynet)
	L4<- layout_as_tree(mynet)
	L5<-layout_in_circle(mynet)
	L6<-layout_nicely(mynet) 
	L7<-layout_on_sphere(mynet) 
	L8<-layout_randomly(mynet) 
	L9<-layout_with_dh(mynet) 
	L10<-layout_with_gem(mynet) 
	L11<-layout_with_graphopt(mynet) 
	L12<-layout_with_lgl(mynet)
	nob<-which(!V(mynet)$name%in%bnames)
	Lev1<-unlist(split(nob,cut(nob,2,labels=F))[[1]])
	#Lev1<-which(V(mynet)$name%in%sugarnames)
	Lev2<-which(V(mynet)$name%in%bnames)
	#Lev3<-which(V(mynet)$name%in%aronames)
	Lev3<-unlist(split(nob,cut(nob,2,labels=F))[[2]])
	L13<-three_layered(L2,Lev1=Lev1,Lev2=Lev2,Lev3=Lev3)
	L14<-three_layered_vert(L2,Lev1=Lev1,Lev2=Lev2,Lev3=Lev3)
	pdf(outgraph)
	par(mar=c(1,2,1,2))
#	for(aaa in paste0("L",seq(1,14)))
	for(aaa in "L14")
	{
	plot(mynet,layout=eval(as.name(aaa)))
	}
	dev.off()
}
correl()

