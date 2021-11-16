# Network-analysis-of-root-and-leaf-transcriptome-integration
This is transcriptomic data analysis by network deconvolution based on correlation inference tool (corto, Giorgi et al., 2020) to find key transcription factors, known as master regulators, responsibles of the symthomps of complex genetic diseases.

For this analysis we used Solanum lycopersicum (tomato) infected by two different strains of the Potato Spindle Tubercule Viroid (PSTVd) as study model.

Some R libraries are need for this analysis:
```R
library(corto)
library(affy)
library(GEOquery)
```
## Obtain data from Gene Expression Omnibus
```R
#Leaves data
getGEOSuppFiles("GSE106912")
#Root data
getGEOSuppFiles("GSE111736")
```
## Generate combined expression matrix
```R
#Leaves data
setwd("C:/Users/octav/OneDrive/Escritorio/Tomato_viroid_project/GSE106912")
untar("GSE106912_RAW.tar", exdir = "data")
cels = list.files("data/", pattern = "cel")
cels = list.files("data/", pattern = "CEL")
sapply(paste("data", cels, sep = "/"), gunzip)
cels = list.files("data/", pattern = "CEL")

#Root data
setwd("C:/Users/octav/OneDrive/Escritorio/Tomato_viroid_project/GSE111736")
untar("GSE111736_RAW.tar", exdir = "data")
cels = list.files("data/", pattern = "cel")
sapply(paste("data", cels, sep = "/"), gunzip)
cels = list.files("data/", pattern = "CEL")
```

Select interest samples and move them to an specific folder (in this case GSE111736/data)
```R
setwd("C:/Users/octav/OneDrive/Escritorio/Tomato_viroid_project/GSE111736/data")
cels = list.files("data/", pattern = "CEL")
#Obtaing raw data
raw.data = ReadAffy(verbose = FALSE, filenames = cels)
```
## Matrix normalization
```R
#As indicated in affy library#
data.rma.norm = rma(raw.data)
inmat = exprs(data.rma.norm)
```

## Assigning centroids for network deconvolution

Use a list of all the TFs in the network you want to analize
```R
Sly_TFs<-read.table("Tomato_TFs.txt", stringsAsFactors = FALSE)
Sly_TFs<-Sly_TFs[[1]]
#Get only the TFs contained in your expression matrix#
centroids_Sly <- Sly_TFs[Sly_TFs %in% rownames(inmat)]
```

## Network deconvolution
```R
#As indicated in the corto library#
regulon_Sly <- corto(inmat = inmat, 
                     centroids = centroids_Sly, 
                     nbootstraps = 100, 
                     p = 1e-7, 
                     nthreads=8)
```


## Identify phenotypes in the matrix

Create an object indicating the phenotype of each sample
```R
pheno <- read.table("Tomato_pheno.csv", header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
```
Reorder the matrix based on the phenotypes
```R
inmat<-inmat[,order(rownames(pheno))] 
```

## Master Regulator Analysis
One for each tretament comparsisson#
```R
#As indicated in corto librery#
#Control vs Severe#
mrs_CvsS <- mra(expmat1 = inmat[, pheno$group=="S"], 
              expmat2 = inmat[, pheno$group=="C"], 
              regulon = regulon_Sly, 
              minsize = 20, 
              nthreads = 6, 
              verbose = TRUE, 
              atacseq = NULL)
```
Obtain corto image as indicated in corto library
```R
pdf(file = "cortoMRS_CvsS.pdf", width = 15, height = 18)
mraplot(mrs_CvsS, mrs = 10)
dev.off()
```
```R
#Control vs Moderate#
mrs_CvsM <- mra(expmat1 = inmat[, pheno$group=="M"], 
              expmat2 = inmat[, pheno$group=="C"], 
              regulon = regulon_Sly, 
              minsize = 20, 
              nthreads = 6, 
              verbose = TRUE, 
```
Obtain corto image as indicated in corto library
```R
pdf(file = "cortoMRS_CvsM.pdf", width = 15, height = 18)
mraplot(mrs_CvsM, mrs = 10)
dev.off()
```
```R
#Severe vs Moderate#
mrs_SvsM <- mra(expmat1 = inmat[, pheno$group=="S"], 
                expmat2 = inmat[, pheno$group=="M"], 
                regulon = regulon_Sly, 
                minsize = 20, 
                nthreads = 6, 
                verbose = TRUE, 
                atacseq = NULL)
```
Obtain corto image as indicated in corto library
```R
pdf(file = "cortoMRS_SvsM.pdf", width = 15, height = 18)
mraplot(mrs_SvsM, mrs = 10)
dev.off()
```
## Get network as a text file
```R
#Function to obtain all regulons automatically#
regulon2sif_Sly <-function(regulon_Sly){
  sif <- data.frame(TF = c(), target = c(), interaction = c())
  for (name in names(regulon_Sly)) {
    source = rep(name, length(regulon_Sly[[name]]$tfmod));
    psif    = cbind(source, names(regulon_Sly[[name]]$tfmod), regulon_Sly[[name]]$tfmod);
    sif <- rbind(sif, psif)}
  colnames(sif) <- c("TF", "target", "int")
  rownames(sif) <- NULL
  return(sif)
}
```
Usage of the previously created function
```R
net_corto<-regulon2sif_Sly(regulon_Sly)
#Create a file to analyze it#
write.table (net_corto, file= "net_corto.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
```
## Obtain top MRA results
One for each treatment comparisson
```R
#Control vs Severe
get_mra <- function(mrs_CvsS, n = 334){names(head(mrs_CvsS$nes[order(mrs_CvsS$pvalue)], n))}
results_CvsS <- get_mra (mrs_CvsS, 334)
#Create a file to analyze it#
write.table (results_CvsS, file= "results C vs S.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

#Control vs Moderate
get_mra2 <- function(mrs_CvsM, n = 334){names(head(mrs_CvsM$nes[order(mrs_CvsM$pvalue)], n))}
results_CvsM <- get_mra2 (mrs_CvsM, 334)
#Create a file to analyze it#
write.table (results_CvsM, file= "results C vs M.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

#Moderate vs Severe
get_mra3 <- function(mrs_SvsM, n = 334){names(head(mrs_SvsM$nes[order(mrs_SvsM$pvalue)], n))}
results_SvsM <- get_mra3 (mrs_SvsM, 334)
#Create a file to analyze it#
write.table (results_SvsM, file= "results S vs M.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
```
## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
