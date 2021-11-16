# Network-analysis-of-root-and-leaf-transcriptome-integration
This is transcriptomic data analysis by network deconvolution based on correlation inference tool (corto, Giorgi et al., 2020) to find key transcription factors, known as master regulators, responsibles of the symthomps of complex genetic diseases, using the results obtained to elucidate the function of certain genes in the symthomps progression.

For this analysis we used Solanum lycopersicum (tomato) infected by two different strains of the Potato Spindle Tubercule Viroid (PSTVd) as study model.

## Obtaining data from Gene Expression Omnibus
Microarray expression mesurments where obtained from Gene Expresion Omnibus (GEO) using the "GEOquery" R library, then, these information was processed (flitrated and normalized) using the "affy" R library (RMA method). The microarray code was then translated to gene ID colapsing a list of corresponce between microarray ID and gene ID. As a result we obtained an expression matrix with samples as columns and gene names as rows.

## Assigning centroids for network deconvolution
For this analysis, we used a list of every TF in tomato known to the date of publication obtained from PlantTFDB, in order to infere their regulons and generate a gene regulation network (GRN). Only the TFs contained in the expression matrix were used.

## Network deconvolution
To obtain a GRN, we employed "corto" R library, using as input the expression matrix and the list of TF descrived above. We mantained "corto" corto function default values.
```R
#As indicated in the corto library#
regulon_Sly <- corto(inmat = inmat, 
                     centroids = centroids_Sly, 
                     nbootstraps = 100, 
                     p = 1e-7, 
                     nthreads=8)
```
## Identify phenotypes in the matrix
For the identificartion of each sample by phenotype, we created an object indicating which phenotype each sample correspond to.

## Master Regulator Analysis
In order to identify de master regulators in the transition from one phenotype to another, we implemented "mra" corto function, using as inputs the GRN and the expression matrix descrived above, indicating which sample corresponded to each phenotype. We mantained "mra" corto function default values. We made one mra for each phenotype comparisson avalailable in the analized data.

#One for each tretament comparsisson#
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
Example image:

![image](https://user-images.githubusercontent.com/94479457/142072824-7523a3d4-a863-4688-8eeb-79134b6503aa.png)

## Networks visualizations
For network visualizations we utilized Cytoscape, generating a visualization for each TF family of interest in which at least one member turned out to be a master regulator. In each visualization the master regulators and their regulon were highlighted.
Example image:

![image](https://user-images.githubusercontent.com/94479457/142077113-50ec6184-d733-4268-9dc8-762be9215f89.png)

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## Finding comunnities in the network
Aditionally to the master regulator analysis, an algorithm to find communities in the infered GRN was aplicated. 
