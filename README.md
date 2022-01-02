# Network-analysis-of-root-and-leaf-transcriptome-integration
A transcriptomic data analysis using "Corto: a tool for network deconvolution based on a correlation inference  (Giorgi et al., 2020). The main objective is to find key transcription factors, known as master regulators, responsible for complex genetic diseases' symptoms. The results obtained help elucidate the function of specific genes in symptom progression.

For this analysis, we used Solanum lycopersicum (tomato) infected by two different Potato Spindle Tubercule Viroid (PSTVd) strains as a study model.

## Obtaining data from Gene Expression Omnibus
Microarray expression measurements were obtained from Gene Expression Omnibus (GEO) using the "GEOquery" R library. Then, this information was processed (filtrated and normalized) using the "affy" R library (RMA method). The microarray code was then translated to gene ID, collapsing a correspondence between microarray ID and gene ID. As a result, we obtained an expression matrix with samples as columns and gene names as rows.

## Assigning centroids for network deconvolution
For this analysis, we used a list of every TF in tomato known to the date of publication obtained from PlantTFDB to infer their regulons and generate a gene regulation network (GRN). Only the TFs contained in the expression matrix was used.

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
To identify each sample by phenotype, we created an object indicating which phenotype each sample corresponds to.

## Master Regulator Analysis
To identify de master regulators in the transition from one phenotype to another, we implemented the "MRA" corto function, using as inputs the GRN and the expression matrix described above, indicating which sample corresponded to each phenotype. We maintained "MRA" Corto function default values. We made one mra for each phenotype comparison available in the analyzed data.

#One for each tretament comparisson #
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
For network visualizations, we utilized Cytoscape, generating a visualization for each TF family of interest in which at least one member turned out to be a master regulator. In each visualization, the master regulators and their regulon were highlighted.
Example image:

![image](https://user-images.githubusercontent.com/94479457/142077113-50ec6184-d733-4268-9dc8-762be9215f89.png)


## Finding communities in the network
Aditionally to the master regulator analysis, an algorithm was applied to find communities in the infered GRN. The code used to generate the communities is "Community detection. ipynb," which can be found on the directory "Programs." 
