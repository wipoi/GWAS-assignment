# GWAS-assignment  
Here are scripts for the Genome-Wide Association Studies (GWAS) assignment of the *BVG-7003* course. Two tools will be used : *rMVP* and *GAPIT*. Example have been made using the genotype data file *geno.hmp.txt* and the phenotype data file *pheno.txt* available in GWAS_data directory.  

## Installing R, Rstudio and Rtools  
> R, Rstudio and Rtools can be downloaded from : https://posit.co/download/rstudio-desktop/.

## rMVP R scripts  
> Scripts for this section are in *rMVP/mvp.R*.  
> Examples of outputs from rMVP section are in *output_rMVP*.  

### Installation and loading  
> *rMVP* package have to be installed and loaded using:  
> ```
> install.packages("rMVP")  
> library(rMVP)  
> ```  

### Set working directory  
> Actual directory can be visualized using:  
> `getwd()`  
> And can be changed using:  
> `setwd("Path/to/GWAS_data/")`  

### Import and adapt genotypic and phenotypic data
> For the example herein, data from *geno.hmp.txt* and *pheno.txt* are imported.  
> *pheno.txt* contains a sole variable named *protein* quantified for all taxa that are also in *geno.hmp.txt*.  
> Data are imported and prepared using *MVP.data()* function.   
> 
> **Note:** A *#* must be placed after *rs* on the fisrt line of *geno.hmp.txt* file to input in the good format for *rMVP*.  
>   
> Six files are created (example in output_rMVP):    
> -*MVP.Data.20221123_211412*  
> -*mvp.hmp.geno.bin*  
> -*mvp.hmp.geno*  
> -*mvp.hmp.geno.ind*  
> -*mvp.hmp.geno.map*  
> -*mvp.hmp.phe*  

### Kinship  
> Kinship can be done using MVP.Data.Kin() function.  
>   
> Two files will be generated (example in output_rMVP):  
> -*mvp.kin.bin*  
> -*mvp.kin*  

### PCA
> Principal component analysis can be generated using *MVP.Data.PC()* function.  
>   
> Two files will be generated (example in output_rMVP):  
> -*mvp.pc.bin*  
> -*mvp.pc*  

### Data input  
> In this section, data are loaded from files generated in the *Import and adapt genotypic and phenotypic data* section.  
>   
> Three data files are loaded:  
> -*genotype*  
> -*phenotype* which contains a column *Taxa* containing sample names and the phenotype column name *protein*  
> -*map* which contains columns *SNP* (SNP IDs), *CHROM* (chromosomes where SNPS are), *POS* (the position of the SNP), *REF* (the reference allele) and *ALT* (the alternative allele).  

### GWAS 
> In this section, GWAS is runt on each phenotype column using the MVP()* function in a *for* loop.  
>   
> The *FarmCPU (Fixed and random model Circulating Probability Unification)*  method is used, but *GLM* and *MLM* methods can also be used.
>     
> Different parameters can my added, for example if covariates are used or if other methods than farmCPU are used.  
>   
> Nine new files are generated (example in output_rMVP):  
> -*MVP.20221123_213323* : Report messages and steps of the GWAS.  
> -*protein.FarmCPU* and *protein.FarmCPU_signals* : containg SNP information and results from the analysis.  
> -*protein.FarmCPU.Circular-Manhattan*, *protein.FarmCPU.Rectangular-Manhattan* : Manhattan plots from the results.  
> -*protein.FarmCPU.QQplot* : a QQplot of the results.  
> -*protein.FarmCPU.SNP-Density* : a graph of SNP densities on each chromosomes.  
> -*protein.PCA_2D* : a 2D PCA.  
> -*protein.Phe_Dist* : a graph of the phenotype distribution.  

## GAPIT R scripts 
> Scripts for this section are in *GAPIT/GAPIT.R*.  
> Examples of outputs from GAPIT section are in *output_gapit*.  

### Installation and loading requiered packages 
> Devtools, GAPIT3 and GAPIT functions need to be installed and loaded using:  
> ```
> install.packages("devtools")  
> devtools::install_github("jiabowang/GAPIT3",force=TRUE)
> library(GAPIT3)
> source("http://zzlab.net/GAPIT/GAPIT.library.R")  
> source("http://zzlab.net/GAPIT/gapit_functions.txt")  
> ```  

### Set working directory  
> Actual directory can be visualized using:  
> `getwd()`  
> And can be changed using:  
> `setwd("Path/to/GWAS_data")`  

### Import genotypic and phenotypic data  
> Data from *geno.hmp.txt* and *pheno.txt* are imported using *read.table()*.  
> **Note**: geno.hmp.txt is not in the same format as the one used for rMVP since the latter scripts needs a # after the rs of the first line, which GAPIT doesn't need.  
 
### Analysing phenotype data  
> Phenotypic data can be analyzed using *str*, *hist*, *mean*, *range* and *sd* function and lines with missing data can be counted using:   
> `which(is.na(pheno$protein))`  
> Example of outputs for these metrics are given in the R script.  

### GWAS
> Here are GWAS analysis using GAPIT and five different models.   
> **Note**: running these sections might produce the error:
```
Error in `[<-`(`*tmp*`, i, 1, value = mean(pieceD, na.rm = T)) : subscript out of bounds
```  
> By retracing the origin of this error using:  
`options(error = recover)`  
> The error seems to came from line 194 of the function GAPIT.Genotype.View, which is:  
`loc[i,1]=mean(pieceD,na.rm=T)`  
> It seems it can't attribute anything to loc[i,1] because i is equal to 1 and loc do not have any lines.  
> Output files are then the ones generated before the error and the empty files named *GAPIT.Genotype.Density_R_sqaure* have not been given in output_gapit.  

### GWAS using *Mixed Linear Model* (*MLM*)   
> In this block, GWAS is performed without compression (*MLM* model) by setting a number of group equal to the number of individual within the population and by regrouping by 1 individual so that all are kept individual.   
> 
> The model produced eigth output files:  
> -*GAPIT.Genotype.Density_R_sqaure*:  graphs of correlations and relations between distance and markers (density).    
> -*GAPIT.Genotype.Kin_Zhang* and *GAPIT.Genotype.Kin_Zhang*: the *.csv* kinship data file and the corresponding kinship plot.   
> -*GAPIT.Genotype.PCA*, *GAPIT.Genotype.PCA_2D* and *GAPIT.Genotype.PCA_3D*: the *.csv* PCA data file and the corresponding PCA graphs in 2D (2 PCs)  and 3D (3PCs).  
> -*GAPIT.Genotype.PCA_eigenvalues* and *GAPIT.Genotype.PCA_eigenValue*: the *.csv* data file associated with the eigenValue of de principal components (PCs) of the PCA and the corresponding graph.  

### GWAS using *compressed MLM model* (*CMLM*)  
> In this block, GWAS is performed using compression (*CMLM* model) by grouping data by 10 as an example.  The grouping number can be changed according to the data, see *"Mixed linear model approach for genome-wide association studies" Zhang et al., Nature Genetics 2010* for more information. 
>  
> The model also produced eigth output files:  
> -*GAPIT.Genotype.Density_R_sqaure*:  graphs of correlations and relations between distance and markers (density).    
> -*GAPIT.Genotype.Kin_Zhang* and *GAPIT.Genotype.Kin_Zhang* which are the *.csv* kinship data file and the corresponding kinship plot.  
> -*GAPIT.Genotype.PCA*, *GAPIT.Genotype.PCA_2D* and *GAPIT.Genotype.PCA_3D* which are the *.csv* PCA data file and the corresponding PCA graphs in 2D (2 PCs)  and 3D (3PCs).  
> -*GAPIT.Genotype.PCA_eigenvalues* and *GAPIT.Genotype.PCA_eigenValue*, the *.csv* data file associated with the eigenValue of the principal components (PCs) of the PCA and the corresponding graph.  
> 
> The differences between compression or not is small herein.  

### GWAS using different kinship clustering methods
> In this block, GWAS is performed using different kinship clustering methods to group individuals according to their kinship. The number of PCs is optimized using the Bayesian information criteria (BIC) based on the trait in *pheno.txt*, which is none herein.
>  
> Again, this model produced eigth output files:  
> -*GAPIT.Genotype.Density_R_sqaure*:  graphs of correlations and relations between distance and markers (density).    
> -*GAPIT.Genotype.Kin_Zhang* and *GAPIT.Genotype.Kin_Zhang* which are the *.csv* kinship data file and the corresponding kinship plot.  
> -*GAPIT.Genotype.PCA*, *GAPIT.Genotype.PCA_2D* and *GAPIT.Genotype.PCA_3D* which are the *.csv* PCA data file and the corresponding PCA graphs in 2D (2 PCs)  and 3D (3PCs).  
> -*GAPIT.Genotype.PCA_eigenvalues* and *GAPIT.Genotype.PCA_eigenValue*, the *.csv* data file associated with the eigenValue of the principal components (PCs) of the PCA and the corresponding graph.  
> 
> The difference obtained herein is that FDR_Adjusted_P-values are smaller.  

### GWAS using *Multiple Locus Mixed Linear Model* (*MLMM*)
> In this block, GWAS is performed using *Multiple Locus Mixed Linear Model* (*MLMM* model) which takes as covariates associated markers.  
>  
> This model produced six output files:  
> -*GAPIT.Genotype.Density_R_sqaure*:  graphs of correlations and relations between distance and markers (density).    
> -*GAPIT.Genotype.PCA*, *GAPIT.Genotype.PCA_2D* and *GAPIT.Genotype.PCA_3D* which are the *.csv* PCA data file and the corresponding PCA graphs in 2D (2 PCs)  and 3D (3PCs).  
> -*GAPIT.Genotype.PCA_eigenvalues* and *GAPIT.Genotype.PCA_eigenValue*, the *.csv* data file associated with the eigenValue of the principal components (PCs) of the PCA and the corresponding graph.    

### GWAS using *FarmCPU* model
> In this block, GWAS is performed using *Fixed and random model Circulating Probability Unification* (*FarmCPU* model). This model aimed to correct for false-positives and cofounding between markers and cofactors. More informations are available in *"Iterative Usage of Fixed and Random Effect Models for
Powerful and Efficient Genome-Wide Association Studies", Liu et al., PLOS Genetics (2016)*. 
>  
> This model also generates six output files:  
> -*GAPIT.Genotype.Density_R_sqaure*: graphs of correlations and relations between distance and markers (density).     
> -*GAPIT.Genotype.PCA*, *GAPIT.Genotype.PCA_2D* and *GAPIT.Genotype.PCA_3D* which are the *.csv* PCA data file and the corresponding PCA graphs in 2D (2 PCs)  and 3D (3PCs).  
> -*GAPIT.Genotype.PCA_eigenvalues* and *GAPIT.Genotype.PCA_eigenValue*, the *.csv* data file associated with the eigenValue of the principal components (PCs) of the PCA and the corresponding graph.  

