# GWAS-assignment  
Here are scripts for the Genome-Wide Association Studies (GWAS) assignment of the *BVG-7003* course. Two software will be used : *rMVP* and *GAPIT*.  

## Installing R, Rstudio and Rtools  
> R, Rstudio and Rtools can be downloaded from : https://posit.co/download/rstudio-desktop/.

## rMVP R scripts  
> Scripts for this section are in *mvp.R*.  
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
> `setwd("Path/to/rMVPdata/for/GWAS/from/actual/directory")`  

### Import and adapt genotypic and phenotypic data
> Data from *geno.hmp.txt* and *pheno.txt* are imported and prepared using *MVP.data()* function.  
>   
> Six files are created:    
> -*MVP.Data.20221123_211412*  
> -*mvp.hmp.geno.bin*  
> -*mvp.hmp.geno*  
> -*mvp.hmp.geno.ind*  
> -*mvp.hmp.geno.map*  
> -*mvp.hmp.phe*  

### Kinship  
> Kinship can be done using MVP.Data.Kin() function.  
>   
> Two files will be generated:  
> -*mvp.kin.bin*  
> -*mvp.kin*  

### PCA
> Principal component analysis can be generated using *MVP.Data.PC()* function.  
>   
> Two files will be generated:  
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
> Different parameters can my added, for exampl if covariates are used or if other methods than farmCPU are used.  
>   
> Nine new files are generated (example in output):  
> *MVP.20221123_213323* : Report messages and steps of the GWAS.  
> *protein.FarmCPU* and *protein.FarmCPU_signals* : containg SNP information and results from the analysis.  
> *protein.FarmCPU.Circular-Manhattan*, *protein.FarmCPU.Rectangular-Manhattan* : Manhattan plots from the results.  
> *protein.FarmCPU.QQplot* : a QQplot of the results.  
> *protein.FarmCPU.SNP-Density* : a graph of SNP densities on each chromosomes.  
> *protein.PCA_2D* : a 2D PCA.  
> *protein.Phe_Dist* : a graph of the phenotype distribution.  

## GAPIT R scripts 
> Scripts for this section are in *GAPIT.R*.  
> Examples of outputs from GAPIT section are in *output_gapit*.  

### Installation and loading requiered packages 
> Devtools, GAPIT3 and GAPIT functions need to be installed using:  
> ```
> install.packages("devtools")  
> devtools::install_github("jiabowang/GAPIT3",force=TRUE)  
> source("http://www.zzlab.net/GAPIT/GAPIT.library.R")  
> source("http://www.zzlab.net/GAPIT/gapit_functions.txt")  
> ```  

### Set working directory  
> Actual directory can be visualized using:  
> `getwd()`  
> And can be changed using:  
> `setwd("Path/to/GAPITdata/for/GWAS/from/actual/directory")`  

### Import and adapt genotypic and phenotypic data  
> Data from *geno.hmp.txt* and *pheno.txt* are imported using *read.table()*.  
> **Note**: geno.hmp.txt is not the same as the one used for rMVP since the latter scripts needs a # after the rs of the first line, which GAPIT doesn't need.  
 
### Analysing phenotype data  
> Phenotypic data can be analyzed using *str*, *hist*, *mean*, *range* and *sd* function and lnes with missing data can be counted using:   
`which(is.na(pheno$protein))`  


