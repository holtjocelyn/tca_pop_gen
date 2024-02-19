#Code following Gunwald lab through Github: https://grunwaldlab.github.io/Population_Genetics_in_R/First_Steps.html

#make sure the following packages are installed. https://grunwaldlab.github.io/Population_Genetics_in_R/Getting_ready_to_use_R.html
#install.packages("poppr")
#install.packages("adegenet")
#install.packages("ade4")
#install.packages("ape")
#install.packages("ggplot2")
#install.packages("mmod")
#install.packages("magrittr")
#install.packages("dplyr")
#install.packages("treemap")
#install.packages("pegas")


#for access to R cheat sheets go here: https://www.rstudio.com/resources/cheatsheets/
  
#load the libraries to be able to analyze data. note that if you load other libraries they may overide content from other packages
library("poppr") #will override the pegas package, need to run pegas separtely for refined anlalysis. automatically loads adegenet and ade4 package 
library("adegenet") #requires ade4 package
library("ade4")
library("ape")
library("ggplot2")
library("mmod")
library("dplyr")
library("treemap")

#library("magrittr") #allows stringing together commands left to right using pipes (%>%)

#Get Session Info 
sessionInfo()

#avoid loading these libraries simultaneously or they will override content
#library("pegas") #loading these other libraries simultaneously will mask features from ape and ade4
#library("ade4") #loading this package will mask 
 
#get more info on other data file formats that poppr takes
??read.genalex


#set your working directory
setwd("/Users/Jocelyn/Desktop/tca_data/poppr_all_tca_26sept2023")


####load info of packages and versions used for workflow
#https://popgen.nescent.org/StartSNP.html
options(width = 100)
devtools::session_info()


########LOAD FILE AND ASSIGN STRATA########## https://grunwaldlab.github.io/Population_Genetics_in_R/Population_Strata.html
#upload the desired file
#example of reading in a structure file and automatically converting it to a genind object
#read.structure(file, n.ind = NULL, n.loc = NULL, onerowperind = NULL, col.lab = NULL, col.pop = NULL, col.others = NULL, row.marknames = NULL, NA.char = "-9", pop = NULL, sep = NULL, ask = TRUE, quiet = FALSE)

#duplicate this file and save formatted as a .tsv. then open it in numbers or excel to view the columns; close file without saving. samples are in column 1, and location info is in column 2
#col.lab: an integer giving the index of the column containing labels of genotypes. '0' if absent.
#col.pop: an integer giving the index of the column containing population to which genotypes belong. '0' if absent.
#if get a dimnames error: (1) go back into .stru file with a text edit and make sure there is not an extra space before the first marker name. (2) make sure that the heading of the file output is deleted, which contains the date and timestamp of the structure output. (3) delete the last blank row of the file
#if the file cannot open or be found, double check in the file properties that the extension is .stru rather than .stru.tsv. manually edit to .stru if this is the case

#make sure have the correct number of loci, this is the number with lesser columns from the structure output
#if there is an error with opening or finding the file, double check that the file extension is correct and not hidden
#use the complete file path
obj <- read.structure("/Users/Jocelyn/Desktop/tca_data/poppr_all_tca_26sept2023/structure_format_sa_us_22sept2023.stru", n.ind = 96, n.loc = 4749, onerowperind = FALSE, col.lab = 0, col.pop = 2, col.others = NULL, row.marknames = 1, NA.char = "-9", pop = NULL, sep = "\t", ask = TRUE, quiet = TRUE)


#which other optional columns should be read? press press enter


#Warning message: In df2genind(X = X, pop = pop, ploidy = 2, sep = sep, ncode = ncode) : duplicate labels detected for some individuals; using generic labels


#call the object of the formatted genind object from the structure file
obj


# make a file with sample names and population names. OR use the same population map from STACKS and STRUCTURE
#Make sure that there are no extra spaces between collection or location names

#upload this .txt file to R called: tca_pop_map_sa_us_22sept2023.txt
#Make sure that the pop map has a header with sample for column 1 and location for column 2
strata.txt <- file.choose()

#now convert this text file into a table that the program can read, and mark that headers are true
strata <- read.table(file = strata.txt, header = TRUE, fill = TRUE)

#view the new table that was created from the text file which is now called strata
View(strata)

#edit the strata in your gen file to now have that information 
#if get Error: Number of rows in data frame are not equal to the number of individuals in object, either add a header, remove a header, or remove the last trailing line
strata(obj) <- strata

#check to see that this information was incorporated into the gen file
View(strata(obj))


#####################
#######
#determine allele frequencies and any missing data, make a locus table. https://grunwaldlab.github.io/Population_Genetics_in_R/Locus_Stats.html
#data("obj")
objlt <- locus_table(obj)

#show the number of alleles at each locus and eveness
objlt

#now plot the missing data. ******Takes ~ 5 minutes to run
#if have loci shared among all individuals, then missing data will be zero
info_table(obj, plot = TRUE) # type = missing



#####Calculate and Plot Genotypic Diversity############
#plot genotypic diversity.*****Takes ~ 2-5 minutes to run. https://grunwaldlab.github.io/Population_Genetics_in_R/Genotypic_EvenRichDiv.html
poppr(obj)

#calculate genotype diversity or MLG for each population. https://grunwaldlab.github.io/Population_Genetics_in_R/Genotypic_EvenRichDiv.html
p.tab <- mlg.table(obj)

#call genotypic diversity among populations
p.tab

#calculate diversity statistics on a matrix with multilocus genotypes per pop
diversity_stats(p.tab, H = TRUE, G = TRUE, lambda = TRUE, E5 = TRUE)


#make a pca of the data, quick visualization
objpca <- indpca(obj)

#now plot the graphic of the pca
plot(objpca, cex = 0.7)



###########Find potentially uninformative loci
nLoc(obj)

iobj <- informloci(obj)

nLoc(iobj)




######CALCULATE GENETIC DISTANCE####### https://grunwaldlab.github.io/poppr/reference/genetic_distance.html
#GST =0 indicates no differentiation, whereas GST =1 indicates that populations are segregating for differing alleles
objdist <- provesti.dist(obj)

#call the genetic distance estimates
objdist





##########Do a DAPC Cross Validation to detect how many PCs there should be in the DAPC###########
#https://grunwaldlab.github.io/Population_Genetics_in_R/DAPC.html

set.seed(999)
objx <- xvalDapc(tab(obj, NA.method = "mean"), pop(obj))
#appears to level off around 10 to 15

#now from the range of PCs do iterations to see which is the best value
set.seed(999)
system.time(objx <- xvalDapc(tab(obj, NA.method = "mean"), pop(obj), n.pca = 10:15, n.rep = 1000, parallel = "multicore", ncpus = 3L))

#give cross validation results
names(objx)
#output the mean success of each PC and the one with the highest success
objx[-1]

#10 or 13 PCAs is the optimal number




########RUN DAPC ON DATA#####https://grunwaldlab.github.io/Population_Genetics_in_R/DAPC.html
#for scattered DAPC use this code: https://grunwaldlab.github.io/Population_Genetics_in_R/gbs_analysis.html
#https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/inst/doc/adegenet-dapc.pdf?revision=887&root=adegenet&pathrev=887


#use data: obj

#load library
library("adegenet")

#how to create a dapc in adegenet
#https://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf
#examples of DAPC graphs: https://grunwaldlab.github.io/Population_Genetics_in_R/DAPC.html
#adegenet documentation: https://www.rdocumentation.org/packages/adegenet/versions/1.2-7/topics/dapc
#further documentation on scatter plot specs: https://www.rdocumentation.org/packages/ade4/versions/1.7-19/topics/add.scatter



#load the data
###read structure data as a genind object
#objgenind <- import2genind("/path/file_name.stru", n.ind = 96, n.loc = 4749, onerowperind = FALSE, col.lab = 0, col.pop = 2, col.others = NULL, row.marknames = 1, NA.char = "-9", pop = NULL, sep = "\t", ask = TRUE, quiet = TRUE)

objgenind <- import2genind("/Users/Jocelyn/Desktop/tca_data/poppr_all_tca_26sept2023/structure_format_sa_us_dapc_labels_22sept2023.stru", n.ind = 96, n.loc = 4749, onerowperind = FALSE, col.lab = 0, col.pop = 2, col.others = NULL, row.marknames = 1, NA.char = "-9", pop = NULL, sep = "\t", ask = TRUE, quiet = TRUE)
#press enter


#import the mapping file to go with this dataset
#upload this .txt file to R called: tca_pop_map_sa_us_dapc_labels_22sept2023.txt
strata.txt <- file.choose()

#now convert this text file into a table that the program can read
strata <- read.table(file = strata.txt, header = TRUE, fill = TRUE)

#view the new table that was created from the text file which is now called strata
View(strata)

#edit the strata in your gen file to now have that information
strata(objgenind) <- strata

#check to see that this information was incorporated into the gen file
View(strata(objgenind))

#set the population to state or location
pop(objgenind) <- strata$location

#check to make sure this worked by calling this argument
pop(objgenind)

#set to the PCs from the dapc cross validation
#this should be equal to or less than the actual pops have. if go too high in number than the discriminant analysis eigenvalues will not be as clear and the data points will not be visible on the scatter plot
#dapc.objgenind <- dapc(objgenind, var.contrib = TRUE, scale = FALSE, n.pca = 10, n.da = NULL)


dapc.objgenind <- dapc(objgenind, var.contrib = TRUE, scale = FALSE, n.pca = NULL)
#this should reflect the pcas identified from cross validation

#dapc.objgenind <- dapc(objgenind, pop = strata$location, var.contrib = TRUE, scale = FALSE, n.pca = 96, n.da = NULL) #n.da means how many discriminant analysis components, if null you will manually type this in and press enter

80#choose the number of PCs to retain

7 #choose the number of discriminant functions to retain based upon the bar plots in the discriminant analysis eigenvalues - this represents the number of pops 

#?pch to find out more about the different symbols associated with the pch number code


#input a color palate to use instead of the default colors, which could overlap and be difficult to tell apart
#myCol <- c("darkblue", "purple", "green", "orange", "brown", "black", "grey", "red", "turquoise")
myCol <- c("darkblue", "purple", "green", "orange", "brown", "black", "grey", "red", "turquoise", "yellowgreen", "pink", "#000033", "#CCFFFF", "#999900", "#660066", "#CC3300", "#CC9999", "#00FF00", "#990066")


#Sample bubbles that are semi-transparent, need to make them a bit smaller or more spread out and change the color palate
#scatter(dapc.objgenind, bg = "white", scree.da = FALSE, pch = 20, solid = 0.5, cell = 0, cstar = 0, cex = 1.7,  clab = 0, col = myCol, legend = TRUE, posi.leg = "bottomleft", cleg = 0.85, xax = 1, yax = 2) #clab = 0 is what gets rid of the sample labels

scatter(dapc.objgenind, bg = "white", scree.da = FALSE, pch = 18:20, solid = 0.5, cell = 0, cstar = 0, cex = 1.7,  clab = 0, col = myCol, legend = TRUE, posi.leg = "left", cleg = 0.85, xax = 1, yax = 2) #clab = 0 is what gets rid of the sample labels
#to have the pops different colors and shapes do 17:22   #17 are triangles, 20 are circles
#to label the pops input this after legend = TRUE: txt.leg=lgn, 

#xlim=c(0,10), ylim=c(0,5)


###CALCULATE Fst and Fis VALUES######## https://search.r-project.org/CRAN/refmans/genepop/html/genedivFis.html
## gene diversity and Fis

#load libraries
library("genepop")
#library("pegas")


#set your working directory
setwd("/Users/Jocelyn/Desktop/tca_data/poppr_all_tca_26sept2023")

#instruction link to genepop program
#https://cran.r-project.org/web/packages/genepop/genepop.pdf
 


##Run Fis values for all collections
#??genedivFis {genepop}
objFis <- genedivFis("snps_in_genepop_format_sa_us_22sept2023.txt", sizes = FALSE, outputFile = "tca_genepop_fis_output.txt", dataType = "Diploid", verbose = interactive())

objFis

## run Fst analysis using genepop file from STACKS (can also create input file using PGDSpider)
#??Fst{genepop}

## run a global Fst analysis using genepop - one runs comparisons for all populations
tca_snps_genepop_global_Fst <- Fst("snps_in_genepop_format_sa_us_22sept2023.txt", sizes = FALSE, pairs = FALSE, outputFile = "tca_fst_global_output.txt", dataType = "Diploid", verbose = interactive())

## run a pairwise Fst for all populations using genepop analysis
tca_snps_genepop_pairwise_Fst <- Fst("snps_in_genepop_format_sa_us_22sept2023.txt", sizes = FALSE, pairs = TRUE, outputFile = "tca_fst_pairwise_output.txt", dataType = "Diploid", verbose = interactive())


###test of genic and genotypic differentiation
#test_diff(
#inputFile,
#genic = TRUE,
#pairs = FALSE,
#outputFile = "",
#settingsFile = "",
#dememorization = 10000,
#batches = 100,
#iterations = 5000,
#verbose = interactive()
#)

###test of genic and genotypic differentiation; will take ~15 mins or less to run
test_diff("snps_in_genepop_format_sa_us_22sept2023.txt", genic = TRUE, pairs = FALSE, outputFile = "tca_genotypic_differentiation_test.txt", settingsFile = "", dememorization = 10000, batches = 100, iterations = 5000, verbose = interactive())




#******Test with current parameters takes ~3 hours to run **********
## run a pairwise differentiation test using genepop
#test parameters to run dememorization 1,000, batches = 100, iterations = 5,000
#then set to 10,000 demorization, 100 batches, and 5,000 iterations
#??differentiation{genepop}
tca_snps_genepop_test_diff <- test_diff("snps_in_genepop_format_sa_us_22sept2023.txt", genic = TRUE, pairs = TRUE, outputFile = "tca_snps_final_differentiation_test.txt", dememorization = 10000, batches = 100, iterations = 5000, verbose = interactive())



####RUN AMOVA######## https://grunwaldlab.github.io/poppr/reference/poppr.amova.html
#https://grunwaldlab.github.io/poppr/reference/poppr-package.html
#example dataset: data("Aeut")
#The implementation of AMOVA in poppr requires two very basic components: (1) A distance matrix derived from the data and (2) a separate table used to partition the data into different stratifications.
#https://grunwaldlab.github.io/poppr/reference/poppr-package.html

#load r packages
library("adegenet")
library("poppr")
#library("ade4")
#library("pegas")
#library("vcfR")

#set working directory
setwd("~/Desktop/tca_data/poppr_all_tca_26sept2023")

#make sure to leave the default heading in the file, otherwise it won't load properly
obj2 <- read.genepop("snps_in_genepop_format_sa_us_22sept2023.gen", ncode = 2L, quiet = FALSE)

#call the object of the formatted object. if structure of genepop it is a genind object. if vcf it is class vcfR. this will also give percent missing data. ~17.58%
obj2

#upload this .txt file to R called: 
#tca_pop_map_sa_us_22sept2023.txt
strata.txt <- file.choose()
#select the updated pop map

#now convert this text file into a table that the program can read
strata <- read.table(file = strata.txt, header = TRUE, fill = TRUE)

#view the new table that was created from the text file which is now called strata
View(strata)

#edit the strata in your gen file to now have that information
strata(obj2) <- strata

#check to see that this information was incorporated into the gen file
View(strata(obj2))



#for interpreting AMOVA data: https://grunwaldlab.github.io/Population_Genetics_in_R/AMOVA.html

#for the simplified equation below for calculating amova: https://grunwaldlab.github.io/poppr/reference/poppr.amova.html

#load library
library("poppr")
library("adegenet")

#for analysis using pegas, run again separately for ade4
pamova <- poppr.amova(
  obj2,
  hier = ~location,
  clonecorrect = FALSE,
  within = TRUE,
  dist = NULL,
  squared = TRUE,
  freq = TRUE,
  correction = "quasieuclid",
  sep = "",
  filter = FALSE,
  threshold = 0,
  algorithm = "farthest_neighbor",
  threads = 1L,
  missing = "loci",
  cutoff = 0.05,
  quiet = FALSE,
  method = c("pegas"),
  nperm = 999
)
#make sure to run the model with pegas soley so dont get potential issues with ade4, can also try to run with ade4 by itself




#to run amova with default settings in poppr
#https://groups.google.com/g/poppr/c/NSag-55d6bs
#poppr.amova(obj2, ~location)


###NOW GET THE RESULTS FROM THE ANOVA ANALYSIS######
pamova


#interpreting amova results
#https://grunwaldlab.github.io/Population_Genetics_in_R/AMOVA.html




#for analysis, run separately with ade4
pamova <- poppr.amova(
  obj2,
  hier = ~location,
  clonecorrect = FALSE,
  within = TRUE,
  dist = NULL,
  squared = TRUE,
  freq = TRUE,
  correction = "quasieuclid",
  sep = "",
  filter = FALSE,
  threshold = 0,
  algorithm = "farthest_neighbor",
  threads = 1L,
  missing = "loci",
  cutoff = 0.05,
  quiet = FALSE,
  method = c("ade4"),
  nperm = 999
)
#make sure to run the model with pegas soley so dont get potential issues with ade4, can also try to run with ade4 by itself

###NOW GET THE RESULTS FROM THE ANOVA ANALYSIS######
pamova



#load package ade4
library("ade4")

library("pegas")

set.seed(1999)
objsig <- randtest(pamova, nrepet = 999) #can do an initial trial of lower nrpeat, then set to 999

objsig

plot(objsig)

write.table(pamova$componentsofcovariance, sep = ",", file = "/Users/Jocelyn/Desktop/tca_data/poppr_all_tca_26setp2023/tca_pamova_ade4_26Sept2023.csv")






#Get Session Info 
sessionInfo()
