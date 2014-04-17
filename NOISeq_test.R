source("http://www.bioconductor.org/biocLite.R")
biocLite("NOISeq")


??NOISeq
library('NOISeq')


# this script checks for differential expression in RNAseq dataset ( post HTSeq/Fluxcapacitor) without replicates
# it loads data in tab delimted format 
# factors must be defined in the myfactors function
# normalisation is carried out   normalisation is recommended, but optional, follow whichever path you think best for your data.
# noiseg()   the replicates'no' invokes noiseq-sim  to simulate technical replicates


# plotting data via plot function


# load samples ( combined HTSeq or Fluxcapacitor output)   e.g
#CBS1954_0h CBS1954_1h CBS1954_6h
#ENSG00000180008          0         22          0
#ENSG00000124209          0         38          0
#ENSG00000064703          0         46          0
#ENSG00000230873          0          0          0


mysamples <-  read.table("~/RNAseq_project//comparisons/genebased_human_NOISeq/CBS6318_genes.txt", 
                  sep="\t", 
                  #col.names=c("id", "name"), 
                  header = TRUE,
                  fill=FALSE, 
                  strip.white=TRUE,
                  row.names=1)

head(mysamples)


# define factors   treated/untreated   liver/ kidney    whatever the conditions are
myfactors = data.frame(Tissue=c("untreated","treated","treated"),
                       TissueRun = c("untreated_1","treated_1","treated_2"))

# combine the two information-sets to one expressiondata set
mydata <- readData(data=mysamples, factors=myfactors)
mydata



# normalisation via trimmed mean of M 
myTMM = tmm(assayData(mydata)$exprs, long = 1000, lc = 0)
myTMMdata <- readData(data=myTMM, factors=myfactors)
mydata
# low count filter is not neccessary

#noiseq-sim needs raw normalized data
myresults <- noiseq(mydata, 
                    factor = "Tissue", 
                    k = NULL, norm = "n", 
                    pnr = 0.2, nss = 5, 
                    v = 0.02, lc = 1, 
                    replicates = "no")

myresults_TMM <- noiseq(myTMMdata, 
                    factor = "Tissue", 
                    k = NULL, norm = "n", 
                    pnr = 0.2, nss = 5, 
                    v = 0.02, lc = 1, 
                    replicates = "no")

#The output myresults@results[[1]]$prob gives the estimated probability of differential expression for each
#feature.

myresults@results[[1]]$prob

DE.plot(myresults, q = 0.9, graphic = "expr", log.scale = TRUE)
DE.plot(myresults_TMM, q = 0.9, graphic = "expr", log.scale = TRUE)

