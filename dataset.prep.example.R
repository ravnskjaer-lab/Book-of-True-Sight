#Example of how to prepare a dataset for the LG/KR bulk RNA-seq shiny app

#------------------------------------
#Load in the data in question 
LA.memory <- read.csv("~/PATH_TO_RAW_COUNTS/Female_pilot_exon_counts_noadj_new.txt", sep = "\t", row.names = 1)
colnames(LA.memory) <- substr(colnames(LA.memory),2,5)  #Cleaning up coloumn names
rownames(LA.memory) <- sub("\\|.*","", LA.memory$nnot) #Raw counts were from HOMER pipeline --> Using the column with gene symbols
LA.memory <- LA.memory[,-c(1,2,3,4,5,6,7)] #Removing unwanted information and retaining only gene counts and samples
LA.memory <- LA.memory[,order(colnames(LA.memory))] #Ordering columns based on seq id number 

#------------------------------------------------------------------------------------------------------------------------
#If required, change ENSEMBL gene ids to gene symbols
biomart.db <- read.csv("PATH_TO_BIOMART_DB/murine.biomart.db.csv", row.names = 1)

biomart.murine.memory <- biomart.db[biomart.db$ensembl_gene_id %in% rownames(LA.memory),] #Subsetting biomart human db
biomart.murine.memory <- biomart.murine.memory[-which(duplicated(biomart.murine.memory$external_gene_name)=="TRUE"), ] #Removing duplicated gene symbols

LA.memory <- LA.memory[rownames(LA.memory) %in% biomart.murine.memory$ensembl_gene_id,] #Subsetting the count matrix

biomart.murine.memory <- biomart.murine.memory[match(rownames(LA.memory), biomart.murine.memory$ensembl_gene_id),] #Correcting the order of biomart db to be that of count matrix

which((rownames(LA.memory) ==  biomart.murine.memory$ensembl_gene_id)=="FALSE") #Checking order of rows are the same in the two dfs 

rownames(LA.memory) <-  biomart.murine.memory$external_gene_name #Inserting gene symbols as rownames

LA.memory <- lx2[rownames(LA.memory) %in% biomart.murine.memory[which(biomart.murine.memory$gene_biotype == "protein_coding"),"external_gene_name"],] #OPTIONAL: Filtering based on protein_coding genes
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------
#Load in and/or prepare metadata
LA.memory.meta <- data.frame(factor(c("WD_CCl4",
                                "WD_CCl4",
                                "WD_CCl4",
                                "WD_CCl4",
                                "WD_CCl4",
                                "Chow",
                                "Chow",
                                "Chow",
                                "Chow",
                                "Chow",
                                "4wREV",
                                "4wREV",
                                "4wREV",
                                "4wREV",
                                "4wREV"), levels = c("Chow","WD_CCl4","4wREV")),row.names = colnames(LA.memory))
colnames(LA.memory.meta) <- "Treatment"

#Do the differential expression 
LA.memory.de <- DESeq(
  DESeqDataSetFromMatrix(countData = round(LA.memory),
                         colData = LA.memory.meta,
                         design = ~Treatment)
)

#Insert ownership information into the DESeq2 object 
mcols(LA.memory.de)$owner <- "LA"

#Remove sizeFactor column and if refitting was computed the replacable column
LA.memory.de$sizeFactor <- NULL

#Save the DESeq2 object and place it on the data folder on github 
saveRDS(LA.memory.de, file = "~/PATH_TO_LOCATION_OF_OBJECT/WD.4week.reversion.rds")
