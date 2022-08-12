# Shiny
#### This repository contains code used for creating the LG/KR bulk RNA-seq shiny app. The app is connected to a directory named "data". Data contains shared bulk RNA-seq data stored as DESeq2 objects.
#### To update the data folder with new datasets, which are applicable to the app.R file follow the below instructions


To be accepted for upload to the data folder on Ucloud, the data **must** comply with the following requirements:
  
- [x] The data must be in a DESeq2 object.

- [x] The DESeq2 object must include owner information using: mcols(DESeq2object)$owner <- "OWNER INITIALS".

- [x] The DESeq2 object must have removed sizeFactor information using: DESeq2object$sizeFactor <- NULL 

- [x] If model was adjusted during intial differential expression analysis, remove Replacable variable using: DESeq2object$Replacable <- NULL


"[](https://tenor.com/view/spiderman-responsibility-gif-4589950)
