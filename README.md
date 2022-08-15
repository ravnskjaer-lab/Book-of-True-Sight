# LG/KR bulk RNA-seq data - Book of True Sight
#### This repository contains code used for creating the LG/KR bulk RNA-seq shiny app. The app is connected to a directory named "data". Data contains shared bulk RNA-seq data stored as DESeq2 objects.
#### To update the data folder with new datasets, which are applicable to the app.R file follow the below instructions


To be accepted for upload to the data folder on Ucloud, the data **must** comply with the following requirements:
  
- [x] The data must be in a DESeq2 object.

- [x] The DESeq2 object must include owner information using: mcols(DESeq2object)$owner <- "OWNER INITIALS".

- [x] The DESeq2 object must have removed sizeFactor information using: DESeq2object$sizeFactor <- NULL 

- [x] If model was adjusted during intial differential expression analysis, remove replaceable variable using: DESeq2object$replaceable <- NULL

- [x] Update the dataset database https://syddanskuni.sharepoint.com.mcas.ms/:x:/r/Sites/Hepatic_fanatics/_layouts/15/Doc.aspx?sourcedoc=%7B9C879D54-4F98-4CC4-A090-50A64DB5B9CD%7D&file=LG.KR_datasets.xlsx&action=default&mobileredirect=true&cid=e05664f4-1030-40ef-bc81-ee665b7f8481

### Running the app on Ucloud: ###
Open Ucloud and start a Rstudio application: https://cloud.sdu.dk/app/jobs/create?app=rstudio&version=4.2.0
Select an appropriate machine type:
- 23 GB RAM is a little better than your own computer 
- 47 GB ram is fast for this application

Select following folders:
- Ravnskjaer -> Software -> R
- Ravnskjaer -> Software -> Shiny

#### Press submit

- Press ppen interface 
- Navigate to the app located: work/R/bookoftruesight

## Enjoy your data surf :) 

<img src="https://images-wixmp-ed30a86b8c4ca887773594c2.wixmp.com/f/da7c68b7-9fd3-4b29-965e-5256b48bab90/d8mr76n-86c50354-e622-4c36-bebd-a87ee3dcc5b3.jpg?token=eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJzdWIiOiJ1cm46YXBwOjdlMGQxODg5ODIyNjQzNzNhNWYwZDQxNWVhMGQyNmUwIiwiaXNzIjoidXJuOmFwcDo3ZTBkMTg4OTgyMjY0MzczYTVmMGQ0MTVlYTBkMjZlMCIsIm9iaiI6W1t7InBhdGgiOiJcL2ZcL2RhN2M2OGI3LTlmZDMtNGIyOS05NjVlLTUyNTZiNDhiYWI5MFwvZDhtcjc2bi04NmM1MDM1NC1lNjIyLTRjMzYtYmViZC1hODdlZTNkY2M1YjMuanBnIn1dXSwiYXVkIjpbInVybjpzZXJ2aWNlOmZpbGUuZG93bmxvYWQiXX0.BVbLm7iPWrEK5ZkSY3rmeVh8dPu_sqz44xUMa6omk0M"/>
