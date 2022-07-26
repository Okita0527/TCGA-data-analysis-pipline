##usage: Rscript download_count_clinic.R -i TCGA-PAAD -r unstranded -n fpkm_unstranded -t rawdata -o results



library(TCGAbiolinks)
library(SummarizedExperiment)
library(getopt)
library(httr)

spec = matrix(c(
  'help'   , 'h', 0, "logical",
  'input'  , 'i', 1, "character",
  'raw_type'  , 'r', 1, "character",
  'normalization_type'  , 'n', 1, "character",
  'rawDataDir'  , 't', 1, "character",
  'outputDir'  ,'o', 1, "character"), byrow=TRUE, ncol=4)
opt = getopt(spec)

cancer_type = opt$input            
raw_count_type = opt$raw_type
normalization_count_type = opt$normalization_type
raw_data_path = opt$rawDataDir
output_path = opt$outputDir


###########  download project directory of GDC
download_porjects <- function(raw_data_path){
    TCGA_projects<-data.frame(id = TCGAbiolinks:::getGDCprojects()$project_id, name = TCGAbiolinks:::getGDCprojects()$name)
    write.csv(TCGA_projects,file = paste(raw_data_path,"/TCGA_projects.csv",sep=""),quote = FALSE) 
    project_directory  = TCGA_projects[,1]
  return(project_directory)
}


###########  Harmonized data options
    ###########  STAR - Counts
download_star_count <- function(cancer_type,count_type){
    query_star_count <- GDCquery(
        project = cancer_type,
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification", 
        workflow.type = "STAR - Counts")
  
    GDCdownload(query_star_count,method = "api",directory = raw_data_path,files.per.chunk = 10)  
    star_count_prepare <- GDCprepare(query_star_count,directory = raw_data_path)
    count_data = switch(  
      count_type,  
      "unstranded"= assay(star_count_prepare,1),  
      "stranded_first"= assay(star_count_prepare,2),  
      "stranded_second"= assay(star_count_prepare,3),  
      "tpm_unstranded"= assay(star_count_prepare,4),
      "fpkm_unstranded"= assay(star_count_prepare,5),
      "fpkm_uq_unstranded"= assay(star_count_prepare,6)
    )      
  return(count_data)
}



###########  clinic INFO
download_clinic_data <- function(cancer_type){
    query_clinic <- GDCquery(
        project = cancer_type,
        data.category = "Clinical",
        data.type = "Clinical",
        data.format = "BCR Biotab")
    GDCdownload(query_clinic,method = "api",directory = raw_data_path,files.per.chunk = 10)                                                     #下载数据，默认为当前路径下GDCdata文件夹中
    clinical_prepare <- GDCprepare(query_clinic,directory = raw_data_path)
    clinical_data <- clinical_prepare$clinical_patient_hnsc
  
  return(clinical_data)
}



###################m    main  #####################
if (dir.exists(raw_data_path) == FALSE){
  dir.create(raw_data_path)
}
all_porject_ID <- download_porjects(raw_data_path)
if(cancer_type %in% all_porject_ID){
  if (dir.exists(output_path) == FALSE){
      dir.create(output_path)
  }

  #####raw_count输出为csv文件#######
  raw_count_data <- download_star_count(cancer_type,raw_count_type)
  rawcount_outputfilename = paste(output_path,"/",cancer_type,".raw.count.csv",sep = "")
  write.csv(raw_count_data,file = rawcount_outputfilename,quote = FALSE) 
  print("RNA-seq raw count has been downloaded")
  
  Sys.sleep(10)
  #####normalization_count输出为csv文件#######
  normalization_count_data<-download_star_count(cancer_type,normalization_count_type)
  normalizationcount_outputfilename = paste(output_path,"/",cancer_type,".normalization.count.csv",sep = "")
  write.csv(normalization_count_data,file = normalizationcount_outputfilename,quote = FALSE)  
  print("RNA-seq normalization count has been downloaded")
  
  Sys.sleep(10)
  #####TCGA clinic 输出为XLS文件#######    
  #clinical_data <- download_clinic_data(cancer_type)
  clinical_data <- GDCquery_clinic(cancer_type, type = "clinical",save.csv = FALSE )
  clinical_outputfilename = paste(output_path,"/",cancer_type,"_clinic.csv",sep = "")
  #write.table(clinical_data,file = clinical_outputfilename,quote = FALSE,sep = "\t",row.names = FALSE)
  write.csv(clinical_data,file = clinical_outputfilename,quote = FALSE) 
  print("clinic data has been downloaded")
  
  
}else{
  print(paste(cancer_type,"is not a correct id in TCGA, 
               Please refer to the document 'rawdata/TCGA_projects.csv' 
               and try again"))
}


















