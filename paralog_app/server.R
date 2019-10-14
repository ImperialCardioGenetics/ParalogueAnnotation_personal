library(shiny)
library(DT)
library(shinythemes)

#library(tidyverse)

shinyServer(function(input, output){
  
  get_paralog<-function(){
    
    #input<-data.frame(chr="1",pos="114713907",ref="T",alt="A")
    
    if(input$format=='pick'){
      req(input$chr)
      req(input$pos)
      req(input$ref)
      req(input$alt)
      var<-paste(input$chr,input$pos,input$ref,input$alt,sep = ":")
      var=var[nzchar(x=var)]
      input_data<-data.frame(mutation=var)
      result<-predict_output(raw_data,input_data)
  }else{
    if(input$format=='paste'){
      req(input$var)
      var<-unlist(strsplit(input$var,split="\\s+"))
      var=var[nzchar(x=var)]
      input_data<-data.frame(mutation=var)
      colnames(input_data)<-"mutation"
      result<-predict_output(raw_data,input_data)
      }
    }
  return(result)
  }
  
  output$paralog<-renderDataTable(DT::datatable(get_paralog(), options = list(paging = FALSE),rownames = F)) # set options for table eg. per page lines
  
  output$download <- downloadHandler(
    filename = function() {
      paste("paralog_annotation", ".tsv",sep="") # need to give specific name?
    },
    content = function(file) {
      write.table(get_paralog(), file, row.names = FALSE,quote = F,sep="\t")
    }
  )
  
})


predict_output<-function(output,input_data){
  
  # I opened formated and outputted again the RData odject from Nick to edit rownames and NAs 
  # This can be done here within a function
  # load the RData or rds file
  # load("data/paralog_tmp.RData")
  
  load("data/place_holder_results_datatable.RData")
  
  # exported RData with hardcoded dataframe name, change to raw_data
  raw_data<-place_holder_results_datatable
  
  # write new column chr:pos:ref:alt to look up
  raw_data$var<-paste(raw_data$CHROM.x,raw_data$POS.x,raw_data$REF.x,raw_data$ALT.x,sep=":")
  
  # select dataframe columns if not formated
  # paralog_tmp<-subset(raw_data,select=c(raw_data$var,raw_data$CHROM.y ,raw_data$POS.y,raw_data$REF.y,raw_data$ALT.y,raw_data$ID.y ,raw_data$SYMBOL,raw_data$Protein_position.y,raw_data$REF_Amino_acids.y,raw_data$ALT_Amino_acids.y ,raw_data$Codons.y,raw_data$Para_Z_score.y))
  raw_data<-subset(raw_data,select=c(var,CHROM.y ,POS.y,REF.y,ALT.y,ID.y ,SYMBOL,Protein_position.y,
                                     REF_Amino_acids.y,ALT_Amino_acids.y ,Codons.y,Para_Z_score.y))
  # rename dataframe columns for webpage
  colnames(raw_data)<-c("Variant_ID","chr","pos","REF","ALT","ClinVar_ID","gene","protein","reference_AA", "alt_AA","codons","para_Z_score" )
  
  # select the vars ## this can be done first to reduce filtering time if final dataset is huge
  output<- raw_data[raw_data$Variant_ID %in% input_data$mutation,]
  
  return(output)
}




