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
  
  output$paralog<-renderDataTable(DT::datatable(get_paralog(),
                                                options = list(paging = FALSE),# set options for table eg. per page lines
                                                rownames = F, 
                                                caption = htmltools::tags$caption(style = 'caption-side: bottom; text-align: center;','Table 1 : ', htmltools::em('Paralogous Variants'))
                                                ) %>%
                                                formatStyle(c("Variant_ID","Query_Gene"),  color = 'black', backgroundColor = 'lightgrey', fontWeight = 'bold')
                                  )
  
  output$download <- downloadHandler(
    filename = function() {
      paste("paralog_annotation", ".tsv",sep="") # need to give specific name?
    },
    content = function(file) {
      write.table(get_paralog(), file, row.names = FALSE,quote = F,sep="\t")
    }
  )
  
})






