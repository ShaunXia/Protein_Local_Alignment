
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(DT)
library(ggplot2)

shinyServer(function(input, output) {
  source("localAlignment.R")
  obClick<-0
  obSub<<-0
  test_str<-localAlignment("KRSLWRSRAPG","SRLWRSG",list("sc_mat"=1,
                                                   "mt_gap"=-10,"mt_ext"=-1))
  
  output$dp_table = DT::renderDataTable(
    test_str$dp_table,class = 'cell-border stripe',
    server = TRUE,
    selection = list(target = 'cell',selected = test_str$optimal),
    
    options = list(dom = 't',bSort=FALSE,pageLength = 100)
  )
  
  output_mat<-matrix(0,1,6)
  opt_res<-data.frame()
  colnames(output_mat)<-c("Score","Start","End","Length","String1","String2")
  output_mat<<-as.data.frame(output_mat)
  output$res_table = DT::renderDataTable(
    output_mat,class = 'cell-border stripe',
    server = TRUE,
    selection = list(target = 'row',mode='single'),
    options = list(pageLength = 10)
  )
  
  select_path<-list()
  
  proxy = dataTableProxy('dp_table')
  
  
  # observeEvent(input$rows, {
  #   print(input$rows)
  #   print(Sys.time())
  #   proxy %>% selectCells(as.numeric(c(2,2)))
  #   #proxy %>% selectRows(1)
  #   
  # })
  
  observeEvent(input$align_but, {
    
    obClick<<-0
    
    if(input$score_matric=="BLOSUM62")
      sc_mat<-1
    else
      sc_mat<-2
    mat_val<-list("sc_mat"=sc_mat,
                  "mt_gap"=input$s_gap,"mt_ext"=input$s_ext)

    test_str<<-localAlignment(input$str1,input$str2,mat_val)
    print(test_str)
    output$dp_table = DT::renderDataTable(
      test_str$dp_table,class = 'cell-border stripe',
      server = TRUE,
      selection = list(target = 'cell',selected = test_str$optimal),
      options = list(dom = 't',bSort=FALSE,pageLength = 100)
    )
    mat<-matrix(0,1,6)
    colnames(mat)<-c("Score","Start","End","Length","String1","String2")
    
    
    
  })
 
  #MFRTKRSALVRRLWRSRAPGGEDEEEGAGGGGGGGELRGEGATDSRAHGAGGGGPGRAGCCLGKAVRGAKGHHHPHPPAAGAGAAGGAEADLKALTHSVLKKLKERQLELLLQAVESRGGTRTACLLLPGRLDCRLGPGAPAGAQPAQPPSSYSLPLLLCKVFRWPDLRHSSEVKRLCCCESYGKINPELVCCNPHHLSRLCELESPPPPYSRYPMDFLKPTADCPDAVPSSAETGGTNYLAPGGLSDSQLLLEPGDRSHWCVVAYWEEKTRVGRLYCVQEPSLDIFYDLPQGNGFCLGQLNSDNKSQLVQKVRSKIGCGIQLTREVDGVWVYNRSSYPIFIKSATLDNPDSRTLLVHKVFPGFSIKAFDYEKAYSLQRPNDHEFMQQPWTGFTVQISFVKGWGQCYTRQFISSCPCWLEVIFNSR
  observeEvent(input$get_sub, {
    
    obSub<<-1
    threshold<-input$s_threshold
    opt_res<<-getSubOptimal(threshold,test_str$dp_table,test_str$path_table)
    opt_res_plot<<-getSubOptimalForPlot(threshold,test_str$dp_table,test_str$path_table)
    
    #output_mat<<-matrix(0,0,6)
    #colnames(output_mat)<<-c("Score","Start","End","Length","String1","String2")
    #output_mat<<-as.data.frame(output_mat)
    #output_mat<<-transform(output_mat,Score = as.numeric(Score))
    
    #tp<-opt_res[[2]]


    output$res_table = DT::renderDataTable(
      opt_res,class = 'cell-border stripe',
      server = TRUE,
      selection = list(target = 'row',mode='single'),
      options = list(pageLength = 10),rownames=NULL
    )
    
    istr1<-input$str1
    istr2<-input$str2
    
    str1_len<-nchar(input$str1)
    str2_len<-nchar(input$str2)
    
    brks=c(1:str2_len)
    labs=strsplit(istr2,"")[[1]]
    
    output$dotplot<-renderPlot({
      ggplot(opt_res_plot,aes(Y,X,group=GP,color="red"))+
        geom_point()+
        geom_line()+ 
        scale_y_reverse()
    })
    #print(output_mat)
    
    #print(output_mat)
    #print(opt_res)
    #output$sub_UI <-renderUI({
    #  selectInput('sub_list', 'Suboptiomal', opt_res, multiple=TRUE, selectize=FALSE)
    # })

    #obSub<<-0
  })
  
  
  observeEvent(input$dp_table_cells_selected, {
    
    if(obClick==0)
    {
      obClick<<-1
      return(NULL)
    }

    last_select<-tail(input$dp_table_cells_selected,n=1)
    
    if (length(last_select)&&(last_select[2]!=0)) 
      {
      print(last_select)
      nowSel<-tail(input$dp_table_cells_selected,n=1)
      select_path<-getPath(test_str$path_table,nowSel[1],nowSel[2])

    mm<-matrix(0,0,2)

    for (j in 1:length(select_path)) {
      mm<-rbind(mm,select_path[[j]])
    }
    proxy %>% selectCells(NULL)
    proxy %>% selectCells(mm)

    
    }
    
  })
  
  
  #SubOptimal Table
  
  observeEvent(input$res_table_rows_selected, {
    
    if(obSub==0)
      return(NULL)
    last_select<-tail(input$res_table_rows_selected,n=1)
    print(last_select)
    tp_row<-as.character(opt_res[last_select,"end"])
    #print(as.character(tp_row))
    end_loc<-strsplit(tp_row,",")[[1]]
    end_loc<-as.numeric(end_loc)+1 # for title
    #print(end_loc[1])
    
    if (length(last_select)) 
    {
      select_path<-getPath(test_str$path_table,end_loc[1],end_loc[2])
      
      mm<-matrix(0,0,2)
      
      for (j in 1:length(select_path)) {
        mm<-rbind(mm,select_path[[j]])
      }
      proxy %>% selectCells(NULL)
      proxy %>% selectCells(mm)

    }
    
  })
  
  
  # print the selected indices
  output$output_area = renderPrint({
    
    last_select<-tail(input$dp_table_cells_selected,n=1)
    
    if (length(last_select)&&(last_select[2]!=0)) {
      
      nowSel<-tail(input$dp_table_cells_selected,n=1)
      select_path<-getPath(test_str$path_table,nowSel[1],nowSel[2])
      cat('Selected Path:\n')
      
      path_str<-getPathStr(select_path,test_str$path_table,test_str$dp_table)
      cat(sprintf("%5d  %s  %5d\n", path_str$start_loc[1],path_str$t1,path_str$end_loc[1]))
      cat(sprintf("%5d  %s  %5d\n", path_str$start_loc[2],path_str$t2,path_str$end_loc[2]))
      cat("\nScore:",test_str$dp_table[nowSel[1],nowSel[2]],"\n")
      
      cat("\n=================\n")
      print(test_str$setting)
      
      
    }
  })
  
  
})
