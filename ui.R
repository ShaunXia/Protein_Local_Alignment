
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(ggplot2)
fluidPage(

  sidebarLayout(
    sidebarPanel(
      textInput('str1', 'String 1',"KRSLWRSRAPG"),
      textInput('str2', 'String 2',"SRLWRSG"),
      
      hr(),
      fluidRow(
      column(5,numericInput('s_gap', 'Gap',-10)),
      column(5,numericInput('s_ext', 'Extension',-1)),
      column(9,selectInput('score_matric', 'Scoring Matrices', c("BLOSUM62","PAM250"), multiple=FALSE, selectize=FALSE))
      ),
      actionButton('align_but', 'Alignment'),

      hr(),
      h2("Filter"),
      
      hr(),
      numericInput('s_threshold', 'Threshold(Score>=)',5),
      actionButton('get_sub', 'Get Suboptiomal'),
      br(),
      hr(),
      #selectInput('sub_list', 'Suboptiomal', state.name, multiple=TRUE, selectize=FALSE),
      htmlOutput("sub_UI")
      #verbatimTextOutput('out3')
      
    ),
    
    mainPanel(
      tabsetPanel(type = "tabs", 
                  tabPanel("Table", 
                           DT::dataTableOutput('dp_table'),
                           hr(),
                           h4("Suboptimal Result"),
                           hr(),
                           DT::dataTableOutput('res_table'),
                           hr(),
                           verbatimTextOutput('output_area')
                           ),
                  tabPanel("Plot", plotOutput("dotplot"))
                  #tabPanel("Summary", verbatimTextOutput("summary"))
                  
      )
    )
    
  )
)