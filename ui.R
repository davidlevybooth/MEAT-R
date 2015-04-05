## Meat-R UI 
## Associated files: 
##  server.R
##  vegan_functions.R

library(shiny)


shinyUI(fluidPage(

  sidebarLayout(
    sidebarPanel(
      br(),
      img(src = "meatr.png"),
      br(),
      br(),
      p("Welcome to MEAT-R, a web-based implementation of the", a("VEGAN package", 
          href = "http://vegan.r-forge.r-project.org/"), " for multivariate ordination [In Progress]."),
      br(),
      fileInput('file1', 'Choose CSV File',
                accept=c('text/csv', 
                         'text/comma-separated-values,text/plain', 
                         '.csv')),
      p("CSV file must contain variable names in the first row (header) and the first column should contain sample labels (groups)."),
      
      tags$hr(),
      
      selectInput("anal_type", label = h4("Select Analysis:"), 
                  choices = list("PCA", 
                                 "RDA",
                                 "PCoA",
                                 "db-RDA"), selected = 1)
      
    ),
    
    # Show a tabset that includes a plot and table view
    mainPanel(
      tabsetPanel(type = "tabs", 
          tabPanel("Plot",  
            plotOutput("plot"), 
              fluidRow(
                column(4, 
                  htmlOutput("primary_slider"),
                  htmlOutput("distance_prompt"),
                  htmlOutput("group_prompt"),
                  htmlOutput("biplot_prompt"),
                  htmlOutput("constrained_slider"),
                  htmlOutput("second_prompt"),
                  htmlOutput("second_slider") 
                  ),
                  htmlOutput("pal_info"),
                column(2, htmlOutput("save_plot"), p(" ", style = "font-si2pt"), htmlOutput("save_eps"))          
              )
          ),
          tabPanel("Data", tableOutput("contents"))
          )
    )
  )
))
