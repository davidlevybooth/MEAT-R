#######################################################
#                                                     #
# Meat-R, a web appliction for multivariate           #
# ordination in vegan for R                           #
#                                                     #
# by David Levy-Booth, UBC Microbiology               #
#                                                     #
# version 0.04                                        #
#   started: 03/29/2015                               #
#   updated: 04/04/2015                               #
#   change log:                                       #
#     added PCA module 03/29/2015                     #
#     added colouring by groups 03/30/2015            #
#     added secondary overlay 04/02/2015              #
#     added colour palette selection 04/02/2015       #
#     added RDA module 04/03/2015                     #
#     added R2-adj to RDA 04/04/2015                  #
#     added PCoA and db-RDA modules 04/04/2015        #
#     added legends 04/04/2015                        #
#                                                     #
# associated files:                                   #
#   ui.R                                              #
#   vegan_functions.R                                 #
#                                                     #
#######################################################


## Load libraries and dependencies
library(shiny)
library(vegan)
library(lattice)

source("vegan_functions.R")

## Start Server Functions
shinyServer(function(input, output) {
## Plot ordination using anal function,
## Removing cases from dataframe
## Store in plotInput

### First produce slider so we can get the input values for plotting,
### Other dynamically rendered UI elements at end of document
  
  output$primary_slider <- renderUI ({ 
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    data1<-read.csv(inFile$datapath, header=TRUE) 
    if(typeof(input$file1) == "list")
      return(sliderInput("primary_range", round = TRUE, step = 1, label = strong("Select data columns:"), min = 1, max = ncol(data1), value = c(2, ncol(data1))))
  })
  
######################   Plots!

## Load the file
  plotInput <- function(){ #1
    inFile <- input$file1
        
## Get data start and end points
    if (is.null(input$primary_range))
      return(NULL)
    if (input$primary_range[1] < 2) {
    data_start = 2 }
    if (input$primary_range[1] >= 2) {
      data_start = input$primary_range[1] }
  
    if (input$primary_range[1] == input$primary_range[2]) {
    data_end = input$primary_range[2]+1
    }
    else { data_end = input$primary_range[2] }


## Get secondary start and end points  
  if (is.null(input$secondary_range) == FALSE ) {
  second_end = input$secondary_range[2]
  if (input$secondary_range[1] <= input$primary_range[2]) {
    second_start = input$primary_range[2]+1
  }
  else { second_start = input$secondary_range[1] } }


## Get constrained start and end points  
if (is.null(input$constrained_range) == FALSE && input$anal_type == "RDA" || is.null(input$constrained_range) == FALSE && input$anal_type == "db-RDA" ) {
  constrained_end = input$constrained_range[2]
  if (input$constrained_range[1] == input$primary_range[2]) {
    constrained_start = input$primary_range[2]+1
  }
  else { constrained_start = input$constrained_range[1] } } 


## Load the data into data frame
data1<-read.csv(inFile$datapath, header=TRUE)

  ##Set the analysis type 
  if(input$anal_type == "PCA") {
    anal_type = "PCA" }
  if(input$anal_type == "PCoA") {
    anal_type = "PCoA" }
  if(input$anal_type == "RDA" && is.null(input$constrained_range) == FALSE) {
    anal_type = "RDA" }
  if(input$anal_type == "RDA" && is.null(input$constrained_range) == TRUE) {
    anal_type = "PCA" }
  if(input$anal_type == "PCA" && input$primary_range[2] == ncol(data1)) {
    anal_type = "PCA" }
  if(input$anal_type == "PCoA" && input$primary_range[2] == ncol(data1)) {
  anal_type = "PCoA" }
  if(input$anal_type == "db-RDA" && is.null(input$constrained_range) == FALSE) {
   anal_type = "db-RDA" }
  if(input$anal_type == "db-RDA" && is.null(input$constrained_range) == TRUE) {
    anal_type = "PCoA" }
    

  ## Colour palette options
  if(is.null(input$pal_choice) == FALSE && input$pal_choice == 1 ) {
    pal_choice = c("#4EA1E8", "#A1FF96","#F963FF", "#FFFA6D", "#E89C4E", "#E85541") }
  if(is.null(input$pal_choice) == FALSE && input$pal_choice == 2 ) {
    pal_choice = c("#B64926", "#64C492", "#FFF0A5", "#FFB03B", "#33634A", "#69C8DF") }
  if(is.null(input$pal_choice) == FALSE && input$pal_choice == 3 ) {
    pal_choice = c("#222222", "#BBBBBB", "#DDDDDD", "#555555", "#000000", "#999999") }
 if(is.null(input$pal_choice) == TRUE) {
    pal_choice = c("#4EA1E8", "#A1FF96","#F963FF", "#FFFA6D", "#E89C4E", "#E85541") }

  
  
## Plot options
  
### PCA
    if(length(data1[data_start:data_end]) > 0 && anal_type == "PCA") {      
      
    if(input$overlay == TRUE && data_end<ncol(data1) && is.null(input$secondary_range) == FALSE) {
    if(input$use_groups == TRUE && is.null(pal_choice) == FALSE && is.null(input$legend) == FALSE) {
      
    anal_group_overlay(anal_type, input$biplot, data1[,1], pal_choice, input$legend, data1[second_start:second_end], data1[data_start:data_end], NULL) } 
    else {  
    anal_overlay(anal_type, input$biplot, data1[second_start:second_end], data1[data_start:data_end], NULL) } } #}
      
    else {
    if(input$use_groups == TRUE && is.null(pal_choice) == FALSE && is.null(input$legend) == FALSE) { ## Groups##d
    anal_group(anal_type, input$biplot, data1[,1], pal_choice, input$legend, data1[data_start:data_end], NULL) }
    else {
    anal(anal_type, input$biplot, data1[data_start:data_end], NULL)}}
    }

### RDA
  if(length(data1[data_start:data_end]) > 0 && is.null(input$constrained_range) == FALSE && anal_type == "RDA")
    if(input$overlay == TRUE && is.null(input$secondary_range) == FALSE && input$primary_range[2] < ncol(data1)) {
      if(input$use_groups == TRUE && is.null(pal_choice) == FALSE && is.null(input$legend) == FALSE && is.null(input$constrained_range) == FALSE) {
        
        anal_group_overlay(anal_type, input$biplot, data1[,1], pal_choice, input$legend, data1[input$secondary_range[1]:input$secondary_range[2]], data1[data_start:data_end], data1[constrained_start:constrained_end]) } 
      else {  
        anal_overlay(anal_type, input$biplot, data1[input$secondary_range[1]:input$secondary_range[2]], data1[data_start:data_end], data1[constrained_start:constrained_end]) } } #}

  else {
    if(input$use_groups == TRUE && is.null(pal_choice) == FALSE && is.null(input$legend) == FALSE && is.null(input$constrained_range) == FALSE) { ## Groups##d
      anal_group(anal_type, input$biplot, data1[,1], pal_choice, input$legend, data1[data_start:data_end], data1[constrained_start:constrained_end]) }
    else {
      anal(anal_type, input$biplot, data1[data_start:data_end], data1[constrained_start:constrained_end])}}

### PCoA
if(length(data1[data_start:data_end]) > 0 && anal_type == "PCoA") {      
  
  if(input$overlay == TRUE && data_end<ncol(data1) && is.null(input$secondary_range) == FALSE) {
    if(input$use_groups == TRUE && is.null(pal_choice) == FALSE && is.null(input$legend) == FALSE) {
      
      anal_group_overlay(anal_type, input$biplot, data1[,1], pal_choice, input$legend, data1[second_start:second_end], data1[data_start:data_end], NULL) } 
    else {  
      anal_overlay(anal_type, input$biplot, data1[second_start:second_end], data1[data_start:data_end], NULL) } } #}
  
  else {
    if(input$use_groups == TRUE && is.null(pal_choice) == FALSE && is.null(input$legend) == FALSE) { ## Groups##d
      anal_group(anal_type, input$biplot, data1[,1], pal_choice, input$legend, data1[data_start:data_end], NULL) }
    else {
      anal(anal_type, input$biplot, data1[data_start:data_end], NULL)}} }

### db-RDA
if(length(data1[data_start:data_end]) > 0 && is.null(input$constrained_range) == FALSE && anal_type == "db-RDA")
  if(input$overlay == TRUE && is.null(input$secondary_range) == FALSE && input$primary_range[2] < ncol(data1)) {
    if(input$use_groups == TRUE && is.null(pal_choice) == FALSE && is.null(input$legend) == FALSE && is.null(input$constrained_range) == FALSE) {
      
      anal_group_overlay(anal_type, input$biplot, data1[,1], pal_choice, input$legend, data1[input$secondary_range[1]:input$secondary_range[2]], data1[data_start:data_end], data1[constrained_start:constrained_end]) } 
    else {  
      anal_overlay(anal_type, input$biplot, data1[input$secondary_range[1]:input$secondary_range[2]], data1[data_start:data_end], data1[constrained_start:constrained_end]) } } #}

else {
  if(input$use_groups == TRUE && is.null(pal_choice) == FALSE && is.null(input$legend) == FALSE && is.null(input$constrained_range) == FALSE) { ## Groups##d
    anal_group(anal_type, input$biplot, data1[,1], pal_choice, input$legend, data1[data_start:data_end], data1[constrained_start:constrained_end]) }
  else {
    anal(anal_type, input$biplot, data1[data_start:data_end], data1[constrained_start:constrained_end])}}

}


## Output plot from plotInput
output$plot <- renderPlot({
  
  print(plotInput())
})


###################### Tables!

## Produce data table from inputed data

  output$contents <- renderTable({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)   
    read.csv(inFile$datapath, header=TRUE)   
  })
  

################### Saving!

## Produce PNG file for download button
  output$download_png <- downloadHandler(
    filename <- function() {
      paste(input$anal_type,"_", Sys.Date(),'.png',sep='') },
    content = function(file) {
      png(file, width = 700, height = 400, units = "px", pointsize = 12,
          bg = "white", res = NA)
      plotInput()
      dev.off()
    }) 


## Only display download button when inputed data is a list
## (Instead of "NULL")
  output$save_plot <- renderUI ({
    if(typeof(input$file1) == "list")
      return(downloadButton('download_png','Download as PNG'))
  })

## Produce  EPS file for download button
output$download_eps<- downloadHandler(
  filename <- function() {
    paste(input$anal_type,"_", Sys.Date(),'.eps',sep='') },
  content = function(file) {
    postscript(file)
    plotInput()
    dev.off()
  }) 


## Only display download button when inputed data is a list
## (Instead of "NULL")
output$save_eps <- renderUI ({
  if(typeof(input$file1) == "list")
    return(downloadButton('download_eps','Download as EPS'))
})


######## Other dynamically rendered UI

## Grouping prompt
output$group_prompt <- renderUI ({
  if(typeof(input$file1) == "list")
return(checkboxInput("use_groups", label = "Colour by group label", value = FALSE)) })

## Biplot prompt
output$biplot_prompt <- renderUI ({
  if(typeof(input$file1) == "list")
     return(checkboxInput("biplot", label = "Display biplot vectors", value = FALSE)) })

## Overlay secondary variables prompt
  output$second_prompt <- renderUI ({
  if(typeof(input$file1) == "list")
    return(checkboxInput("overlay", label = "Overlay secondary variables", value = FALSE)) })

## Secondary overlay range slider
output$second_slider <- renderUI ({ 
  inFile <- input$file1
  if (is.null(inFile))
    return(NULL)
  data1<-read.csv(inFile$datapath, header=TRUE)
  if (is.null(input$overlay))
    return (NULL)
  if (input$overlay == TRUE && input$primary_range[2] < ncol(data1)) {
    if(input$anal_type == "PCA" || input$anal_type == "PCoA")
    return(sliderInput("secondary_range", round = TRUE, step = 1, label = strong("Select overlay columns:"), min = 1, max = ncol(data1), value = c((input$primary_range[2]+1), ncol(data1))))}
  if (input$overlay == TRUE && is.null(input$constrained_range) == FALSE && input$constrained_range[2] < ncol(data1)) {
    if(input$anal_type == "RDA" || input$anal_type == "db-RDA")
    return(sliderInput("secondary_range", round = TRUE, step = 1, label = strong("Select overlay columns:"), min = 1, max = ncol(data1), value = c((input$constrained_range[2]+1), ncol(data1))))}
  
})

### Colour Pallet selection
output$pal_info <-  renderUI ({
  if(typeof(input$file1) == "list")
    if(is.null(input$use_groups) == FALSE && input$use_groups == FALSE)
      return(column(4, p(" ")))
    if(is.null(input$use_groups) == FALSE && input$use_groups == TRUE)
  return(column(4, 
         strong("Color Palettes:"), 
         p(" ", style = "font-si2pt"),
         img(src = "pal_both.png"), 
         selectInput("pal_choice", 
                      label = p(" "), 
                      choices = list("West End" = 1, 
                                     "North Isle" = 2,
                                     "Point Grey" = 3
                                     
                                     ),
                      selected = 1),
         checkboxInput("legend", label = "Plot Legend", value = FALSE)
          ))

})


############  RDA/db-RDA Contrained variable slider

output$constrained_slider <- renderUI ({ 
  inFile <- input$file1
  if (is.null(inFile))
    return(NULL)
  data1<-read.csv(inFile$datapath, header=TRUE)
  if (input$anal_type == "RDA" || input$anal_type == "db-RDA")
    if(is.null(input$primary_range) == FALSE && input$primary_range[2] < ncol(data1)) {
    return(sliderInput("constrained_range", round = TRUE, step = 1, label = strong("Constraining variables:"), min = 1, max = ncol(data1), value = c((input$primary_range[2]+1), ncol(data1))))}
})


##########  Distance matrix prompt calculation for PCOA/NDMS 

output$distance_prompt <- renderUI ({
  inFile <- input$file1
  if (is.null(inFile))
    return(NULL)
  data1<-read.csv(inFile$datapath, header=TRUE)
  if (input$anal_type == "PCoA" || input$anal_type == "db-RDA")
    return(selectInput("dist_choice", 
                  label = p("Dissimilarity Index:"), 
                  choices = list("Bray-Curtis" = "bray"
                  ),
                  selected = "bray"))  })

})
