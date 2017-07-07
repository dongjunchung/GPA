
shinyGPA <- function(out=NULL){

  smat <- out$pTestPval
  pmat <- out$pmat
  if (is.null(rownames(pmat))) {
    rownames(pmat) <- paste("SNP_",1:nrow(pmat),sep="")
  }
  if (is.null(colnames(pmat))) {
    colnames(pmat) <- paste("GWAS_",1:ncol(pmat),sep="")
    rownames(smat) <- colnames(smat) <- colnames(pmat)
  }
  
  
  
  ui <- fluidPage(
    
    headerPanel( 'ShinyGPA'),
    
    sidebarPanel(
      wellPanel(
        
        textInput("title", "Plot Title", "", width="100%"),
        
        
       # fileInput("file", label ="File input"),
      #  fluidRow(column(12, verbatimTextOutput("value"))),
        
        
        downloadButton("downloadPlot", label = "Download Plot"),
        
        hr(),
        
        conditionalPanel( "output.fileUploaded",
                          checkboxGroupInput(inputId="checkGroup", 
                                             label="Choose phenotypes",
                                             choices = list("1" = 1, "2" = 2, "3" = 3),
                                             inline=TRUE,
                                             selected = 1
                          )
        )
        ,
        
        
        numericInput(inputId = "clusters", 
                     label="Number of Clusters",
                     1,
                     min=0, max=1),
        
        
        sliderInput(inputId="lambda",
                    label="Choose a lambda value",
                    value=0, min=-1, max=1, step=0.05, #for simulation
                    width="100%"), 
        
        plotOutput("flexible", height = 200),
        
        
        numericInput(inputId = "fontsize", 
                     label="Label Font Size",
                     5,
                     min=1, max=100),
        
        numericInput(inputId = "axisfontsize", 
                     label="Axis Font Size",
                     15,
                     min=1, max=100),
        
        
        
        bsCollapse(id="clusterInfo", open="Clustering Details",
                   
                   bsCollapsePanel("Clustering Details",
                                   
                                   selectInput(inputId="clustertype",
                                               label="Type of clustering",
                                               choices= list("K-means", "Hierarchical"),
                                               selected=1),
                                   
                                   uiOutput("clusterdeets")
                                   
                   )),
        
        bsCollapse(id="lamInfo", open=NULL,
                   
                   bsCollapsePanel("Lambda Details",
                                   
                                   numericInput('lowLam', 'Lower Lambda Value', -1, min = -100, max = 100),
                                   
                                   numericInput('highLam', 'Upper Lambda Value', 1, min = -100, max = 100),
                                   
                                   numericInput('setLam', 'Default Lambda Value', 0, min = -100, max = 100),
                                   
                                   numericInput('stepLam', 'Lambda Value Increments', 0.01, min = -100, max = 100)
                   ))
        
        
        
        
      )
    ),
    
    
    
    mainPanel(
      tabsetPanel(
        tabPanel("Plot", plotOutput(outputId = "plot", 
                                    height=700, width=700,
                                    brush = brushOpts(id = "plot_brush")
        )),
        
        tabPanel("Info",
                 wellPanel(
                   conditionalPanel( 
                     "output.fileUploaded",
                     radioButtons(inputId="checklist1", 
                                  label="Choose Phenotype 1",
                                  choices = list("1" = 1, "2" = 2, "3" = 3),
                                  inline=TRUE,
                                  selected = 1)
                     ,
                     
                     radioButtons(inputId="checklist2", 
                                  label="Choose Phenotype 2",
                                  choices = list("1" = 1, "2" = 2, "3" = 3),
                                  inline=TRUE,
                                  selected = 1)),
                   
                   numericInput(inputId = "fdr", 
                                label="Set false discovery rate (FDR)",
                                value= 0.1,
                                min=0, max=1,
                                step = 0.01,
                                width= "20%")
                   
                   
                 ),
                 tableOutput(outputId = "table")
                 )
      )
    ))
  
  
  
  server <- function(input, output, session) {
    
    
    
    
    
    output$clusterdeets <- renderUI({
      
      if (is.null(input$clustertype))
        return()
      
      switch(input$clustertype,
             "K-means" = list(selectInput("algorK", label="Algorithm",
                                          choices = c("Hartigan-Wong", "Lloyd",
                                                      "MacQueen"),
                                          selected = 1),
                              
                              numericInput("iter", label="Number of Iterations",
                                           value = 1000)),
             
             "Hierarchical"= selectInput("algorH", label="Algorithm",
                                         choices = c( "ward.D", "ward.D2", "single", 
                                                      "complete", "average", 
                                                      "Mcquitty"),
                                         selected = "complete") )
    })
    
    #uploading smat 
    
    
    data <- reactive({
      #req(input$file)
      
      #inFile <- input$file
      
      
      
      
      #dataname <- load(inFile$datapath)
      #smat <- get(dataname[1])
      if ( is.null(colnames(smat)) ) {
        colnames(smat) <- paste("V",1:ncol(smat),sep="")
      }
      
      low<- input$lowLam
      high<-input$highLam
      set<-input$setLam
      incr<-input$stepLam
      
      updateSliderInput(session, "lambda", 
                        value=set, 
                        min=low,
                        max=high,
                        step=incr)
      
      updateCheckboxGroupInput(session, "checkGroup",
                               choices= colnames(smat),
                               inline=TRUE,
                               selected = colnames(smat)
      )
      
      updateNumericInput(session, inputId = "clusters", "Number of Clusters",
                         3, min=0, max=ncol(smat)-1) 
      
      updateRadioButtons(session, "checklist2",
                         choices= colnames(smat),
                         inline=TRUE,
                         selected = ""
      )
      
      updateRadioButtons(session, "checklist1",
                         choices= colnames(smat),
                         inline=TRUE,
                         selected = ""
      )
      
      
      #return(smat)
    })
    
    output$fileUploaded <- reactive({
      return(!is.null(data()))
    })
    outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)
    
    
    
    
    
    
    
    
    #Making distance transformation plot
    
    output$flexible <- renderPlot({
      
      #lambda flexible plot 
      values <- seq(0,1,by=0.001)
      
      
      lambda2 = input$lambda
      if ( lambda2 != 0 ) {
        dd <- ( values^lambda2 -1 ) / lambda2
      } else {
        dd <- log10(values)
      }
      
      ddmin <- min( dd[ dd > -Inf ] )
      dd[ dd == -Inf ] <- 1.5 * ddmin
      dd <- dd - 2 * ddmin
      
      plot( values, dd, type="l", 
            main="",
            xlab="raw p-values",
            ylab="transformed distance")
      lines(range(values), range(dd), col="lawngreen")
      
      
    })
    
    
    plotInput <- function(){
      
      palette(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"))
      
      #gpa.sel <- data()[ match( input$checkGroup, colnames(data()) ), 
      #                   match( input$checkGroup, colnames(data()) ) ]
      gpa.sel <- smat[ match( input$checkGroup, colnames(smat) ), 
                         match( input$checkGroup, colnames(smat) ) ]
      
      #boxcox distance transformation
      lambda <- input$lambda
      if ( lambda != 0 ) {
        dd <- ( gpa.sel^lambda -1 ) / lambda
      } else {
        dd <- log10(gpa.sel)
      }
      ddmin <- min( dd[ dd > -Inf ] )
      dd[ dd == -Inf ] <- 1.5 * ddmin
      dd <- dd - 2 * ddmin
      
      a <- as.list(1:ncol(gpa.sel)) 
      
      
      #dimensional scaling
      fit.isomap <- isomap(vegdist(dd), ndim=2, epsilon=0.45)
      plotdata <- as.data.frame(fit.isomap$points)
      
      
      
      colnames(plotdata) <- c("PC1","PC2")
      
      #clustering options
      clusterCols <- reactive({
        
        switch(input$clustertype,
               "K-means" = kmeans(dist(plotdata), input$clusters, nstart=100,
                                  iter.max =input$iter, 
                                  algorithm =input$algorK)$cluster,
               "Hierarchical"= cutree(hclust(dist(plotdata), method=input$algorH), 
                                      k=input$clusters)
               
               
        )}) 
      
      
      #output main plot
      ggplot(plotdata) +
        geom_point(aes(PC1, PC2), 
                   color = clusterCols(), 
                   size=4, shape=(clusterCols()+14)) +
        geom_text_repel(aes(PC1, PC2), label = rownames(plotdata), size=input$fontsize) +
        theme_bw() +
        theme(panel.grid.major= element_blank(),
              panel.grid.minor=element_blank()
        )+
        labs( x = "Coordinate 1", y = "Coordinate 2",
              caption = paste("Lambda =", input$lambda)
              #subtitle= paste("Lambda =", input$lambda)
        )+ 
        ggtitle(input$title)   +
        theme(plot.title = element_text(size = 30, face = "bold"),
              axis.title = element_text(size=input$axisfontsize),
              plot.caption=element_text(margin=margin(t=15),
                                        face="italic", size=16))
    }
    
    
    output$plot <- renderPlot({
      print(plotInput())
    })
    
    
    
    #highlight info
    output$plot_brushinfo <- renderPrint({
      cat("input$plot_brush:\n")
      str(input$plot_brush)
    })
    
    output$plot_brushed_points <- DT::renderDataTable({
      
      gpa.sel <- data()[ match( input$checkGroup, colnames(data()) ), 
                         match( input$checkGroup, colnames(data()) ) ]
      
      
      #boxcox distance transformation
      lambda <- input$lambda
      if ( lambda != 0 ) {
        dd <- ( gpa.sel^lambda -1 ) / lambda
      } else {
        dd <- log10(gpa.sel)
      }
      ddmin <- min( dd[ dd > -Inf ] )
      dd[ dd == -Inf ] <- 1.5 * ddmin
      dd <- dd - 2 * ddmin
      
      a <- as.list(1:ncol(gpa.sel)) 
      
      
      #dimensional scaling
      fit.isomap <- isomap(vegdist(dd), ndim=2, epsilon=0.45)
      plotdata <- as.data.frame(fit.isomap$points)
      
      
      
      colnames(plotdata) <- c("PC1","PC2")
      
      # dat <- as.data.frame(cbind( plotdata, gpa.sel ))
      
      res <- brushedPoints(plotdata, input$plot_brush , xvar="PC1", yvar="PC2", allRows=TRUE)
      #problem with the x and y variables
      
      res
    })
    
    output$downloadPlot <- downloadHandler(
      filename = function() {
        if( input$title != "" ) {
          paste(input$title, '-Shinyplot','.pdf', sep='')
        } else {
          'data-Shinyplot.pdf'
        }
      },
      content = function(file) {
        pdf(file)
        print(plotInput())
        dev.off()
      })
    
    
    #table for the info tab
    
     tableInput <- function(){
        pid <- sort(c(input$checklist1, input$checklist2)) 
        fit <-which(out$combs[,1] == which( colnames(smat) == pid[1] ) & out$combs[,2] == which( colnames(smat) == pid[2]) )
        tableInfo <- assoc(out$fitGPA[[fit]], FDR=input$fdr, pattern="11")
        
        if (!is.null(rownames(pmat))) {
             snpnames <- rownames(pmat)
           } else {
             snpnames <- paste("SNP_",1:nrow(pmat),sep="")
           }
           
           snptable <- data.frame(snpnames, fdr(out$fitGPA[[fit]], pattern="11"))[tableInfo==1,]
           snptable <- snptable[ order(snptable[,2]), ]
           colnames(snptable) <- c("SNP ID", "local FDR")
           
           return(snptable)
        
     }
    
     output$table <- renderTable({
       print(tableInput())
     })
    
  }
  
  shinyApp(ui=ui, server=server)

}

