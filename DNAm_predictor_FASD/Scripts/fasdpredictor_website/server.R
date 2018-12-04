#load packages
library(shiny)
library(shinyjs)
library(caret)
library(pROC)
library(plotROC)
library(ggplot2)
library(reshape)
library(gbm)
library(e1071)
library(rsconnect)

#load data
#setwd("~/Desktop/fasdpredictorapp/")
load("gbmFit.648.FASD.RData")

#runApp("~/Desktop/fasdpredictorapp/")


#training.data <- read.csv("training.beta.noPAE.csv", row.names=1)
#training.meta <- read.csv("training.meta.noPAE.csv", row.names =1)

#shinyApp(ui = ui, server = server)

server <- function(input, output) {
  
   output$CpG.input <-  downloadHandler(
    filename <- function() {
      paste("cpgs2input", "txt", sep=".")
    },
    
    content <- function(file) {
      write.table(colnames(gbmFit.648.FASD$trainingData)[1:648], file, row.names=F, quote =F,col.names=F)
    }
  )
   
   output$example.data <-  downloadHandler(
     filename <- function() {
       paste("example", "csv", sep=".")
     },
     
     content <- function(file) {
       data.2.write <- gbmFit.648.FASD$trainingData[1:20,1:648]
       rownames(data.2.write) <- paste("Sample_", 1:20, sep ="")
       write.csv(data.2.write, file, quote=F)
     }
   )
   
   output$example.meta <-  downloadHandler(
     filename <- function() {
       paste("meta.example", "csv", sep=".")
     },
     
     content <- function(file) {
       write.csv(gbmFit.648.FASD$trainingData[1:20,649], file, quote=F, row.names=F)
     }
   )
  
  output$Predicted_class <- renderTable({
    inFile <- input$file1
    known <- input$file2
    if (is.null(inFile))
      return(NULL)
    
    if(is.null(known)) {
      predict.results <- data.frame(Sample = rownames(read.csv(inFile$datapath, row.names=1)),
                                    Score = predict(gbmFit.648.FASD, read.csv(inFile$datapath, 
                                                                              row.names =1), 
                                                    type = "prob"),
                                    Predicted = predict(gbmFit.648.FASD, read.csv(inFile$datapath, 
                                                                                   row.names = 1))
                                    )
      print(predict.results)}
    else{
      predict.results <- data.frame(Sample = rownames(read.csv(inFile$datapath, row.names=1)),
                                    Score = predict(gbmFit.648.FASD, read.csv(inFile$datapath, 
                                                                              row.names =1), 
                                                    type = "prob"),
                                    Real = read.csv(known$datapath)[,1],
                                    Predicted = predict(gbmFit.648.FASD, read.csv(inFile$datapath, 
                                                                                   row.names = 1))
                                    )
      print(predict.results)
    }
    
  })
  
  output$results.down <-  downloadHandler(
    
    
    filename <- function() {
      paste("fasd_predictor_results",Sys.Date(), "csv", sep=".")
    },
    
    
    content <- function(file) {
      inFile <- input$file1
      known <- input$file2
      if (is.null(inFile))
        return(NULL)
      
      if(is.null(known)) {
        predict.results <- data.frame(Sample = rownames(read.csv(inFile$datapath, row.names=1)),
                                      Score = predict(gbmFit.648.FASD, read.csv(inFile$datapath, 
                                                                                row.names =1), 
                                                      type = "prob"),
                                      Predicted = predict(gbmFit.648.FASD, read.csv(inFile$datapath, 
                                                                                    row.names = 1))
        )}
      else{
        if(length(grep("fasd", known)) >0){
          known[which(known == "fasd")] <- "FASD"
        }
        
        if(length(grep("control", known)) >0){
          known[which(known == "control")] <- "Control"
        }
        predict.results <- data.frame(Sample = rownames(read.csv(inFile$datapath, row.names=1)),
                                      Score = predict(gbmFit.648.FASD, read.csv(inFile$datapath, 
                                                                                row.names =1), 
                                                      type = "prob"),
                                      Real = read.csv(known$datapath)[,1],
                                      Predicted = predict(gbmFit.648.FASD, read.csv(inFile$datapath, 
                                                                                    row.names = 1))
        )
      }
      
      write.csv(predict.results, file, quote=F, row.names=F)
    }
  )
  
  
  output$Confusion_matrix <- renderPrint({
    inFile <- input$file1
    known <- input$file2
    if (is.null(known))
      return(invisible(NULL))
    
    if(length(grep("fasd", known)) >0){
      known[which(known == "fasd")] <- "FASD"
    }
    
    if(length(grep("control", known)) >0){
      known[which(known == "control")] <- "Control"
    }
    
    confusion.matrix <- confusionMatrix(predict(gbmFit.648.FASD, read.csv(inFile$datapath, 
                                                                          row.names=1)), 
                                        factor(read.csv(known$datapath)[,1], levels = c("FASD","Control")))
    print(confusion.matrix)
  })
  
  output$hist <- renderPlot({
    inFile <- input$file1
    known <- input$file2
    if (is.null(known))
      return(NULL)
    
    if(length(grep("fasd", known)) >0){
      known[which(known == "fasd")] <- "FASD"
    }
    
    if(length(grep("control", known)) >0){
      known[which(known == "control")] <- "Control"
    }
    
    votes <- predict(gbmFit.648.FASD, read.csv(inFile$datapath,  row.names =1), type = "prob")
    histogram(~votes$Control|factor(read.csv(known$datapath)[,1], levels = c("FASD","Control")), 
              xlab="Probability of Poor Segmentation")
    
    
  })
  
  output$roc <- renderPlot({
    inFile <- input$file1
    known <- input$file2
    if (is.null(known))
      return(NULL)
    
    if(length(grep("fasd", known)) >0){
      known[which(known == "fasd")] <- "FASD"
    }
    
    if(length(grep("control", known)) >0){
      known[which(known == "control")] <- "Control"
    }
    
    votes <- predict(gbmFit.648.FASD, read.csv(inFile$datapath,  row.names =1), type = "prob")
    gbm.ROC.plot <- data.frame(predictor= votes$Control,
                               response=factor(read.csv(known$datapath)[,1], levels = c("FASD","Control")))
    
    ggplot(gbm.ROC.plot, aes(d = response, m = predictor))+
      geom_roc(n.cuts = 0, labels = F)+
      #geom_rocci(labels =F, size=0, sig.level = 0.05)+
      ggtitle("ROC of input data")+
      style_roc(xlab = "1 - Specificity")+
      coord_flip()+
      theme(axis.title.y = element_text(size=25, face="bold", vjust=1.2),
            axis.title.x = element_text(size=25, face="bold", vjust=1.2),
            axis.text.y = element_text(size=18, color="black"),
            axis.text.x = element_text(size=18, color="black"),
            plot.title = element_text(face="bold", size=30))
  })
 
  output$varImpPlot <- renderPlot({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    cpg.names <- colnames(read.csv(inFile$datapath, row.names = 1))

    var.importance <- data.frame(importance = varImp(gbmFit.648.FASD)$importance[order(varImp(gbmFit.648.FASD)$importance, decreasing = T)[1:183],])
    var.importance$cpg <- rownames(varImp(gbmFit.648.FASD)$importance)[order(varImp(gbmFit.648.FASD)$importance, decreasing = T)[1:183]]
    melt.var.importance <- melt(var.importance)
    melt.var.importance$cpg <- factor(melt.var.importance$cpg, levels = melt.var.importance$cpg[1:183])
    melt.var.importance$var.col <- rep("Absent",183)
    melt.var.importance$var.col[which(melt.var.importance$cpg %in% cpg.names)] <- "Present"
    melt.var.importance$var.col <- factor(melt.var.importance$var.col, levels = c("Present","Absent"))
    col.cpgs <- rep("red",183)
    col.cpgs[which(melt.var.importance$var.col == "Present")] <- "white"
    
    ggplot(melt.var.importance, aes(x= cpg, y= value, fill = var.col))+
      geom_bar(stat="identity")+
      ggtitle("Weighted impact of top 183 most important CpGs (missing in red)")+
      scale_fill_manual(values = c("darkgrey","red"))+
      xlab("Missing CpGs in red")+
      ylab("Weighted importance of CpGs")+
      theme_classic()+
      guides(fill = guide_legend(title = "CpG in dataset"))+
      theme(axis.text.x= element_text(angle = 45, hjust = 1, size= 6, face = "bold", colour= col.cpgs))
   
      })
  
  output$missing.CpGs <- renderTable({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    cpg.names <- colnames(read.csv(inFile$datapath, row.names = 1))
    
    var.importance <- data.frame(importance = varImp(gbmFit.648.FASD)$importance[order(varImp(gbmFit.648.FASD)$importance, decreasing = T)[1:183],])
    var.importance$cpg <- rownames(varImp(gbmFit.648.FASD)$importance)[order(varImp(gbmFit.648.FASD)$importance, decreasing = T)[1:183]]
    melt.var.importance <- melt(var.importance)
    melt.var.importance$cpg <- factor(melt.var.importance$cpg, levels = melt.var.importance$cpg[1:183])
    
    all.cpgs <- colnames(gbmFit.648.FASD$trainingData)[1:648]
    if(length(which(all.cpgs %in%cpg.names)) == 648){return(print("No CpGs are missing from your dataset!"))}
    
    missing.cpgs <- data.frame(cpg = all.cpgs[-which(all.cpgs %in% cpg.names)], value = 0)
    
    missing.var <- melt.var.importance[which(melt.var.importance$cpg %in% missing.cpgs$cpg),c(1,3)]
    missing.var.2 <- rbind(missing.var,missing.cpgs[-which(missing.cpgs$cpg %in% melt.var.importance$cpg),] )
    missing.var.2$value <- round(missing.var.2$value, 5)
    missing.var.2$value[which(missing.var.2$value == 0)] <- "negligible"
    colnames(missing.var.2) <- c("Missing CpG", "Importance")
    
    print(missing.var.2)
    })
  
}

