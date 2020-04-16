library(shiny)
library(datasets)
library(ggplot2)
library(tidyverse)
library(readxl)

is_outlier <- function(x) {
    return(x < quantile(x, 0.25,na.rm =TRUE) - 1.5 * IQR(x,na.rm =TRUE) | x > quantile(x, 0.75,na.rm =TRUE) + 1.5 * IQR(x,na.rm =TRUE))
}

ui <- shinyUI(fluidPage(
    titlePanel("HMPPMP Drug Database"),
    tabsetPanel(
        tabPanel("Upload File",
                 titlePanel("Upload File"),
                 sidebarLayout(
                     sidebarPanel(
                         fileInput('file1', 'Choose File',
                                   accept=c('text/csv', 
                                            'text/comma-separated-values,text/plain', 
                                            '.csv','.xlsx')),
                         
                         # added interface for uploading data from
                         # http://shiny.rstudio.com/gallery/file-upload.html
                         tags$br(),
                         selectInput('library','Library',c("Cambridge","Kinase","FDA")),
                         #checkboxInput('header', 'Header', TRUE),
                         #radioButtons('sep', 'Separator',
                                      #c(Comma=',',
                                        #Semicolon=';',
                                        #Tab='\t'),
                                     # ','),
                         #radioButtons('quote', 'Quote',
                                      #c(None='',
                                        #'Double Quote'='"',
                                        #'Single Quote'="'"),
                                      #'"')
                     ),
                     mainPanel(
                         textOutput('contents'),
                         tableOutput('sampleinfo')
                     )
                 )
        ),
        tabPanel("Drug Z score Distribution",
                 pageWithSidebar(
                     headerPanel('Drug of Interest'),
                     sidebarPanel(
                         
                         # "Empty inputs" - they will be updated after the data is uploaded
                         selectInput('xcol', 'Drug', ""),
                         #selectInput('ycol', 'Y Variable', "", selected = "")
                         
                     ),
                     mainPanel(
                         textOutput('numberoflines'),
                         plotOutput('MyPlot'),
                         tableOutput('contents2'),
                         tableOutput('contents2_1')
                         
                     )
                 )
        ),
        
        tabPanel("Sample Hits",
                 pageWithSidebar(
                     headerPanel('Select the Sample and Drug'),
                     sidebarPanel(
                         
                         # "Empty inputs" - they will be updated after the data is uploaded
                         selectInput('ycol', 'Sample', ""),
                         selectInput('xcol2', 'Drug', "")
                         
                     ),
                     mainPanel(
                         textOutput("numberofhits"),
                         plotOutput("Outlier"),
                         tableOutput('contents3')
                     )
                 )
        ),
        tabPanel("Target View",
                 pageWithSidebar(
                     headerPanel('Select'),
                     sidebarPanel(
                         
                         # "Empty inputs" - they will be updated after the data is uploaded
                         
                         textInput('target', 'Target',""),
                         selectInput('ycol2', 'Sample to highlight', ""),
                         actionButton("submit", label = "Submit")
                         
                     ),
                     mainPanel(
                         textOutput("numberofdrugs"),
                         plotOutput("drugswithsametarget"),
                         tableOutput('contents4')
                     )
                 )
        )
    )
)
)









server <- shinyServer(function(input, output, session) {
    # added "session" because updateSelectInput requires it
    
    
    data <- reactive({ 
        req(input$file1) ## ?req #  require that the input is available
        
        inFile <- input$file1 
        
        # tested with a following dataset: write.csv(mtcars, "mtcars.csv")
        # and                              write.csv(iris, "iris.csv")
        df <- read_excel(inFile$datapath, sheet=input$library)#header = input$header, sep = input$sep,quote = input$quote,
                      
        colnames(df)<-gsub(" ","_",colnames(df))
        colnames(df)<-gsub("-","_",colnames(df))
        infocol <- length(grep("Info",df[1,]))
        # Update inputs (you could create an observer with both updateSel...)
        # You can also constraint your choices. If you wanted select only numeric
        # variables you could set "choices = sapply(df, is.numeric)"
        # It depends on what do you want to do later on.
        
        updateSelectInput(session, inputId = 'xcol', label = 'Drug',
                          choices = df$Drug[-1], selected = df$Drug[-1])
        updateSelectInput(session, inputId = 'xcol2', label = 'Drug',
                          choices = df$Drug[-1], selected = df$Drug[-1])
        updateSelectInput(session, inputId = 'ycol', label = 'Sample',
                          choices = colnames(df)[-c(1:infocol)], selected = colnames(df)[-c(1:infocol)])
        updateSelectInput(session, inputId = 'ycol2', label = 'Sample to Highlight',
                          choices = colnames(df)[-c(1:infocol)], selected = colnames(df)[-c(1:infocol)])
        
        return(df)
    })
    
    
    output$contents <- renderText({
       df<- data()
       infocol <- length(grep("Info",df[1,]))
       num1<-nrow(df)-1
       num2<-ncol(df)-infocol
       paste("There are ", num1, " drugs.","\n","There are ", ncol(df)-infocol, " lines screened.")
    })
    
    output$sampleinfo <- renderTable({
        df<-data()
        samplecol <- length(grep("Info",df[1,]))+1
        df2<-data.frame(Sample=colnames(df)[samplecol:ncol(df)],Type=t(df[1,samplecol:ncol(df)]))
        df2
    })
    
    
    output$MyPlot <- renderPlot({
        df<-data()
        infocol <- length(grep("Info",df[1,]))
        datacol <- infocol+1
        x <- as.data.frame(t(df[which(df$Drug==input$xcol), datacol:ncol(df)]))
        colnames(x)[1]<-"Zscore"
        x$Sample<-rownames(x)
        df2<-data.frame(Sample=colnames(df)[datacol:ncol(df)],Type=t(df[1,datacol:ncol(df)]))
        df2$Type <- as.character(df2$Type)
        highlight<- c("Atypical_Teratoid_Rhabdoid_Tumour_(ATRT)","High_Grade_Glioma","Diffuse_intrinsic_pontine_glioma_(DIPG)")
        for (i in 1:length(df2$Sample)){
            df2$Type[i]<-ifelse(is.na(match(df2$Type[i],highlight)),"Other", df2$Type[i])
        }
        x<-merge(df2,x,by="Sample")
        x$Zscore <- as.numeric(as.character(x$Zscore))
        ggplot(x,aes(x=Zscore))+ 
            geom_density(color="lightblue",fill="lightblue",alpha=0.5)+
            theme(legend.position="bottom")+
            theme_classic()+ 
            geom_vline(xintercept=-1.5, linetype="dashed", color = "red")+
            geom_rug(aes(color=Type))+
            theme(legend.position="bottom",legend.text = element_text( size=12, face="bold"))+
            scale_color_manual(values=c('deeppink3','#E69F00','#56B4E9','#999999'))
        
    })

    
    output$contents2 <- renderTable({
        df<-data()
        infocol <- length(grep("Info",df[1,]))
        mediancol <- infocol -2
        df2<-df[which(df$Drug==input$xcol),1:infocol]
        df2[,infocol] <- as.numeric(as.character(df2[,infocol]))
        df2[,mediancol] <- as.numeric(as.character(df2[,mediancol]))
        df2
    })
    
    output$contents2_1 <- renderTable({
        df<-data()
        datacol <- length(grep("Info",df[1,]))+1
        df2<-df[which(df$Drug==input$xcol), ]
        df3<-data.frame(Sample =colnames(df)[datacol:ncol(df)], Type= unlist(df[1,datacol:ncol(df)]), Zscore=t(df2[1,datacol:ncol(df)]))
        df3$Zscore <- as.numeric(as.character(df3$Zscore))
        df3$Sample <- as.character(df3$Sample)
        df3$Type <- as.character(df3$Type)
        df3<-df3[df3$Zscore<=-1.5,]
        df3<-df3[order(df3$Zscore),]
        df3
    })
    
    output$numberoflines <- renderText({
        df<-data()
        datacol <- length(grep("Info",df[1,]))+1
        df2<-df[which(df$Drug==input$xcol), ]
        df3<-data.frame(Info =colnames(df)[datacol:ncol(df)], Zscore=t(df2[1,datacol:ncol(df)]))
        df3$Zscore <- as.numeric(as.character(df3$Zscore))
        df3$Info <- as.character(df3$Info)
        df3<-df3[df3$Zscore<=-1.5,]
        Lines <- ncol(df)-datacol+1
        percentage<-as.integer(nrow(df3)*100/Lines)
        paste("Percentage of inhibition is ", percentage, "% of", Lines, "lines.")
    })
    
    output$contents3 <- renderTable({
        df<-data()
        infocol <- length(grep("Info",df[1,]))
        datacol <- infocol+1
        sampleinfo <- df[1,]
        mediancol <- infocol -2
        outlier<- gather(df[-1,c(1,datacol:ncol(df))],key="Sample",value="Zscore",-c("Drug"))
        outlier$Zscore <- as.numeric(as.character(outlier$Zscore))
        outlier <- outlier %>% group_by(Drug) %>% mutate(Outlier = ifelse(is_outlier(Zscore), "Yes", "No"))
        df2<-cbind(df[-1,1:infocol],df[-1,colnames(df)==input$ycol],Outlier=outlier$Outlier[outlier$Sample==input$ycol])
        df2[,datacol] <- as.numeric(as.character(df2[,datacol]))
        df2<- df2[which(df2[,datacol]<=-1.5),]
        df2<-df2[order(df2[,datacol]),]
        df2[,infocol] <- as.numeric(as.character(df2[,infocol]))
        df2[,mediancol] <- as.numeric(as.character(df2[,mediancol]))
        df2
    })
    
    output$numberofhits <- renderText({
        df<-data()
        infocol <- length(grep("Info",df[1,]))
        df2<-cbind(df[-1,1:infocol],df[-1,colnames(df)==input$ycol])
        datacol <- infocol+1
        df2[,datacol]<-as.numeric(as.character(df2[,datacol]))
        df2<- df2[which(df2[,datacol]<=-1.5),]
        df2<-df2[order(df2[,datacol]),]
        type<-df[1,colnames(df)==input$ycol]
        paste0("This line is ", type,". There are ", nrow(df2), " hits.")
    })
    
    output$Outlier <- renderPlot({
        df<-data()
        infocol <- length(grep("Info",df[1,]))
        datacol <- infocol+1
        sampleinfo <- df[1,]
        outlier<- gather(df[-1,c(1,datacol:ncol(df))],key="Sample",value="Zscore",-c("Drug"))
        outlier$Zscore <- as.numeric(as.character(outlier$Zscore))
        outlier <- outlier %>% group_by(Drug) %>% mutate(Outlier = ifelse(is_outlier(Zscore), Sample, NA))
        df2<-outlier[outlier$Drug==input$xcol2,]
        df2$Type <- t(df[1,datacol:ncol(df)])
        df2<-as.data.frame(df2)
        matchtype<-as.character(df[1,colnames(df)==input$ycol])
        df2$Type[df2$Type != matchtype] <- "Other"
        df2$Type[df2$Sample == input$ycol] <- input$ycol
        df2$Type <- factor(df2$Type, levels = c("Other", matchtype, input$ycol))
        ggplot(df2, aes(x=Drug, y=Zscore))+
            geom_boxplot(color="gray24",width=0.3) +
            geom_dotplot(aes(fill=Type,color=Type),binaxis='y', stackratio=1.2,stackdir='center', dotsize=1,alpha=0.8) +
            scale_fill_manual(values=c("gray30","chartreuse1","firebrick2"))+scale_color_manual(values=c("gray30","chartreuse1","firebrick2"))+
            theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold"),legend.key.size = unit(1, "cm"),legend.title=element_text(size=12), legend.text=element_text(size=12),axis.text.x = element_text(angle = 45, hjust = 1))
    })
    
    
    # observe event for updating the reactiveValues
    text_reactive <- eventReactive( input$submit, {
        input$target})
    
    
    
    output$numberofdrugs<- renderText({
        df<-data()
        type<-text_reactive()
        df2<-df[grep(type,df$Target),]
        paste0("There are ", nrow(df2)," drugs targeting ", type, ".  Note: * indicate outliers.")
    })
    output$drugswithsametarget<-renderPlot({
        df<-data()
        type<-text_reactive()
        infocol <- length(grep("Info",df[1,]))
        datacol = infocol+1
        mediancol = infocol-2
        df2<-df[grep(type,df$Target),c(1,mediancol,datacol:ncol(df))]
        df3<-gather(df2,key="Sample",value="Zscore",-c("Drug",colnames(df)[mediancol]))
        df3<-as.data.frame(df3)
        df3$Zscore<-as.numeric(unlist(df3$Zscore))
        for (i in 1:length(df3$Sample)){
            df3$Type[i] <- df[1,match(df3$Sample[i],colnames(df))]
}
        df3$Type<-unlist(df3$Type)
        matchtype<-as.character(df[1,colnames(df)==input$ycol2])
        df3<-as.data.frame(df3)
        df3$Type[df3$Type != matchtype] <- "Other"
        df3$Type[df3$Sample == input$ycol2] <- input$ycol2
        df3$Type <- factor(df3$Type, levels = c("Other", matchtype, input$ycol2))
        df3$Type[df3$Type == "Other"]<-NA
        ggplot(df3,aes(x=Drug, y=Zscore))+
            geom_boxplot(color="gray24",width=0.3,outlier.colour="gray30", outlier.shape=8, outlier.size=2) +
            geom_dotplot(aes(fill=Type,color=Type),binaxis='y', stackratio=1.2,stackdir='center', dotsize=1,alpha=0.8) +
            scale_fill_manual(values=c("chartreuse1","firebrick2"))+scale_color_manual(values=c("chartreuse1","firebrick2"))+
            theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold"),legend.key.size = unit(1, "cm"),legend.title=element_text(size=12), legend.text=element_text(size=12),axis.text.x = element_text(angle = 45, hjust = 1),legend.position="bottom")
        
    })
    output$contents4<-renderTable({
        df<-data()
        type<-text_reactive()
        infocol <- length(grep("Info",df[1,]))
        mediancol = infocol-2
        df2<-df[grep(type,df$Target),1:infocol]
        df2[,infocol] <- round(as.numeric(unlist(df2[,infocol])),2)
        df2[,mediancol] <- round(as.numeric(unlist(df2[,mediancol])),2)
        df2
    })
})

shinyApp(ui, server)
