library(shiny)
library(shinythemes)
library(shinyjs)

ui <-  fluidPage(
  fluidRow(column(width = 12, align="center",h1(strong("FASD prediction by DNA methylation")))),
  #fluidRow(column(width = 12, align="center",em(h4("for research purposes only*")))),
  fluidRow(column(width = 12, align = "center", img(src="banner4.png", align = "center", height="40"))),
  fluidRow(column(width = 12, align="center",h4(em("***This tool is for research purposes only and is not suitable for clinical applications***")))),
  fluidRow(column(width = 12,h2(" "))),

  #  titlePanel("FASD prediction by DNA methylation"),
  theme = shinytheme("cosmo"),
  sidebarLayout(
    sidebarPanel(
      strong(em(h4("Input data here", align = "center"))),
      
      fileInput("file1", "Choose CSV file of beta values
                (rows = sample, columns = CpG)",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")),
      
      fileInput("file2","Optional: Choose CSV file of class (FASD/Control)",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")),
      tags$hr(),
      strong(em(h4("Examples", align = "center"))),
      fluidRow(downloadButton("CpG.input", label = "CpGs to input", class = "btn-primary"),
               downloadButton("example.data", label = "Data format",class = "btn-primary"),
               downloadButton("example.meta", label = "Metadata format",class = "btn-primary"))
      ),
  
  
  mainPanel(
    tabsetPanel(type = "tabs",
                tabPanel("Welcome",strong("Welcome to our DNA methylation predictor of FASD!"),
                         br(),br(),
                         "Fetal alcohol spectrum disorder (FASD) is a developmental disorder that manifests through a range of cognitive, adaptive, physiological, 
                         and neurobiological deficits resulting from prenatal alcohol exposure. Although the North American prevalence is currently estimated at 
                         2-5%, FASD has proven difficult to identify in the absence of the overt physical features characteristic of more severe cases of fetal 
                         alcohol syndrome. Building on our previous work identifying a distinct DNA methylation signature in children and adolescents with FASD, 
                         we have developed a preliminary epigenetic predictor of FASD using DNA methylation data for 648 CpG sites across the genome.",
                         br(),br(),
                         "The predictor was developed using buccal epithelial cell (BEC) DNA methylation data from the Illumina 450K array. As such, results using 
                         data from other sources or tissue types may vary. Please note that this tool may also be confounded for ethnicity or other factors, as 
                         limited by the clinical cohorts described in the publications below, and may yield different results in other populations.",
                         br(),br(),
                         "We would appreciate any feedback on the tool and its performance in your hands would be greatly appreciated. We also welcome any inquiries 
                         for collaborations and suggestions to improve the predictor. See below for contact information."),
                
                tabPanel("Results", 
                         downloadButton("results.down", label = "Download results", class = "btn-primary", style="float:right"),
                         "Please input your data in CSV format in the panel to the left. Data input examples can be found in the bottom lefthand panel.",
                         "Note that all 648 CpGs are not necessary to run the predictor, but missing CpGs may impact the quality of the results.",
                         tableOutput("Predicted_class")),
                
                tabPanel("Missing CpGs", "Of the CpGs used in this tool, 183 are most important for accurate prediction. Missing variables can impact the performance of the predictor, some more strongly than others. Please use all 648 CpGs for best results.",
                         br(),
                         plotOutput("varImpPlot"),
                         br(),
                         tableOutput("missing.CpGs")),
                tabPanel("Histogram", "A histogram of segmentation errors will be output if you provide known sample groups (metadata).",
                         br(),
                         plotOutput("hist")),
                
                tabPanel("Confusion matrix", "A confusion matrix will be output if you provide known sample groups (metadata).",
                         br(),
                         verbatimTextOutput("Confusion_matrix")),
                
                tabPanel("ROC curve", "A receiver operating characteristic (ROC) curve will be output if you provide known sample groups (metadata).", 
                         br(),
                         plotOutput("roc")),
                tabPanel("More information",
                         em("A brief overview of our methods and results can be found here. For more information, please refer to the manuscripts cited at the bottom of the page."), 
                         br(),br(),
                         
                         fluidRow(column(width = 6, align = "left", 
                                         "To assess whether DNA methylation data could be used to predict FASD status, we created a predictive algorithm of FASD using 
                                         machine learning approaches. First, we selected the normalized DNA methylation data (beta-values) of 179 samples from the NeuroDevNet cohort 
                                         (NDN - GSE80261) in the 648 initial probes that were also found in the Kid's Brain Health Network data (KBHN - GSE109042). 
                                         Our strategy was to build the predictor using an initial training cohort (NDN), followed by subsequent evaluation in the test cohort (KBHN).",
                                         br(),br(),
                                         "In parallel, buccal epithelial cell samples from an independent published autism spectrum disorder (ASD) cohort were used to assess the specificity of the model 
                                         in the FASD cohorts. To this end, we used a publically available dataset of 450K array data from the BECs of 48 individuals with 
                                         ASD and 48 typically developing controls from the gene expression omnibus (GSE50759).",
                                         br(), br(),
                                         "The predictive model of FASD status was created using DNA methylation data and the caret package in R. First, a predictive model was created using 
                                         stochastic gradient boosting on the NDN cohort (83 FASD cases, 96 controls) using the beta-values of the differentially methylated probes identified 
                                         in the NDN study (648 probes). The parameters of the model were optimized for area under the receiver operating characteristic (ROC) curve by 
                                         grid tuning for repeated cross-validation (number of trees 50–1500; 1, 5, or 9 interaction depth; 0.1 shrinkage). The optimal model for predicting 
                                         clinical FASD status using 648 probes was 550 trees, 1 of interaction depth, and 20 minimum observations per node. "),
                                  column(width = 6, align = "left", img(src="figure3.png", align = "left", height="500"))),
                         br(),br(),
                         strong("Description of FASD cohorts"),
                         br(),"Genome-wide DNA methylation patterns were analyzed using the Illumina HumanMethylation450 array in buccal epithelial cells. 
                         Data from FASD cases and typically-developing controls from the below cohorts were used in the development of our epigenetic predictor. 
                         Importantly, there was a discrepancy in self-declared ethnicity between FASD cases and controls. Although we partially accounted for
                         ethnicity-related effects by selecting CpGs with less ethnic bias (Portales-Camasar, Lussier, et al.), this particularity could potentially 
                         influence the results of the predictor in other populations.",
                         br(),br(),
                         fluidRow(column(width = 6, align = "left", em("NeuroDevNet cohort (training dataset)"), br(),
                                         img(src="table1_ndn.png",align="left", width = "350")),
                                  column(width = 6, align = "left", em("Kid's Brain Health Network cohort (test dataset)"), br(),
                                         img(src="table1.png",align = "left",width = "350"))),
                         
                         br(),br(),
                         strong("Prediction results on KBHN and ASD cohorts"),br(),
                         fluidRow(column(width = 6, align ="left","The predicted sensitivity and specificity for the training cohort were 0.879 and 0.944, respectively, for an area under the curve of 0.977 
                         (95% confidence intervals, 0.972–0.982; Fig. 4a). The performance of the predictor on the training data indicated that DNA methylation could 
                         be used to distinguish FASD cases and controls, although these results will need to be carefully assessed in independent test sets or clinical settings.",
                         br(),br(),
                         "To get a better understanding of the utility of this tool, we next assessed the predictive model using the normalized, batch-corrected 
                         DNA methylation data of the KBHN cohort as a test set. Of note, these data were not corrected for any covariates or surrogate variables 
                         other than batch correction. As expected for analysis of an independent test set, the model performed at a lower level in this cohort, 
                         displaying 0.917 sensitivity, 0.75 specificity, and 0.920 area under the ROC curve (Table 3; Fig. 4a). The balanced accuracy of the model 
                         in this cohort was 0.833 (95% CI 0.698–0.925), and the ROC curve was not significantly different from the one obtained in the training cohort (p = 0.192).",
                         br(),br(),
                         "Given the discrepancies in ethnic backgrounds between FASD and control groups, the misclassified samples were assessed for differences
                         in self-reported ethnicity, caregiver status, age, and buccal cell-type proportions in the classification. We did not identify any skew 
                         of these data in the misclassified controls. Although every misclassified individual with FASD had a previous diagnosis of ARND, a category 
                         that was present in high proportion within this cohort, no other patterns emerged between the correctly and incorrectly classified individuals with FASD.
                         Taken together, these findings suggested that differences in these demographic variables between the groups did not drive their classification.",
                         br(),br(),
                         "BEC samples from an independent published autism spectrum disorder (ASD) cohort were used to assess the specificity of the model in the FASD cohorts. 
                         To this end, we used a publically available dataset of 450K array data from the BECs of 48 individuals with ASD and 48 typically developing controls 
                         from the gene expression omnibus (GSE50759). Using processed GEO data from this cohort, the predictor correctly identified the vast majority of individuals 
                         in the cohort as non-FASD. The model only misclassified 1 individual as FASD, for a specificity of 0.990 (95% CI 0.943–0.9997), higher than the predicted 
                         specificity in the training set. This sample, a 3-year-old female with ASD (51% African ancestry, 41% European ancestry) did not have any particular d
                         istinguishing features compared to the correctly classified samples, suggesting that the predictive model was not biased for ASD, sex, age, or African 
                         ancestry in this independent cohort."),
                         column(width = 6, align = "left", 
                                img(src="table3.png", align="left",height = "450"),
                                br(),br(),br(),
                                img(src="figure4.png", align = "left", height = "450")))
                         
                         
    )
    )
  )
),

fluidRow(column(width = 12, h4("Please cite:"))),
fluidRow(column(width = 10, offset = 1, h4("Portales-Casamar E*, Lussier AA*, Jones MJ, MacIsaac JL, Edgar RD, Mah SM,
                                          Barhdadi A, Provost S, Lemieux-Perreault LP, Cynader MS, Chudley AE, Dubé MP,
                                          Reynolds JN, Pavlidis P, Kobor MS.", tags$a(href="https://doi.org/10.1186/s13072-016-0074-4", "DNA methylation signature of human fetal
                                          alcohol spectrum disorder."), "Epigenetics Chromatin. 2016 Jun 29;9:25. doi:
                                          10.1186/s13072-016-0074-4."))),

fluidRow(column(width = 10, offset = 1,h4("Lussier AA, Morin AM, MacIsaac JL, Salmon J, Weinberg J, Reynolds JN, Pavlidis
                                          P, Chudley AE, Kobor MS.", tags$a(href="https://doi.org/10.1186/s13148-018-0439-6", "DNA methylation as a predictor of fetal alcohol spectrum
                                          disorder."), 
                                          "Clin Epigenetics. 2018 Jan 12;10:5. doi: 10.1186/s13148-018-0439-6."))),

fluidRow(column(width = 12, align = "center", em(h4("For more information or collaboration opportunities, please contact Dr. Alexandre Lussier, PhD - fasdpredictor[at]gmail.com"))))

)

