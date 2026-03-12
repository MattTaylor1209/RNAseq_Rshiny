#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

required_packages <- c("limma",
                       
                       "edgeR",
                       
                       "Glimma",
                       
                       "org.Dm.eg.db",
                       
                       "org.Mm.eg.db",
                       
                       "org.Rn.eg.db",
                       
                       "org.Hs.eg.db",
                       
                       "org.Dr.eg.db",
                       
                       "gplots",
                       
                       "RColorBrewer",
                       
                       "NMF",
                       
                       "BiasedUrn",
                       
                       "GO.db",
                       
                       "qusage",
                       
                       "Rsamtools",
                       
                       "Rsubread",
                       
                       "GenomicFeatures",
                       
                       "Rfastp",
                       
                       "biomaRt",
                       
                       "DESeq2",
                       
                       "IHW",
                       
                       "topGO",
                       
                       "apeglm",
                       
                       "clusterProfiler",
                       
                       "ashr",
                       
                       "goseq",
                       
                       "KEGGREST",
                       
                       "msigdbr",
                       
                       "GSEABase",
                       
                       "enrichplot",
                       
                       "gage",
                       
                       "gageData",
                       
                       "pathview",
                       
                       "sva",
                       
                       "RUVSeq",
                       
                       "vsn",
                       
                       "biomaRt",
                       
                       "DBI",
                       
                       "tidyverse",
                       
                       "readxl",
                       
                       "multcomp",
                       
                       "ggthemes",
                       
                       "ggpubr",
                       
                       "ggsignif",
                       
                       "ggrepel",
                       
                       "lsmeans",
                       
                       "rstatix",
                       
                       "ggtext",
                       
                       "RColorBrewer",
                       
                       "ggsci",
                       
                       "ggprism",
                       
                       "patchwork",
                       
                       "minpack.lm",
                       
                       "markdown",
                       
                       "treemap",
                       
                       "VennDiagram",
                       
                       "grid",
                       
                       "ggVennDiagram",
                       
                       "eulerr",
                       
                       "pheatmap",
                       
                       "plotly",
                       
                       "htmlwidgets",
                       
                       "shiny",
                       
                       "shinyFiles",
                       
                       "DT",
                       
                       "tidyr",
                       
                       "viridisLite")



# Install BiocManager if needed

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  
  install.packages("BiocManager")
  
}


# Try BiocManager first, then fallback to install.packages

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    message("Trying BiocManager::install('", pkg, "')")
    
    tryCatch({
      
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
      
    }, error = function(e_bioc) {
      
      message("BiocManager failed. Trying install.packages('", pkg, "')")
      
      tryCatch({
        
        install.packages(pkg, dependencies = TRUE)
        
      }, error = function(e_cran) {
        
        message("Failed to install '", pkg, "' via both methods.")
        
      })
      
    })
    
  }
  
}


# Packages loading


invisible(lapply(required_packages, library, character.only = TRUE))



# Increase upload size 

options(shiny.maxRequestSize = 50 * 1024^2)  # 50 MB


ui <- fluidPage(
  titlePanel("RNA-seq Analysis"),
  
  sidebarLayout(
    sidebarPanel(width = 3,
                 fileInput("countsFile", "Upload Count Data (.rds, .csv, .tsv, .txt)",
                           accept = c(".rds", ".csv", ".tsv", ".txt")),
                 conditionalPanel(
                   condition = "output.isPlainCountTable",
                   uiOutput("geneLengthStatusUI"),
                   actionButton("fetchBiomartBtn", HTML("Fetch gene lengths from biomaRt<br><small>(check organism first)</small>"),
                                icon = icon("download")),
                   helpText("Fetches median transcript lengths per gene from Ensembl (requires internet)."),
                   fileInput("geneLengthFile",
                             "Or upload Gene Lengths (.csv/.tsv/.txt — two columns: GeneID, Length)",
                             accept = c(".csv", ".tsv", ".txt")),
                   helpText("Gene lengths enable length-bias correction in GO analysis.",
                            "Uploaded file will override biomaRt-fetched or auto-detected lengths.")
                 ),
                 radioButtons("sampleInfoMode", "Sample info source:",
                              choices = c("Upload file" = "upload", "Define groups manually (must know column order)" = "manual"),
                              selected = "upload", inline = TRUE),
                 conditionalPanel(
                   condition = "input.sampleInfoMode == 'upload'",
                   fileInput("sampleInfoFile", "Upload SampleInfo.txt", accept = c(".txt", ".tsv"))
                 ),
                 conditionalPanel(
                   condition = "input.sampleInfoMode == 'manual'",
                   textInput("manualGroupName", "Group name:", placeholder = "e.g. Control"),
                   numericInput("manualGroupN", "Number of replicates (n):", value = 3, min = 1, step = 1),
                   actionButton("addGroupBtn", "Add group", icon = icon("plus")),
                   br(), br(),
                   tableOutput("manualGroupsTable"),
                   fluidRow(
                     column(6, actionButton("clearGroupsBtn", "Clear all", icon = icon("trash"))),
                     column(6, actionButton("confirmGroupsBtn", "Confirm groups", icon = icon("check"),
                                            class = "btn-success"))
                   ),
                   br(),
                   uiOutput("manualSampleCountMsg")
                 ),
                 conditionalPanel(
                   condition = "output.inputsReady",
                   tags$hr(),
                   h5("Sample filter"),
                   helpText("Deselect samples to exclude them from all analyses."),
                   uiOutput("sampleFilterUI"),
                   tags$hr()
                 ),
                 selectInput("organism", "Select organism:",
                             choices = c(
                               "Rat (Rattus norvegicus)" = "org.Rn.eg.db",
                               "Mouse (Mus musculus)" = "org.Mm.eg.db",
                               "Human (Homo sapiens)" = "org.Hs.eg.db",
                               "Drosophila (Drosophila melanogaster)" = "org.Dm.eg.db",
                               "Zebrafish (Danio rerio)" = "org.Dr.eg.db"
                             ),
                             selected = "org.Dm.eg.db"
                 ),
                 numericInput("lfcThreshold", "log2 Fold Change threshold", value = 0.5, min = 0),
                 numericInput("padjThreshold", "Adjusted p-value threshold (FDR)", value = 0.05, min = 0, max = 1),
                 checkboxInput("filterlow", "Automatically filter low-expression genes?", value = TRUE),
                 checkboxInput("deseqfilter", "Use DESeq2 independent filtering?", value = TRUE),
                 uiOutput("groupOrderUI"),
                 uiOutput("contrastSelectUI"),
                 uiOutput("convertToSymbolUI"),
                 checkboxInput("filterLncRNA", "Filter out lncRNAs (keep only protein-coding)", value = FALSE),
                 actionButton("analyzeBtn", "Run Analysis"),
                 conditionalPanel(
                   condition = "output.inputsReady",
                   checkboxInput("viewPreview", "Preview input files", value = TRUE)
                 )
                 
    ),
    
    mainPanel(
      conditionalPanel(
        condition = "output.inputsReady",
        tabsetPanel(
          tabPanel("Input Preview", 
                   tableOutput("previewTable"),
                   h4("Log"),
                   verbatimTextOutput("log", )),
          tabPanel("Sample Info", tableOutput("sampleInfo")),
          tabPanel("PCA", 
                   column(2,
                          sliderInput("pointsize", "Point size", min = 0, max = 20, step = 0.5, 
                                      value = 10)
                   ),
                   column(2,
                          sliderInput("labelsize", "Label font size", min = 0, max = 20, step = 0.5, 
                                      value = 10)
                   ),
                   column(2,
                          sliderInput("legendtitlesize", "Legend title size", min = 0, max = 20, step = 0.5, 
                                      value = 10)
                   ),
                   column(2,
                          sliderInput("axistitlesize", "Axis title size", min = 0, max = 30, step = 0.5, 
                                      value = 14)
                   ),
                   column(2,
                          sliderInput("axistextsize", "Axis text size", min = 0, max = 20, step = 0.5, 
                                      value = 10),
                   ),
                   column(2,
                          sliderInput("legendtextsize", "Legend text size", min = 0, max = 20, step = 0.5, 
                                      value = 10),
                   ),
                   plotOutput("pcaPlot")),
          tabPanel("DE Results", dataTableOutput("deTable"),
                   conditionalPanel(
                     condition = "output.analysisReady",
                     downloadButton("dl_res", "Download DE table")
                   )
          ),
          tabPanel("Heatmap",
                   fluidRow(
                     column(3,
                            radioButtons("heatmapGeneSource", "Gene selection:",
                                         choices = c("Top DE genes" = "de",
                                                     "Most variable genes (all samples)" = "variance"),
                                         selected = "de")
                     ),
                     column(3,
                            numericInput("heatmapTopN", "Number of genes", 
                                         value = 50, min = 10, max = 500, step = 10)
                     ),
                     column(3,
                            conditionalPanel(
                              condition = "input.heatmapGeneSource == 'de'",
                              selectInput("heatmapSortBy", "Sort genes by:",
                                          choices = c("Adjusted p-value" = "padj", 
                                                      "log2 Fold Change (absolute)" = "abs_lfc",
                                                      "Test statistic (absolute)" = "abs_stat"),
                                          selected = "padj")
                            )
                     ),
                     column(3,
                            selectInput("heatmapScale", "Scale data:",
                                        choices = c("Row (gene)" = "row",
                                                    "Column (sample)" = "column", 
                                                    "None" = "none"),
                                        selected = "row")
                     )
                   ),
                   fluidRow(
                     column(3,
                            checkboxInput("heatmapClusterRows", "Cluster rows (genes)", value = TRUE)
                     ),
                     column(3,
                            checkboxInput("heatmapClusterCols", "Cluster columns (samples)", value = TRUE)
                     ),
                     column(3,
                            selectInput("heatmapColorScheme", "Color scheme:",
                                        choices = c("Blue-White-Red" = "RdBu",
                                                    "Red-White-Blue" = "RdYlBu", 
                                                    "Viridis" = "viridis",
                                                    "Plasma" = "plasma"),
                                        selected = "RdBu")
                     ),
                     column(3,
                            checkboxInput("heatmapShowRowNames", "Show gene names", value = TRUE)
                     ),
                     column(3,
                            checkboxInput("heatmapShowColNames", "Show sample names", value = TRUE)
                     )
                   ),
                   fluidRow(
                     column(6,
                            selectInput("heatmapGroups", "Groups to include:",
                                        choices = NULL,
                                        selected = NULL,
                                        multiple = TRUE)
                     )
                   ),
                   plotOutput("heatmapPlot", height = "600px")
          ),
          tabPanel("GO",
                   fluidRow(
                     column(4,
                            actionButton("runGO", "Run GO"),
                            selectInput("GOontology", "Ontology", c("BP","MF","CC"), selected = "BP"),
                            numericInput("GOnumber", "No. of categories to show", value = 10, min = 1, step = 1)
                     ),
                     column(8,
                            plotOutput("GOBarplot"),
                            DT::dataTableOutput("goTable"),
                            downloadButton("dl_go", "Download GO table")
                     )
                   )
          ),
          tabPanel("GSEA",
                   fluidRow(
                     column(4, 
                            selectInput("gsea_type", "Pathways to analyse:",
                                        choices = c(
                                          "GO",
                                          "Wiki Pathways",
                                          "KEGG"
                                        ),
                                        selected = "GO"
                            ),
                            conditionalPanel(
                              condition = "input.gsea_type == 'Wiki Pathways'",
                              selectInput("wikiorg", "Select Wiki Pathways Organism",
                                          choices = c("Homo sapiens", 
                                                      "Mus musculus",
                                                      "Drosophila melanogaster",
                                                      "Rattus norvegicus",
                                                      "Danio rerio"),
                                          selected = "Drosophila melanogaster")
                            ),
                            conditionalPanel(
                              condition = "input.gsea_type == 'KEGG'",
                              selectInput("keggorg", "Select KEGG Organism",
                                          choices = c(
                                            "Human (Homo sapiens)" = "hsa",
                                            "Mouse (Mus musculus)" = "mmu",
                                            "Rat (Rattus norvegicus)" = "rno",
                                            "Drosophila (Drosophila melanogaster)" = "dme",
                                            "Zebrafish (Danio rerio)" = "dre"
                                          ),
                                          selected = "dme")
                            ),
                            conditionalPanel(
                              condition = "input.gsea_type == 'GO'",
                              selectInput("gseaOnt", "Ontology", c("BP","MF","CC"), selected = "BP")
                            ),
                            selectInput("gseaMetric", "Rank by", c("stat","log2FoldChange"), selected = "stat"),
                            numericInput("gseaP", "p-value cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
                            numericInput("gseaMin", "Min gene set size", value = 10, min = 5, step = 5),
                            numericInput("gseaMax", "Max gene set size", value = 500, min = 50, step = 50),
                            numericInput("numcategories", "No. of categories to show", value = 10, min = 1, step = 1),
                            actionButton("runGSEA", "Run GSEA"),
                     ),
                     column(8,
                            plotOutput("gseaDotplot"),
                            selectizeInput("gseaTerm", "Show enrichment plot for:", choices = NULL, multiple = FALSE),
                            plotOutput("gseaEnrichPlot"),
                            DT::dataTableOutput("gseaTable"),
                            downloadButton("dl_gsea", "Download GSEA table")
                     )
                   )
                   
          ),
          tabPanel("GAGE (KEGG)",
                   fluidRow(
                     column(4,
                            actionButton("runGAGE", "Run GAGE (KEGG)"),
                            numericInput("gageP", "q-value cutoff", value = 0.1, min = 0, max = 1, step = 0.01),
                            numericInput("gageMin", "Min geneset size", value = 10, min = 2, max = 1000, step = 1),
                            numericInput("gageMax", "Max geneset size", value = 1000, min = 2, max = 5000, step = 10),
                            selectInput("gageMetric", "Use metric",
                                        choices = c("log2FoldChange","stat"), selected = "stat"),
                            checkboxInput("gageSameDir", "Same direction (one-sided)", value = TRUE)
                     ),
                     column(8,
                            h4("Significant pathways (Upregulated)"),
                            DT::dataTableOutput("gageUpTable"),
                            downloadButton("dl_gage_up", "Download GAGE up table"),
                            h4("Significant pathways (Downregulated)"),
                            DT::dataTableOutput("gageDownTable"),
                            downloadButton("dl_gage_down", "Download GAGE down table")
                     )
                   ),
                   hr(),
                   fluidRow(
                     column(4,
                            selectInput("keggPathway", "Pathway to render (KEGG ID):", choices = NULL),
                            downloadButton("dl_pathview", "Download pathway PNG")
                     ),
                     column(8,
                            imageOutput("keggPathview", height = "600px")
                     )
                   )
          ),
          tabPanel("Standard volcano plot",
                   column(3,
                          sliderInput("topgenes", "How many top genes to label?",min = 1, max = 100, step = 1, 
                                      value = 50),
                   ),
                   column(3,
                          sliderInput("volcpointsize", "Point size", min = 0, max = 20, step = 0.5, 
                                      value = 2)
                   ),
                   column(3,
                          sliderInput("volclabelsize", "Label font size", min = 0, max = 20, step = 0.5, 
                                      value = 3)
                   ),
                   column(3,
                          sliderInput("volclegendtitlesize", "Legend title size", min = 0, max = 20, step = 0.5, 
                                      value = 14)
                   ),
                   column(3,
                          sliderInput("volclegendtextsize", "Legend text size", min = 0, max = 20, step = 0.5, 
                                      value = 12)
                   ),
                   column(3,
                          sliderInput("volcaxistitlesize", "Axis title size", min = 0, max = 30, step = 0.5, 
                                      value = 16)
                   ),
                   column(3,
                          sliderInput("volcaxistextsize", "Axis text size", min = 0, max = 20, step = 0.5, 
                                      value = 12)
                   ),
                   column(3,
                          checkboxInput("pvalueselect", "Use adjusted p-value?", value = TRUE)
                   ),
                   plotOutput("standardvolcanoplot",
                              width = "90%", height = "500px")
                   
          ),
          tabPanel("Interactive Volcano Plot", plotlyOutput("volcanoPlot"),
                   plotOutput("geneBoxplot")),
          
          # ---- Explore Gene Sets Tab ----
          tabPanel("Explore Gene Sets",
                   fluidRow(
                     column(12,
                            h4("Explore Genes in a GO / GSEA / KEGG Term"),
                            p("Select a term from a completed analysis, then view its genes as a heatmap, boxplots, or highlighted on a volcano plot.")
                     )
                   ),
                   hr(),
                   fluidRow(
                     column(4,
                            selectInput("exploreSource", "Source analysis:",
                                        choices = c("GO" = "go", "GSEA" = "gsea", "KEGG (GAGE)" = "kegg")),
                            uiOutput("exploreTermUI"),
                            actionButton("loadTermGenesBtn", "Load genes for this term",
                                         icon = icon("search")),
                            br(), br(),
                            uiOutput("exploreGeneCountMsg"),
                            hr(),
                            h5("Heatmap options"),
                            selectInput("exploreHeatScale", "Scale:",
                                        choices = c("Row (gene)" = "row", "Column (sample)" = "column", "None" = "none"),
                                        selected = "row"),
                            selectInput("exploreHeatColor", "Color scheme:",
                                        choices = c("Blue-White-Red" = "RdBu", "Viridis" = "viridis"),
                                        selected = "RdBu"),
                            hr(),
                            h5("Volcano options"),
                            sliderInput("exploreVolcLabelSize", "Label size", 1, 8, 3, step = 0.5),
                            sliderInput("exploreVolcPtSize", "Point size", 0.5, 6, 2, step = 0.5),
                            sliderInput("exploreVolcAlpha", "Background point opacity", 0.01, 1, 0.1, step = 0.01)
                     ),
                     column(8,
                            tabsetPanel(id = "exploreSubTabs",
                                        tabPanel("Gene Table",
                                                 DT::dataTableOutput("exploreGeneTable"),
                                                 downloadButton("dl_explore_genes", "Download gene list")
                                        ),
                                        tabPanel("Heatmap",
                                                 plotOutput("exploreHeatmap", height = "550px")
                                        ),
                                        tabPanel("Boxplots",
                                                 uiOutput("exploreBoxplotGeneUI"),
                                                 plotOutput("exploreBoxplot", height = "400px")
                                        ),
                                        tabPanel("Volcano",
                                                 plotOutput("exploreVolcano", height = "500px")
                                        )
                            )
                     )
                   )
          ),
          
          # ---- Compare Contrasts Tab ----
          tabPanel("Compare Contrasts",
                   fluidRow(
                     column(12,
                            h4("Compare DEGs Across Multiple Contrasts"),
                            p("Run additional contrasts from the same DESeq2 dataset and compare their DEGs via Venn/Euler diagrams, 
                   volcano plots, and downstream analyses (GO, GSEA, GAGE/pathview) on shared or unique gene sets.")
                     )
                   ),
                   hr(),
                   
                   # --- Contrast builder ---
                   fluidRow(
                     column(3,
                            h5("Add a contrast to compare"),
                            uiOutput("cmpNumeratorUI"),
                            uiOutput("cmpDenominatorUI"),
                            numericInput("cmpLFC", "LFC threshold", value = 0.5, min = 0),
                            numericInput("cmpPadj", "Adj. p-value threshold", value = 0.05, min = 0, max = 1, step = 0.01),
                            actionButton("addContrastBtn", "Add Contrast", icon = icon("plus")),
                            br(), br(),
                            actionButton("clearContrastsBtn", "Clear All Contrasts", icon = icon("trash"))
                     ),
                     column(9,
                            h5("Contrasts queued (including primary contrast from sidebar):"),
                            tableOutput("contrastListTable")
                     )
                   ),
                   hr(),
                   
                   # --- Venn / Euler ---
                   fluidRow(
                     column(12, h4("Venn / Euler Diagram"))
                   ),
                   fluidRow(
                     column(3,
                            selectInput("vennType", "Diagram type", 
                                        choices = c("Venn (ggVennDiagram)" = "venn",
                                                    "Euler (eulerr)" = "euler"),
                                        selected = "venn"),
                            selectInput("vennDEGdirection", "DEG direction",
                                        choices = c("Both (|LFC| threshold)" = "both",
                                                    "Upregulated only" = "up",
                                                    "Downregulated only" = "down"),
                                        selected = "both"),
                            actionButton("runVennBtn", "Find DEGs and Draw Diagram", icon = icon("circle")),
                            br(), br(),
                            downloadButton("dl_venn", "Download PNG")
                     ),
                     column(9,
                            plotOutput("vennPlot")
                     )
                   ),
                   hr(),
                   
                   # --- Gene tables from Venn regions ---
                   fluidRow(
                     column(12, h4("Explore Gene Sets from Diagram"))
                   ),
                   fluidRow(
                     column(4,
                            uiOutput("vennSetSelectorUI"),
                            actionButton("loadVennGenesBtn", "Load genes for this set", icon = icon("table"))
                     ),
                     column(8,
                            DT::dataTableOutput("vennGenesTable"),
                            downloadButton("dl_venn_genes", "Download gene list")
                     )
                   ),
                   hr(),
                   
                   # --- Volcano plots coloured by Venn membership ---
                   fluidRow(
                     column(12, h4("Volcano Plots Coloured by Overlap"))
                   ),
                   fluidRow(
                     column(3,
                            uiOutput("vennVolcContrastUI"),
                            uiOutput("vennVolcContextUI"),
                            sliderInput("cmpVolcTop", "Top genes to label", 1, 80, 20, step = 1),
                            sliderInput("cmpVolcPtSz", "Point size", 0.5, 6, 2, step = 0.5),
                            actionButton("runCmpVolcBtn", "Draw Volcano", icon = icon("chart-line"))
                     ),
                     column(9,
                            plotOutput("cmpVolcanoPlot", height = "450px")
                     )
                   ),
                   hr(),
                   
                   # --- Downstream analyses on a selected gene set ---
                   fluidRow(
                     column(12, h4("Downstream Analysis on Selected Gene Set"))
                   ),
                   fluidRow(
                     column(3,
                            uiOutput("dsSetSelectorUI"),
                            hr(),
                            h5("GO analysis"),
                            selectInput("cmpGOont", "Ontology", c("BP","MF","CC"), selected = "BP"),
                            numericInput("cmpGOnum", "Categories to show", value = 10, min = 1, step = 1),
                            actionButton("runCmpGOBtn", "Run GO", icon = icon("dna")),
                            hr(),
                            h5("GSEA"),
                            selectInput("cmpGseaType", "Pathways to analyse:",
                                        choices = c("GO", "Wiki Pathways", "KEGG"),
                                        selected = "GO"),
                            conditionalPanel(
                              condition = "input.cmpGseaType == 'Wiki Pathways'",
                              selectInput("cmpWikiorg", "Select Wiki Pathways Organism",
                                          choices = c("Homo sapiens",
                                                      "Mus musculus",
                                                      "Drosophila melanogaster",
                                                      "Rattus norvegicus",
                                                      "Danio rerio"),
                                          selected = "Drosophila melanogaster")
                            ),
                            conditionalPanel(
                              condition = "input.cmpGseaType == 'KEGG'",
                              selectInput("cmpKeggorg", "Select KEGG Organism",
                                          choices = c(
                                            "Human (Homo sapiens)" = "hsa",
                                            "Mouse (Mus musculus)" = "mmu",
                                            "Rat (Rattus norvegicus)" = "rno",
                                            "Drosophila (Drosophila melanogaster)" = "dme",
                                            "Zebrafish (Danio rerio)" = "dre"
                                          ),
                                          selected = "dme")
                            ),
                            conditionalPanel(
                              condition = "input.cmpGseaType == 'GO'",
                              selectInput("cmpGseaOnt", "Ontology", c("BP","MF","CC"), selected = "BP")
                            ),
                            selectInput("cmpGseaMetric", "Rank by", c("stat","log2FoldChange"), selected = "stat"),
                            numericInput("cmpGseaP", "p-value cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
                            numericInput("cmpGseaMin", "Min gene set size", value = 10, min = 5, step = 5),
                            numericInput("cmpGseaMax", "Max gene set size", value = 500, min = 50, step = 50),
                            numericInput("cmpGseaNum", "Categories to show", value = 10, min = 1, step = 1),
                            actionButton("runCmpGseaBtn", "Run GSEA", icon = icon("dna")),
                            hr(),
                            h5("GAGE / KEGG"),
                            numericInput("cmpGageP", "q-value cutoff", value = 0.1, min = 0, max = 1, step = 0.01),
                            numericInput("cmpGageMin", "Min geneset size", value = 10, min = 2, step = 1),
                            numericInput("cmpGageMax", "Max geneset size", value = 1000, min = 2, step = 10),
                            selectInput("cmpGageMetric", "Metric", c("log2FoldChange","stat"), selected = "stat"),
                            checkboxInput("cmpGageSameDir", "Same direction (one-sided)", value = TRUE),
                            actionButton("runCmpGageBtn", "Run GAGE", icon = icon("dna"))
                     ),
                     column(9,
                            tabsetPanel(id = "cmpDownstreamTabs",
                                        tabPanel("GO Barplot",
                                                 plotOutput("cmpGOplot", height = "450px"),
                                                 DT::dataTableOutput("cmpGOtable"),
                                                 downloadButton("dl_cmp_go", "Download GO table")
                                        ),
                                        tabPanel("GSEA",
                                                 plotOutput("cmpGseaDotplot", height = "350px"),
                                                 selectizeInput("cmpGseaTerm", "Show enrichment plot for:", choices = NULL),
                                                 plotOutput("cmpGseaEnrichPlot", height = "300px"),
                                                 DT::dataTableOutput("cmpGseaTable"),
                                                 downloadButton("dl_cmp_gsea", "Download GSEA table")
                                        ),
                                        tabPanel("GAGE / Pathview",
                                                 fluidRow(
                                                   
                                                   h5("Upregulated pathways"),
                                                   DT::dataTableOutput("cmpGageUpTable")
                                                 ),
                                                 fluidRow(
                                                   h5("Downregulated pathways"),
                                                   DT::dataTableOutput("cmpGageDownTable")
                                                   
                                                 ),
                                                 hr(),
                                                 fluidRow(
                                                   column(4,
                                                          selectInput("cmpKeggPathway", "Pathway to render:", choices = NULL),
                                                          actionButton("drawCmpPathviewBtn", "Draw Pathview"),
                                                          downloadButton("dl_cmp_pathview", "Download PNG")
                                                   ),
                                                   column(8,
                                                          imageOutput("cmpPathviewImg", height = "500px")
                                                   )
                                                 )
                                        )
                            )
                     )
                   )
          ) # end Compare Contrasts tabPanel
        )
      ),
    )
  )
)

# Define server logic 
server <- function(input, output, session) {
  
  logText <- reactiveVal("")
  
  appendLog <- function(msg) {
    isolate({
      current <- logText()
      logText(paste(current, msg, sep = "\n"))
    })
    # Trigger UI update
    logText(logText())
  }
  
  output$log <- renderText({
    logText()
  })
  
  # --- Helper: auto-detect gene ID type from a sample of IDs ---
  # Returns list(orgdb = "ENSEMBL"/"SYMBOL"/"FLYBASE", biomart = "ensembl_gene_id"/...)
  detect_gene_id_type <- function(gene_ids) {
    ids <- head(gene_ids[!is.na(gene_ids) & nzchar(gene_ids)], 200)
    
    # Check what fraction match known Ensembl patterns
    ens_pattern <- "^ENS[A-Z]*G\\d+"
    fb_pattern  <- "^FBgn\\d+"
    
    frac_ens <- mean(grepl(ens_pattern, ids))
    frac_fb  <- mean(grepl(fb_pattern, ids))
    
    if (frac_ens > 0.5) {
      return(list(orgdb = "ENSEMBL", biomart = "ensembl_gene_id"))
    } else if (frac_fb > 0.5) {
      return(list(orgdb = "FLYBASE", biomart = "flybase_gene_id"))
    } else {
      return(list(orgdb = "SYMBOL", biomart = "external_gene_name"))
    }
  }
  
  # Expose whether the uploaded file is a plain count table (not an RDS)
  # Used by the conditionalPanel to show/hide the gene length file input
  output$isPlainCountTable <- reactive({
    req(input$countsFile)
    tolower(tools::file_ext(input$countsFile$name)) != "rds"
  })
  outputOptions(output, "isPlainCountTable", suspendWhenHidden = FALSE)
  
  # Optional gene length table reactive (only relevant for plain count tables)
  gene_lengths_table <- reactive({
    req(input$geneLengthFile)
    path <- input$geneLengthFile$datapath
    ext  <- tolower(tools::file_ext(input$geneLengthFile$name))
    sep  <- if (ext == "csv") "," else "\t"
    tbl  <- read.delim(path, sep = sep, header = TRUE,
                       stringsAsFactors = FALSE, check.names = FALSE)
    # Accept flexible column names: first column = GeneID, second = Length
    # but also honour explicit "GeneID" / "Length" headers
    if (!all(c("GeneID", "Length") %in% colnames(tbl))) {
      colnames(tbl)[1] <- "GeneID"
      colnames(tbl)[2] <- "Length"
    }
    tbl <- tbl[, c("GeneID", "Length")]
    tbl$Length <- as.numeric(tbl$Length)
    validate(need(nrow(tbl) > 0, "Gene length file appears to be empty."))
    tbl
  })
  
  # Load uploaded count data — supports featureCounts .rds OR plain count tables (.csv/.tsv/.txt)
  counts_data <- reactive({
    req(input$countsFile)
    path <- input$countsFile$datapath
    ext  <- tolower(tools::file_ext(input$countsFile$name))
    
    if (ext == "rds") {
      # Original featureCounts RDS path — return as-is
      obj <- readRDS(path)
      # Validate it looks like a featureCounts object
      validate(
        need(!is.null(obj$counts),      "RDS file does not contain a $counts matrix."),
        need(!is.null(obj$annotation),  "RDS file does not contain an $annotation data frame.")
      )
      return(obj)
    }
    
    # --- Plain count table or featureCounts text output (.csv / .tsv / .txt) ---
    sep <- if (ext == "csv") "," else "\t"
    
    # Read the full table, skipping comment lines (featureCounts text output
    # starts with a '# Program:featureCounts ...' header line).
    raw <- read.delim(path, sep = sep, header = TRUE,
                      comment.char = "#", check.names = FALSE,
                      stringsAsFactors = FALSE)
    
    validate(need(nrow(raw) > 0 && ncol(raw) > 0,
                  "Count table appears to be empty. Check file format and delimiter."))
    
    # --- Detect featureCounts text output ---
    # The standard featureCounts summary has columns:
    #   Geneid, Chr, Start, End, Strand, Length, <sample1>, <sample2>, ...
    # With --extraAttributes it may also include gene_id, gene_type, etc.
    # In all cases the first column is 'Geneid' and there are non-numeric
    # annotation columns before the integer count columns.
    #
    # Strategy: identify which columns are numeric (counts) vs character
    # (annotation metadata) and split them apart.
    
    # Known featureCounts annotation column names (case-insensitive match)
    fc_annot_names <- c("geneid", "chr", "start", "end", "strand", "length",
                        "gene_id", "gene_type", "gene_name", "gene_biotype")
    
    # Use the first column as gene IDs regardless of its name
    gene_ids <- raw[[1]]
    
    # Classify remaining columns as annotation (character/non-numeric) or count (numeric)
    remaining <- raw[, -1, drop = FALSE]
    
    is_count_col <- vapply(remaining, function(col) {
      # A count column should be coercible to integer with no (or very few) NAs
      suppressWarnings({
        nums <- as.numeric(col)
      })
      # If >90% parse as finite numbers, treat as a count column
      sum(is.finite(nums)) / length(nums) > 0.9
    }, logical(1))
    
    annot_cols <- remaining[, !is_count_col, drop = FALSE]
    count_cols <- remaining[,  is_count_col, drop = FALSE]
    
    validate(need(ncol(count_cols) > 0,
                  "No numeric count columns detected. Check file format."))
    
    # Build the count matrix
    mat <- as.matrix(count_cols)
    mode(mat) <- "integer"
    rownames(mat) <- gene_ids
    
    # Log what was detected
    if (ncol(annot_cols) > 0) {
      showNotification(paste0("Detected featureCounts text format. ",
                              "Stripped ", ncol(annot_cols), " annotation column(s): ",
                              paste(names(annot_cols), collapse = ", "), "."),
                       type = "message", duration = 8)
    }
    
    # Handle duplicate gene IDs (can happen when -g gene_name is used)
    if (anyDuplicated(rownames(mat))) {
      n_dup <- sum(duplicated(rownames(mat)))
      showNotification(paste0(n_dup, " duplicate gene IDs found — summing counts."),
                       type = "warning", duration = 8)
      mat <- rowsum(mat, rownames(mat))
      gene_ids <- rownames(mat)
    }
    
    # --- Build annotation data frame ---
    annotation_df <- data.frame(
      GeneID = rownames(mat),
      stringsAsFactors = FALSE
    )
    
    # Extract gene length from the featureCounts 'Length' column if present
    fc_length_extracted <- FALSE
    if ("Length" %in% names(annot_cols)) {
      len_vec <- as.numeric(annot_cols$Length)
      # If there were duplicate gene IDs, average the lengths per unique ID
      if (length(len_vec) != nrow(annotation_df)) {
        len_df <- data.frame(GeneID = gene_ids, Length = len_vec,
                             stringsAsFactors = FALSE)
        len_df <- aggregate(Length ~ GeneID, data = len_df, FUN = mean, na.rm = TRUE)
        annotation_df <- dplyr::left_join(annotation_df, len_df, by = "GeneID")
      } else {
        annotation_df$Length <- len_vec
      }
      fc_length_extracted <- TRUE
      showNotification("Extracted gene lengths from featureCounts 'Length' column.",
                       type = "message", duration = 6)
    }
    
    # If user uploaded a separate gene-length file, it takes priority
    if (!is.null(input$geneLengthFile)) {
      len_tbl <- gene_lengths_table()
      # Drop any existing Length column to replace with the user-supplied one
      annotation_df$Length <- NULL
      annotation_df <- dplyr::left_join(annotation_df, len_tbl, by = "GeneID")
      n_missing <- sum(is.na(annotation_df$Length))
      if (n_missing > 0) {
        showNotification(
          paste0(n_missing, " gene(s) in the count table had no matching entry in the ",
                 "gene length file and will have NA lengths."),
          type = "warning", duration = 8
        )
      }
    } else if (!fc_length_extracted) {
      annotation_df$Length <- NA_real_
      # Note: gene length status is shown via geneLengthStatusUI in the sidebar
    }
    
    list(counts = mat, annotation = annotation_df)
  })
  
  # Gene length status indicator
  output$geneLengthStatusUI <- renderUI({
    req(counts_data_with_lengths())
    lengths <- counts_data_with_lengths()$annotation$Length
    has_lengths <- !all(is.na(lengths))
    n_with <- sum(!is.na(lengths))
    n_total <- length(lengths)
    
    if (has_lengths) {
      tags$div(
        style = "padding: 6px 10px; margin-bottom: 8px; border-radius: 4px; background-color: #d4edda; border: 1px solid #c3e6cb;",
        icon("check-circle", style = "color: #28a745;"),
        tags$b(paste0("Gene lengths detected: ", n_with, "/", n_total, " genes.")),
        tags$br(),
        tags$small("Length-bias correction will be applied in GO analysis.")
      )
    } else {
      tags$div(
        style = "padding: 6px 10px; margin-bottom: 8px; border-radius: 4px; background-color: #fff3cd; border: 1px solid #ffc107;",
        icon("exclamation-triangle", style = "color: #856404;"),
        tags$b("No gene lengths found in count file."),
        tags$br(),
        tags$small("GO analysis will run without length-bias correction.",
                   "Fetch lengths via biomaRt or upload a gene length file.")
      )
    }
  })
  
  # ---- biomaRt gene length fetching ----
  biomartLengths <- reactiveVal(NULL)
  
  # Wrapper reactive: overlays biomaRt-fetched lengths onto counts_data() when
  # the file itself has no lengths and no user-uploaded length file
  counts_data_with_lengths <- reactive({
    cd <- counts_data()
    bm <- biomartLengths()
    
    # Only apply biomaRt lengths if the annotation currently has no lengths
    # (user-uploaded file or featureCounts Length column take priority)
    if (!is.null(bm) && all(is.na(cd$annotation$Length))) {
      ann <- cd$annotation[, "GeneID", drop = FALSE]
      ann <- dplyr::left_join(ann, bm, by = "GeneID")
      # Fill any unmatched genes with NA
      if (!"Length" %in% names(ann)) ann$Length <- NA_real_
      cd$annotation <- ann
    }
    cd
  })
  
  observeEvent(input$fetchBiomartBtn, {
    req(counts_data())
    
    withProgress(message = "Connecting to Ensembl biomaRt...", value = 0, {
      incProgress(0.05)
      
      # Map organism selection to biomaRt dataset
      orgdb_name <- isolate(input$organism)
      bm_dataset <- switch(orgdb_name,
                           "org.Hs.eg.db" = "hsapiens_gene_ensembl",
                           "org.Mm.eg.db" = "mmusculus_gene_ensembl",
                           "org.Rn.eg.db" = "rnorvegicus_gene_ensembl",
                           "org.Dm.eg.db" = "dmelanogaster_gene_ensembl",
                           "org.Dr.eg.db" = "drerio_gene_ensembl"
      )
      
      # Determine which biomaRt attribute matches the gene IDs in the count table
      gene_ids <- unique(counts_data()$annotation$GeneID)
      gene_ids <- gene_ids[!is.na(gene_ids) & nzchar(gene_ids)]
      id_type <- detect_gene_id_type(gene_ids)
      bm_gene_attr <- id_type$biomart
      
      tryCatch({
        setProgress(value = 0.05, message = "Connecting to Ensembl...")
        ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = bm_dataset)
        
        # --- Chunked query with progress updates ---
        chunk_size <- 5000
        chunks <- split(gene_ids, ceiling(seq_along(gene_ids) / chunk_size))
        n_chunks <- length(chunks)
        bm_results_list <- vector("list", n_chunks)
        
        for (i in seq_along(chunks)) {
          setProgress(
            value   = 0.1 + 0.8 * (i / n_chunks),
            message = paste0("Querying Ensembl... batch ", i, "/", n_chunks,
                             " (", length(chunks[[i]]), " genes)")
          )
          bm_results_list[[i]] <- biomaRt::getBM(
            attributes = c(bm_gene_attr, "transcript_length"),
            filters    = bm_gene_attr,
            values     = chunks[[i]],
            mart       = ensembl
          )
        }
        
        bm_results <- do.call(rbind, bm_results_list)
        
        setProgress(value = 0.92, message = "Computing median lengths per gene...")
        
        if (is.null(bm_results) || nrow(bm_results) == 0) {
          showNotification("biomaRt returned no results. Check organism selection.",
                           type = "error", duration = 8)
          return()
        }
        
        # Use median transcript length per gene as the length estimate
        colnames(bm_results) <- c("GeneID", "transcript_length")
        gene_lengths <- aggregate(transcript_length ~ GeneID, data = bm_results, FUN = median)
        colnames(gene_lengths) <- c("GeneID", "Length")
        gene_lengths$Length <- round(gene_lengths$Length)
        
        n_matched <- sum(gene_ids %in% gene_lengths$GeneID)
        n_total   <- length(gene_ids)
        
        biomartLengths(gene_lengths)
        
        setProgress(value = 1, message = "Done!")
        showNotification(
          paste0("biomaRt: fetched lengths for ", nrow(gene_lengths),
                 " genes (", n_matched, "/", n_total, " matched to count data)."),
          type = "message", duration = 8
        )
        
      }, error = function(e) {
        showNotification(
          paste0("biomaRt fetch failed: ", conditionMessage(e)),
          type = "error", duration = 10
        )
      })
    })
  })
  
  # ---- Manual group builder ----
  manualGroups <- reactiveVal(data.frame(
    Group = character(0), n = integer(0), stringsAsFactors = FALSE
  ))
  manualConfirmed <- reactiveVal(FALSE)
  
  observeEvent(input$addGroupBtn, {
    req(input$manualGroupName)
    gname <- trimws(input$manualGroupName)
    if (nchar(gname) == 0) {
      showNotification("Please enter a group name.", type = "error")
      return()
    }
    current <- manualGroups()
    if (gname %in% current$Group) {
      showNotification(paste("Group", gname, "already added."), type = "warning")
      return()
    }
    manualGroups(rbind(current, data.frame(Group = gname, n = input$manualGroupN,
                                           stringsAsFactors = FALSE)))
    manualConfirmed(FALSE)
    updateTextInput(session, "manualGroupName", value = "")
  })
  
  observeEvent(input$clearGroupsBtn, {
    manualGroups(data.frame(Group = character(0), n = integer(0),
                            stringsAsFactors = FALSE))
    manualConfirmed(FALSE)
  })
  
  observeEvent(input$confirmGroupsBtn, {
    mg <- manualGroups()
    if (nrow(mg) == 0) {
      showNotification("Add at least one group first.", type = "error")
      return()
    }
    total_n <- sum(mg$n)
    n_cols <- tryCatch(ncol(counts_data()$counts), error = function(e) NA)
    if (!is.na(n_cols) && total_n != n_cols) {
      showNotification(
        paste0("Total samples (", total_n, ") does not match the number of ",
               "count columns (", n_cols, "). Please adjust."),
        type = "error", duration = 8
      )
      return()
    }
    manualConfirmed(TRUE)
    showNotification(paste0("Groups confirmed: ", total_n, " samples across ",
                            nrow(mg), " group(s)."), type = "message")
  })
  
  output$manualGroupsTable <- renderTable({
    mg <- manualGroups()
    if (nrow(mg) == 0) return(NULL)
    # Show preview of sample names
    mg$`Sample names` <- vapply(seq_len(nrow(mg)), function(i) {
      paste0(mg$Group[i], seq_len(mg$n[i]), collapse = ", ")
    }, character(1))
    mg
  })
  
  output$manualSampleCountMsg <- renderUI({
    mg <- manualGroups()
    if (nrow(mg) == 0) return(NULL)
    total_n <- sum(mg$n)
    n_cols <- tryCatch(ncol(counts_data()$counts), error = function(e) NA)
    if (!is.na(n_cols)) {
      col <- if (total_n == n_cols) "green" else "red"
      tags$p(style = paste0("color:", col, "; font-weight:bold;"),
             paste0("Total samples: ", total_n, " / ", n_cols, " count columns"))
    } else {
      tags$p(paste0("Total samples: ", total_n, " (upload count data to verify)"))
    }
  })
  
  # Build sample_info data frame from manual groups
  manual_sample_info <- reactive({
    req(manualConfirmed())
    mg <- manualGroups()
    req(nrow(mg) > 0)
    
    rows <- lapply(seq_len(nrow(mg)), function(i) {
      data.frame(
        SampleName = paste0(mg$Group[i], seq_len(mg$n[i])),
        Group      = rep(mg$Group[i], mg$n[i]),
        stringsAsFactors = TRUE
      )
    })
    do.call(rbind, rows)
  })
  
  # Unified sample_info reactive — works in both modes
  sample_info <- reactive({
    if (input$sampleInfoMode == "upload") {
      req(input$sampleInfoFile)
      
      path <- input$sampleInfoFile$datapath
      
      # Try to read the file, catching malformed input
      si <- tryCatch({
        read.delim(path, stringsAsFactors = TRUE)
      }, error = function(e) {
        # Fall back: try comma-separated in case user uploaded .csv-style
        tryCatch({
          read.csv(path, stringsAsFactors = TRUE)
        }, error = function(e2) {
          showNotification(
            paste("Could not read sample info file:", e$message,
                  "\nExpected a tab- or comma-delimited file with a header row."),
            type = "error", duration = NULL
          )
          NULL
        })
      })
      
      validate(need(!is.null(si), "Sample info file could not be parsed. Check the format."))
      validate(need(ncol(si) >= 2,
                    paste0("Sample info has only ", ncol(si), " column(s). ",
                           "Expected at least 2 (e.g. SampleName, Group). ",
                           "Check the delimiter — the file should be tab-separated.")))
      validate(need("SampleName" %in% colnames(si),
                    paste0("Column 'SampleName' not found. Found columns: ",
                           paste(colnames(si), collapse = ", "),
                           ". Check header names and delimiter.")))
      validate(need("Group" %in% colnames(si),
                    paste0("Column 'Group' not found. Found columns: ",
                           paste(colnames(si), collapse = ", "),
                           ". Check header names and delimiter.")))
      
      si
      
    } else {
      req(manualConfirmed())
      manual_sample_info()
    }
  })
  
  # ---- Sample filter UI + filtered reactive ----
  output$sampleFilterUI <- renderUI({
    req(sample_info())
    sn <- as.character(sample_info()$SampleName)
    checkboxGroupInput("selectedSamples", "Include samples:",
                       choices = sn, selected = sn)
  })
  
  filtered_sample_info <- reactive({
    si <- sample_info()
    req(si)
    sel <- input$selectedSamples
    if (is.null(sel)) return(si)
    out <- si[as.character(si$SampleName) %in% sel, , drop = FALSE]
    validate(need(nrow(out) >= 2, "Please keep at least 2 samples selected."))
    out
  })
  
  # Tell UI when both files are ready
  output$inputsReady <- reactive({
    has_counts <- !is.null(input$countsFile)
    has_info   <- if (input$sampleInfoMode == "upload") {
      !is.null(input$sampleInfoFile)
    } else {
      isTRUE(manualConfirmed())
    }
    has_counts && has_info
  })
  outputOptions(output, "inputsReady", suspendWhenHidden = FALSE)
  
  # Reactive to get unique group names
  group_levels <- reactive({
    req(filtered_sample_info())
    unique(as.character(filtered_sample_info()$Group))
  })
  
  observeEvent(filtered_sample_info(), {
    updateSelectInput(session, "groupOrder",
                      choices = group_levels(),
                      selected = group_levels())
  })
  
  # UI to specify the order of those groups
  output$groupOrderUI <- renderUI({
    req(group_levels())
    selectInput("groupOrder", "Specify group order:",
                choices = group_levels(),
                selected = group_levels(),
                multiple = TRUE)
  })
  
  # Populate contrast options dynamically
  
  output$contrastSelectUI <- renderUI({
    req(input$groupOrder)
    
    # Ensure at least 2 levels to contrast
    if (length(input$groupOrder) < 2) {
      return(h4("Please select at least 2 groups."))
    }
    
    tagList(
      selectInput("contrastNumerator", "Compare (numerator):", choices = input$groupOrder),
      selectInput("contrastDenominator", "Against (denominator):", choices = input$groupOrder)
    )
  })
  
  observeEvent(input$groupOrder, {
    updateSelectInput(session, "contrastNumerator",
                      choices = input$groupOrder,
                      selected = tail(input$groupOrder, 1))  # or whatever default
    updateSelectInput(session, "contrastDenominator",
                      choices = input$groupOrder,
                      selected = head(input$groupOrder, 1))
  })
  
  # Preview count data input
  output$previewTable <- renderTable({
    counts <- counts_data()$counts
    # For featureCounts RDS, columns are BAM paths — rename to sample names.
    # For plain count tables, columns are already sample names; only rename if lengths match.
    if (ncol(counts) == nrow(sample_info())) {
      colnames(counts) <- sample_info()$SampleName
    }
    # Subset to selected samples for preview
    fsi <- filtered_sample_info()
    keep_cols <- colnames(counts) %in% as.character(fsi$SampleName)
    if (any(keep_cols)) counts <- counts[, keep_cols, drop = FALSE]
    head_df <- head(counts)
    head_df <- cbind(GeneID = rownames(head_df), head_df)
    head_df
  })
  
  # Preview sample info
  output$sampleInfo <- renderTable({
    filtered_sample_info()
  })
  
  output$convertToSymbolUI <- renderUI({
    req(counts_data())
    gene_ids <- counts_data()$annotation$GeneID
    id_type <- detect_gene_id_type(gene_ids)
    
    if (id_type$orgdb != "SYMBOL") {
      id_label <- id_type$orgdb  # "ENSEMBL" or "FLYBASE"
      tagList(
        checkboxInput("convertToSymbol",
                      paste0("Convert ", id_label, " IDs to gene symbols?"),
                      value = TRUE),
        helpText(paste0("Detected ", id_label, " identifiers. Ticking this will map them to gene symbols for readability."))
      )
    }
  })
  
  
  analysisResults <- eventReactive(input$analyzeBtn, {
    withProgress(message = "Running RNA-seq analysis", value = 0, {
      orgdb_name <- isolate(input$organism)
      orgdb <- isolate(get(orgdb_name))
      
      fc <- counts_data_with_lengths()
      
      # Decide what your GeneID column represents
      id_type <- detect_gene_id_type(fc$annotation$GeneID)
      from_type <- id_type$orgdb
      
      # Map -> ENTREZ, keep it clean and 1:1 on the input id
      symbol_to_entrez <- clusterProfiler::bitr(
        fc$annotation$GeneID,
        fromType = from_type,
        toType   = "ENTREZID",
        OrgDb    = orgdb
      ) |>
        dplyr::filter(!is.na(ENTREZID)) |>
        dplyr::distinct(.data[[from_type]], .keep_all = TRUE)
      
      # Join using the correct key name for the mapping table
      fc$annotation <- dplyr::left_join(
        fc$annotation,
        symbol_to_entrez,
        by = c("GeneID" = from_type)
      )
      
      # Proceed
      counts <- counts_data()$counts
      # For featureCounts RDS, column names are BAM paths — rename to sample names.
      # For plain count tables, columns are already sample names; only rename if lengths match.
      if (ncol(counts) == nrow(sample_info())) {
        colnames(counts) <- sample_info()$SampleName
      }
      col_data <- filtered_sample_info()
      
      # Subset counts to only selected samples
      keep_samples <- colnames(counts) %in% as.character(col_data$SampleName)
      counts <- counts[, keep_samples, drop = FALSE]
      
      # set group levels
      col_data$Group <- factor(col_data$Group, levels = isolate(input$groupOrder))
      
      appendLog("Running DESeq2 analysis...")
      showNotification("Running DESeq2 analysis...", type="message")
      incProgress(0.1)
      
      
      # Create DESEq2 data set
      appendLog("Creating DESeqDataSet...")
      showNotification("Creating DESeqDataSet...", type="message")
      dds <- DESeqDataSetFromMatrix(
        countData = counts,
        colData = col_data,
        design = ~ Group
      )
      
      # Convert non-symbol identifiers to gene symbols if checkbox ticked
      if (isTRUE(isolate(input$convertToSymbol)) && from_type != "SYMBOL") {
        appendLog(paste0("Converting ", from_type, " IDs to gene symbols..."))
        showNotification(paste0("Converting ", from_type, " → gene symbols..."), type = "message")
        
        gene_symbols <- mapIds(
          orgdb,
          keys = rownames(dds),
          column = "SYMBOL",
          keytype = from_type,
          multiVals = "first"
        )
        
        # Drop any rows with NA gene symbols
        keep_idx <- !is.na(gene_symbols)
        n_dropped <- sum(!keep_idx)
        dds <- dds[keep_idx, ]
        rownames(dds) <- gene_symbols[keep_idx]
        
        if (n_dropped > 0) {
          appendLog(paste0("Dropped ", n_dropped, " genes with no symbol mapping."))
          showNotification(paste0("Dropped ", n_dropped, " genes with no symbol mapping."),
                           type = "warning", duration = 6)
        }
        appendLog(paste0("Converted to symbols. ", nrow(dds), " genes remaining."))
      }
      
      if (isTRUE(isolate(input$filterLncRNA))) {
        appendLog("Filtering out lncRNAs (keeping protein-coding only)...")
        showNotification("Filtering out lncRNAs...", type="message")
        
        # Determine the current ID format from the dds rownames
        # (may have changed after FlyBase → Symbol conversion above)
        current_keytype <- detect_gene_id_type(rownames(dds))$orgdb
        
        # Map rownames to their Gene Type using the selected OrgDb
        gene_types <- mapIds(
          orgdb,
          keys = rownames(dds),
          column = "GENETYPE",
          keytype = current_keytype,
          multiVals = "first"
        )
        
        # Keep ONLY genes explicitly annotated as protein-coding
        keep_coding <- which(gene_types == "protein-coding")
        
        if (length(keep_coding) > 0) {
          dds <- dds[keep_coding, ]
          appendLog(paste("Kept", length(keep_coding), "protein-coding genes."))
        } else {
          appendLog("Warning: No protein-coding genes found. Skipping lncRNA filter.")
          showNotification("Warning: No protein-coding genes found. Skipping filter.", type="warning")
        }
      }
      incProgress(0.1)
      
      filt <- isTRUE(input$filterlow)
      
      if (filt) {
        # filter our low counts
        appendLog("Filtering low-expression genes...")
        showNotification("Filtering low-expression genes...", type="message")
        keep <- edgeR::filterByExpr(counts(dds), group = dds$Group)
        
        dds <- dds[keep,]
      }
      
      incProgress(0.1)
      
      appendLog("Running DESeq()...")
      showNotification("Running DESeq()...", type="message")
      dds <- DESeq(dds)
      
      incProgress(0.2)
      
      appendLog("Getting results...")
      showNotification("Getting results...", type="message")
      
      res <- results(
        dds,
        lfcThreshold = isolate(input$lfcThreshold),
        altHypothesis = "greaterAbs",
        alpha = isolate(input$padjThreshold),
        independentFiltering = isolate(input$deseqfilter),
        contrast = c("Group", isolate(input$contrastNumerator), isolate(input$contrastDenominator))
      )
      
      message("filterThreshold = ", metadata(res)$filterThreshold)
      
      # Match your UI thresholds in the summary
      sum_txt <- paste(
        capture.output(summary(res, alpha = isolate(input$padjThreshold))),
        collapse = "\n"
      )
      appendLog("DESeq2 results summary:")
      appendLog(sum_txt)  
      
      incProgress(0.1)
      
      appendLog("Computing vst...")
      showNotification("Computing vst...", type="message")
      
      vsd <- vst(dds, blind = TRUE)
      
      incProgress(0.1)
      
      # Create res object with rownames being ENTREZIDs
      
      # Get gene symbols from DESeq2 results (res)
      uni_gene_symbols <- rownames(res)  # Get gene symbols from DESeq2 results
      
      # Map gene symbols to Entrez IDs using bitr
      
      # Map gene IDs from DESeq2 results to Entrez IDs
      # Detect from rownames(res) since these may differ from original GeneIDs
      # (e.g. after FlyBase → Symbol conversion)
      from_type <- detect_gene_id_type(rownames(res))$orgdb
      
      
      uni_entrez_ids <- clusterProfiler::bitr(
        uni_gene_symbols,
        fromType = from_type,
        toType   = "ENTREZID",
        OrgDb    = orgdb
      )
      
      # Merge the DESeq2 results (res) with the Entrez IDs 
      res_entrez <- merge(as.data.frame(res), uni_entrez_ids, by.x = "row.names", by.y = from_type)
      
      # rename the first column to "GeneIDs"
      res_entrez <- res_entrez %>% rename(GeneIDs = Row.names)
      
      # Add gene length information
      
      m <- match(res_entrez$ENTREZID, fc$annotation$ENTREZID)
      gene_lengths <- as.numeric(fc$annotation$Length[m])
      
      # Step 2: Create a named vector of gene lengths, where the names are the ENTREZ IDs
      names(gene_lengths) <- res_entrez$ENTREZID
      
      
      # Step 3: Add the gene_length column back to your 'res_entrez' data frame if needed
      res_entrez <- cbind(res_entrez, gene_lengths)
      
      appendLog("PCA analysis...")
      showNotification("PCA analysis...", type="message")
      
      # PCA on vst
      mat <- assay(vsd)
      pcDat <- prcomp(t(mat), center = TRUE, scale. = FALSE)
      percentVar <- round(100 * (pcDat$sdev^2 / sum(pcDat$sdev^2)), 1)
      
      # Convert PCA results to a dataframe
      pca_df <- as.data.frame(pcDat$x)
      pca_df$SampleName <- filtered_sample_info()$SampleName  # Add sample names
      pca_df$Group <- filtered_sample_info()$Group  # Add group info
      
      
      incProgress(0.2)
      
      
      appendLog("Analysis complete.")
      showNotification("Analysis complete.", type="message")
      
    })
    # Output from reactive expression
    list(dds = dds, res = res, vsd = vsd, pca_df = pca_df, percentVar = percentVar,
         res_entrez = res_entrez)
    
  })
  
  
  # PCA plot
  output$pcaPlot <- renderPlot({
    req(analysisResults())
    # Recall the variables from the reactive analysisResults event
    pca_df <- analysisResults()$pca_df
    percentVar <- analysisResults()$percentVar
    
    ggplot(pca_df, aes(x = PC1, y = PC2, fill = Group, shape = Group)) +
      geom_point(size = input$pointsize) +
      geom_text_repel(aes(label = SampleName), size = input$labelsize,
                      box.padding = 1,      # Increases space around labels
                      point.padding = 1,    # Increases distance from points
                      min.segment.length = 0,
                      segment.size = if (input$labelsize == 0) 0 else 0.5) +
      scale_shape_manual(values = c(21:25, 21:25)) + # Change this depending on how many shapes you want, 21-25 are decent
      guides(fill = guide_legend(override.aes = list(shape = 22))) +
      theme_minimal() +
      theme(
        legend.title = element_text(face = "bold", size = input$legendtitlesize),
        legend.text = element_text(size = input$legendtextsize),
        axis.text = element_text(face = "bold", size = input$axistextsize),
        axis.title = element_text(face = "bold", size = input$axistitlesize)
      )+
      labs(
        title = "PCA Plot of RNA-seq Samples",
        x = paste0("PC1 (", percentVar[1], "% variance)"),
        y = paste0("PC2 (", percentVar[2], "% variance)")
      )
  })
  
  # DE results table
  output$deTable <- renderDataTable({
    req(analysisResults())
    datatable(as.data.frame(analysisResults()$res))
  })
  
  # Update heatmap group choices when analysis results are available
  observeEvent(analysisResults(), {
    req(analysisResults())
    available_groups <- unique(as.character(colData(analysisResults()$dds)$Group))
    updateSelectInput(session, "heatmapGroups",
                      choices = available_groups,
                      selected = available_groups)  # Default: all groups selected
  })
  
  # Heatmap plot
  output$heatmapPlot <- renderPlot({
    req(analysisResults())
    
    res <- analysisResults()$res
    vsd <- analysisResults()$vsd
    
    # Get expression data
    mat <- assay(vsd)
    
    # Select genes based on chosen source
    if (input$heatmapGeneSource == "de") {
      # --- Top DE genes (contrast-dependent) ---
      res_df <- as.data.frame(res)
      res_df$abs_lfc <- abs(res_df$log2FoldChange)
      res_df$abs_stat <- abs(res_df$stat)
      
      top_genes <- switch(input$heatmapSortBy,
                          "padj"     = head(rownames(res_df[order(res_df$padj, na.last = TRUE), ]), input$heatmapTopN),
                          "abs_lfc"  = head(rownames(res_df[order(res_df$abs_lfc, decreasing = TRUE, na.last = TRUE), ]), input$heatmapTopN),
                          "abs_stat" = head(rownames(res_df[order(res_df$abs_stat, decreasing = TRUE, na.last = TRUE), ]), input$heatmapTopN)
      )
      
      heatmap_title <- paste("Top", length(top_genes[!is.na(top_genes)]), "DE genes —",
                             switch(input$heatmapSortBy,
                                    "padj"     = "sorted by adjusted p-value",
                                    "abs_lfc"  = "sorted by absolute log2 fold change",
                                    "abs_stat" = "sorted by absolute test statistic"))
      
    } else {
      # --- Most variable genes (contrast-independent) ---
      row_vars <- apply(mat, 1, var)
      row_vars <- sort(row_vars, decreasing = TRUE)
      n_genes  <- min(input$heatmapTopN, length(row_vars))
      top_genes <- names(row_vars)[seq_len(n_genes)]
      
      heatmap_title <- paste("Top", n_genes, "most variable genes (by variance across all samples)")
    }
    
    # Remove any NA gene names
    top_genes <- top_genes[!is.na(top_genes)]
    
    # Subset expression matrix
    heatmap_mat <- mat[top_genes, , drop = FALSE]
    
    # Filter samples by selected groups
    if (!is.null(input$heatmapGroups) && length(input$heatmapGroups) > 0) {
      selected_samples <- colData(vsd)$Group %in% input$heatmapGroups
      heatmap_mat <- heatmap_mat[, selected_samples, drop = FALSE]
      
      # Update sample annotation to match filtered samples
      sample_annotation <- data.frame(
        Group = colData(vsd)$Group[selected_samples],
        row.names = colnames(heatmap_mat)
      )
    } else {
      # If no groups selected, show message
      validate(need(length(input$heatmapGroups) > 0, "Please select at least one group to display."))
    }
    
    
    # Set up colors
    if (input$heatmapColorScheme %in% c("viridis", "plasma")) {
      colors <- if (input$heatmapColorScheme == "viridis") {
        viridisLite::viridis(100)
      } else {
        viridisLite::plasma(100)
      }
    } else {
      colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, input$heatmapColorScheme)))(100)
      
    }
    
    # Create annotation colors - handle any number of groups
    n_groups <- length(unique(sample_annotation$Group))
    if (n_groups <= 8) {
      group_colors <- RColorBrewer::brewer.pal(max(n_groups, 3), "Set2")[seq_len(n_groups)]
    } else {
      group_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_groups)
    }
    names(group_colors) <- unique(sample_annotation$Group)
    annotation_colors <- list(Group = group_colors)
    
    validate(need(nrow(heatmap_mat) >= 2 && ncol(heatmap_mat) >= 2,
                  "Need at least 2 genes and 2 samples to draw the heatmap."))
    
    # Generate heatmap
    pheatmap::pheatmap(
      heatmap_mat,
      scale = input$heatmapScale,
      clustering_distance_rows = "correlation",
      clustering_distance_cols = "correlation",
      cluster_rows = input$heatmapClusterRows,
      cluster_cols = input$heatmapClusterCols,
      show_rownames = input$heatmapShowRowNames,
      show_colnames = input$heatmapShowColNames,
      annotation_col = sample_annotation,
      annotation_colors = annotation_colors,
      color = colors,
      main = heatmap_title,
      fontsize = 10,
      fontsize_row = if (input$heatmapShowRowNames && length(top_genes) > 50) 6 else 8,
      fontsize_col = if (input$heatmapShowColNames) 8 else 8
    )
  }, res = 96)
  
  
  ###---GSEA---###
  
  gseaResults <- eventReactive(input$runGSEA, {
    withProgress(message = "Running GSEA analysis", value = 0, {
      appendLog("Running GSEA analysis...")
      showNotification("Running GSEA analysis...", type="message")
      incProgress(0.2)
      req(analysisResults())
      ba <- analysisResults()
      
      
      orgdb_name <- isolate(input$organism)
      orgdb <- isolate(get(orgdb_name))
      
      # 1) Choose the ranking metric (use unshrunk 'stat' if available)
      metric_col <- switch(input$gseaMetric,
                           stat = "stat",
                           log2FoldChange = "log2FoldChange")
      
      validate(need(metric_col %in% names(ba$res_entrez),
                    sprintf("Column '%s' not found in DE results.", metric_col)))
      
      m <- ba$res_entrez[[metric_col]]
      names(m) <- ba$res_entrez$ENTREZID
      
      # 2) Clean: drop NA/Inf, collapse duplicates (keep max), sort
      ok <- is.finite(m)
      m <- m[ok]
      
      # tapply returns an array; coerce back to a plain named numeric vector
      m <- tapply(m, names(m), max)                  # collapse dup ENTREZ to max
      m <- sort(m, decreasing = TRUE)
      m <- setNames(as.vector(m), names(m))          # <-- ensure class "numeric" w/ names
      
      validate(need(length(m) >= input$gseaMin,
                    "Not enough ranked genes to run GSEA."))
      
      # Run the appropriate gse function — each branch wrapped in tryCatch
      # so the eventReactive ALWAYS returns a value (result or NULL).
      
      gse <- NULL
      
      if (input$gsea_type == "GO") {
        gse <- tryCatch({
          clusterProfiler::gseGO(
            geneList     = m,
            OrgDb        = orgdb,
            keyType      = "ENTREZID",
            ont          = input$gseaOnt,
            minGSSize    = input$gseaMin,
            maxGSSize    = input$gseaMax,
            pvalueCutoff = input$gseaP,
            verbose      = TRUE
          )
        }, error = function(e) {
          showNotification(paste("GO GSEA failed:", e$message),
                           type = "error", duration = NULL)
          NULL
        })
        
      } else if (input$gsea_type == "Wiki Pathways") {
        validate(need(nzchar(input$wikiorg),
                      "Please select an organism for WikiPathways."))
        gse <- tryCatch({
          clusterProfiler::gseWP(
            geneList     = m,
            organism     = input$wikiorg,
            minGSSize    = input$gseaMin,
            maxGSSize    = input$gseaMax,
            pvalueCutoff = input$gseaP,
            verbose      = TRUE
          )
        }, error = function(e) {
          showNotification(paste("WikiPathways GSEA failed:", e$message),
                           type = "error", duration = NULL)
          NULL
        })
        
      } else if (input$gsea_type == "KEGG") {
        gse <- tryCatch({
          clusterProfiler::gseKEGG(
            geneList     = m,
            organism     = input$keggorg,
            keyType      = "ncbi-geneid",
            minGSSize    = input$gseaMin,
            maxGSSize    = input$gseaMax,
            pvalueCutoff = input$gseaP,
            verbose      = TRUE
          )
        }, error = function(e) {
          showNotification(paste("KEGG GSEA failed:", e$message),
                           type = "error", duration = NULL)
          NULL
        })
      }
      
      incProgress(0.8)
      
      if (!is.null(gse)) {
        appendLog("GSEA analysis complete.")
        showNotification("GSEA analysis complete.", type = "message")
      } else {
        appendLog("GSEA analysis failed — see notification for details.")
      }
      
      # Always return gse (possibly NULL) so the reactive has a value
      gse
    })
  })
  
  observeEvent(gseaResults(), {
    gse <- gseaResults()
    if (is.null(gse)) return()
    df <- as.data.frame(gse)
    updateSelectizeInput(
      session, "gseaTerm",
      choices  = df$ID,
      selected = head(df$ID, 1),
      server   = TRUE
    )
  })
  
  output$gseaDotplot <- renderPlot({
    gseres <- gseaResults()
    validate(
      need(!is.null(gseres), "Run GSEA first."),
      need(inherits(gseres, c("gseaResult","enrichResult","compareClusterResult")),
           paste("dotplot() needs a gseaResult/enrichResult. Got:", paste(class(gseres), collapse=", ")))
    )
    enrichplot::dotplot(gseres, showCategory = input$numcategories)  # <- enrichplot
  })
  
  
  
  output$gseaEnrichPlot <- renderPlot({
    gseres <- gseaResults()
    validate(
      need(!is.null(gseres), "No GSEA results to display."),
      need(input$gseaTerm, "Select a term to display.")
    )
    df_check <- as.data.frame(gseres)
    validate(need(input$gseaTerm %in% df_check$ID, "Selected term not found in results."))
    enrichplot::gseaplot2(gseres, geneSetID = input$gseaTerm, title = input$gseaTerm)
  })
  
  output$gseaTable <- DT::renderDataTable({
    gseres <- gseaResults()
    validate(need(!is.null(gseres), "No GSEA results to display."))
    DT::datatable(as.data.frame(gseres),
                  extensions = "Buttons",
                  options = list(dom = "Bfrtip", buttons = c("copy","csv")))
  })
  
  output$dl_gsea <- downloadHandler(
    filename = function() "gsea_results.csv",
    content  = function(f) {
      gse <- gseaResults()
      if (is.null(gse)) return()
      readr::write_csv(as.data.frame(gse), f)
    }
  )
  
  ###---GO---###
  
  goSpeciesCode <- reactive({
    org <- isolate(input$organism)
    if (org == "org.Hs.eg.db") return("Hs")
    if (org == "org.Mm.eg.db") return("Mm")
    if (org == "org.Rn.eg.db") return("Rn")
    if (org == "org.Dm.eg.db") return("Dm")
    if (org == "org.Dr.eg.db") return("Dr")
    "hsa"
  })
  
  goResults <- eventReactive(input$runGO, {
    withProgress(message = "Running GO analysis", value = 0, {
      appendLog("Running GO analysis...")
      showNotification("Running GO analysis...", type="message")
      incProgress(0.2)
      req(analysisResults())
      ba <- analysisResults()
      
      orgdb_name <- isolate(input$organism)
      orgdb <- isolate(get(orgdb_name))
      
      # First, ensure we have a clean data frame with no length mismatches
      res_clean <- ba$res_entrez[complete.cases(ba$res_entrez[c("padj", "log2FoldChange", "ENTREZID")]), ]
      
      # Now apply your filters on the cleaned data
      sig_genes <- res_clean[which(res_clean$padj < isolate(input$padjThreshold) &
                                     abs(res_clean$log2FoldChange) >= isolate(input$lfcThreshold)), ]
      
      # Extract vectors
      entrez_id_vector <- sig_genes$ENTREZID
      universe_entrez_ids <- res_clean$ENTREZID
      
      # Fix gene lengths
      gene_lengths <- as.numeric(res_clean$gene_lengths)
      names(gene_lengths) <- res_clean$ENTREZID
      
      # If lengths are unavailable (plain count table), set covariate to NULL
      covariate_arg <- if (all(is.na(gene_lengths))) NULL else gene_lengths
      
      if (is.null(covariate_arg)) {
        appendLog("Note: Gene lengths not available — GO running without length-bias correction.")
        showNotification("GO running without length-bias correction (no gene lengths available).",
                         type = "warning", duration = 6)
      } else {
        appendLog("Gene lengths available — applying length-bias correction in GO.")
      }
      
      # Read species code 
      sp <- isolate(goSpeciesCode())
      
      # Run GO enrichment analysis with goana
      go_results <- limma::goana(de = entrez_id_vector, species = sp, 
                                 universe = universe_entrez_ids, 
                                 covariate = covariate_arg)
      
      incProgress(0.8)
      appendLog("GO analysis complete.")
      showNotification("GO analysis complete.", type="message")
      
      list(go_results = go_results)
    })
  })
  
  output$GOBarplot <- renderPlot({
    req(goResults(), input$GOontology)
    
    go_results <- goResults()$go_results
    topgo <- limma::topGO(
      go_results,
      ontology = input$GOontology,
      number   = input$GOnumber
    )
    
    ggplot(data = topgo, aes(x = reorder(Term, -log10(P.DE)), y = -log10(P.DE))) +
      geom_bar(stat = "identity") +
      coord_flip() +
      theme_minimal(base_size = 13) +
      labs(x = "GO Terms", y = "-log10(p-value)", title = paste("Top GO Terms", 
                                                                isolate(input$GOontology), sep = ": ")
      )
  })
  
  output$goTable <- DT::renderDataTable({
    req(goResults(), input$GOontology)
    go_results <- goResults()$go_results
    topgo <- limma::topGO(go_results, ontology = input$GOontology, number = Inf)
    topgo$GOID <- rownames(topgo)
    topgo <- topgo[, c("GOID", setdiff(names(topgo), "GOID"))]
    DT::datatable(topgo, options = list(pageLength = 10, scrollX = TRUE),
                  rownames = FALSE)
  })
  
  output$dl_go <- downloadHandler(
    filename = function() paste0("GO_", input$GOontology, "_results.csv"),
    content  = function(f) {
      go_results <- goResults()$go_results
      topgo <- limma::topGO(go_results, ontology = input$GOontology, number = input$GOnumber)
      topgo$GOID <- rownames(topgo)
      readr::write_csv(topgo, f)
    }
  )
  
  # ============================================================
  # EXPLORE GENE SETS – server logic
  # ============================================================
  
  exploreGenes <- reactiveVal(NULL)  # data.frame with at least GeneSymbol column
  exploreTermLabel <- reactiveVal("")
  
  # --- Dynamic term selector ---
  output$exploreTermUI <- renderUI({
    src <- input$exploreSource
    
    if (src == "go") {
      req(goResults())
      go_df <- goResults()$go_results
      # Show top terms across all ontologies, ordered by P.DE
      go_df <- go_df[order(go_df$P.DE), ]
      term_choices <- setNames(rownames(go_df),
                               paste0(rownames(go_df), ": ", go_df$Term, " (p=", signif(go_df$P.DE, 3), ")"))
      #term_choices <- head(term_choices, 200) # commented out so that all terms are available (slower)
      selectizeInput("exploreTerm", "Select GO term:", choices = term_choices, selected = term_choices[1])
      
    } else if (src == "gsea") {
      req(gseaResults())
      gsea_df <- as.data.frame(gseaResults())
      term_choices <- setNames(gsea_df$ID,
                               paste0(gsea_df$ID, ": ", gsea_df$Description, " (p.adj=", signif(gsea_df$p.adjust, 3), ")"))
      selectizeInput("exploreTerm", "Select GSEA term:", choices = term_choices, selected = term_choices[1])
      
    } else if (src == "kegg") {
      req(gageRes())
      gr <- gageRes()
      up_ids   <- if (NROW(gr$up))   gr$up$PathwayID   else character(0)
      down_ids <- if (NROW(gr$down)) gr$down$PathwayID else character(0)
      up_desc   <- if (NROW(gr$up))   gr$up$Description   else character(0)
      down_desc <- if (NROW(gr$down)) gr$down$Description else character(0)
      
      all_ids   <- c(up_ids, down_ids)
      all_descs <- c(up_desc, down_desc)
      if (length(all_ids) == 0) return(helpText("No significant KEGG pathways found. Run GAGE first."))
      
      term_choices <- setNames(all_ids, paste0(all_ids, ": ", all_descs))
      selectizeInput("exploreTerm", "Select KEGG pathway:", choices = term_choices, selected = term_choices[1])
    }
  })
  
  # --- Load genes for the selected term ---
  observeEvent(input$loadTermGenesBtn, {
    req(input$exploreTerm, analysisResults())
    
    src <- input$exploreSource
    term_id <- input$exploreTerm
    orgdb_name <- isolate(input$organism)
    orgdb <- get(orgdb_name)
    ba <- analysisResults()
    res_df <- as.data.frame(ba$res)
    res_df$GeneSymbol <- rownames(res_df)
    res_entrez <- ba$res_entrez
    
    gene_symbols <- character(0)
    
    withProgress(message = "Looking up genes in term...", value = 0.3, {
      
      if (src == "go") {
        # Get ENTREZID for this GO term from the OrgDb
        tryCatch({
          entrez_in_term <- AnnotationDbi::select(orgdb, keys = term_id,
                                                  columns = "ENTREZID",
                                                  keytype = "GOALL")$ENTREZID
          entrez_in_term <- unique(entrez_in_term[!is.na(entrez_in_term)])
          
          # Map to gene symbols via res_entrez
          matched <- res_entrez[res_entrez$ENTREZID %in% entrez_in_term, ]
          gene_symbols <- unique(matched$GeneIDs)
          
          # Store the term label
          go_df <- goResults()$go_results
          if (term_id %in% rownames(go_df)) {
            exploreTermLabel(paste0(term_id, ": ", go_df[term_id, "Term"]))
          } else {
            exploreTermLabel(term_id)
          }
        }, error = function(e) {
          showNotification(paste("GO lookup failed:", conditionMessage(e)), type = "error")
        })
        
      } else if (src == "gsea") {
        gsea_df <- as.data.frame(gseaResults())
        row <- gsea_df[gsea_df$ID == term_id, ]
        if (nrow(row) == 0) {
          showNotification("Term not found in GSEA results.", type = "error")
          return()
        }
        # core_enrichment contains slash-separated ENTREZID
        core_entrez <- unlist(strsplit(row$core_enrichment, "/"))
        matched <- res_entrez[res_entrez$ENTREZID %in% core_entrez, ]
        gene_symbols <- unique(matched$GeneIDs)
        exploreTermLabel(paste0(term_id, ": ", row$Description))
        
      } else if (src == "kegg") {
        gr <- gageRes()
        sp <- gr$species
        
        tryCatch({
          kegg_entrez <- character(0)
          
          # Use the same kegg.gsets that GAGE used — guaranteed Entrez IDs
          kg <- gage::kegg.gsets(species = sp, id.type = "entrez", check.new = FALSE)
          matching_idx <- grep(paste0("^", term_id), names(kg$kg.sets))
          
          if (length(matching_idx) > 0) {
            kegg_entrez <- unique(kg$kg.sets[[matching_idx[1]]])
          }
          
          kegg_entrez <- kegg_entrez[!is.na(kegg_entrez) & nzchar(kegg_entrez)]
          
          # Force both to character for safe matching
          kegg_entrez_chr <- as.character(kegg_entrez)
          entrez_col_chr  <- as.character(res_entrez$ENTREZID)
          
          matched <- res_entrez[entrez_col_chr %in% kegg_entrez_chr, ]
          gene_symbols <- unique(matched$GeneIDs)
          
          # Label
          all_paths <- rbind(
            if (NROW(gr$up)) gr$up else NULL,
            if (NROW(gr$down)) gr$down else NULL
          )
          desc_row <- all_paths[all_paths$PathwayID == term_id, ]
          if (NROW(desc_row) > 0) {
            exploreTermLabel(paste0(term_id, ": ", desc_row$Description[1]))
          } else {
            exploreTermLabel(term_id)
          }
        }, error = function(e) {
          showNotification(paste("KEGG lookup failed:", conditionMessage(e)), type = "error")
        })
      }
      
      incProgress(0.7)
    })
    
    # Filter to genes that exist in the DESeq2 results
    gene_symbols <- gene_symbols[gene_symbols %in% rownames(res_df)]
    
    if (length(gene_symbols) == 0) {
      showNotification("No genes from this term found in the DE results.", type = "warning")
      exploreGenes(NULL)
      return()
    }
    
    # Build output table
    out_df <- res_df[gene_symbols, c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"), drop = FALSE]
    out_df$GeneSymbol <- rownames(out_df)
    out_df <- out_df[, c("GeneSymbol", setdiff(names(out_df), "GeneSymbol"))]
    out_df <- out_df[order(out_df$padj), ]
    rownames(out_df) <- NULL
    
    exploreGenes(out_df)
    showNotification(paste(nrow(out_df), "genes loaded from", exploreTermLabel()), type = "message")
  })
  
  output$exploreGeneCountMsg <- renderUI({
    eg <- exploreGenes()
    if (is.null(eg)) return(NULL)
    tags$p(style = "font-weight:bold;",
           paste0(nrow(eg), " genes in: ", exploreTermLabel()))
  })
  
  # --- Gene table ---
  output$exploreGeneTable <- DT::renderDataTable({
    req(exploreGenes())
    DT::datatable(exploreGenes(), options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE)
  })
  
  output$dl_explore_genes <- downloadHandler(
    filename = function() paste0("geneset_", gsub("[^A-Za-z0-9]", "_", exploreTermLabel()), ".csv"),
    content  = function(f) readr::write_csv(req(exploreGenes()), f)
  )
  
  # --- Heatmap ---
  output$exploreHeatmap <- renderPlot({
    req(exploreGenes(), analysisResults())
    
    genes <- exploreGenes()$GeneSymbol
    vsd   <- analysisResults()$vsd
    mat   <- assay(vsd)
    
    genes_in_mat <- genes[genes %in% rownames(mat)]
    validate(need(length(genes_in_mat) >= 2, "Need at least 2 genes to draw a heatmap."))
    
    hm_mat <- mat[genes_in_mat, , drop = FALSE]
    
    # Annotation
    sample_annotation <- data.frame(
      Group = colData(vsd)$Group,
      row.names = colnames(hm_mat)
    )
    n_groups <- length(unique(sample_annotation$Group))
    if (n_groups <= 8) {
      group_colors <- RColorBrewer::brewer.pal(max(n_groups, 3), "Set2")[seq_len(n_groups)]
    } else {
      group_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_groups)
    }
    names(group_colors) <- unique(sample_annotation$Group)
    
    # Colors
    if (input$exploreHeatColor == "viridis") {
      colors <- viridisLite::viridis(100)
    } else {
      colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, input$exploreHeatColor)))(100)
    }
    
    pheatmap::pheatmap(
      hm_mat,
      scale = input$exploreHeatScale,
      clustering_distance_rows = "correlation",
      clustering_distance_cols = "correlation",
      show_rownames = (length(genes_in_mat) <= 80),
      annotation_col = sample_annotation,
      annotation_colors = list(Group = group_colors),
      color = colors,
      main = exploreTermLabel(),
      fontsize = 10,
      fontsize_row = if (length(genes_in_mat) > 50) 6 else 8
    )
  }, res = 96)
  
  # --- Boxplots ---
  output$exploreBoxplotGeneUI <- renderUI({
    req(exploreGenes())
    genes <- exploreGenes()$GeneSymbol
    selectizeInput("exploreBoxGene", "Select gene(s) to plot:",
                   choices = genes, selected = head(genes, min(6, length(genes))),
                   multiple = TRUE,
                   options = list(maxItems = 20))
  })
  
  output$exploreBoxplot <- renderPlot({
    req(input$exploreBoxGene, analysisResults())
    
    vsd <- analysisResults()$vsd
    mat <- assay(vsd)
    meta <- as.data.frame(colData(vsd))
    
    genes_to_plot <- input$exploreBoxGene
    genes_to_plot <- genes_to_plot[genes_to_plot %in% rownames(mat)]
    validate(need(length(genes_to_plot) >= 1, "Select at least one gene."))
    
    # Build long-form data frame
    plot_list <- lapply(genes_to_plot, function(g) {
      data.frame(
        Gene = g,
        Expression = mat[g, ],
        Sample = colnames(mat),
        Group = meta$Group,
        stringsAsFactors = FALSE
      )
    })
    plot_df <- do.call(rbind, plot_list)
    
    ggplot(plot_df, aes(x = Group, y = Expression)) +
      geom_boxplot(aes(fill = Group), outlier.shape = NA, alpha = 0.5) +
      geom_jitter(aes(colour = Group), width = 0.2, size = 2) +
      facet_wrap(~ Gene, scales = "free_y") +
      theme_minimal(base_size = 13) +
      labs(title = exploreTermLabel(),
           y = "VST-normalised expression", x = "Group") +
      theme(legend.position = "bottom")
  })
  
  # --- Volcano with highlighted genes ---
  output$exploreVolcano <- renderPlot({
    req(exploreGenes(), analysisResults())
    
    res <- analysisResults()$res
    res_df <- as.data.frame(res)
    res_df$GeneSymbol <- rownames(res_df)
    
    term_genes <- exploreGenes()$GeneSymbol
    res_df$InTerm <- res_df$GeneSymbol %in% term_genes
    
    # Significance for colouring the term genes
    lfc_cutoff  <- input$lfcThreshold
    padj_cutoff <- input$padjThreshold
    
    res_df$TermStatus <- "Background"
    res_df$TermStatus[res_df$InTerm & !is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange >= lfc_cutoff]  <- "In term (Up)"
    res_df$TermStatus[res_df$InTerm & !is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange <= -lfc_cutoff] <- "In term (Down)"
    res_df$TermStatus[res_df$InTerm & res_df$TermStatus == "Background"] <- "In term (NS)"
    
    res_df$TermStatus <- factor(res_df$TermStatus,
                                levels = c("Background", "In term (NS)", "In term (Up)", "In term (Down)"))
    
    # Labels for term genes (top by padj)
    term_df <- res_df[res_df$InTerm, ]
    term_df <- term_df[order(term_df$padj), ]
    top_labels <- head(term_df$GeneSymbol, 30)
    res_df$label <- ifelse(res_df$GeneSymbol %in% top_labels, res_df$GeneSymbol, "")
    
    # Alpha: background dim, term genes full
    res_df$pt_alpha <- ifelse(res_df$InTerm, 1, input$exploreVolcAlpha)
    
    # Size: term genes slightly larger
    res_df$pt_size <- ifelse(res_df$InTerm, input$exploreVolcPtSize * 1.5, input$exploreVolcPtSize)
    
    term_colors <- c(
      "Background"     = "grey80",
      "In term (NS)"   = "#f4a736",
      "In term (Up)"   = "#d62728",
      "In term (Down)"  = "#1f77b4"
    )
    
    # Plot background first, then term genes on top
    ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(data = res_df[!res_df$InTerm, ],
                 colour = "grey80", size = input$exploreVolcPtSize,
                 alpha = input$exploreVolcAlpha) +
      geom_point(data = res_df[res_df$InTerm, ],
                 aes(colour = TermStatus),
                 size = input$exploreVolcPtSize * 1.5, alpha = 1) +
      ggrepel::geom_text_repel(
        data = res_df[res_df$label != "", ],
        aes(label = label),
        size = input$exploreVolcLabelSize,
        max.overlaps = 25, show.legend = FALSE,
        segment.size = 0.3, segment.color = "grey40"
      ) +
      scale_colour_manual(values = term_colors, name = "Gene status") +
      geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", colour = "grey40") +
      geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", colour = "grey40") +
      theme_minimal(base_size = 13) +
      labs(title = exploreTermLabel(),
           x = "log2 Fold Change", y = "-log10(adjusted p-value)")
  }, res = 96)
  
  ###---Interactive volcano plot output---###
  
  output$volcanoPlot <- renderPlotly({
    req(analysisResults())
    
    res <- analysisResults()$res
    res_df <- as.data.frame(res)
    
    res_df$log10padj <- -log10(res_df$padj)
    res_df$log10padj[is.infinite(res_df$log10padj)] <- NA
    
    # Thresholding
    lfc_cutoff <- input$lfcThreshold
    padj_cutoff <- input$padjThreshold
    
    res_df$significance <- "Not significant"
    res_df$significance[res_df$padj < padj_cutoff & res_df$log2FoldChange >= lfc_cutoff] <- "Upregulated"
    res_df$significance[res_df$padj < padj_cutoff & res_df$log2FoldChange <= -lfc_cutoff] <- "Downregulated"
    
    # Set levels so even missing categories are recognized
    res_df$significance <- factor(res_df$significance,
                                  levels = c("Downregulated", "Not significant", "Upregulated"))
    
    volcano_colors <- c(
      "Downregulated" = "blue",
      "Not significant" = "grey",
      "Upregulated" = "red"
    )
    
    volcano_data <- res_df
    volcano_data$Gene <- rownames(volcano_data)
    
    p <- plotly::plot_ly(
      data = volcano_data,
      x = ~log2FoldChange,
      y = ~log10padj,
      type = "scatter",
      mode = "markers",
      text = ~paste("Gene: ", Gene, "<br>Log2FC: ", signif(log2FoldChange, 3), "<br>padj: ", signif(padj, 4)),
      color = ~significance,
      colors = volcano_colors,
      key = ~Gene,
      source = "volcano",
      marker = list(size = 6)
    )
    
    p <- plotly::event_register(p, "plotly_click")
    
    plotly::layout(
      p,
      title = "Interactive Volcano Plot",
      xaxis = list(title = "log2 Fold Change"),
      yaxis = list(title = "-log10 Adjusted p-value"),
      shapes = list(
        list(type = "line", x0 = -lfc_cutoff, x1 = -lfc_cutoff, y0 = 0, y1 = max(volcano_data$log10padj, na.rm = TRUE), line = list(dash = "dash")),
        list(type = "line", x0 =  lfc_cutoff, x1 =  lfc_cutoff, y0 = 0, y1 = max(volcano_data$log10padj, na.rm = TRUE), line = list(dash = "dash")),
        list(type = "line", y0 = -log10(padj_cutoff), y1 = -log10(padj_cutoff),
             x0 = min(volcano_data$log2FoldChange), x1 = max(volcano_data$log2FoldChange), line = list(dash = "dash"))
      )
    )
    
  })
  
  
  output$geneBoxplot <- renderPlot({
    click <- plotly::event_data("plotly_click", source = "volcano")
    
    req(click$key)
    clicked_gene <- click$key
    
    res <- analysisResults()$res
    vst_mat <- assay(vst(analysisResults()$dds), blind = TRUE)
    sample_meta <- as.data.frame(colData(analysisResults()$dds))
    
    
    df <- data.frame(
      Expression = vst_mat[clicked_gene, ],
      Sample = colnames(vst_mat),
      Group = sample_meta$Group
    )
    
    ggplot(df, aes(x = Group, y = Expression)) +
      geom_boxplot(aes(fill = Group), outlier.shape = NA, alpha = 0.5) +
      geom_jitter(aes(color = Group), width = 0.2, size = 2) +
      theme_minimal() +
      labs(
        title = paste("Expression of", clicked_gene),
        y = "vst-normalized expression",
        x = "Group"
      )
  })
  
  # Normal volcano plot
  output$standardvolcanoplot <- renderPlot({
    req(analysisResults())
    
    res <- analysisResults()$res
    res_df <- as.data.frame(res)
    
    # Basic sanity checks (DESeq2 results should have these)
    validate(
      need("log2FoldChange" %in% names(res_df), "Missing log2FoldChange in results."),
      need("padj" %in% names(res_df), "Missing padj in results."),
      need("pvalue" %in% names(res_df), "Missing pvalue in results.")
    )
    
    # Compute -log10 metrics, handling padj/pvalue == 0
    min_nonzero_padj <- min(res_df$padj[res_df$padj > 0], na.rm = TRUE)
    min_nonzero_p    <- min(res_df$pvalue[res_df$pvalue > 0], na.rm = TRUE)
    
    res_df$log10padj <- -log10(ifelse(res_df$padj == 0, min_nonzero_padj, res_df$padj))
    res_df$log10p    <- -log10(ifelse(res_df$pvalue == 0, min_nonzero_p, res_df$pvalue))
    
    # Cap at 99.9th percentile to prevent axis stretching
    y_col_name <- if (isTRUE(input$pvalueselect)) "log10padj" else "log10p"
    y_cap <- quantile(res_df[[y_col_name]], 0.999, na.rm = TRUE)
    res_df$capped <- !is.na(res_df[[y_col_name]]) & res_df[[y_col_name]] > y_cap
    res_df$log10padj <- pmin(res_df$log10padj, quantile(res_df$log10padj, 0.999, na.rm = TRUE))
    res_df$log10p    <- pmin(res_df$log10p, quantile(res_df$log10p, 0.999, na.rm = TRUE))
    
    # Decide which p column is active
    use_padj <- isTRUE(input$pvalueselect)  # TRUE -> use padj; FALSE -> use pvalue
    p_col    <- if (use_padj) "padj" else "pvalue"
    y_col    <- if (use_padj) "log10padj" else "log10p"
    y_lab    <- if (use_padj) "-Log10 Adjusted P-value" else "-Log10 P-value"
    
    # Thresholds
    lfc_cutoff <- input$lfcThreshold
    p_cutoff   <- input$padjThreshold  # keep your existing input name as the cutoff value
    
    # Significance (based on chosen p metric)
    res_df$significance <- "Not significant"
    sig_ok <- !is.na(res_df[[p_col]]) & (res_df[[p_col]] < p_cutoff)
    
    res_df$significance[sig_ok & res_df$log2FoldChange >=  lfc_cutoff] <- "Upregulated"
    res_df$significance[sig_ok & res_df$log2FoldChange <= -lfc_cutoff] <- "Downregulated"
    
    res_df$significance <- factor(
      res_df$significance,
      levels = c("Downregulated", "Not significant", "Upregulated")
    )
    
    # Label top genes by chosen p metric (and avoid NA ordering issues)
    n_top <- input$topgenes
    
    res_filt <- res_df %>%
      filter(significance %in% c("Upregulated", "Downregulated"))
    
    ord <- order(res_filt[[p_col]], na.last = NA) # removes NA p-values from ranking
    
    top_genes <- rownames(res_filt)[head(ord, n_top)]
    res_df$label <- ifelse(rownames(res_df) %in% top_genes, rownames(res_df), NA)
    
    volcano_colors <- c(
      "Downregulated"   = "blue",
      "Not significant" = "grey",
      "Upregulated"     = "red"
    )
    
    ggplot(res_df, aes(x = log2FoldChange, y = .data[[y_col]])) +
      geom_point(aes(color = significance, shape = capped), size = input$volcpointsize) +
      scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17),
                         labels = c("FALSE" = "Within range", "TRUE" = "Capped"),
                         guide = guide_legend(title = NULL)) +
      scale_color_manual(values = volcano_colors) +
      theme_minimal() +
      theme(
        legend.title = element_text(face = "bold", size = input$volclegendtitlesize),
        legend.text  = element_text(size = input$volclegendtextsize),
        axis.text    = element_text(face = "bold", size = input$volcaxistextsize),
        axis.title   = element_text(face = "bold", size = input$volcaxistitlesize)
      ) +
      labs(
        title = "Volcano Plot",
        x = "Log2 Fold Change",
        y = y_lab
      ) +
      geom_hline(
        yintercept = -log10(p_cutoff),
        linetype = "dashed",
        color = "black"
      ) +
      geom_vline(
        xintercept = c(-lfc_cutoff, lfc_cutoff),
        linetype = "dashed",
        color = "black"
      ) +
      ggrepel::geom_text_repel(
        aes(label = label),
        size = input$volclabelsize,
        max.overlaps = 10
      )
  }, res = 96)
  
  output$dl_res <- downloadHandler(
    filename = function() "deseq2_results.csv",
    content  = function(f) utils::write.csv(as.data.frame(analysisResults()$res),
                                            f, row.names = TRUE)
  )
  
  #---GAGE---#
  
  gageSpeciesCode <- reactive({
    org <- isolate(input$organism)
    if (org == "org.Hs.eg.db") return("hsa")
    if (org == "org.Mm.eg.db") return("mmu")
    if (org == "org.Rn.eg.db") return("rno")
    if (org == "org.Dm.eg.db") return("dme")
    if (org == "org.Dr.eg.db") return("dre")
    "hsa"
  })
  
  gageRes     <- reactiveVal(NULL)
  gageRunning <- reactiveVal(FALSE)
  
  # SINGLE handler for the button
  observeEvent(input$runGAGE, ignoreInit = TRUE, {
    if (isTRUE(gageRunning())) return(invisible(NULL))
    gageRunning(TRUE); on.exit(gageRunning(FALSE), add = TRUE)
    
    withProgress(message = "Running GAGE (KEGG)", value = 0, {
      incProgress(0.05)
      appendLog(sprintf("runGAGE value: %s", isolate(input$runGAGE)))
      appendLog("Running GAGE (KEGG)...")
      
      # Read all inputs once
      sp      <- isolate(gageSpeciesCode())
      min_sz  <- isolate(input$gageMin)
      max_sz  <- isolate(input$gageMax)
      sameDir <- isTRUE(isolate(input$gageSameDir))
      qcut    <- isolate(input$gageP)
      metric  <- isolate(input$gageMetric)
      
      # Build ENTREZ-named vector from DE results
      ba <- isolate(analysisResults())
      req(!is.null(ba), "Run the DE analysis first.")
      req("ENTREZID" %in% names(ba$res_entrez), "ENTREZID not found in DE results.")
      
      metric_col <- if (metric %in% names(ba$res_entrez)) metric else "log2FoldChange"
      v   <- ba$res_entrez[[metric_col]]
      ids <- as.character(ba$res_entrez[["ENTREZID"]])
      keep <- is.finite(v) & !is.na(ids) & nzchar(ids)
      v <- v[keep]; ids <- ids[keep]
      
      agg <- tapply(v, ids, function(z) z[which.max(abs(z))])
      fc  <- as.numeric(agg); names(fc) <- as.character(names(agg))
      
      # Drop NA/empty names and enforce uniqueness
      ok_names <- !is.na(names(fc)) & nzchar(names(fc))
      fc <- fc[ok_names]
      fc <- fc[!duplicated(names(fc))]
      
      appendLog(sprintf("fc length (unique ENTREZ, cleaned): %d", length(fc)))
      
      # KEGG sets (all pathways), size filter
      gs <- gage::kegg.gsets(species = sp, id.type = "entrez", check.new = FALSE)
      kegg.gs <- lapply(gs$kg.sets, unique)
      kegg.gs <- kegg.gs[vapply(kegg.gs, function(g) length(g) >= min_sz && length(g) <= max_sz, logical(1))]
      
      universe <- unique(unlist(kegg.gs, use.names = FALSE))
      ov <- length(intersect(names(fc), universe))
      appendLog(sprintf("KEGG universe: %d genes; overlap with fc: %d", length(universe), ov))
      
      incProgress(0.5)
      
      # Run GAGE
      gr <- gage::gage(fc, gsets = kegg.gs, same.dir = sameDir, set.size = c(min_sz, max_sz))
      
      # Post-process like your script
      sig_up   <- as.data.frame(gr$greater, stringsAsFactors = FALSE)
      sig_down <- as.data.frame(gr$less,    stringsAsFactors = FALSE)
      
      # filter (prefer q.val; fallback p.val); don't drop rows for other NA columns
      if ("q.val" %in% names(sig_up))   sig_up   <- sig_up[!is.na(sig_up$q.val)   & sig_up$q.val   <= qcut, , drop = FALSE]
      else if ("p.val" %in% names(sig_up))   sig_up   <- sig_up[!is.na(sig_up$p.val)   & sig_up$p.val   <= qcut, , drop = FALSE]
      
      if ("q.val" %in% names(sig_down)) sig_down <- sig_down[!is.na(sig_down$q.val) & sig_down$q.val <= qcut, , drop = FALSE]
      else if ("p.val" %in% names(sig_down)) sig_down <- sig_down[!is.na(sig_down$p.val) & sig_down$p.val <= qcut, , drop = FALSE]
      
      # Add KEGG ID column (first 8 chars of rownames)
      rn_up   <- rownames(sig_up);   sig_up$PathwayID   <- if (!is.null(rn_up))   substr(rn_up,   1, 8) else character(0)
      rn_down <- rownames(sig_down); sig_down$PathwayID <- if (!is.null(rn_down)) substr(rn_down, 1, 8) else character(0)
      sig_up$Direction   <- if (nrow(sig_up))   rep("Up",   nrow(sig_up))   else character(0)
      sig_down$Direction <- if (nrow(sig_down)) rep("Down", nrow(sig_down)) else character(0)
      
      # human-readable pathway name at the end of the row
      sig_up$Description   <- if (!is.null(rn_up))   trimws(sub("^\\S+\\s*", "", rn_up))   else character(0)
      sig_down$Description <- if (!is.null(rn_down)) trimws(sub("^\\S+\\s*", "", rn_down)) else character(0)
      
      # Order for display
      if ("q.val" %in% names(sig_up))   sig_up   <- sig_up[order(sig_up$q.val),   , drop = FALSE]
      if ("q.val" %in% names(sig_down)) sig_down <- sig_down[order(sig_down$q.val), , drop = FALSE]
      if (!"q.val" %in% names(sig_up)   && "p.val" %in% names(sig_up))   sig_up   <- sig_up[order(sig_up$p.val),   , drop = FALSE]
      if (!"q.val" %in% names(sig_down) && "p.val" %in% names(sig_down)) sig_down <- sig_down[order(sig_down$p.val), , drop = FALSE]
      
      appendLog(sprintf("GAGE significant: up=%d, down=%d (cutoff=%.3f)", nrow(sig_up), nrow(sig_down), qcut))
      
      # Publish results ONCE for the rest of the app
      gageRes(list(up = sig_up, down = sig_down, fc = fc, species = sp))
      
      incProgress(0.98)
      showNotification("GAGE analysis complete.", type = "message")
      appendLog("GAGE analysis complete.")
    })
  })
  
  # ---- Tables ----
  output$gageUpTable <- DT::renderDataTable({
    gr <- req(gageRes()); df <- gr$up
    if (!NROW(df)) return(DT::datatable(data.frame(Message = "No significant upregulated pathways."),
                                        options = list(dom = 't'), rownames = FALSE))
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    rownames(df) <- NULL
    # ensure PathwayID present (you already add it upstream)
    if (!"PathwayID" %in% names(df)) {
      rn <- gr$up |> rownames()
      df$PathwayID <- if (!is.null(rn)) substr(rn, 1, 8) else ""
    }
    # ensure Description present (you added upstream)
    if (!"Description" %in% names(df)) {
      rn <- gr$up |> rownames()
      df$Description <- if (!is.null(rn)) trimws(sub("^\\S+\\s*", "", rn)) else ""
    }
    # (Optional) keep Description as the last column
    last_cols <- c(setdiff(names(df), "Description"), "Description")
    df <- df[, last_cols, drop = FALSE]
    
    DT::datatable(df, options = list(pageLength = 10), selection = "single", rownames = FALSE)
  })
  
  output$gageDownTable <- DT::renderDataTable({
    gr <- req(gageRes()); df <- gr$down
    if (!NROW(df)) return(DT::datatable(data.frame(Message = "No significant downregulated pathways."),
                                        options = list(dom = 't'), rownames = FALSE))
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    rownames(df) <- NULL
    if (!"PathwayID" %in% names(df)) {
      rn <- gr$down |> rownames()
      df$PathwayID <- if (!is.null(rn)) substr(rn, 1, 8) else ""
    }
    if (!"Description" %in% names(df)) {
      rn <- gr$down |> rownames()
      df$Description <- if (!is.null(rn)) trimws(sub("^\\S+\\s*", "", rn)) else ""
    }
    last_cols <- c(setdiff(names(df), "Description"), "Description")
    df <- df[, last_cols, drop = FALSE]
    
    DT::datatable(df, options = list(pageLength = 10), selection = "single", rownames = FALSE)
  })
  
  
  output$dl_gage_up <- downloadHandler(
    filename = function() "GAGE_results_UP.csv",
    content  = function(f) readr::write_csv(as.data.frame(gageRes()$up), f)
  )
  
  output$dl_gage_down <- downloadHandler(
    filename = function() "GAGE_results_DOWN.csv",
    content  = function(f) readr::write_csv(as.data.frame(gageRes()$down), f)
  )
  
  
  # ---- Update KEGG ID dropdown when results arrive ----
  observeEvent(gageRes(), ignoreInit = TRUE, {
    gr <- req(gageRes())
    get_ids <- function(df) {
      if (is.null(df) || !NROW(df)) return(character(0))
      if ("PathwayID" %in% names(df)) return(df$PathwayID)
      rn <- rownames(df); raw <- if (!is.null(rn)) rn else names(df)
      if (is.null(raw)) return(character(0))
      substr(raw, 1, 8)
    }
    ids <- unique(stats::na.omit(c(get_ids(gr$up), get_ids(gr$down))))
    updateSelectInput(session, "keggPathway",
                      choices  = ids,
                      selected = if (length(ids)) ids[[1]] else NULL)
  })
  
  # ---- Pathview image ----
  output$keggPathview <- renderImage({
    gr <- req(gageRes())
    fc <- gr$fc
    sp <- gr$species
    
    # Require a valid selection
    pid_raw <- input$keggPathway
    validate(need(!is.null(pid_raw) && length(pid_raw) == 1 && nzchar(pid_raw),
                  "Select a pathway to render."))
    pid <- substr(pid_raw, 1, 8)  # e.g. "dme01200"
    
    # Nothing to plot if fc is empty
    validate(need(length(fc) > 0, "No gene data to render."))
    
    tmpdir <- tempfile("pathview_"); dir.create(tmpdir)
    oldwd <- getwd(); setwd(tmpdir); on.exit(setwd(oldwd), add = TRUE)
    
    # Try native KEGG renderer first; if it errors, retry with graphviz
    pv_ok <- TRUE
    try({
      pathview::pathview(
        gene.data   = fc,
        pathway.id  = pid,
        species     = sp,
        gene.idtype = "entrez",
        kegg.native = TRUE,
        na.col      = "transparent",
        out.suffix  = "shiny"
      )
    }, silent = TRUE) -> try_native
    
    if (inherits(try_native, "try-error")) {
      pv_ok <- FALSE
    }
    
    if (!pv_ok) {
      try({
        pathview::pathview(
          gene.data   = fc,
          pathway.id  = pid,
          species     = sp,
          gene.idtype = "entrez",
          kegg.native = FALSE,      # fallback renderer
          na.col      = "transparent",
          out.suffix  = "shiny-fallback"
        )
      }, silent = TRUE) -> try_fallback
      
      validate(need(!inherits(try_fallback, "try-error"),
                    "Pathview failed to render this pathway. Try another pathway."))
    }
    
    # Prefer the data-coloured overlay(s) that include 'pathview' in the filename
    pngs <- list.files(tmpdir, pattern = "\\.png$", full.names = TRUE)
    overlay <- pngs[grepl(paste0("^", pid, ".*pathview.*\\.png$"), basename(pngs))]
    png_file <- if (length(overlay)) overlay[which.max(file.mtime(overlay))] else pngs[which.max(file.mtime(pngs))]
    validate(need(!is.na(png_file) && file.exists(png_file), "No image produced by pathview."))
    
    list(src = png_file, contentType = "image/png", alt = pid)
  }, deleteFile = FALSE)
  
  # ---- Download currently selected pathway image ----
  output$dl_pathview <- downloadHandler(
    filename = function() paste0("pathview_", substr(req(input$keggPathway), 1, 8), ".png"),
    content  = function(file) {
      gr  <- req(gageRes()); fc <- gr$fc; sp <- gr$species
      pid <- substr(req(input$keggPathway), 1, 8)
      tmpdir <- tempfile("pathview_"); dir.create(tmpdir)
      oldwd <- getwd(); setwd(tmpdir); on.exit(setwd(oldwd), add = TRUE)
      pathview::pathview(gene.data = fc, pathway.id = pid, species = sp,
                         kegg.native = TRUE, out.suffix = "shiny")
      
      # Prefer the *coloured overlay* for download
      pngs <- list.files(tmpdir, pattern = "\\.png$", full.names = TRUE)
      overlay <- pngs[grepl(paste0("^", pid, ".*pathview.*\\.png$"), basename(pngs))]
      png_file <- if (length(overlay)) overlay[which.max(file.mtime(overlay))] else pngs[which.max(file.mtime(pngs))]
      file.copy(png_file, file, overwrite = TRUE)
    }
  )
  
  # ============================================================
  # COMPARE CONTRASTS – server logic
  # ============================================================
  
  # --- Reactive stores ---
  extraContrasts <- reactiveVal(list())   # list of contrast specs added by user
  vennGenesets   <- reactiveVal(NULL)     # named list of gene-symbol vectors per contrast
  vennSets       <- reactiveVal(NULL)     # the Venn set membership breakdown
  cmpGageRes     <- reactiveVal(NULL)     # GAGE results for selected gene set
  cmpGseaRes     <- reactiveVal(NULL)     # GSEA results for selected gene set
  cmpSelectedGenes <- reactiveVal(NULL)   # DEG data.frame for selected Venn region
  
  # --- Dynamic UI: numerator / denominator selects for adding contrasts ---
  output$cmpNumeratorUI <- renderUI({
    req(input$groupOrder)
    selectInput("cmpNumerator", "Numerator (treatment):", choices = input$groupOrder)
  })
  output$cmpDenominatorUI <- renderUI({
    req(input$groupOrder)
    selectInput("cmpDenominator", "Denominator (reference):", choices = input$groupOrder)
  })
  
  # --- Add a contrast ---
  observeEvent(input$addContrastBtn, {
    req(input$cmpNumerator, input$cmpDenominator)
    if (input$cmpNumerator == input$cmpDenominator) {
      showNotification("Numerator and denominator must be different groups.", type = "error")
      return()
    }
    label <- paste0(input$cmpNumerator, "_vs_", input$cmpDenominator)
    new_entry <- list(
      label     = label,
      numerator = input$cmpNumerator,
      denominator = input$cmpDenominator,
      lfc       = input$cmpLFC,
      padj      = input$cmpPadj
    )
    current <- extraContrasts()
    # Avoid duplicate labels
    if (any(sapply(current, `[[`, "label") == label)) {
      showNotification(paste("Contrast", label, "already added."), type = "warning")
      return()
    }
    extraContrasts(c(current, list(new_entry)))
    showNotification(paste("Added contrast:", label), type = "message")
  })
  
  observeEvent(input$clearContrastsBtn, {
    extraContrasts(list())
    vennGenesets(NULL)
    vennSets(NULL)
    showNotification("All extra contrasts cleared.", type = "message")
  })
  
  # --- Show table of contrasts (primary + extras) ---
  output$contrastListTable <- renderTable({
    req(analysisResults())
    primary <- data.frame(
      Label       = paste0(input$contrastNumerator, "_vs_", input$contrastDenominator, " [PRIMARY]"),
      Numerator   = input$contrastNumerator,
      Denominator = input$contrastDenominator,
      LFC         = input$lfcThreshold,
      Padj        = input$padjThreshold,
      stringsAsFactors = FALSE
    )
    extras <- extraContrasts()
    if (length(extras) == 0) return(primary)
    extra_df <- do.call(rbind, lapply(extras, function(x) {
      data.frame(Label = x$label, Numerator = x$numerator,
                 Denominator = x$denominator, LFC = x$lfc, Padj = x$padj,
                 stringsAsFactors = FALSE)
    }))
    rbind(primary, extra_df)
  })
  
  # ---- Helper: run a DESeq2 contrast and return sig gene symbols ----
  run_contrast_degs <- function(dds, numerator, denominator, lfc_thr, padj_thr, direction = "both") {
    res_c <- results(dds,
                     lfcThreshold = lfc_thr,
                     altHypothesis = "greaterAbs",
                     alpha = padj_thr,
                     contrast = c("Group", numerator, denominator))
    res_df <- as.data.frame(res_c)
    res_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]
    
    if (direction == "both") {
      sig <- res_df[res_df$padj < padj_thr & abs(res_df$log2FoldChange) >= lfc_thr, ]
    } else if (direction == "up") {
      sig <- res_df[res_df$padj < padj_thr & res_df$log2FoldChange >= lfc_thr, ]
    } else {
      sig <- res_df[res_df$padj < padj_thr & res_df$log2FoldChange <= -lfc_thr, ]
    }
    list(genes = rownames(sig), res_df = res_df, sig_df = sig)
  }
  
  # ---- Run Venn / Euler ----
  observeEvent(input$runVennBtn, {
    req(analysisResults())
    dds <- analysisResults()$dds
    dir <- input$vennDEGdirection
    
    withProgress(message = "Computing contrasts for Venn...", value = 0, {
      
      # Primary contrast
      primary_label <- paste0(input$contrastNumerator, "_vs_", input$contrastDenominator)
      primary_data  <- run_contrast_degs(dds, 
                                         input$contrastNumerator, input$contrastDenominator,
                                         input$cmpLFC, input$cmpPadj, dir)
      genesets  <- list()
      sig_dfs   <- list()   # store per-contrast significant results with correct stats
      genesets[[primary_label]] <- primary_data$genes
      sig_dfs[[primary_label]]  <- primary_data$sig_df
      
      extras <- extraContrasts()
      incProgress(0.3)
      
      for (i in seq_along(extras)) {
        x <- extras[[i]]
        dat <- run_contrast_degs(dds, x$numerator, x$denominator, x$lfc, x$padj, dir)
        genesets[[x$label]] <- dat$genes
        sig_dfs[[x$label]]  <- dat$sig_df
        incProgress(0.3 / max(length(extras), 1))
      }
      
      vennGenesets(genesets)
      
      # Compute set membership
      all_genes <- unique(unlist(genesets))
      membership_mat <- sapply(genesets, function(gs) all_genes %in% gs)
      rownames(membership_mat) <- all_genes
      vennSets(list(genesets = genesets, sig_dfs = sig_dfs, membership = membership_mat))
      
      incProgress(0.4)
    })
    
    showNotification("Venn diagram ready.", type = "message")
  })
  
  output$vennPlot <- renderPlot({
    req(vennSets())
    gs <- vennSets()$genesets
    validate(need(length(gs) >= 2, "Please add at least one extra contrast and click 'Draw Diagram'."))
    
    if (input$vennType == "venn") {
      ggVennDiagram::ggVennDiagram(gs, label_alpha = 0) +
        scale_fill_gradient(low = "#DDF093", high = "#638475") +
        coord_equal(clip = "off") +
        theme_void() +
        theme(
          legend.position = "none",
          plot.margin = margin(t = 30, r = 80, b = 30, l = 180, unit = "pt")  # <- more padding
        )
    } else {
      fit <- eulerr::euler(gs)
      plot(fit,
           fills = list(fill = RColorBrewer::brewer.pal(min(length(gs), 8), "Set2"), alpha = 0.5),
           edges = list(col = "black", lex = 2),
           quantities = TRUE,
           labels = list(font = 2))
    }
  })
  
  output$dl_venn <- downloadHandler(
    filename = function() paste0("venn_diagram_", Sys.Date(), ".png"),
    content = function(file) {
      req(vennSets())
      gs <- vennSets()$genesets
      png(file, width = 1800, height = 1400, res = 150)
      if (input$vennType == "euler") {
        fit <- eulerr::euler(gs)
        print(plot(fit,
                   fills = list(fill = RColorBrewer::brewer.pal(min(length(gs), 8), "Set2"), alpha = 0.5),
                   edges = list(col = "black", lex = 2),
                   quantities = TRUE,
                   labels = list(font = 2)))
      } else {
        p <- ggVennDiagram::ggVennDiagram(gs, label_alpha = 0) +
          scale_fill_gradient(low = "#DDF093", high = "#638475") +
          theme_void() + theme(legend.position = "none")
        print(p)
      }
      dev.off()
    }
  )
  
  # ---- Venn set selector ----
  # Build human-readable names for all Venn regions
  venn_region_choices <- reactive({
    req(vennSets())
    gs <- vennSets()$genesets
    nms <- names(gs)
    n <- length(nms)
    regions <- list()
    for (i in 1:n) {
      combos <- combn(n, i, simplify = FALSE)
      for (co in combos) {
        label <- paste(nms[co], collapse = " ∩ ")
        if (i < n) label <- paste0(label, " (unique to these)")
        regions[[label]] <- co
      }
    }
    # Add "all shared" shortcut at top if >2
    regions
  })
  
  output$vennSetSelectorUI <- renderUI({
    req(vennSets())
    gs <- vennSets()$genesets
    nms <- names(gs)
    
    tagList(
      radioButtons("vennSetMode", "Gene set type:",
                   choices = c("Shared across ALL contrasts" = "all_shared",
                               "Shared between selected contrasts" = "selected_shared",
                               "Unique to one contrast" = "unique"),
                   selected = "all_shared"),
      conditionalPanel(
        condition = "input.vennSetMode == 'selected_shared'",
        checkboxGroupInput("vennSharedContrasts", "Select contrasts (2 or more):",
                           choices = nms, selected = nms[1:min(2, length(nms))]),
        helpText("Genes must be significant in ALL ticked contrasts and in NONE of the unticked contrasts (exclusive overlap).")
      ),
      conditionalPanel(
        condition = "input.vennSetMode == 'unique'",
        selectInput("vennUniqueContrast", "Unique to:", choices = nms, selected = nms[1])
      ),
      conditionalPanel(
        condition = "input.vennSetMode != 'unique'",
        selectInput("vennConcordance", "Direction filter (shared sets):",
                    choices = c(
                      "All shared (average stats)"     = "all",
                      "Concordant UP in all contrasts" = "up",
                      "Concordant DOWN in all contrasts" = "down",
                      "Discordant only (mixed direction)" = "discordant"
                    ),
                    selected = "all"),
        helpText("Concordant: gene goes same direction in every relevant contrast.")
      )
    )
  })
  
  # Reactive to build the choice string from the new UI inputs
  vennSetChoice <- reactive({
    mode <- input$vennSetMode
    if (is.null(mode)) return(NULL)
    if (mode == "all_shared") return("ALL_SHARED")
    if (mode == "unique") {
      req(input$vennUniqueContrast)
      return(paste0("UNIQUE__", input$vennUniqueContrast))
    }
    if (mode == "selected_shared") {
      sel <- input$vennSharedContrasts
      validate(need(length(sel) >= 2, "Please select at least 2 contrasts."))
      return(paste0("PAIR__", paste(sel, collapse = "|||")))
    }
    NULL
  })
  
  # Helper to extract genes for a region choice
  extract_venn_genes <- function(gs, choice) {
    nms <- names(gs)
    if (choice == "ALL_SHARED") {
      return(Reduce(intersect, gs))
    }
    if (startsWith(choice, "UNIQUE__")) {
      focal <- sub("^UNIQUE__", "", choice)
      others <- setdiff(nms, focal)
      unique_genes <- gs[[focal]]
      for (o in others) unique_genes <- setdiff(unique_genes, gs[[o]])
      return(unique_genes)
    }
    if (startsWith(choice, "PAIR__")) {
      pair_nms <- strsplit(sub("^PAIR__", "", choice), "\\|\\|\\|")[[1]]
      shared <- Reduce(intersect, gs[pair_nms])
      # Exclude genes that also appear in other contrasts (exclusive overlap)
      others <- setdiff(nms, pair_nms)
      for (o in others) shared <- setdiff(shared, gs[[o]])
      return(shared)
    }
    character(0)
  }
  
  # ── Helper: classify concordance for shared gene sets ──────────────────────
  # Returns a data frame (one row per gene) with:
  #   GeneID, Direction_<contrast> columns (UP/DOWN), Concordance (concordant_up/down/discordant)
  #   plus mean log2FoldChange, mean stat, mean padj across contrasts
  build_shared_gene_table <- function(chosen_genes, src_names, sig_dfs) {
    parts <- lapply(src_names, function(nm) {
      df <- sig_dfs[[nm]]
      df <- df[rownames(df) %in% chosen_genes, , drop = FALSE]
      if (nrow(df) == 0) return(NULL)
      df$GeneID   <- rownames(df)
      df$Contrast <- nm
      df[, c("GeneID", "Contrast", "log2FoldChange", "stat", "padj",
             intersect(c("baseMean", "lfcSE", "pvalue"), names(df)))]
    })
    parts <- Filter(Negate(is.null), parts)
    if (length(parts) == 0) return(NULL)
    
    long_df <- do.call(rbind, parts)
    rownames(long_df) <- NULL
    
    # Wide summary: one row per gene
    lfc_wide <- reshape(long_df[, c("GeneID", "Contrast", "log2FoldChange")],
                        idvar = "GeneID", timevar = "Contrast", direction = "wide")
    names(lfc_wide) <- gsub("log2FoldChange.", "LFC_", names(lfc_wide), fixed = TRUE)
    
    padj_wide <- reshape(long_df[, c("GeneID", "Contrast", "padj")],
                         idvar = "GeneID", timevar = "Contrast", direction = "wide")
    names(padj_wide) <- gsub("padj.", "padj_", names(padj_wide), fixed = TRUE)
    
    stat_wide <- reshape(long_df[, c("GeneID", "Contrast", "stat")],
                         idvar = "GeneID", timevar = "Contrast", direction = "wide")
    names(stat_wide) <- gsub("stat.", "stat_", names(stat_wide), fixed = TRUE)
    
    # Direction per contrast
    lfc_cols <- grep("^LFC_", names(lfc_wide), value = TRUE)
    dir_wide  <- lfc_wide
    for (col in lfc_cols) {
      nm  <- sub("^LFC_", "", col)
      dir_wide[[paste0("Direction_", nm)]] <- ifelse(dir_wide[[col]] > 0, "UP ↑", "DOWN ↓")
    }
    
    # Concordance: all UP, all DOWN, or discordant
    dir_cols   <- grep("^Direction_", names(dir_wide), value = TRUE)
    directions <- dir_wide[, dir_cols, drop = FALSE]
    dir_wide$Concordance <- apply(directions, 1, function(row) {
      vals <- unique(row)
      if (length(vals) == 1 && vals == "UP ↑")   return("Concordant UP ↑")
      if (length(vals) == 1 && vals == "DOWN ↓") return("Concordant DOWN ↓")
      return("Discordant ↕")
    })
    
    # Mean stats
    lfc_mat  <- as.matrix(lfc_wide[, lfc_cols, drop = FALSE])
    stat_mat <- as.matrix(stat_wide[, grep("^stat_", names(stat_wide)), drop = FALSE])
    padj_mat <- as.matrix(padj_wide[, grep("^padj_", names(padj_wide)), drop = FALSE])
    dir_wide$mean_log2FoldChange <- rowMeans(lfc_mat,  na.rm = TRUE)
    dir_wide$mean_stat           <- rowMeans(stat_mat, na.rm = TRUE)
    dir_wide$mean_padj           <- rowMeans(padj_mat, na.rm = TRUE)
    
    # Merge padj columns in
    merged <- merge(dir_wide, padj_wide, by = "GeneID")
    merged <- merge(merged,   stat_wide, by = "GeneID")
    
    # Reorder: GeneID, Concordance, Direction_*, LFC_*, padj_*, stat_*, means
    first_cols  <- c("GeneID", "Concordance")
    dir_c  <- grep("^Direction_", names(merged), value = TRUE)
    lfc_c  <- grep("^LFC_",       names(merged), value = TRUE)
    padj_c <- grep("^padj_",      names(merged), value = TRUE)
    stat_c <- grep("^stat_",      names(merged), value = TRUE)
    mean_c <- c("mean_log2FoldChange", "mean_stat", "mean_padj")
    merged <- merged[, c(first_cols, dir_c, lfc_c, padj_c, stat_c, mean_c)]
    merged
  }
  
  # ── Helper: apply concordance filter to a wide gene table ───────────────────
  apply_concordance_filter <- function(gene_df, concordance_choice) {
    if (concordance_choice == "all")        return(gene_df)
    if (concordance_choice == "up")         return(gene_df[gene_df$Concordance == "Concordant UP ↑",   , drop = FALSE])
    if (concordance_choice == "down")       return(gene_df[gene_df$Concordance == "Concordant DOWN ↓", , drop = FALSE])
    if (concordance_choice == "discordant") return(gene_df[gene_df$Concordance == "Discordant ↕",      , drop = FALSE])
    gene_df
  }
  
  # Load genes for selected region
  observeEvent(input$loadVennGenesBtn, {
    req(vennSets(), analysisResults())
    vs       <- vennSets()
    gs       <- vs$genesets
    sig_dfs  <- vs$sig_dfs
    choice   <- vennSetChoice()
    req(choice)
    chosen_genes <- extract_venn_genes(gs, choice)
    
    if (length(chosen_genes) == 0) {
      showNotification("No genes in this set.", type = "warning")
      cmpSelectedGenes(NULL)
      return()
    }
    
    if (startsWith(choice, "UNIQUE__")) {
      # Unique set: one row per gene, stats from that contrast only
      focal   <- sub("^UNIQUE__", "", choice)
      src_df  <- sig_dfs[[focal]]
      gene_df <- src_df[rownames(src_df) %in% chosen_genes, , drop = FALSE]
      gene_df$GeneID   <- rownames(gene_df)
      gene_df$Contrast <- focal
      gene_df <- gene_df[, c("GeneID", "Contrast", setdiff(names(gene_df), c("GeneID", "Contrast")))]
      rownames(gene_df) <- NULL
      
    } else {
      # Shared set: wide table with concordance annotation
      src_names <- if (choice == "ALL_SHARED") names(gs) else strsplit(sub("^PAIR__", "", choice), "\\|\\|\\|")[[1]]
      gene_df   <- build_shared_gene_table(chosen_genes, src_names, sig_dfs)
      
      if (is.null(gene_df) || nrow(gene_df) == 0) {
        showNotification("No matching genes found in contrast results.", type = "warning")
        cmpSelectedGenes(NULL)
        return()
      }
      
      # Apply concordance filter
      conc_choice <- if (!is.null(input$vennConcordance)) input$vennConcordance else "all"
      gene_df <- apply_concordance_filter(gene_df, conc_choice)
      
      if (nrow(gene_df) == 0) {
        showNotification("No genes match the selected direction filter.", type = "warning")
        cmpSelectedGenes(NULL)
        return()
      }
    }
    
    cmpSelectedGenes(gene_df)
    showNotification(paste(nrow(gene_df), "genes loaded."), type = "message")
  })
  
  output$vennGenesTable <- DT::renderDataTable({
    req(cmpSelectedGenes())
    DT::datatable(cmpSelectedGenes(), options = list(pageLength = 10), rownames = FALSE)
  })
  
  output$dl_venn_genes <- downloadHandler(
    filename = function() paste0("venn_genes_", Sys.Date(), ".csv"),
    content  = function(f) readr::write_csv(req(cmpSelectedGenes()), f)
  )
  
  # ---- Volcano plot coloured by Venn membership ----
  output$vennVolcContrastUI <- renderUI({
    req(vennSets())
    selectInput("vennVolcContrast", "Choose contrast to plot:", 
                choices = names(vennSets()$genesets))
  })
  
  # Context selector: which gene set scope to use for shared/unique classification
  output$vennVolcContextUI <- renderUI({
    req(vennSets())
    gs  <- vennSets()$genesets
    nms <- names(gs)
    
    choices <- c("All contrasts (any sharing counts)" = "ALL_OTHER")
    for (nm in nms) choices[paste0("Unique vs: ", nm, " only")] <- paste0("PAIR__", nm)
    if (length(nms) > 1) {
      pairs <- combn(nms, 2, simplify = FALSE)
      for (p in pairs) {
        lbl <- paste0("Shared between: ", paste(p, collapse = " & "))
        choices[lbl] <- paste0("CTX_PAIR__", paste(p, collapse = "|||"))
      }
    }
    choices["Shared across ALL contrasts"] <- "CTX_ALL_SHARED"
    
    tagList(
      selectInput("vennVolcContext", "Classify shared/unique relative to:",
                  choices = choices, selected = "ALL_OTHER"),
      helpText("Match this to your 'Select gene set' choice above for consistency.")
    )
  })
  
  observeEvent(input$runCmpVolcBtn, {
    req(vennSets(), analysisResults())
    output$cmpVolcanoPlot <- renderPlot({
      gs    <- vennSets()$genesets
      dds   <- analysisResults()$dds
      focal <- input$vennVolcContrast
      
      # Get the contrast spec
      all_contrasts <- c(
        list(list(label = paste0(input$contrastNumerator, "_vs_", input$contrastDenominator),
                  numerator = input$contrastNumerator, denominator = input$contrastDenominator,
                  lfc = input$lfcThreshold, padj = input$padjThreshold)),
        extraContrasts()
      )
      spec <- Filter(function(x) x$label == focal, all_contrasts)[[1]]
      
      res_c <- as.data.frame(results(dds,
                                     lfcThreshold = spec$lfc,
                                     altHypothesis = "greaterAbs",
                                     alpha = spec$padj,
                                     contrast = c("Group", spec$numerator, spec$denominator)))
      res_c <- res_c[!is.na(res_c$padj) & !is.na(res_c$log2FoldChange), ]
      res_c$GeneID <- rownames(res_c)
      
      # Compute -log10(padj) with capping for extreme values
      # Replace padj == 0 with the smallest non-zero padj so -log10 is finite
      min_nonzero_padj <- min(res_c$padj[res_c$padj > 0], na.rm = TRUE)
      res_c$padj_safe <- ifelse(res_c$padj == 0, min_nonzero_padj, res_c$padj)
      res_c$neg_log10_padj <- -log10(res_c$padj_safe)
      
      # Cap at a sensible ceiling and flag capped genes with a different shape
      y_cap <- quantile(res_c$neg_log10_padj, 0.999)
      y_cap <- max(y_cap, -log10(spec$padj) + 10)  # ensure cap is above the significance line
      res_c$capped <- res_c$neg_log10_padj > y_cap
      res_c$neg_log10_padj_plot <- pmin(res_c$neg_log10_padj, y_cap)
      
      # Determine shared vs unique membership with concordance for shared genes.
      # The context (input$vennVolcContext) controls which contrasts count as "other":
      #   ALL_OTHER      -> any other contrast that has the gene counts as sharing
      #   CTX_PAIR__A|||B -> only contrasts A and B count (pairwise scope)
      #   CTX_ALL_SHARED  -> only genes in ALL contrasts are "shared"; rest are unique
      #   PAIR__X         -> focal is compared only against X (gene unique if absent from X)
      
      vs_local      <- vennSets()
      sig_dfs_local <- vs_local$sig_dfs
      ctx           <- if (!is.null(input$vennVolcContext)) input$vennVolcContext else "ALL_OTHER"
      
      # Determine which other contrasts to check against, based on context
      all_other_nms <- setdiff(names(gs), focal)
      
      context_nms <- if (ctx == "ALL_OTHER") {
        all_other_nms
      } else if (ctx == "CTX_ALL_SHARED") {
        all_other_nms   # shared = in every other contrast; handled below
      } else if (startsWith(ctx, "CTX_PAIR__")) {
        pair_members <- strsplit(sub("^CTX_PAIR__", "", ctx), "\\|\\|\\|")[[1]]
        setdiff(pair_members, focal)   # only the other member(s) of this pair
      } else if (startsWith(ctx, "PAIR__")) {
        # "Unique vs X only" - compare focal against just X
        nm_x <- sub("^PAIR__", "", ctx)
        if (nm_x == focal) all_other_nms else c(nm_x)
      } else {
        all_other_nms
      }
      
      focal_genes <- gs[[focal]]
      concordant_up <- concordant_down <- discordant <- unique_genes <- character(0)
      
      for (g in focal_genes) {
        focal_lfc <- res_c$log2FoldChange[res_c$GeneID == g]
        if (length(focal_lfc) == 0 || is.na(focal_lfc)) next
        focal_dir <- sign(focal_lfc)
        
        # Check presence/direction in the context contrasts only
        other_dirs <- sapply(context_nms, function(nm) {
          v <- sig_dfs_local[[nm]]$log2FoldChange[rownames(sig_dfs_local[[nm]]) == g]
          if (length(v) == 0) return(NA_real_) else sign(v[1])
        })
        
        if (ctx == "CTX_ALL_SHARED") {
          # "Shared" only if present in ALL context contrasts
          if (any(is.na(other_dirs))) {
            unique_genes <- c(unique_genes, g)
          } else if (all(other_dirs == focal_dir)) {
            if (focal_dir > 0) concordant_up   <- c(concordant_up,   g)
            else               concordant_down <- c(concordant_down, g)
          } else {
            discordant <- c(discordant, g)
          }
        } else {
          # "Shared" if present in at least one context contrast
          other_dirs_present <- other_dirs[!is.na(other_dirs)]
          if (length(other_dirs_present) == 0) {
            unique_genes <- c(unique_genes, g)
          } else if (all(other_dirs_present == focal_dir)) {
            if (focal_dir > 0) concordant_up   <- c(concordant_up,   g)
            else               concordant_down <- c(concordant_down, g)
          } else {
            discordant <- c(discordant, g)
          }
        }
      }
      
      res_c$Category <- "Not significant"
      res_c$Category[res_c$GeneID %in% unique_genes]    <- paste0("Unique to ", focal)
      res_c$Category[res_c$GeneID %in% concordant_up]   <- "Shared concordant UP ↑"
      res_c$Category[res_c$GeneID %in% concordant_down] <- "Shared concordant DOWN ↓"
      res_c$Category[res_c$GeneID %in% discordant]      <- "Shared discordant ↕"
      
      # Top genes to label
      sig_rows <- res_c[res_c$Category != "Not significant", ]
      top_lab  <- head(sig_rows[order(sig_rows$padj), "GeneID"], input$cmpVolcTop)
      res_c$label <- ifelse(res_c$GeneID %in% top_lab, res_c$GeneID, "")
      
      unique_col_name <- paste0("Unique to ", focal)
      cat_colors <- c(
        "Not significant"                 = "grey75",
        "Shared concordant UP ↑"     = "#2166ac",
        "Shared concordant DOWN ↓"   = "#d6604d",
        "Shared discordant ↕"        = "#9970ab"
      )
      cat_colors[unique_col_name] <- "#f4a736"
      
      ggplot(res_c, aes(x = log2FoldChange, y = neg_log10_padj_plot, colour = Category, label = label)) +
        geom_point(aes(shape = capped), size = input$cmpVolcPtSz, alpha = 0.7) +
        scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17),
                           labels = c("FALSE" = "Within range", "TRUE" = "Capped (more significant)"),
                           guide = guide_legend(title = NULL, order = 2)) +
        ggrepel::geom_text_repel(size = 3, max.overlaps = 20, show.legend = FALSE) +
        scale_colour_manual(values = cat_colors) +
        geom_vline(xintercept = c(-spec$lfc, spec$lfc), linetype = "dashed", colour = "grey40") +
        geom_hline(yintercept = -log10(spec$padj), linetype = "dashed", colour = "grey40") +
        theme_minimal(base_size = 13) +
        labs(title = paste("Volcano:", focal),
             x = "log2 Fold Change", y = "-log10(adjusted p-value)",
             colour = "Gene category")
    })
  })
  
  # ---- Downstream selector UI ----
  output$dsSetSelectorUI <- renderUI({
    req(vennSets())
    gs  <- vennSets()$genesets
    nms <- names(gs)
    
    tagList(
      radioButtons("dsSetMode", "Gene set type:",
                   choices = c("Shared across ALL contrasts" = "all_shared",
                               "Shared between selected contrasts" = "selected_shared",
                               "Unique to one contrast" = "unique"),
                   selected = "all_shared"),
      conditionalPanel(
        condition = "input.dsSetMode == 'selected_shared'",
        checkboxGroupInput("dsSharedContrasts", "Select contrasts (2 or more):",
                           choices = nms, selected = nms[1:min(2, length(nms))]),
        helpText("Exclusive overlap: genes in ALL ticked contrasts and NONE of the unticked.")
      ),
      conditionalPanel(
        condition = "input.dsSetMode == 'unique'",
        selectInput("dsUniqueContrast", "Unique to:", choices = nms, selected = nms[1])
      ),
      conditionalPanel(
        condition = "input.dsSetMode != 'unique'",
        selectInput("dsConcordance", "Direction filter for downstream analyses:",
                    choices = c(
                      "All shared (average stats)"        = "all",
                      "Concordant UP in all contrasts"    = "up",
                      "Concordant DOWN in all contrasts"  = "down",
                      "Discordant only (mixed direction)" = "discordant"
                    ),
                    selected = "all"),
        helpText("Filters which genes are passed to GO / GSEA / GAGE based on directional concordance across contrasts.")
      ),
      p(em("Tip: load genes above first to preview concordance before running analyses."))
    )
  })
  
  # Reactive to build the choice string from the downstream UI inputs
  dsSetChoice <- reactive({
    mode <- input$dsSetMode
    if (is.null(mode)) return(NULL)
    if (mode == "all_shared") return("ALL_SHARED")
    if (mode == "unique") {
      req(input$dsUniqueContrast)
      return(paste0("UNIQUE__", input$dsUniqueContrast))
    }
    if (mode == "selected_shared") {
      sel <- input$dsSharedContrasts
      validate(need(length(sel) >= 2, "Please select at least 2 contrasts."))
      return(paste0("PAIR__", paste(sel, collapse = "|||")))
    }
    NULL
  })
  
  # Helper: get ENTREZ IDs for a selected Venn region
  get_entrez_for_set <- function(choice) {
    req(vennSets(), analysisResults())
    vs         <- vennSets()
    gs         <- vs$genesets
    sig_dfs    <- vs$sig_dfs
    symbols    <- extract_venn_genes(gs, choice)
    res_entrez <- analysisResults()$res_entrez  # kept for gene_lengths / universe
    orgdb_name <- isolate(input$organism)
    orgdb      <- get(orgdb_name)
    
    from_type <- detect_gene_id_type(symbols)$orgdb
    
    mapped <- clusterProfiler::bitr(symbols, fromType = from_type, toType = "ENTREZID", OrgDb = orgdb)
    
    # Build a stats subset using the correct contrast sig_dfs
    # For downstream ranking (GSEA/GAGE), pool genes from all relevant contrasts
    # and use mean stat/LFC when a gene appears in multiple contrasts.
    if (startsWith(choice, "UNIQUE__")) {
      focal    <- sub("^UNIQUE__", "", choice)
      src_df   <- sig_dfs[[focal]]
      stats_df <- src_df[rownames(src_df) %in% symbols, c("log2FoldChange", "stat", "padj"), drop = FALSE]
      stats_df$GeneIDs <- rownames(stats_df)
    } else {
      src_names <- if (choice == "ALL_SHARED") names(gs) else strsplit(sub("^PAIR__", "", choice), "\\|\\|\\|")[[1]]
      
      # Build wide concordance table to apply direction filter
      wide_tbl <- build_shared_gene_table(symbols, src_names, sig_dfs)
      
      # Apply concordance filter (uses input$vennConcordance if available)
      conc_choice <- if (!is.null(input$dsConcordance)) input$dsConcordance else "all"
      wide_tbl <- apply_concordance_filter(wide_tbl, conc_choice)
      
      if (is.null(wide_tbl) || nrow(wide_tbl) == 0) {
        validate(need(FALSE, "No genes remain after concordance filter. Try a less restrictive direction filter."))
      }
      
      # Use mean stats for ranking in GSEA/GAGE
      stats_df <- data.frame(
        GeneIDs        = wide_tbl$GeneID,
        log2FoldChange = wide_tbl$mean_log2FoldChange,
        stat           = wide_tbl$mean_stat,
        padj           = wide_tbl$mean_padj,
        stringsAsFactors = FALSE
      )
      rownames(stats_df) <- stats_df$GeneIDs
      # Restrict symbols to filtered set
      symbols <- wide_tbl$GeneID
    }
    
    # Merge with ENTREZID mapping
    res_entrez_subset <- merge(stats_df, mapped, by.x = "GeneIDs", by.y = from_type, all.x = FALSE)
    # Also pull gene_lengths from the master res_entrez
    res_entrez_subset <- merge(res_entrez_subset, 
                               res_entrez[, c("GeneIDs", "gene_lengths")], 
                               by = "GeneIDs", all.x = TRUE)
    
    list(symbols = symbols, entrez = mapped$ENTREZID, mapped = mapped,
         res_entrez_subset = res_entrez_subset)
  }
  
  # ---- Compare GO ----
  observeEvent(input$runCmpGOBtn, {
    req(vennSets(), analysisResults())
    withProgress(message = "Running GO on selected gene set...", value = 0, {
      info <- get_entrez_for_set(dsSetChoice())
      sp   <- isolate(goSpeciesCode())
      
      universe_entrez <- analysisResults()$res_entrez$ENTREZID
      gene_lengths_vec <- as.numeric(analysisResults()$res_entrez$gene_lengths)
      names(gene_lengths_vec) <- analysisResults()$res_entrez$ENTREZID
      
      # If lengths are unavailable (plain count table), set covariate to NULL
      cov_vec <- gene_lengths_vec[universe_entrez]
      covariate_arg <- if (all(is.na(cov_vec))) NULL else cov_vec
      
      if (is.null(covariate_arg)) {
        showNotification("GO running without length-bias correction (no gene lengths available).",
                         type = "warning", duration = 6)
      }
      
      validate(need(length(info$entrez) >= 2, "Fewer than 2 genes with Entrez IDs; cannot run GO."))
      
      go_res <- limma::goana(de = info$entrez, species = sp,
                             universe = universe_entrez,
                             covariate = covariate_arg)
      incProgress(0.9)
      
      output$cmpGOplot <- renderPlot({
        topgo <- limma::topGO(go_res, ontology = input$cmpGOont, number = input$cmpGOnum)
        ggplot(topgo, aes(x = reorder(Term, -log10(P.DE)), y = -log10(P.DE))) +
          geom_bar(stat = "identity", fill = "#638475") +
          coord_flip() + theme_minimal() +
          labs(x = "GO Term", y = "-log10(p-value)",
               title = paste("GO", input$cmpGOont, "–", dsSetChoice()))
      })
      
      output$cmpGOtable <- DT::renderDataTable({
        topgo <- limma::topGO(go_res, ontology = input$cmpGOont, number = Inf)
        DT::datatable(topgo, options = list(pageLength = 10))
      })
      
      output$dl_cmp_go <- downloadHandler(
        filename = function() paste0("cmp_GO_", Sys.Date(), ".csv"),
        content  = function(f) readr::write_csv(as.data.frame(go_res), f)
      )
      showNotification("GO complete.", type = "message")
    })
  })
  
  # ---- Compare GSEA ----
  observeEvent(input$runCmpGseaBtn, {
    req(vennSets(), analysisResults())
    withProgress(message = "Running GSEA on selected gene set...", value = 0, {
      info       <- get_entrez_for_set(dsSetChoice())
      orgdb_name <- isolate(input$organism)
      orgdb      <- get(orgdb_name)
      
      sub_res <- info$res_entrez_subset
      validate(need(nrow(sub_res) >= input$cmpGseaMin, 
                    "Not enough genes in this set for GSEA."))
      
      metric_col <- input$cmpGseaMetric
      m <- sub_res[[metric_col]]
      names(m) <- sub_res$ENTREZID
      m <- m[is.finite(m)]
      m <- tapply(m, names(m), max)
      m <- sort(setNames(as.vector(m), names(m)), decreasing = TRUE)
      
      validate(need(length(m) >= input$cmpGseaMin,
                    "Not enough ranked genes to run GSEA."))
      
      # Run the appropriate GSEA — each branch wrapped in tryCatch
      gse <- NULL
      
      if (input$cmpGseaType == "GO") {
        gse <- tryCatch({
          clusterProfiler::gseGO(
            geneList     = m,
            OrgDb        = orgdb,
            keyType      = "ENTREZID",
            ont          = input$cmpGseaOnt,
            minGSSize    = input$cmpGseaMin,
            maxGSSize    = input$cmpGseaMax,
            pvalueCutoff = input$cmpGseaP,
            verbose      = FALSE
          )
        }, error = function(e) {
          showNotification(paste("GO GSEA failed:", e$message),
                           type = "error", duration = NULL)
          NULL
        })
        
      } else if (input$cmpGseaType == "Wiki Pathways") {
        validate(need(nzchar(input$cmpWikiorg),
                      "Please select an organism for WikiPathways."))
        gse <- tryCatch({
          clusterProfiler::gseWP(
            geneList     = m,
            organism     = input$cmpWikiorg,
            minGSSize    = input$cmpGseaMin,
            maxGSSize    = input$cmpGseaMax,
            pvalueCutoff = input$cmpGseaP,
            verbose      = FALSE
          )
        }, error = function(e) {
          showNotification(paste("WikiPathways GSEA failed:", e$message),
                           type = "error", duration = NULL)
          NULL
        })
        
      } else if (input$cmpGseaType == "KEGG") {
        gse <- tryCatch({
          clusterProfiler::gseKEGG(
            geneList     = m,
            organism     = input$cmpKeggorg,
            keyType      = "ncbi-geneid",
            minGSSize    = input$cmpGseaMin,
            maxGSSize    = input$cmpGseaMax,
            pvalueCutoff = input$cmpGseaP,
            verbose      = FALSE
          )
        }, error = function(e) {
          showNotification(paste("KEGG GSEA failed:", e$message),
                           type = "error", duration = NULL)
          NULL
        })
      }
      
      # Store result (possibly NULL) in reactiveVal
      cmpGseaRes(gse)
      
      # Update term selector only if we have results
      if (!is.null(gse)) {
        df <- as.data.frame(gse)
        updateSelectizeInput(session, "cmpGseaTerm",
                             choices = df$ID, selected = head(df$ID, 1), server = TRUE)
        showNotification("GSEA complete.", type = "message")
      } else {
        updateSelectizeInput(session, "cmpGseaTerm", choices = NULL)
      }
      
      output$cmpGseaDotplot <- renderPlot({
        gse_local <- cmpGseaRes()
        validate(need(!is.null(gse_local), "No GSEA results to display."))
        df_check <- as.data.frame(gse_local)
        validate(need(nrow(df_check) > 0, "No significant GSEA terms found."))
        n_show <- min(input$cmpGseaNum, nrow(df_check))
        enrichplot::dotplot(gse_local, showCategory = n_show)
      })
      output$cmpGseaEnrichPlot <- renderPlot({
        gse_local <- cmpGseaRes()
        validate(
          need(!is.null(gse_local), "No GSEA results to display."),
          need(input$cmpGseaTerm, "Select a term to display.")
        )
        df_check <- as.data.frame(gse_local)
        validate(need(nrow(df_check) > 0, "No significant GSEA terms found."))
        validate(need(input$cmpGseaTerm %in% df_check$ID, "Selected term not found in results."))
        enrichplot::gseaplot2(gse_local, geneSetID = input$cmpGseaTerm, title = input$cmpGseaTerm)
      })
      output$cmpGseaTable <- DT::renderDataTable({
        gse_local <- cmpGseaRes()
        validate(need(!is.null(gse_local), "No GSEA results to display."))
        DT::datatable(as.data.frame(gse_local), options = list(pageLength = 10))
      })
      output$dl_cmp_gsea <- downloadHandler(
        filename = function() paste0("cmp_GSEA_", Sys.Date(), ".csv"),
        content  = function(f) {
          gse_local <- cmpGseaRes()
          if (is.null(gse_local)) return()
          readr::write_csv(as.data.frame(gse_local), f)
        }
      )
      incProgress(0.9)
    })
  })
  
  # ---- Compare GAGE ----
  observeEvent(input$runCmpGageBtn, {
    req(vennSets(), analysisResults())
    withProgress(message = "Running GAGE on selected gene set...", value = 0, {
      info       <- get_entrez_for_set(dsSetChoice())
      orgdb_name <- isolate(input$organism)
      
      # Choose KEGG prefix
      sp_prefix <- switch(orgdb_name,
                          "org.Rn.eg.db" = "rno",
                          "org.Mm.eg.db" = "mmu",
                          "org.Hs.eg.db" = "hsa",
                          "org.Dm.eg.db" = "dme",
                          "org.Dr.eg.db" = "dre",
                          "hsa"
      )
      
      kg <- gage::kegg.gsets(sp_prefix, id.type = "entrez")
      gsets <- kg$kg.sets[kg$sigmet.idx]
      
      sub_res    <- info$res_entrez_subset
      metric_col <- input$cmpGageMetric
      fc_vec     <- sub_res[[metric_col]]
      names(fc_vec) <- sub_res$ENTREZID
      fc_vec <- fc_vec[is.finite(fc_vec)]
      
      validate(need(length(fc_vec) >= input$cmpGageMin, 
                    "Not enough genes with Entrez IDs for GAGE."))
      
      qcut <- input$cmpGageP
      gr <- gage::gage(fc_vec, gsets = gsets,
                       same.dir   = input$cmpGageSameDir,
                       set.size   = c(input$cmpGageMin, input$cmpGageMax))
      
      sig_up   <- as.data.frame(gr$greater)
      sig_down <- as.data.frame(gr$less)
      
      if ("q.val" %in% names(sig_up))   sig_up   <- sig_up[!is.na(sig_up$q.val)   & sig_up$q.val   <= qcut, , drop = FALSE]
      if ("q.val" %in% names(sig_down)) sig_down <- sig_down[!is.na(sig_down$q.val) & sig_down$q.val <= qcut, , drop = FALSE]
      
      # Safe column assignment - guard against 0-row data frames
      if (nrow(sig_up) > 0) {
        rn_up <- rownames(sig_up)
        sig_up$PathwayID   <- substr(rn_up, 1, 8)
        sig_up$Direction   <- "Up"
        sig_up$Description <- trimws(sub("^\\S+\\s*", "", rn_up))
      }
      if (nrow(sig_down) > 0) {
        rn_down <- rownames(sig_down)
        sig_down$PathwayID   <- substr(rn_down, 1, 8)
        sig_down$Direction   <- "Down"
        sig_down$Description <- trimws(sub("^\\S+\\s*", "", rn_down))
      }
      cmpGageRes(list(up = sig_up, down = sig_down, fc = fc_vec, species = sp_prefix))
      
      all_ids <- unique(c(sig_up$PathwayID, sig_down$PathwayID))
      updateSelectInput(session, "cmpKeggPathway",
                        choices  = all_ids,
                        selected = if (length(all_ids)) all_ids[[1]] else NULL)
      
      output$cmpGageUpTable <- DT::renderDataTable({
        if (!NROW(sig_up)) return(DT::datatable(data.frame(Message = "None"), options = list(dom = "t"), rownames = FALSE))
        DT::datatable(sig_up, options = list(pageLength = 10), rownames = FALSE)
      })
      output$cmpGageDownTable <- DT::renderDataTable({
        if (!NROW(sig_down)) return(DT::datatable(data.frame(Message = "None"), options = list(dom = "t"), rownames = FALSE))
        DT::datatable(sig_down, options = list(pageLength = 10), rownames = FALSE)
      })
      
      incProgress(0.9)
      showNotification("GAGE complete.", type = "message")
    })
  })
  
  # ---- Pathview for Compare tab ----
  observeEvent(input$drawCmpPathviewBtn, {
    req(cmpGageRes())
    output$cmpPathviewImg <- renderImage({
      gr  <- cmpGageRes()
      fc  <- gr$fc
      sp  <- gr$species
      pid <- substr(req(input$cmpKeggPathway), 1, 8)
      validate(need(nzchar(pid), "Select a pathway."))
      validate(need(length(fc) > 0, "No fold change data."))
      
      tmpdir <- tempfile("cmp_pv_"); dir.create(tmpdir)
      oldwd  <- getwd(); setwd(tmpdir); on.exit(setwd(oldwd), add = TRUE)
      
      try(pathview::pathview(gene.data = fc, pathway.id = pid, species = sp,
                             gene.idtype = "entrez", kegg.native = TRUE,
                             out.suffix = "cmp_shiny"), silent = TRUE)
      
      pngs    <- list.files(tmpdir, pattern = "\\.png$", full.names = TRUE)
      overlay <- pngs[grepl(paste0("^", pid, ".*pathview.*\\.png$"), basename(pngs))]
      png_f   <- if (length(overlay)) overlay[which.max(file.mtime(overlay))] else pngs[which.max(file.mtime(pngs))]
      validate(need(length(png_f) && file.exists(png_f), "Pathview failed for this pathway."))
      list(src = png_f, contentType = "image/png", alt = pid)
    }, deleteFile = FALSE)
  })
  
  output$dl_cmp_pathview <- downloadHandler(
    filename = function() paste0("cmp_pathview_", substr(req(input$cmpKeggPathway), 1, 8), ".png"),
    content  = function(file) {
      gr  <- req(cmpGageRes()); fc <- gr$fc; sp <- gr$species
      pid <- substr(req(input$cmpKeggPathway), 1, 8)
      tmpdir <- tempfile("cmp_pv_dl_"); dir.create(tmpdir)
      oldwd  <- getwd(); setwd(tmpdir); on.exit(setwd(oldwd), add = TRUE)
      pathview::pathview(gene.data = fc, pathway.id = pid, species = sp,
                         kegg.native = TRUE, out.suffix = "cmp_dl")
      pngs    <- list.files(tmpdir, pattern = "\\.png$", full.names = TRUE)
      overlay <- pngs[grepl(paste0("^", pid, ".*pathview.*\\.png$"), basename(pngs))]
      png_f   <- if (length(overlay)) overlay[which.max(file.mtime(overlay))] else pngs[which.max(file.mtime(pngs))]
      file.copy(png_f, file, overwrite = TRUE)
    }
  )
  
  # ============================================================
  # END Compare Contrasts
  # ============================================================
  
  output$analysisReady <- reactive({ !is.null(analysisResults()) })
  outputOptions(output, "analysisReady", suspendWhenHidden = FALSE)
  outputOptions(output, "heatmapPlot", suspendWhenHidden = FALSE)
  outputOptions(output, "pcaPlot", suspendWhenHidden = FALSE)
  outputOptions(output, "volcanoPlot", suspendWhenHidden = FALSE)
  outputOptions(output, "deTable", suspendWhenHidden = FALSE)
  outputOptions(output, "geneBoxplot", suspendWhenHidden = TRUE)
  outputOptions(output, "GOBarplot", suspendWhenHidden = FALSE)
  outputOptions(output, "gseaEnrichPlot", suspendWhenHidden = FALSE)
  
}


# Run the application 
shinyApp(ui = ui, server = server)