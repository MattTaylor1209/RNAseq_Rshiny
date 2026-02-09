#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

# Required packages

library(apeglm)
library(ashr)
library(BiasedUrn)
library(biomaRt)
library(clusterProfiler)
library(DBI)
library(DESeq2)
library(dplyr)
library(DT)
library(edgeR)
library(enrichplot)
library(eulerr)
library(gage)
library(gageData)
library(GenomicFeatures)
library(ggprism)
library(ggpubr)
library(ggrepel)
library(ggsci)
library(ggsignif)
library(ggtext)
library(ggthemes)
library(ggVennDiagram)
library(GO.db)
library(goseq)
library(gplots)
library(grid)
library(GSEABase)
library(htmlwidgets)
library(KEGGREST)
library(limma)
library(lsmeans)
library(markdown)
library(minpack.lm)
library(msigdbr)
library(multcomp)
library(NMF)
library(org.Dm.eg.db)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(patchwork)
library(pathview)
library(pheatmap)
library(plotly)
library(qusage)
library(RColorBrewer)
library(readr)
library(readxl)
library(rstatix)
library(Rsubread)
library(shiny)
library(shinyFiles)
library(stats)
library(sva)
library(tidyr)
library(tidyverse)
library(treemap)
library(utils)
library(VennDiagram)
library(viridisLite)
library(vsn)




# Increase upload size 

options(shiny.maxRequestSize = 50 * 1024^2)  # 50 MB


ui <- fluidPage(
  titlePanel("RNA-seq Analysis from featureCounts (.rds)"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("rdsFile", "Upload featureCounts .rds File", accept = ".rds"),
      fileInput("sampleInfoFile", "Upload SampleInfo.txt", accept = c(".txt", ".tsv")),
      selectInput("organism", "Select organism:",
                  choices = c(
                    "Rat (Rattus norvegicus)" = "org.Rn.eg.db",
                    "Mouse (Mus musculus)" = "org.Mm.eg.db",
                    "Human (Homo sapiens)" = "org.Hs.eg.db",
                    "Drosophila (Drosophila melanogaster)" = "org.Dm.eg.db"
                  ),
                  selected = "org.Rn.eg.db"
      ),
      numericInput("lfcThreshold", "log2 Fold Change threshold", value = 1, min = 0),
      numericInput("padjThreshold", "Adjusted p-value threshold (FDR)", value = 0.1, min = 0, max = 1),
      uiOutput("groupOrderUI"),
      uiOutput("contrastSelectUI"),
      uiOutput("flybaseCheckboxUI"),
      actionButton("analyzeBtn", "Run Analysis"),
      conditionalPanel(
        condition = "output.inputsReady",
        checkboxInput("viewPreview", "Preview input files", value = TRUE)
      ),
      conditionalPanel(
        condition = "output.analysisReady",
        downloadButton("dl_res", "Download DE table")
      ),
      tags$hr(),
      h4("GSEA"),
      selectInput("gseaOnt", "Ontology", c("BP","MF","CC"), selected = "BP"),
      selectInput("gseaMetric", "Rank by", c("stat","log2FoldChange"), selected = "stat"),
      numericInput("gseaP", "p-value cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
      numericInput("gseaMin", "Min gene set size", value = 10, min = 5, step = 5),
      numericInput("gseaMax", "Max gene set size", value = 500, min = 50, step = 50),
      numericInput("numcategories", "No. of categories to show", value = 10, min = 1, step = 1),
      actionButton("runGSEA", "Run GSEA"),
      tags$hr(),
      h4("GO"),
      actionButton("runGO", "Run GO")
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
          tabPanel("DE Results", dataTableOutput("deTable")),
          tabPanel("Heatmap",
                   fluidRow(
                     column(3,
                            numericInput("heatmapTopN", "Number of top DE genes", 
                                         value = 50, min = 10, max = 500, step = 10)
                     ),
                     column(3,
                            selectInput("heatmapSortBy", "Sort genes by:",
                                        choices = c("Adjusted p-value" = "padj", 
                                                    "log2 Fold Change (absolute)" = "abs_lfc",
                                                    "Test statistic (absolute)" = "abs_stat"),
                                        selected = "padj")
                     ),
                     column(3,
                            selectInput("heatmapScale", "Scale data:",
                                        choices = c("Row (gene)" = "row",
                                                    "Column (sample)" = "column", 
                                                    "None" = "none"),
                                        selected = "row")
                     ),
                     column(3,
                            checkboxInput("heatmapClusterRows", "Cluster rows (genes)", value = TRUE)
                     )
                   ),
                   fluidRow(
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
          tabPanel("GSEA",
                   plotOutput("gseaDotplot"),
                   selectizeInput("gseaTerm", "Show enrichment plot for:", choices = NULL, multiple = FALSE),
                   plotOutput("gseaEnrichPlot"),
                   DT::dataTableOutput("gseaTable"),
                   downloadButton("dl_gsea", "Download GSEA table")
          ),
          tabPanel("GO",
                   selectInput("GOontology", "Ontology", c("BP","MF","CC"), selected = "BP"),
                   numericInput("GOnumber", "No. of categories to show", value = 10, min = 1, step = 1),
                   plotOutput("GOBarplot")
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
                            h4("Significant pathways (Downregulated)"),
                            DT::dataTableOutput("gageDownTable")
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
                   plotOutput("standardvolcanoplot",
                              width = "90%", height = "500px")
                   
          ),
          tabPanel("Interactive Volcano Plot", plotlyOutput("volcanoPlot"),
                   plotOutput("geneBoxplot"))
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
  
  # Load uploaded featureCounts data
  counts_data <- reactive({
    req(input$rdsFile)
    readRDS(input$rdsFile$datapath)
  })
  
  # Load uploaded sample metadata
  sample_info <- reactive({
    req(input$sampleInfoFile)
    read.delim(input$sampleInfoFile$datapath, stringsAsFactors = TRUE)
  })
  
  # Tell UI when both files are ready
  output$inputsReady <- reactive({
    !is.null(input$rdsFile) && !is.null(input$sampleInfoFile)
  })
  outputOptions(output, "inputsReady", suspendWhenHidden = FALSE)
  
  # Reactive to get unique group names
  group_levels <- reactive({
    req(sample_info())
    unique(as.character(sample_info()$Group))
  })
  
  observeEvent(sample_info(), {
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
  
  # Preview featureCounts input
  output$previewTable <- renderTable({
    counts <- counts_data()$counts
    colnames(counts) <- sample_info()$SampleName
    head_df <- head(counts)
    head_df <- cbind(GeneID = rownames(head_df), head_df)
    head_df
  })
  
  # Preview sample info
  output$sampleInfo <- renderTable({
    sample_info()
  })
  
  output$flybaseCheckboxUI <- renderUI({
    req(input$organism)
    if (input$organism == "org.Dm.eg.db") {
      checkboxInput("convertFlybase", "Convert FlyBase IDs?", value = TRUE)
    }
  })
  
  
  analysisResults <- eventReactive(input$analyzeBtn, {
    withProgress(message = "Running RNA-seq analysis", value = 0, {
      orgdb_name <- isolate(input$organism)
      orgdb <- isolate(get(orgdb_name))
      
      fc <- counts_data()
      
      # Decide what your GeneID column represents
      from_type <- if (orgdb_name == "org.Dm.eg.db" && isTRUE(input$convertFlybase)) {
        "FLYBASE"
      } else {
        "SYMBOL"
      }
      
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
      colnames(counts) <- sample_info()$SampleName
      col_data <- sample_info()
      
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
      
      # Convert FLYBASE identifiers to symbols if checkbox ticked:
      
      if (orgdb_name == "org.Dm.eg.db") {
        if (isTRUE(isolate(input$convertFlybase))) {
          # Get the mapped gene symbols
          gene_symbols <- mapIds(
            orgdb,
            keys = rownames(dds),
            column = "SYMBOL",
            keytype = "FLYBASE",
            multiVals = "first"
          )
          
          # Drop any rows with NA gene symbols
          keep_idx <- !is.na(gene_symbols)
          dds <- dds[keep_idx, ]
          rownames(dds) <- gene_symbols[keep_idx]
        }
      }
      incProgress(0.1)
      
      # filter our low counts
      appendLog("Filtering low-expression genes...")
      showNotification("Filtering low-expression genes...", type="message")
      keep <- edgeR::filterByExpr(counts(dds), group = dds$Group)
      
      dds <- dds[keep,]
      incProgress(0.1)
      
      appendLog("Running DESeq()...")
      showNotification("Running DESeq()...", type="message")
      dds <- DESeq(dds)
      
      incProgress(0.2)
      
      appendLog("Getting results...")
      showNotification("Getting results...", type="message")
      res <- results(dds, lfcThreshold = isolate(input$lfcThreshold), alpha =
                       isolate(input$padjThreshold), 
                     contrast = c("Group", isolate(input$contrastNumerator), isolate(input$contrastDenominator)))
      
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
      
      from_type <- if (orgdb_name == "org.Dm.eg.db") {
        if (isTRUE(isolate(input$convertFlybase))) {
          "SYMBOL"
        } else {
          "FLYBASE"
        }
      } else {
        "SYMBOL"
      }
      
      
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
      pca_df$SampleName <- sample_info()$SampleName  # Add sample names
      pca_df$Group <- sample_info()$Group  # Add group info
      
      
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
                      min.segment.length = 0) +  # Add sample labels
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
    
    # Convert results to dataframe and add sorting columns
    res_df <- as.data.frame(res)
    res_df$abs_lfc <- abs(res_df$log2FoldChange)
    res_df$abs_stat <- abs(res_df$stat)
    
    # Sort by selected metric and get top N genes
    top_genes <- switch(input$heatmapSortBy,
                        "padj" = head(rownames(res_df[order(res_df$padj, na.last = TRUE), ]), input$heatmapTopN),
                        "abs_lfc" = head(rownames(res_df[order(res_df$abs_lfc, decreasing = TRUE, na.last = TRUE), ]), input$heatmapTopN),
                        "abs_stat" = head(rownames(res_df[order(res_df$abs_stat, decreasing = TRUE, na.last = TRUE), ]), input$heatmapTopN)
    )
    
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
    
    # Create annotation colors - handle cases with <3 groups
    n_groups <- length(unique(sample_annotation$Group))
    if (n_groups == 1) {
      group_colors <- "#66C2A5"  # Single color
    } else if (n_groups == 2) {
      group_colors <- c("#66C2A5", "#FC8D62")  # Two colors from Set2
    } else {
      group_colors <- RColorBrewer::brewer.pal(min(n_groups, 8), "Set2")
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
      main = paste("Top", length(top_genes), "DE genes -", 
                   switch(input$heatmapSortBy,
                          "padj" = "sorted by adjusted p-value",
                          "abs_lfc" = "sorted by absolute log2 fold change", 
                          "abs_stat" = "sorted by absolute test statistic")),
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
      
      # 3) Call gseGO; add eps=0 if your clusterProfiler supports it
      gse <- clusterProfiler::gseGO(
        geneList     = m,
        OrgDb        = orgdb,
        keyType      = "ENTREZID",
        ont          = input$gseaOnt,
        minGSSize    = input$gseaMin,
        maxGSSize    = input$gseaMax,
        pvalueCutoff = input$gseaP,
        verbose      = TRUE
      )
      
      incProgress(0.8)
      appendLog("GSEA analysis complete.")
      showNotification("GSEA analysis complete.", type="message")
      
      
      gse
    })
  })
  
  observeEvent(gseaResults(), {
    df <- as.data.frame(gseaResults())
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
    req(gseaResults(), input$gseaTerm)
    enrichplot::gseaplot2(gseaResults(), geneSetID = input$gseaTerm, title = input$gseaTerm)
  })
  
  output$gseaTable <- DT::renderDataTable({
    req(gseaResults())
    DT::datatable(as.data.frame(gseaResults()),
                  extensions = "Buttons",
                  options = list(dom = "Bfrtip", buttons = c("copy","csv")))
  })
  
  output$dl_gsea <- downloadHandler(
    filename = function() "gsea_results.csv",
    content  = function(f) readr::write_csv(as.data.frame(gseaResults()), f)
  )
  
  ###---GO---###
  
  goSpeciesCode <- reactive({
    org <- isolate(input$organism)
    if (org == "org.Hs.eg.db") return("Hs")
    if (org == "org.Mm.eg.db") return("Mm")
    if (org == "org.Rn.eg.db") return("Rn")
    if (org == "org.Dm.eg.db") return("Dm")
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
      
      # Read species code 
      sp <- isolate(goSpeciesCode())
      
      # Run GO enrichment analysis with goana
      go_results <- limma::goana(de = entrez_id_vector, species = sp, 
                          universe = universe_entrez_ids, 
                          covariate = gene_lengths)
      
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
      theme_minimal() +
      labs(x = "GO Terms", y = "-log10(p-value)", title = paste("Top GO Terms", 
                                                                isolate(input$GOontology), sep = ": ")
           )
  })
  
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
        list(type = "line", x0 = -lfc_cutoff, x1 = -lfc_cutoff, y0 = 0, y1 = max(-log10(volcano_data$padj)), line = list(dash = "dash")),
        list(type = "line", x0 =  lfc_cutoff, x1 =  lfc_cutoff, y0 = 0, y1 = max(-log10(volcano_data$padj)), line = list(dash = "dash")),
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
    
    # Label for top genes
    # Select the top n most significant genes based on padj
    top_genes <- rownames(res_df[order(res_df$padj), ])[1:input$topgenes]
    
    # Add a 'label' column to label only the top n genes
    res_df$label <- ifelse(rownames(res_df) %in% top_genes, rownames(res_df), NA)
    
    volcano_colors <- c(
      "Downregulated" = "blue",
      "Not significant" = "grey",
      "Upregulated" = "red"
    )
    
    volcano_data <- res_df
    volcano_data$Gene <- rownames(volcano_data)
    
    ggplot(as.data.frame(res_df), aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(aes(color = significance), size = input$volcpointsize) +
      scale_color_manual(values = volcano_colors) +  
      theme_minimal() +
      theme(
        legend.title = element_text(face = "bold", size = input$volclegendtitlesize),
        legend.text = element_text(size = input$volclegendtextsize),
        axis.text = element_text(face = "bold", size = input$volcaxistextsize),
        axis.title = element_text(face = "bold", size = input$volcaxistitlesize)
      )+
      labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
      geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "black") +  # Horizontal line at p-value threshold
      geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "black") +  # Vertical lines at LFC = -1 and 1
      geom_text_repel(aes(label = label), size = input$volclabelsize, max.overlaps = 10)  # Add gene labels
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