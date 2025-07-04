#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

# Required packages

biocmanager_packages <- c(
  "limma",
  "edgeR",
  "Glimma",
  "org.Dm.eg.db",
  "org.Mm.eg.db",
  "org.Rn.eg.db",
  "org.Hs.eg.db",
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
  "gage",
  "gageData",
  "pathview",
  "sva",
  "RUVSeq",
  "vsn",
  "biomaRt"
)

other_packages <- c(
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
  "DT"
)

# Install packages not yet installed

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

installed_biocpackages <- biocmanager_packages %in% rownames(installed.packages())
if (any(installed_biocpackages == FALSE)) {
  BiocManager::install(biocmanager_packages[!installed_biocpackages], dependencies=TRUE)
}

installed_otherpackages <- other_packages %in% rownames(installed.packages())
if (any(installed_otherpackages == FALSE)) {
  install.packages(other_packages[!installed_otherpackages], dependencies=TRUE)
}

# Packages loading

invisible(lapply(biocmanager_packages, library, character.only = TRUE))
invisible(lapply(other_packages, library, character.only = TRUE))



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
      numericInput("padjThreshold", "Adjusted p-value threshold", value = 0.05, min = 0, max = 1),
      uiOutput("groupOrderUI"),
      uiOutput("contrastSelectUI"),
      uiOutput("flybaseCheckboxUI"),
      actionButton("analyzeBtn", "Run Analysis"),
      verbatimTextOutput("log"),
      tags$hr(),
      conditionalPanel(
        condition = "output.inputsReady",
        checkboxInput("viewPreview", "Preview input files", value = TRUE)
      )
    ),
    
    mainPanel(
      conditionalPanel(
        condition = "output.inputsReady",
        tabsetPanel(
          tabPanel("Input Preview", tableOutput("previewTable")),
          tabPanel("Sample Info", tableOutput("sampleInfo")),
          tabPanel("PCA", plotOutput("pcaPlot")),
          tabPanel("DE Results", dataTableOutput("deTable")),
          tabPanel("Pathways", plotOutput("pathwayPlot")),
          tabPanel("Volcano Plot", plotlyOutput("volcanoPlot"),
                   plotOutput("geneBoxplot"))
        )
      )
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
      keep <- filterByExpr(counts(dds), group = dds$Group)
      
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
      
      incProgress(0.1)
      
      appendLog("Computing vst...")
      showNotification("Computing vst...", type="message")
      
      vsd <- vst(dds)
      
      incProgress(0.1)
      
      # Create res object with rownames being ENTREZIDs
      
      # Step 1: Get gene symbols from DESeq2 results (res)
      uni_gene_symbols <- rownames(res)  # Get gene symbols from DESeq2 results
      
      # Step 2: Map gene symbols to Entrez IDs using bitr
      
      from_type <- if (orgdb_name == "org.Dm.eg.db") {
        if (isTRUE(isolate(input$convertFlybase))) {
          "SYMBOL"
        } else {
          "FLYBASE"
        }
      } else {
        "SYMBOL"
      }
      
      
      uni_entrez_ids <- bitr(uni_gene_symbols, fromType = from_type, toType = "ENTREZID", OrgDb = orgdb)
      
      # Step 3: Merge the DESeq2 results (res) with the Entrez IDs 
      res_entrez <- merge(as.data.frame(res), uni_entrez_ids, by.x = "row.names", by.y = from_type)
      
      # Step 4: rename the first column to "GeneIDs"
      res_entrez <- res_entrez %>% rename(GeneIDs = Row.names)
      
      
      appendLog("PCA analysis...")
      showNotification("PCA analysis...", type="message")
      
      rlogcounts <- rlog(counts)
      # run PCA
      pcDat <- prcomp(t(rlogcounts))
      # Calculate variance explained
      percentVar <- round(100 * (pcDat$sdev^2 / sum(pcDat$sdev^2)), 1)
      
      # Convert PCA results to a dataframe
      pca_df <- as.data.frame(pcDat$x)
      pca_df$SampleName <- sample_info()$SampleName  # Add sample names
      pca_df$Group <- sample_info()$Group  # Add group info
      
      incProgress(0.1)
      
      # GSEGO
      
      appendLog("GSE analysis...")
      showNotification("GSE analysis...", type="message")
      # Rank the genes by their log2 fold change
      ranked_genes <- res_entrez$log2FoldChange
      names(ranked_genes) <- res_entrez$ENTREZID
      
      # Remove any NA values from the ranked list
      ranked_genes <- na.omit(ranked_genes)
      
      # Sort the list by the ranking metric (log2FoldChange in this case)
      ranked_genes <- sort(ranked_genes, decreasing = TRUE)
      
      # Perform GSEA for GO Biological Process (BP) terms
      gsego_results <- gseGO(geneList = ranked_genes,
                             OrgDb = orgdb,    # Use the appropriate organism database
                             ont = "BP",             # Biological Process ontology
                             keyType = "ENTREZID",   # The gene identifiers used (ENTREZ IDs)
                             pvalueCutoff = 0.05,    # p-value threshold
                             verbose = FALSE,
                             eps = 1e-300)
      incProgress(0.1)
      
      
      appendLog("Analysis complete.")
      showNotification("Analysis complete.", type="message")
      
    })
    # Output from reactive expression
    list(dds = dds, res = res, vsd = vsd, pca_df = pca_df, percentVar = percentVar,
         res_entrez = res_entrez, gsego_results = gsego_results)
    
  })
  
  # PCA plot
  output$pcaPlot <- renderPlot({
    req(analysisResults())
    # Recall the variables from the reactive analysisResults event
    pca_df <- analysisResults()$pca_df
    percentVar <- analysisResults()$percentVar
    
    ggplot(pca_df, aes(x = PC1, y = PC2, fill = Group, shape = Group)) +
      geom_point(size = 5) +
      geom_text_repel(aes(label = SampleName), size = 4,
                      box.padding = 0.5,      # Increases space around labels
                      point.padding = 0.3,    # Increases distance from points
                      min.segment.length = 0) +  # Add sample labels
      scale_shape_manual(values = c(21:25, 21:25)) + # Change this depending on how many shapes you want, 21-25 are decent
      guides(fill = guide_legend(override.aes = list(shape = 22))) +
      theme_minimal() +
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
  
  output$pathwayPlot <- renderPlot({
    req(analysisResults())
    
    gsego_results <- analysisResults()$gsego_results
    
    dotplot(gsego_results, showCategory = 10)
  })
  
  # Volcano plot output
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
    
    plot_ly(
      data = volcano_data,
      x = ~log2FoldChange,
      y = ~log10padj,
      type = "scatter",
      mode = "markers",
      text = ~paste("Gene: ", Gene, "<br>Log2FC: ", signif(log2FoldChange, 3), "<br>padj: ", signif(padj, 4)),
      color = ~significance,
      colors = volcano_colors,  # named mapping
      source = "volcano",
      marker = list(size = 6)
    ) %>%
      layout(
        title = "Interactive Volcano Plot",
        xaxis = list(title = "log2 Fold Change"),
        yaxis = list(title = "-log10 Adjusted p-value"),
        shapes = list(
          list(type = "line", x0 = -lfc_cutoff, x1 = -lfc_cutoff, y0 = 0, y1 = max(-log10(volcano_data$padj)), line = list(dash = "dash")),
          list(type = "line", x0 = lfc_cutoff, x1 = lfc_cutoff, y0 = 0, y1 = max(-log10(volcano_data$padj)), line = list(dash = "dash")),
          list(type = "line", y0 = -log10(padj_cutoff), y1 = -log10(padj_cutoff), x0 = min(volcano_data$log2FoldChange), x1 = max(volcano_data$log2FoldChange), line = list(dash = "dash"))
        )
      )
  })
  
  
  output$geneBoxplot <- renderPlot({
    click <- event_data("plotly_click", source = "volcano")
    req(click)
    
    res <- analysisResults()$res
    rlog_mat <- assay(rlog(analysisResults()$dds))
    sample_meta <- as.data.frame(colData(analysisResults()$dds))
    
    clicked_gene <- rownames(res)[which(res$log2FoldChange == click$x & -log10(res$padj) == click$y)]
    req(length(clicked_gene) == 1)
    
    df <- data.frame(
      Expression = rlog_mat[clicked_gene, ],
      Sample = colnames(rlog_mat),
      Group = sample_meta$Group
    )
    
    ggplot(df, aes(x = Group, y = Expression)) +
      geom_boxplot(aes(fill = Group), outlier.shape = NA, alpha = 0.5) +
      geom_jitter(aes(color = Group), width = 0.2, size = 2) +
      theme_minimal() +
      labs(
        title = paste("Expression of", clicked_gene),
        y = "rlog-normalized expression",
        x = "Group"
      )
  })
  
  outputOptions(output, "pcaPlot", suspendWhenHidden = FALSE)
  outputOptions(output, "volcanoPlot", suspendWhenHidden = FALSE)
  outputOptions(output, "deTable", suspendWhenHidden = FALSE)
  outputOptions(output, "geneBoxplot", suspendWhenHidden = FALSE)
  outputOptions(output, "pathwayPlot", suspendWhenHidden = FALSE)
  
}


# Run the application 
shinyApp(ui = ui, server = server)
