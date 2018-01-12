install.packages("shiny")
library(shiny)
install.packages("ggplot2")
library(ggplot2)
install.packages("stringr")
library(stringr)
install.packages("rmarkdown")
library(rmarkdown)
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("GenomicRanges")
library(GenomicRanges)
biocLite("AnnotationHub")
library(AnnotationHub)

############################################################################################################################
#READ IN THE ENSEMBLE GENE ANNOTATIONS THROUGH ANNOTATIONHUB, GRCh37, release 75, SNAPSHOT DATE 2017-10-27
############################################################################################################################
hub = AnnotationHub()
Ensembl = hub[["AH10684"]]
GRexons = Ensembl[Ensembl$type == "CDS", ]  #subset to include coding DNA sequences (CDS)
seqlevels(GRexons) = paste("chr", seqlevels(GRexons), sep ="")  #change format of chromosome name to UCSC format

############################################################################################################################
#READ IN THE COSMIC DATASET, COSMICMUTANTEXPORT.TSV, VERSION 83, GRCh37
#MANUALLY DOWNLOADED THROUGH SFTP SITE (SEE http://cancer.sanger.ac.uk/cosmic/download FOR INSTRUCTIONS)
############################################################################################################################
cosmic = read.delim("cosmic.tsv", header = TRUE)
cosmic = cosmic[!is.na(cosmic$GRCh), ]  #exclude variants with no genomic coordinates to allow interval operations

#separate the genomic coordiante from "chr:start-end" format into "chrZ", "start", and "end" columns (UCSC format)
coordinate = str_split(cosmic$Mutation.genome.position, ":")
cosmic$chr = sapply(coordinate, function(x){x[1]})
cosmic$chr = paste("chr", cosmic$chr, sep = "")
position = sapply(coordinate, function(x){x[2]})
splitPosition = str_split(position, "-")
cosmic$start = sapply(splitPosition, function(x){x[1]})
cosmic$end = sapply(splitPosition, function(x){x[2]})

#change chromosome names to UCSC format
cosmic[cosmic$chr == "chr23", "chr"] = "chrX"
cosmic[cosmic$chr == "chr24", "chr"] = "chrY"
cosmic[cosmic$chr == "chr25", "chr"] = "chrMT"

GRcosmic = makeGRangesFromDataFrame(cosmic, keep.extra.columns = TRUE)  #convert to GRanges object
GRcosmic = sort(GRcosmic)

############################################################################################################################
#SHINY USER INTERFACE
############################################################################################################################
ui = fluidPage(
    titlePanel(
      title = em("NGS Panel Analyzer v1.0", style = "color:blue", align = "center"), 
      windowTitle = "NGS Panel Analyzer"
    ),
    
    
    tabsetPanel(
      tabPanel(title = h4("Analyzer"),
               br(),
               fileInput(
                 inputId = "bed",
                 label = "Choose BED file (csv format)", 
                 accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"),
                 multiple = FALSE), 
               
               checkboxInput(                     
                 inputId = "header", 
                 label = "My BED file has header", 
                 value = TRUE),
               
               radioButtons('format', 'Document format', c('PDF', 'HTML', 'Word'),
                            inline = TRUE),
               
               downloadButton('downloadReport'),
               
               textOutput("summary"),
               plotOutput(outputId = "insertLength"),
               plotOutput(outputId = "amplicons"),
               plotOutput(outputId = "coding"),
               plotOutput(outputId = "IDs"),
               plotOutput(outputId = "counts")
      ),
               
      tabPanel(title = h4("About"),
               br(),
               h4(p("PURPOSE:  The NGS Panel Analyzer is a Shiny App that provides summary statistics for targeted amplicon-based next generation sequencing (NGS) panels intendend to detect somatic mutations in cancer:")),
               h5(p("1.  Number of genes and amplicons in the panel")),
               h5(p("2.  Insert length for each amplicon")),
               h5(p("3.  Number of amplicons for each gene")),
               h5(p("4.  Percent of coding sequence (CDS) covered for each gene")),
               h5(p("5.  Percent of unique COSMIC IDs covered for each gene")),
               h5(p("6.  Percent of COSMIC counts covered for each gene")),
               br(),
               h4(p("REQUIREMENTS:  The Analyzer requires a BED file in CSV format with the following columns in order:")),
               h5(p("1.  chr (chromosome name)")),
               h5(p("2.  start (genomic start coordinate)")),
               h5(p("3.  end (genomic end coordinate)")),
               h5(p("4.  gene (HUGO gene symbol)")),
               h5(p("NOTE1:  Chromosome names must be in UCSC format")),
               h5(p("NOTE2:  Genomic coordinates must be relative to the hg19 build")),
               br(),
               h4(p("LIMITATIONS:")),
               h5(p("1.  Analyses are limited to amplicon-based sequencing panels.")),
               h5(p("2.  Only COSMIC variants with genomic coordinates are used.")),
               h5(p("3.  Only COSMIC variants entirely within coding regions are considered (i.e., promoter and intronic variants are excluded)")),
               h5(p("4.  Analyses are limited to COSMIC-preferred transcripts.")),
               h5(p("5.  If the COSMIC-preferred transcript is a RefSeq transcript (rare), that gene will not be analyzed.")),
               br(),
               h4(p("WORK IN PROGRESS:")),
               h5(p("1.  Allow similar analyses of capture-based sequencing panels.")),
               h5(p("2.  Preclude the need to include a gene column in the BED file.")),
               h5(p("3.  Allow subsetting of COSMIC variants based on tumor site/histology, somatic/germline status, and FATHMM prediction")),
               h5(p("4.  Allow BED files with genomic coordinates relative to other builds of the genome.")),
               h5(p("5.  Allow COSMIC-preferred RefSeq transcripts.")),
               h5(p("6.  Allow user-preferred transcripts.")),
               h5(p("7.  Preclude need for user to convert chrX, chrY, and chrM to chr23, chr24, and chr25, respectively.")),
               h5(p("8.  Compare multiple NGS panels."))
      )
    )
)

############################################################################################################################
#SHINY SERVER
############################################################################################################################
server = function(input, output) {
    bed = reactive({
    bed = input$bed
    if (is.null(bed))
      return(NULL)
    bed = read.csv(bed$datapath, header = input$header)
    bed$insert_length = bed$end - bed$start
    bed
  })
  
  dfamplicons = reactive({
    bed = bed()
    dfamplicons = as.data.frame(table(bed$gene))
    colnames(dfamplicons) = c("gene", "Number of amplicons")
    dfamplicons = dfamplicons[order(dfamplicons$`Number of amplicons`, decreasing = TRUE), ]
    dfamplicons$gene = factor(dfamplicons$gene, levels = dfamplicons$gene)
    dfamplicons
  })
  

  dfcoding = reactive({
    bed = bed()
    GRbed = makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE)
    flatGRbed = reduce(GRbed)  #merge bed file
    
    #Restore metadata (i.e., gene symbols) to flatGbed
    ov = findOverlaps(flatGRbed, GRbed)
    qhits = queryHits(ov)
    DupflatGRbed = flatGRbed[qhits]
    DupflatGRbed$gene = GRbed$gene
    flatGRbed = unique(DupflatGRbed)
    
    geneList = unique(bed$gene)  #obtain gene list from BED file
    GRcosmic = GRcosmic[GRcosmic$Gene.name %in% geneList, ]  #subset COSMIC to include genes in the gene list
    transcriptList = unique(GRcosmic$Accession.Number)  #extract COSMIC-preferred transcripts
    GRexons = GRexons[GRexons$transcript_id %in% transcriptList, ]  #subset to include transcripts within the transcript list
    
    GRintersect = intersect(GRexons, flatGRbed, ignore.strand = TRUE)  #Intersect merged bed file with exon intervals
    
    #Restore metadata (i.e., gene symbols) to GRintersect
    ov = findOverlaps(GRintersect, flatGRbed)
    subhits = subjectHits(ov)
    GRintersect$gene = flatGRbed[subhits]$gene
    
    #calculate basepairs of coding sequence (CDS) per gene
    CDS = tapply(width(GRexons), GRexons$gene_name, sum)
    df1 = data.frame(gene = names(CDS), CDS = CDS)
    
    #calculate basepairs of coding sequence covered by panel per gene
    coverage = tapply(width(GRintersect), GRintersect$gene, sum)
    df2 = data.frame(gene = names(coverage), coverage = coverage)
    
    #create data frame with coding sequence (CDS) and covered sequence
    dfcoding = merge(x = df1, y = df2, by = "gene")
    
    #calculate/plot percent of coding sequence covered by panel per gene and plot
    dfcoding$percent_coding_coverage = dfcoding$coverage/dfcoding$CDS * 100
    dfcoding = dfcoding[order(dfcoding$percent_coding_coverage, decreasing = TRUE), ]
    dfcoding$gene = factor(dfcoding$gene, levels = dfcoding$gene)
    dfcoding
  })
  
  dfIDs = reactive({
    bed = bed()
    GRbed = makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE)
    flatGRbed = reduce(GRbed)  #merge bed file
    
    #Restore metadata (i.e., gene symbols) to flatGbed
    ov = findOverlaps(flatGRbed, GRbed)
    qhits = queryHits(ov)
    DupflatGRbed = flatGRbed[qhits]
    DupflatGRbed$gene = GRbed$gene
    flatGRbed = unique(DupflatGRbed)
    
    geneList = unique(bed$gene)  #obtain gene list from BED file
    GRcosmic = GRcosmic[GRcosmic$Gene.name %in% geneList, ]  #subset COSMIC to include genes in the gene list
    transcriptList = unique(GRcosmic$Accession.Number)  #extract COSMIC-preferred transcripts
    GRexons = GRexons[GRexons$transcript_id %in% transcriptList, ]  #subset to include transcripts within the transcript list
    
    GRintersect = intersect(GRexons, flatGRbed, ignore.strand = TRUE)  #Intersect merged bed file with exon intervals
    
    #Restore metadata (i.e., gene symbols) to GRintersect
    ov = findOverlaps(GRintersect, flatGRbed)
    subhits = subjectHits(ov)
    GRintersect$gene = flatGRbed[subhits]$gene
    
    ov = findOverlaps(query = GRcosmic, subject = GRexons, type = "within")
    qhits = queryHits(ov)
    x = data.frame(GRcosmic[qhits])
    x$Gene.name = factor(x$Gene.name, levels = unique(x$Gene.name))
    codingIDs = tapply(x$Mutation.ID, x$Gene.name, function(x){length(unique(x))})
    codingIDs[is.na(codingIDs)] = 0
    df1 = data.frame("gene" = names(codingIDs), codingIDs = codingIDs)
    
    #Find unique COSMIC IDs that overlap exactly within the coverage region
    ov = findOverlaps(query = GRcosmic, subject = GRintersect, type = "within")
    qhits = queryHits(ov)
    x = data.frame(GRcosmic[qhits])
    x$Gene.name = factor(x$Gene.name, levels = unique(x$Gene.name))
    detectedIDs = tapply(x$Mutation.ID, x$Gene.name, function(x){length(unique(x))})
    detectedIDs[is.na(detectedIDs)] = 0
    df2 = data.frame("gene" = names(detectedIDs), detectedIDs = detectedIDs)
    
    
    dfIDs = merge(df1, df2, by = "gene")
    
    #Calculate & plot percent of COSMIC IDs covered by the panel
    dfIDs$Percent_ID_coverage = dfIDs$detectedIDs/dfIDs$codingIDs * 100
    dfIDs = dfIDs[order(dfIDs$Percent_ID_coverage, decreasing = TRUE), ]
    dfIDs$gene = factor(dfIDs$gene, levels = dfIDs$gene)
    dfIDs
  })
  
  dfcounts = reactive({
    bed = bed()
    GRbed = makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE)
    flatGRbed = reduce(GRbed)  #merge bed file
    
    #Restore metadata (i.e., gene symbols) to flatGbed
    ov = findOverlaps(flatGRbed, GRbed)
    qhits = queryHits(ov)
    DupflatGRbed = flatGRbed[qhits]
    DupflatGRbed$gene = GRbed$gene
    flatGRbed = unique(DupflatGRbed)
    
    geneList = unique(bed$gene)  #obtain gene list from BED file
    GRcosmic = GRcosmic[GRcosmic$Gene.name %in% geneList, ]  #subset COSMIC to include genes in the gene list
    transcriptList = unique(GRcosmic$Accession.Number)  #extract COSMIC-preferred transcripts
    GRexons = GRexons[GRexons$transcript_id %in% transcriptList, ]  #subset to include transcripts within the transcript list
    
    GRintersect = intersect(GRexons, flatGRbed, ignore.strand = TRUE)  #Intersect merged bed file with exon intervals
    
    #Restore metadata (i.e., gene symbols) to GRintersect
    ov = findOverlaps(GRintersect, flatGRbed)
    subhits = subjectHits(ov)
    GRintersect$gene = flatGRbed[subhits]$gene
    
    #Find counts that overlap exactly with exons
    ov = findOverlaps(query = GRcosmic, subject = GRexons, type = "within")
    qhits = queryHits(ov)
    codingCounts = as.data.frame(table(GRcosmic[qhits]$Gene.name))
    y = subset.data.frame(x = codingCounts, subset = codingCounts$Var1 %in% geneList)
    
    #Find unique COSMIC IDs that overlap exactly within the coverage region
    ov = findOverlaps(query = GRcosmic, subject = GRintersect, type = "within")
    qhits = queryHits(ov)
    detectedCounts = as.data.frame(table(GRcosmic[qhits]$Gene.name))
    x = subset.data.frame(x = detectedCounts, subset = detectedCounts$Var1 %in% geneList)
    
    Percent_count_coverage = x$Freq/y$Freq * 100
    names(Percent_count_coverage) = x$Var1
    dfcounts = data.frame("gene" = names(Percent_count_coverage), "Percent_count_coverage" = Percent_count_coverage)
    dfcounts = dfcounts[order(dfcounts$Percent_count_coverage, decreasing = TRUE), ]
    dfcounts$gene = factor(dfcounts$gene, levels = dfcounts$gene)
    dfcounts
  })

  output$summary = renderText({
    bed = input$bed
    if (is.null(bed))
      return(NULL)
    paste("your BED file includes", length(unique(bed()$gene)), "genes and", nrow(bed()), "amplicons.")
  })
  
  output$insertLength = renderPlot({
    bed = input$bed
    if (is.null(bed))
      return(NULL)
    ggplot(data = bed(), mapping = aes(y = insert_length, x = factor(0))) + 
      geom_jitter(col = "cornflowerblue") + 
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5)) + 
      ylab("Insert length (bp)") + 
      xlab(label = "") + 
      ggtitle("Length distribution of DNA inserts") + 
      theme(plot.title = element_text(hjust = 0.5))
  })
  
  output$amplicons = renderPlot({
    bed = input$bed
    if (is.null(bed))
      return(NULL)
    ggplot(data = dfamplicons(), mapping = aes(x = gene, y = `Number of amplicons`)) + 
      geom_col(fill = "cornflowerblue") + 
      theme(legend.position = "right", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), plot.title = element_text(hjust = 0.5)) + 
      ggtitle(label = "Amplicons per gene") + 
      ylab("Number of amplicons") + 
      xlab("Gene")
  })
  
  output$coding = renderPlot({
    bed = input$bed
    if (is.null(bed))
      return(NULL)
    dfcoding = dfcoding()
    ggplot(data = dfcoding, mapping = aes(x = gene, y = percent_coding_coverage)) + 
      geom_col(fill = "cornflowerblue") + 
      theme(legend.position = "right", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), plot.title = element_text(hjust = 0.5)) + 
      ggtitle(label = "Coding sequence coverage") + 
      ylab("Percent coverage (%)") + 
      xlab("Gene")
  })
  
  output$IDs = renderPlot({
    bed = input$bed
    if (is.null(bed))
      return(NULL)
    dfIDs = dfIDs()
    ggplot(data = dfIDs, mapping = aes(x = gene, y = Percent_ID_coverage)) + 
      geom_col(fill = "cornflowerblue") + 
      theme(legend.position = "right", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), plot.title = element_text(hjust = 0.5)) + 
      ggtitle(label = "COSMIC ID coverage") + 
      ylab("Percent coverage (%)") + 
      xlab("Gene")
  })
  
  output$counts = renderPlot({
    bed = input$bed
    if (is.null(bed))
      return(NULL)
    dfcounts = dfcounts()
    ggplot(data = dfcounts, mapping = aes(x = gene, y = Percent_count_coverage)) + 
      geom_col(fill = "cornflowerblue") + 
      theme(legend.position = "right", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), plot.title = element_text(hjust = 0.5)) + 
      ggtitle(label = "COSMIC counts coverage") + 
      ylab("Percent coverage (%)") + 
      xlab("Gene")
  })
    
}
#bed = read.csv(file = "example_BED_file.csv", header = TRUE)


############################################################################################################################
#SHINY APP EXECUTE
############################################################################################################################
shinyApp(ui = ui, server = server)

flatGRbed = reactive({
  bed = bed()
  GRbed = makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE)
  flatGRbed = reduce(GRbed)  #merge bed file
  
  #Restore metadata (i.e., gene symbols) to flatGbed
  ov = findOverlaps(flatGRbed, GRbed)
  qhits = queryHits(ov)
  DupflatGRbed = flatGRbed[qhits]
  DupflatGRbed$gene = GRbed$gene
  flatGRbed = unique(DupflatGRbed)
  flatGRbed
})

GRexons = reactive({
  bed = bed()
  geneList = unique(bed$gene)  #obtain gene list from BED file
  GRcosmic = GRcosmic[GRcosmic$Gene.name %in% geneList, ]  #subset COSMIC to include genes in the gene list
  transcriptList = unique(GRcosmic$Accession.Number)  #extract COSMIC-preferred transcripts
  GRexons = GRexons[GRexons$transcript_id %in% transcriptList, ]  #subset to include transcripts within the transcript list
  GRexons
})

GRintersect = reactive({
  flatGRbed = flatGRbed()
  GRexons = GRexons()
  GRintersect = intersect(GRexons, flatGRbed, ignore.strand = TRUE)  #Intersect merged bed file with exon intervals
  
  #Restore metadata (i.e., gene symbols) to GRintersect
  ov = findOverlaps(GRintersect, flatGRbed)
  subhits = subjectHits(ov)
  GRintersect$gene = flatGRbed[subhits]$gene
  GRintersect
})





