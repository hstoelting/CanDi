####loading packages ####
print("loading packages")
print(Sys.time())
library(shiny)
library(bslib)
library(plotly)
library(HelensTools)
library(ggplot2)
library(scales)
library(egg)
library(gridBase)
library(plyr)
library(dplyr)
library(tidyr)
library(ggpubr)
library(rPraat)
library(data.table)
library(magrittr)
library(stringr)
library(ggfortify)
library(clipr)
library(rstatix)
library(onewaytests)
library(dunn.test)
library(MESS)
library(dr4pl)
library(HelensTools)
library(eulerr)
library(ggrepel)
library(limma)
library(edgeR)
library(jsonlite)
#library(randomForest)
library(qs2)
library(DT)
library(shinyjs)
library(colourpicker)
library(shinycssloaders)
library(waiter)
library(clusterProfiler)
library(msigdbr)
library(renv)
source("R/HelperFunctions.R")
#renv::snapshot()
#### loading vectors and files #### 
print("loading vectors etc ")
print(Sys.time())
rm(list = ls())
# mainDir <- "C:/Users/hsto0009/OneDrive - Monash University/Experiments/HS17 - 1x10^5 C. albicans in BKS.db for omics/R_HS17"
# #mainDir <- "C:/Users/helen/OneDrive - Monash University/Experiments/HS17 - 1x10^5 C. albicans in BKS.db for omics/R_HS17"
# setwd(mainDir)

Groups <- c("Con_PBS", "T2D_PBS", "Con_Can",  "T2D_Can")

coloursT <- c("PBS" = "#AFCBFF", "Candida" = "#006DAE")
coloursP <- c("Control" = "white", "Diabetic" = "grey50")

Genotypes <- c("+/+" ,"Het", "Homo")
Sexes <- c("Female", "Male")
Labels1 <- c("+/+", "+/db", "db/db")
Phenotypes <- c("Control", "Diabetic")
Treatments <- c("PBS", "Candida")
Labels2 <- c(expression("PBS"), expression(italic("C. albicans")))
Labels <- c("Con_PBS" = expression("Control - PBS"),
            "T2D_PBS" = expression("Diabetic - PBS"),
            "Con_Can" = expression(paste("Control - ", italic("C. albicans"))),
            "T2D_Can" = expression(paste("Diabetic - ", italic("C. albicans"))))
Labs <- c(expression("Control - PBS"), expression("Diabetic - PBS"),
          expression(paste("Control - ", italic("C. albicans"))), expression(paste("Diabetic - ", italic("C. albicans"))))
Labs2 <- c("Con_PBS" = "Control/PBS",
           "T2D_PBS" = "Diabetic/PBS",
           "Con_Can" = "Control/Candida",
           "T2D_Can" = "Diabetic/Candida")
LabelsT <- c("PBS" = expression("PBS"), "Candida" = expression(italic("C. albicans")))
LabelsP <- c("Control", "Diabetic")

metadata <- read.csv(file = file.path("data", "HS17_metadata.csv"), stringsAsFactors = FALSE) %>%
  mutate(Genotype = factor(Genotype, levels = Genotypes), 
         Phenotype = factor(Phenotype, levels = Phenotypes), 
         Treatment = factor(Treatment, levels = Treatments),
         Group = factor(Group, levels = Groups))

## TRANSCRIPTOMICS DATA ##
# trans.plotdata <- read.csv(file = "data/HS17_transcriptomics_plotdata.csv", stringsAsFactors = FALSE) %>%
#   mutate(Genotype = factor(Genotype, levels = Genotypes),
#          Phenotype = factor(Phenotype, levels = Phenotypes),
#          Treatment = factor(Treatment, levels = Treatments),
#          Group = factor(Group, levels = Groups, labels = Labs2))
#
# qs_save(object = trans.plotdata, file = "data/HS17_transcriptomics_plotdata.qs2")
# trans.plotdata <- qs_read(file = "data/HS17_transcriptomics_plotdata.qs2")
# all_genes <- trans.plotdata$gene %>% unique() %>% sort()
# trans.statdata <- read.csv(file = "data/HS17_transcriptomics_statdata.csv", stringsAsFactors = FALSE) %>%
#   mutate(Group1 = factor(Group1, levels = Groups, labels = Labs2),
#          Group2 = factor(Group2, levels = Groups, labels = Labs2))
#
# qs_save(object = trans.statdata, file = "data/HS17_transcriptomics_statdata.qs2")
# trans.statdata <- qs_read(file = "data/HS17_transcriptomics_statdata.qs2")
###
# PROTEOMICS DATA ##
# prot.plotdata <- read.csv(file = "data/HS17_proteomics_plotdata.csv", stringsAsFactors = FALSE) %>%
#     mutate(Genotype = factor(Genotype, levels = Genotypes),
#            Phenotype = factor(Phenotype, levels = Phenotypes),
#            Treatment = factor(Treatment, levels = Treatments),
#            Group = factor(Group, levels = Groups, labels = Labs2))
# 
# qs_save(object = prot.plotdata, file = "data/HS17_proteomics_plotdata.qs2")
# prot.plotdata <- qs_read(file = "data/HS17_proteomics_plotdata.qs2") %>%
#   subset(PG.Organisms == "Mus musculus")
# all_prot_genes <- prot.plotdata$gene %>% unique()
# all_proteins <- prot.plotdata$ProteinID %>% unique() %>% sort()
# prot.statdata <- read.csv(file = "data/HS17_proteomics_statdata.csv", stringsAsFactors = FALSE) %>%
#   mutate(Group1 = factor(Group1, levels = Groups, labels = Labs2),
#          Group2 = factor(Group2, levels = Groups, labels = Labs2))
# qs_save(object = prot.statdata, file = "data/HS17_proteomics_statdata.qs2")
# 
# prot.statdata <- qs_read(file = file.path(mainDir, "data/HS17_proteomics_statdata.qs2")) %>%
#   subset(PG.Organisms == "Mus musculus")


# prot.statdata <- qs_read(file = "data/HS17_proteomics_statdata.qs2") %>%
#   subset(PG.Organisms == "Mus musculus")
# 
# prot.lookup <- setNames(prot.statdata$ProteinID %>% unique(), prot.statdata$gene %>% unique())




# metstd <- read.csv("data/HS17_Metabolomics_Metabolite-standards.csv", stringsAsFactors = FALSE)
# names <- read.csv("data/HS17_Metabolomics_Metabolite-names.csv", stringsAsFactors = FALSE) %>%
#   merge.data.frame(metstd, all.x = TRUE) %>%
#   mutate(Standard = if_else(is.na(Standard), FALSE, TRUE)) %>%
#   mutate(Label_Conf = paste0(Label, " (conf: ", Confidence, ")"))
# 
# write.csv(x = names, file = "data/HS17_Metabolomics_Metabolite-names.csv", row.names = FALSE)
names <- read.csv(file = file.path("data", "HS17_Metabolomics_Metabolite-names.csv"), stringsAsFactors = FALSE)
met.lookup <- setNames(names$Label_Conf, names$Metabolite)

# met.data <- read.csv(file = "C:/Users/hsto0009/OneDrive - Monash University/Experiments/HS17 - 1x10^5 C. albicans in BKS.db for omics/PAG_Metabolomics and proteomics/Metabolomics files/00_From Facility/P25_0932_Exp2_all-features.csv") %>%
#   pivot_longer(cols = 3:ncol(.), names_to = "Metabolite", values_to = "Intensity") %>%
#   mutate(Mouse = substr(sample, 11, 12) %>% as.numeric()) %>%
#   dplyr:: select(-label) %>%
#   merge.data.frame(names) %>%
#   merge.data.frame(y = subset(., Intensity > 0) %>% #creating a data frame that finds the
#                      group_by(Metabolite) %>% #lowest value above zero for each metabolite
#                      summarise_at(vars(Intensity), list(min = min)) %>% # and then creates a Limit of Detection (LoD)
#                      mutate(LoD = 0.2*min) %>% #column that's 1/5 of the lowest positive value
#                      dplyr::select(-min), all.x = TRUE) %>%
#   mutate(Intensity = if_else(Intensity > 0, Intensity, LoD)) %>% #replacing zeros with the LoD
#   dplyr::select(-LoD)  #dropping the LoD column because it's no longer needed
# 
# qs_save(object = met.data, file = "data/HS17_metabolomics_data.qs2")
# met.data <- qs_read(file = "data/HS17_metabolomics_data.qs2")


#BiocManager::install("org.Mm.eg.db", site_repository = "https://bioconductor.org/packages/3.21")
# library(org.Mm.eg.db )
# library(AnnotationDbi)
# 
# gene_aliases <- select(org.Mm.eg.db,
#                        keys = union(all_genes, all_prot_genes),
#                        keytype = "SYMBOL",
#                        columns = c("ALIAS")) %>%
#   group_by(SYMBOL) %>%
#   summarise(aliases = paste(ALIAS, collapse = ", ")) %>%
#   dplyr::rename(value = SYMBOL) %>%
#   mutate(label = if_else(aliases == "NA", paste(value), paste(value, " (", aliases,")", sep = ""))) %>%
#   dplyr::select(-aliases)
# 
# qs_save(object = gene_aliases, file = "data/HS17_omics_gene-aliases.qs2")
gene_aliases <- qs_read(file = file.path("data", "HS17_omics_gene-aliases.qs2"))
# gene_aliases_trans <- subset(gene_aliases, value %in% all_genes)
# qs_save(object = gene_aliases_trans, file = "data/HS17_omics_gene-aliases_trans.qs2")
gene_aliases_trans <- qs_read(file = file.path("data", "HS17_omics_gene-aliases_trans.qs2"))
# gene_aliases_prot <- subset(gene_aliases, value %in% all_prot_genes)
# qs_save(object = gene_aliases_prot, file = "data/HS17_omics_gene-aliases_prot.qs2")
gene_aliases_prot <- qs_read(file = file.path("data", "HS17_omics_gene-aliases_prot.qs2"))
#downloading pathway lists 
# library(msigdbr)
# msig.MH <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "MH")
# pathwaylist.MH <- split(msig.MH$gene_symbol, msig.MH$gs_name)
# pathwaylist.MH[["HALLMARK_ADIPOGENESIS"]]
# msig.REACTOME <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "M2", subcollection = "CP:REACTOME")
# pathwaylist.REACTOME <- split(msig.REACTOME$gene_symbol, msig.REACTOME$gs_name)
# msig.GOBP <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "M5", subcollection = "GO:BP")
# pathwaylist.GOBP <- split(msig.GOBP$gene_symbol, msig.GOBP$gs_name)
# 
# pathwaylist.all <- c(pathwaylist.MH, pathwaylist.REACTOME, pathwaylist.GOBP)
# qs_save(object = pathwaylist.all, file = "data/HS17_omics_pathwaylist.all.qs2")
# pathwaylist.all <- qs_read(file = "data/HS17_omics_pathwaylist.all.qs2")

# pathlist <- data.frame(value = names(pathwaylist.all)) %>%
#   mutate(Collection = if_else(str_detect(value, "HALLMARK_"), "HALLMARK",  
#                               if_else(str_detect(value, "REACTOME_"), "REACTOME", 
#                                       if_else(str_detect(value, "GOBP_"), "GO:BP", NA))),
#          Pathway = value %>% str_remove("HALLMARK_") %>% str_remove("REACTOME_") %>% str_remove("GOBP_") %>%
#            str_replace_all("_", " ") %>% str_to_title(), 
#          Label = paste0(Collection, ": ", Pathway))
# qs_save(object = pathlist, file = "data/HS17_omics_pathlist.qs2")
pathlist <- qs_read(file = file.path("data", "HS17_omics_pathlist.qs2"))

# pathway.options <-  msigdbr_collections(db_species = "MM") %>%
#   mutate(db_species = "Mouse", dbs = "MM")%>%
#   rbind(msigdbr_collections(db_species = "HS") %>% 
#           mutate(db_species = "Human", dbs = "HS")) %>% 
#   mutate(Label = paste0(gs_collection_name, " (", db_species, ", ", num_genesets, " gene sets)"), 
#          id = 1:nrow(.))
# qs_save(object = pathway.options, file = "data/HS17_omics_pathway.options.qs2")
pathway.options <- qs_read(file = file.path("data", "HS17_omics_pathway.options.qs2"))
# 
# all_features <- data.frame(data = "trans", Category = "Transcriptomics", 
#                            Feature = all_genes, label = all_genes) %>%
#   rbind(data.frame(data = "prot", Category = "Proteomics", Feature = all_proteins, label = all_proteins)) %>%
#   rbind(data.frame(data = "met", Category = "Metabolomics", Feature = names$Metabolite, label = names$Label)) %>% 
#   mutate(value = paste(data, Feature, sep = "_"), 
#          name = paste(Category, label, sep = ": "))
# 
# misc.data <- read.csv(file = "data/HS17_data_Shiny.csv", stringsAsFactors = FALSE) %>%
#   select(c(Mouse, KIF:Weightloss_abs)) %>%
#   pivot_longer(KIF:Weightloss_abs) %>%
#   pivot_wider(names_prefix = "Misc_")%>%
#   mutate(across(2:14, log2))
# qs_save(object = misc.data, file = "data/HS17_omics_misc.data.qs2")
# misc.data.names <- read.csv(file = "data/HS17_misc data_names.csv", stringsAsFactors = FALSE)
# qs_save(object = misc.data.names, file = "data/HS17_omics_misc.data.names.qs2")
misc.data.names <- read.csv(file = file.path("data", "HS17_misc data_names.csv"), stringsAsFactors = FALSE)
misc.lookup <- setNames(misc.data.names$Name, misc.data.names$Measurement)

# UI ------
print("starting UI")
print(Sys.time())
ui <- page_navbar(

  header = tagList(
    use_waiter(), 
    waiter_show_on_load(
      spin_fading_circles(), 
      html = h3("Loading omics dashboard...")
    )
  ),


  theme = bs_theme(
    version = 5, 
    primary = "#006DAE",
    secondary = "#862164"
  ), 
  #CSS styles ----
  tags$head(
    tags$style(
      HTML("body {font-size: 12px; }
    h5 { font-size: 12px; font-weight: bold; }
    h4 { font-size: 14px; font-weight: bold; margin-top: 0px; margin-bottom: 0px; }
    h3 { font-size: 16px; font-weight: bold; }
    h2 { font-size: 18px;  font-weight: bold; }
    h1 { font-size: 20px;  font-weight: bold; }
    t14 { font-size: 14px; }
    p.msg {font-size: 14px; line-height: 2; margin-bottom: 14px; font-weight: 500;   }
    .shiny-input-container input {font-size: 12px; }
    .shiny-text-output {font-size: 12px; text-align: center; font-weight: bold; }
    .selectize-input {font-size: 12px; }
    .selectize-dropdown {font-size: 12px; }
    .control-label {font-size: 13px; font-weight: 600; }
    .form-check-label{font-size: 13px; font-weight: 600; }
    .card-header{text-align: center; display: flex; justify-content: center; align-items: center; font-weight: bold; font-size: 14px; }
    .card-body{text-align: center; display: flex; flex-direction: column; justify-content: center}
    .navbar .nav-link{color: black; font-weight: 400; }
    .navbar .nav-link:hover {color: #862164; }
    .navbar .nav-link.active {color: #862164; font-weight: 700 !important; }
    .navcard-custom .card {background-color: white ; border-radius: 0px; padding: 0px ; border: 1px solid #F6F6F6; box-shadow: 2px 5px 2px #F6F6F6;  }
    .navcard-custom .nav-link{color: black; font-weight: 400; }
    .navcard-custom .nav-link:hover {color: #006dae; }
    .navcard-custom .nav-link.active {color: #006dae; font-weight: 700 !important; }
    .card-custom {background-color: #F6F6F6; border-radius: 0px; padding: 5px ; border: 1px solid #bebebe; box-shadow: 2px 5px 2px #bebebe;  }
    .card-custom .card-body {color: black; font-size: 14px; text-align: center; display: flex; flex-direction: column; justify-content: center}
    .card-custom .card-header {text-align: center; display: flex; justify-content: center; align-items: center; }
    .card-content {background-color: white ; border-radius: 0px; padding: 0px ; border: 1px solid #F6F6F6; box-shadow: 2px 5px 2px #F6F6F6; 
                  min-height: 400px; height: auto;  }
    .card-content .card-body {color: black; font-size: 14px; text-align: center; padding: 5px;  }
    .card-content .card-header {text-align: center; }
    .card-content-sm {background-color: white ; border-radius: 0px; padding: 0px ; border: 1px solid #F6F6F6; box-shadow: 2px 5px 2px #F6F6F6; 
                  min-height: 200px; height: auto;  }
    .card-content-sm .card-body {color: black; font-size: 14px; text-align: center; padding: 5px;  }
    .card-content-sm .card-header {text-align: center; }
    .button-grey {background-color: #006dae; font-size: 14px; text-align: center; color: white; border: none; 
                  border-radius: 1px; font-weight: 600; height: 60px; }
    .button-grey-sm {background-color: #006dae; font-size: 12px; text-align: center; color: white; border: none; 
                  border-radius: 1px; font-weight: 600; height: 30px; margin-bottom: 16px;  }
    .dt-button {background-color: #006dae; font-size: 12px; text-align: center; color: white; border: none; 
                  border-radius: 1px; font-weight: 600; height: 30px; margin-bottom: 16px;  }
    .dataTables_wrapper .dt-buttons {display: inline-block; margin-right: 15px; }
    .dataTables_filter {display: inline-block; float: right; }
    .flip-card {background-color: #F6F6F6; perspective: 1000px; cursor: pointer; width: 100%; max-width: 700px; min-height: 350px; 
                border-radius: 0px; padding: 5px ; border: 1px solid #bebebe; box-shadow: 2px 5px 2px #bebebe;  }
    .flip-card-inner {position: relative; width: 100%; height: 100%; transition: transform 0.6s; transform-style: preserve-3d; }
    .flip-card.flipped .flip-card-inner {transform: rotateY(180deg); }
    .flip-card-front, .flip-card-back {position: absolute; width:  100%; height: 100%; backface-visibility: hidden; 
                  background-color: white ; border-radius: 0px; padding: 15px ; border: 1px solid #F6F6F6; box-shadow: 2px 5px 2px #F6F6F6; 
                  box-sizing: border-box; 
                  display: flex; align-items: center; justify-contents: center; }
    .flip-card-back {transform: rotateY(180deg); flex-direction: column; }
    .flip-img {width: 100%; height: auto; max-height: 100%; object-fit: contain; display: block; }

         ")
      ),#tags$style
    tags$script(
      HTML("
     $(document).on('click', '.flip-card', function() {
      $(this).toggleClass('flipped');
    });
         ")#HTML
    )#tags$script
  ), #$tags$head
  

  # App title ----
  title = "Omics dashboard",
  id = "main_nav",
  

  #tab 0: landing page ----
  nav_panel(
    title = "Home",
    fluidRow(
      layout_column_wrap(
        width = 1, 
        card(
          style = "padding-left: 20px; padding-right: 20px; padding-top: 10px; padding-bottom: 10px; ", 
          class = "card-custom", 
          card_header(
            style = "background-color: #006dae; color: white; font-size: 14px; font-weight: bold; border-radius: 0px;",
            "Welcome"),
          height = 120,
          p("Choose from the options below to explore the multi-omics data sets of Stölting et al.")
          )#card
      ), #column 
      layout_column_wrap(
        width = "250px",
        gap = "15px", 
        card(
          card_header(
            actionButton(
            "goabout",
            "About the omics dashboard",
            class = "button-grey",
            width = "100%", 
            icon = icon("circle-info")
            )#end of actionButton
            ),
          class = "card-custom",
          height = 150,
          card_body(
            p("Learn about the experimental set-up, find details about the data analysis and explanations on how to use the Omics dashboard")
          )#end of card_body
        ), #end of card
        card(
        card_header(
          actionButton(
            "gotransprot",
            "Transcriptomics | Proteomics",
            class = "button-grey",
            width = "100%", 
            icon = icon("dna")
          )#end of actionButton
        ),
        class = "card-custom",
        height = 150,
        card_body(
          p("Explore transcriptomics and proteomics side-by-side with volcano plots and boxplots")
        )#end of card_body
      ), #end of card
        card(
          card_header(
            actionButton(
            "gomet",
            "Metabolomics",
            class = "button-grey",
            width = "100%",
            icon = icon("atom")
          )#end of actionButton),
          ), 
          class = "card-custom", 
          card_body(
                height = 150,
          p("Explore metabolomics data as volcano plots and boxplots")
          
          )#card_body
          ),
        card(
          card_header(
            actionButton(
              "godimred", 
              "Dimensionality reduction", 
              class = "button-grey",, 
              width = "100%", 
              icon = icon("diagram-project")
            )
            ), 
          height = 150, 
          class = "card-custom", 
          card_body(
            p("Display the data with the help of principal component analysis (PCA) or multi-dimensional scaling (MDS) plots")
           
          )
        ), #card
        card(
          card_header(
            actionButton(
              "gopathway", 
              "Pathway analysis", 
              class = "button-grey",
              width = "100%", 
              icon = icon("arrows-turn-to-dots")
            )
          ), 
          height = 150, 
          class = "card-custom", 
          card_body(
            p("Perform gene set enrichment analysis (GSEA) or over-representating analysis (ORA) on transcriptomics and proteomics data")
            
          )
        ), #card
        card(
          card_header(
            actionButton(
              "gocorr", 
              "Correlations", 
              class = "button-grey",
              width = "100%", 
              icon = icon("chart-line")
            )
          ), 
          height = 150, 
          class = "card-custom", 
          card_body(
            p("Explore correlations between the different features of the three datasets")
            
          )
        ), #card
        
        card(
          card_header(
            actionButton(
              "goeuler",
              "Euler diagrams",
              class = "button-grey",
              width = "100%", 
              icon = icon("circle-half-stroke")
            )#end of actionButton
            ),
          height = 150,
          class = "card-custom", 
          card_body(
             p("Sharing is caring: Find shared differential expression and create Euler diagrams")
        
          )#card_body
         
        ), #card
      card(
        card_header(
          actionButton(
            "gofungal",
            "Fungal Omics",
            class = "button-grey",
            width = "100%", 
            icon = icon("bugs")
          )#end of actionButton
        ),
        height = 150,
        class = "card-custom", 
        card_body(
          p("Explore differences in fungal transcripts and proteins between phenotypes")
          
        )#card_body
        
      ), #end of card#end of card
        card(
          card_header(
            actionButton(
              "gomixomics2",
              "mixOmics",
              class = "button-grey",
              width = "100%", 
              icon = icon("network-wired")
            )#end of actionButton
          ),
          height = 150,
          class = "card-custom", 
          card_body(
            p("Perform multi-omics integration with the help of mixOmics")
            
          )#card_body
          
        )#end of card#end of card
      ) #end of column

    )#end of fluidRow
  ), #end of tabPanel
  #About panel----
  nav_panel(
    title = "About", 
    value = "about",
    navset_pill_list(
      widths = c(3, 9), 
      nav_panel(
        title = "Experiment", 
        value = "aboutexp", 
        layout_column_wrap(
          width = "400px", 
          min_height = "500px", 
          gap = "10px", 
          card(
            class = "flip-card", 
            
            div(
              class = "flip-card-inner", 
              #FRONT
              div(
                class = "flip-card-front", 
                img(src = "ExpSchem.svg?v=2", 
                    class = "flip-img")
              ), #div front
              #BACK
              div(
                style = "max-height: 100%; max-width: 100%; display: flex; align-items: left; justify-content: center; text-align: left; padding: 5px; ",
                class = "flip-card-back", 
                p(HTML("Male BKS.Cg-Dock7<sup>m</sup> +/+ Lepr<sup>db</sup>/J (diabetic Lepr<sup>db/db</sup> and Lepr<sup>+/db</sup> or Lepr<sup>+/+</sup> 
                   controls) were infected with 1&times;10<sup>5</sup>&nbsp;CFU&nbsp;<i>C. albicans</i> SC5314 via tail vein injection at 9-11 weeks of age. 
                   Animals were culled and tissues harvested 24 hours after infection."), class = "msg")
              )#div back 
            )#div card inner
          ), #div flip card 
          card(
            class = "flip-card", 
            
            div(
              class = "flip-card-inner", 
              #FRONT
              div(
                class = "flip-card-front", 
                img(src = "TissueHarvest.svg", 
                    class = "flip-img")
              ), #div front
              #BACK
              div(
                style = "max-height: 100%; max-width: 100%; display: flex; align-items: left; justify-content: center; text-align: left; padding: 5px; 
                          font-size: clamp(14px, 1.5vw, 14px); ",
                class = "flip-card-back", 
                p(HTML("<ol>
<li>Mice were euthanised by cervical dislocation, immediately decapitated, trunk blood collected and blood glucose levels measured.&nbsp;</li>
<li>Right kidneys&nbsp;were diseccted, capsules removed and the tissues were washed in ice-cold saline and blotted dry on Whatman paper. 
Kidneys were added to centrifuge tube filters with a 0.45 &micro;m-pore cellulose acetate membrane and spun for 10 minutes at 8000&times;
<i>g</i> and 4&deg;C. Recovered kidney interstitial fluids were transferred to safe-lock centrifuge tubes and snap-frozen on dry ice.&nbsp;</li>
<li>Simultaneously, left kidneys were removed, added to cryovials and immediately snap-frozen in liquid nitrogen.&nbsp;</li>
</ol>"), class = "msg")
              )#div back 
            )#div card inner
          )#div flip card 
        )
        
      ), #nav_panel
      nav_panel(
        title = "Transcriptomics", 
        value = "abouttrans"
      ), #nav_panel
      nav_panel(
        title = "Proteomics", 
        value = "aboutprot"
      ), #nav_panel
      nav_panel(
        title = "Metabolomics", 
        value = "aboutmet"
      ), #nav_panel
      nav_panel(
        title = "Miscellaneous data", 
        value = "aboutmisc", 
        card(
          class = "card-custom", 
          div(
            style = "
            display:grid; 
            grid-template-columns: 150px 1fr 150px;
            gap: 15px; 
            align-items: end; 
            margin-bottom: 15px; 
            ", 
            actionButton(
              "aboutmiscprev", 
              label = "Previous", 
              icon = icon("backward"), 
              class = "button-grey-sm"
            ), 
            selectizeInput(
              inputId = "selected_misc",
              label = "Measurement", 
              choices = setNames(c("", misc.data.names$Measurement), c("Type to search...", misc.data.names$Name)),
              selected = "",
              width = "100%", 
              options = list(
                placeholder = "Type to search..."
              )
            ), #selectizeInput
            actionButton(
              "aboutmiscnext", 
              label = "Next", 
              icon = icon("forward"), 
              class = "button-grey-sm"
            )
          ), #div
          plotOutput(outputId = "miscplot", height = "350px")
        )#card
      ), #nav_panel
      nav_panel(
        title = "Helpful hints", 
        value = "abouthints"
      ) #nav_panel
    )#navset_pill_list
  ), 

  #Transcriptomics | Proteomics panel----
  nav_panel(
    title = "Transcriptomics | Proteomics",
    value = "transprot",
    layout_sidebar(
      sidebar = sidebar(
        title = "Settings", 
        width = 250,
        position = "left", 
        h4("Volcano plot settings"),
        radioButtons(
          inputId = "volccomp",
          label = "Choose comparison",
          choices = list(
            "PBS: Diabetic vs. Control" = "PBS",
            "Candida: Diabetic vs. Control" = "Can",
            "Control: Candida vs. PBS" = "Con",
            "Diabetic: Candida vs. PBS" = "T2D"),
          selected = "Con"),
        numericInput(
          inputId = "volcfc",
          label = "Minimum fold change",
          value = 2
        ),
        numericInput(
          inputId = "volcp",
          label = "FDR cut-off",
          value = 0.05
        ),
        actionButton(
          inputId = "volcreset",
          label = "Reset",
          class = "btn btn-secondary btn-secondary btn-sm"
        ),
        input_switch(
          id = "volcpath",
          label = "Colour by pathway",
          value = FALSE),

        conditionalPanel(
          condition = "input.volcpath == true",
          selectizeInput(
            inputId = "volcpathsel",
            label = "Pathway",
            choices = NULL# setNames(c("", pathlist$value), c("", pathlist$Label)), 
            #options = list(placeholder = 'Type to search...', create = FALSE)
          )
        ),

        h4("Boxplot settings"),

        selectizeInput(
          label = "Feature",
          inputId = "tpfeature",
          choices = NULL
          ),



        input_switch(
          id = "showstats",
          label = "Show stats",
          value = FALSE),

        conditionalPanel(
          condition = "input.showstats == true",
          radioButtons(
            inputId = "statrep",
            label = "Stats representation",
            choices = list("numerical (p = 2e-3)" = 1, "asterisks (***)" = 2),
            selected = 2))

      ), #end of sidebar Panel
      div(
        style = "display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 10px; ", 
        card(
            class = "card-content",
            card_body(
              div(
                style = "max-width: 600px;",
                plotlyOutput(outputId = "transvolcano", height = "340px", width = "100%")%>% withSpinner(type = 7)
                )
              )

          ),
          card(
            class = "card-content",
            card_body(
              div(
                style = "max-width: 600px;",
                plotlyOutput(outputId = "protvolcano", height = "340px") %>% withSpinner(type = 7)
                )
              )

          ),
          card(
            class = "card-content",
            card_body( uiOutput("transint_content")%>% withSpinner(type = 7))
          ), #card
          card(
            class = "card-content",
            card_body( uiOutput("protint_content")%>% withSpinner(type = 7))
          ) #card
      )#div 
    )# end of sidebarLayout
  ), #end of TabPanel
  #Metabolomics panel----
  nav_panel(
    title = "Metabolomics",
    value = "met",
    layout_sidebar(
      sidebar = sidebar(
        title = "Settings", 
        width = 250,
        position = "left", 
        h4("Confidence filtering"),
        input_switch(
          id = "metfiltstd",
          label = "Only metabolites validated by standards",
          value = TRUE
        ),
        conditionalPanel(
          condition = "input.metfiltstd == false",
          sliderInput(
            inputId = "metfiltconf",
            label = "Confidence cut-off",
            min = 2.5,
            max = 10,
            value = 8,
            step = 0.5

          )#sliderInput
        ), ##conditionalPanel
        actionButton(
          inputId = "recalcmet", 
          label = "Apply changes", 
          class = "btn btn-secondary btn-sm", 
          icon = icon("arrows-rotate")
        ), 
        h4("Volcano plot settings"),
        radioButtons(
          inputId = "metvolccomp",
          label = "Choose comparison",
          choices = list(
            "PBS: Diabetic vs. Control" = "PBS",
            "Candida: Diabetic vs. Control" = "Can",
            "Control: Candida vs. PBS" = "Con",
            "Diabetic: Candida vs. PBS" = "T2D"),
          selected = "PBS"),
        numericInput(
          inputId = "metvolcfc",
          label = "Minimum fold change",
          value = 1.5
        ),
        numericInput(
          inputId = "metvolcp",
          label = "FDR cut-off",
          value = 0.05
        ),
        actionButton(
          inputId = "metvolcreset",
          label = "Reset",
          class = "btn btn-secondary btn-sm"
        ), 
        h4("Boxplot settings"),
        
        selectizeInput(
          label = "Feature",
          inputId = "mfeature",
          choices = NULL),
        
        
        
        input_switch(
          id = "metshowstats",
          label = "Show stats",
          value = FALSE),
        
        conditionalPanel(
          condition = "input.metshowstats == true",
          radioButtons(
            inputId = "metstatrep",
            label = "Stats representation",
            choices = list("numerical (p = 2e-3)" = 1, "asterisks (***)" = 2),
            selected = 2))
      ), #sidebarPanel
        layout_column_wrap(
          width = 1,
          fixed_width = FALSE, 
          div(
            class = "navcard-custom", 
            navset_card_tab(
              nav_panel(title = "Plots",
                        layout_column_wrap(
                          width = "400px", 
                          gap = "10px", 
                          max_height = "400px", 
                          max_width = "600px", 
                          fixed_width = FALSE, 
                          heights_equal = "row", 
                          card(
                            class = "card-content", 
                            
                            card_body(
                              plotlyOutput("metvolcano", height = "340px")
                            )
                          ), 
                          card(
                            class = "card-content", 
                            card_body(
                              uiOutput("metint_content")
                            )
                          )#card 
                        )#layout_column_wrap
                        
              ), 
              nav_panel(
                title = "Data (normalised)", 
                div(
                  style = "width: 100%; overflow: visible; ", 
                  DTOutput(outputId = "mettable")
                  )#div
                ),#nav_panel 
              nav_panel(
                title = "Stats results", 
                div(
                  style = "width: 100%; overflow: visible; ", 
                  DTOutput(outputId = "metstats")
                )#div 
              )#nav_panel
            )#navsetcardtab
          )#div 
          )#layout column wrap 
          
        
    )#sidebarLayout
  ), #tabPanel
  #Dimensionality reduction panel ----
  nav_panel(
    title = "Dimensionality reduction", 
    value = "dimred", 
    layout_sidebar(
      sidebar = sidebar(
        title = "Settings", 
        width = 250,
        position = "left", 
        checkboxGroupInput(
          inputId = "dimreddata", 
          label = "Data set(s)", 
          choices = c("Transcriptomics" = "trans", "Proteomics" = "prot", "Metabolomics" = "met")
        ),
        div(
          style = "margin-bottom: 20px;",
          conditionalPanel(
          condition = "input.dimreddata.includes('met')", 
              actionButton(
                inputId = "gomet2", 
                label = "Adjust filters for metabolomics data",
                class = "btn-outline-secondary btn-sm"
                ) #actionButton
              )#conditional Panel
        ), #div
        #       radioButtons(
        #   inputId = "dimredweigh", 
        #   label = "Data set weight", 
        #   choices = c("None" = 0, 
        #               "Weight by number of features" = 1, 
        #               "Equalise variance per dataset" = 2, 
        #               "Use correlation distance" = 3), 
        #   selected = 0
        # ), 
        radioButtons(
          inputId = "dimredtype", 
          label = "Approach", 
          choices = c("Principal Component Analysis" = "pca", 
                      "Multi-dimensional scaling (correlation)" = "mds"
                      ), 
          selected = "pca"
        ) ,#radioButtons
        checkboxGroupInput(
          inputId = "dimredflip", 
          label = "Flip axes", 
          choices = c("x-axis" = "x", 
                      "y-axis" = "y")
        )
  
      ), #sidebarPanel
        div(
          card( 
            class = "card-custom", 
            style = "margin-bottom: 20px; ",
            #height = "100px",
            tags$span(
              HTML(" 
               <h4>Use different approaches to explore different questions</h4> <br> 
               <i><b>What dominates variance overall?</b></i> → Principal Component Analysis ("), 
              shiny::actionLink("selectPCA", "PCA"), 
              HTML(")
               <br><br> 
               <i><b>Which samples behave similarly?</b></i> → Multi-dimensional scaling of correlation distances ("), 
              shiny::actionLink("selectMDS", "MDS-corr"), 
              HTML(")
               <br><br> 
               <i><b>What is shared across omics layers?</b></i> → 
               "),#HTML 
              shiny::actionLink("gomixomics", "MixOmics") 
            )
          ), #card
          card(
            class = "card-content", 
            card_body(uiOutput("dimred_content"))
          )#card 
        )
          
    )#Sidebarlayout
  ),#tabpanel
  #Pathway Analysis panel ---- 
  nav_panel(
    title = "Pathway Analysis", 
    value = "pathway", 
    layout_sidebar(
      sidebar = sidebar(
        title = "Settings", 
        width = 250,
        position = "left", 
        h3("Settings"),
        radioButtons(
          inputId = "pathwaymethod", 
          label = "Select method", 
          choices = c("Gene Set Enrichment Analysis (GSEA)" = "GSEA", "Over-Representation Analysis (ORA)" = "ORA")
        ), 
        selectInput(
          inputId = "pathwayid", 
          label = "Pathway collection", 
          choices = setNames(pathway.options$id, pathway.options$Label), 
          selected = pathway.options[pathway.options$gs_collection == "MH", "id"]
        ),
        selectInput(
          inputId = "pathwaycomp",
          label = "Comparison",
          choices = list(
            "PBS: Diabetic vs. Control" = "PBS",
            "Candida: Diabetic vs. Control" = "Can",
            "Control: Candida vs. PBS" = "Con",
            "Diabetic: Candida vs. PBS" = "T2D"),
          selected = "Con"
        ), #end of selectInput
        
        conditionalPanel(
          condition = "input.pathwaymethod == 'ORA'", 
          numericInput(
            inputId = "pathwayfc",
            label = "Minimum fold change",
            value = 2
          ),
          numericInput(
            inputId = "pathwayp",
            label = "FDR cut-off",
            value = 0.05
          ),
          selectInput(
            inputId = "pathwaydir",
            label = "Direction of change",
            choices = list("up" = "up", "down" = "down", "either" = "both"),
            selected = "both"
          )  #end of selectInput
        ), 
        actionButton(
          inputId = "pathwaycalc", 
          label = "Analyse", 
          icon = icon("calculator"), 
          class = "btn btn-secondary btn-sm"
        )
      ), #sidebarPanel
      div(
        class = "navcard-custom", 
        navset_card_tab(
          id = "pathwaytabs", 
          nav_panel(
            title = "Plot", 
            value = "pathwayplots", 
            div(
              style = "max-width: 1000px; ", 
              uiOutput(outputId = "PAplotcont") %>% withSpinner(type = 7)
              #plotOutput(outputId = "PAplot", height = "500px") %>% withSpinner(type = 7)
            )
            
          ), 
          nav_panel(
            title = "Table", 
            value = "pathwaytable", 
            DTOutput(outputId = "PAtable") %>% withSpinner(type = 7)
            
            
          ), 
        )#navset_card_tab
      )#div 
      
    )#sidebarLayout
  ), #tabPanel
  #Correlations panel ----
  nav_panel(
    title = "Correlations",
    value = "corr",
    layout_sidebar(
      sidebar = sidebar(
        title = "Settings", 
        width = 250,
        position = "left", 
        h4("Feature 1 (x-axis)"),
        selectizeInput(
          inputId = "corrxtype", 
          label = "Data set", 
          choices = c("", "Transcriptomics", "Proteomics", "Metabolomics", "Miscellaneous"), 
          options = list(
            placeholder = "Select one..."
          )), 
          conditionalPanel(
            condition = "input.corrxtype == 'Transcriptomics'",
            selectizeInput(
              label = "Feature", 
              inputId = "corrx_trans", 
              choices = NULL), 
            ), 
        conditionalPanel(
          condition = "input.corrxtype == 'Proteomics'",
          selectizeInput(
            label = "Feature", 
            inputId = "corrx_prot", 
            choices = NULL)
        
        ), 
        conditionalPanel(
          condition = "input.corrxtype == 'Metabolomics'", 
          selectizeInput(
            label = "Feature", 
            inputId = "corrx_met", 
              choices = NULL
          )
        ), #conditional panel
        conditionalPanel(
          condition = "input.corrxtype == 'Miscellaneous'", 
          selectizeInput(
            label = "Feature", 
            inputId = "corrx_misc", 
            choices = NULL
          )
        ), #conditional panel
        h4("Feature 2 (y-axis)"),
        selectizeInput(
          inputId = "corrytype", 
          label = "Data set", 
          choices = c("", "Transcriptomics", "Proteomics", "Metabolomics", "Miscellaneous"), 
          options = list(
            placeholder = "Select one..."
          )), 
        conditionalPanel(
          condition = "input.corrytype == 'Transcriptomics'",
          selectizeInput(
            label = "Feature", 
            inputId = "corry_trans", 
            choices = NULL
           ) 
        ), 
        conditionalPanel(
          condition = "input.corrytype == 'Proteomics'",
          selectizeInput(
            label = "Feature", 
            inputId = "corry_prot", 
            choices = NULL
            )
          
        ), 
        conditionalPanel(
          condition = "input.corrytype == 'Miscellaneous'", 
          selectizeInput(
            label = "Feature", 
            inputId = "corry_misc",
            choices = NULL
          )
        ), #conditional panel
        conditionalPanel(
          condition = "input.corrytype == 'Metabolomics'", 
          selectizeInput(
            label = "Feature", 
            inputId = "corry_met",
            choices = NULL
          )
        ), #conditional panel
        actionButton(
          inputId = "corrplotgo", 
          label = "Plot correlation", 
          icon = icon("circle-right"), 
          class = "btn btn-secondary btn-sm"
        )

      ), #sidebarPanel
        layout_column_wrap(
          width = 1,  
          div(
            class = "navcard-custom", 
            navset_card_tab(
              id = "corrtabs", 
              nav_panel(
                title = "Plot", 
                value = "corrplot", 
                #plotOutput(outputId = "corrplot", width = "500px")
                uiOutput(outputId = "corrplotcont") %>% withSpinner(type = 7)
                
              ), 
              nav_panel(
                title = "Corr. with Feature 1", 
                value = "xcorr", 
                uiOutput(outputId = "corr_x_title"), 
                DTOutput(outputId = "corrdatax") %>% withSpinner(type = 7)
              ), 
              nav_panel(
                title = "Corr. with Feature 2", 
                value = "ycorr", 
                uiOutput(outputId = "corr_y_title"), 
                DTOutput(outputId = "corrdatay") %>% withSpinner(type = 7)
              )#navpanel
            )#navsetcardtab 
          )#div
         
        )#layout column wrap 
    )#sidebarLayout
  ), #tabPanel

  #Euler panel ----
  nav_panel(
    title = "Euler diagrams",
    value = "eulers",
    layout_sidebar(
      sidebar = sidebar(
        title = "Settings", 
        width = 250,
        position = "left", 
        h4("First set"),
        selectizeInput(
          inputId = "eulerdata1",
          label = "Data set",
          choices = c("", "Transcriptomics" = "trans", "Proteomics" = "prot", "Metabolomics" = "met"),
          options = list(
            placeholder = "Select one..."
          )),  #end of selectInput
        selectizeInput(
          inputId = "eulercomp1",
          label = "Comparison",
          choices = c(
            "",
            "PBS: Diabetic vs. Control" = "PBS",
            "Candida: Diabetic vs. Control" = "Can",
            "Control: Candida vs. PBS" = "Con",
            "Diabetic: Candida vs. PBS" = "T2D"),
          selected = "",
          options = list(
            placeholder = "Select one..."
          )
        ), #end of selectInput
        selectizeInput(
          inputId = "eulerdir1",
          label = "Direction of change",
          choices = c("", "up" = "up", "down" = "down", "either" = "both"),
          options = list(placeholder = "Select one...")
        ),  #end of selectInput
        numericInput(
          inputId = "eulerfc1",
          label = "Minimum fold change",
          value = 2
        ),
        numericInput(
          inputId = "eulerp1",
          label = "FDR cut-off",
          value = 0.05
        ),
        h4("Second set"),
        selectizeInput(
          inputId = "eulerdata2",
          label = "Data set",
          choices = c("", "Transcriptomics" = "trans", "Proteomics" = "prot", "Metabolomics" = "met"),
          options = list(
            placeholder = "Select one..."
          )),  #end of selectInput
        selectizeInput(
          inputId = "eulercomp2",
          label = "Comparison",
          choices = c(
            "",
            "PBS: Diabetic vs. Control" = "PBS",
            "Candida: Diabetic vs. Control" = "Can",
            "Control: Candida vs. PBS" = "Con",
            "Diabetic: Candida vs. PBS" = "T2D"),
          selected = "",
          options = list(
            placeholder = "Select one..."
          )
        ), #end of selectInput
        selectizeInput(
          inputId = "eulerdir2",
          label = "Direction of change",
          choices = c("", "up" = "up", "down" = "down", "either" = "both"),
          options = list(placeholder = "Select one...")
        ),  #end of selectInput
        numericInput(
          inputId = "eulerfc2",
          label = "Minimum fold change",
          value = 2
        ),
        numericInput(
          inputId = "eulerp2",
          label = "FDR cut-off",
          value = 0.05
        )
      ),  #end of sidebarPanel
      div(
        uiOutput(outputId = "eulercontent") %>% withSpinner(type= 7)
      )
      # div(
      #   #style = "display: grid; grid-template-columns: auto auto; grid-template-rows: auto auto; gap: 10px; width: 100%; ", 
      #   div(
      #     style = "width: 100%; ", 
      #     card(
      #       class = "card-content", 
      #       card_header(h3("Comparison results", class = "text-center", style = "margin: 0;")),
      #       # h4("First list"),
      #       navset_tab(
      #         nav_panel("First list", DTOutput(outputId = "eulertable1")),
      #         nav_panel("Second list", DTOutput(outputId = "eulertable2")),
      #         nav_panel("Shared", DTOutput(outputId = "eulershared")),
      #         nav_panel("Unique to first list", DTOutput(outputId = "eulerunique1")),
      #         nav_panel("Unique to second list", DTOutput(outputId = "eulerunique2"))
      #       ) #end of navset_tab
      #     )  #end of card
      #   ), #table div 
      #     card(
      #       style = "width: fit-content; ", 
      #       max_height = "300px", 
      #       class = "card-content-sm", 
      #       card_header(h3("Euler diagram", class = "text-center", style = "margin: 0;")),
      #       div(
      #         style = "display: grid; grid-template-columns: auto auto; gap: 10px; align-items: center;",
      #         div(
      #           style = "width: 320px; height: 250px; display: flex; align-items: center; justify-content: center; ", 
      #           plotOutput(outputId = "eulerdiagram" , width = "300px", height = "230px")
      #         ),
      #         div(
      #           style = "max-width: 200px; display: flex; flex-direction: column; align-items: center;  ", 
      #           h4("Euler settings", style = "margin-bottom: 10px; " ), 
      #           colourInput(
      #             inputId = "eulercol1",
      #             label = "Colour for first list",
      #             value = "#106107AF",
      #             allowTransparent = TRUE,
      #             closeOnClick = FALSE,
      #             width = "150px"
      #           ),
      #           colourInput(
      #             inputId = "eulercol2",
      #             label = "Colour for second list",
      #             value = "#1B7DBF9B",
      #             allowTransparent = TRUE,
      #             closeOnClick = FALSE,
      #             width = "150px"
      #           ), 
      #           actionButton(
      #             inputId = "reseteulercolours", 
      #             label = "Reset colours",
      #             class = "btn-outline-secondary btn-sm"
      #           )
      #         ) #colour settings div 
      #       )#whole card div 
      #       )# end of card 
      #   )#end of largediv 
    )#end of sidebarLayout
  ), #end of tabPanel
      #Fungal Omics Panel ----
  nav_panel(
        title = "Fungal Omics",
        value = "fungalomics",
        layout_sidebar(
          sidebar = sidebar(
            title = "Settings", 
            width = 250,
            position = "left", 
            h4("Transcriptomics filtering"), 
            numericInput(
              inputId = "fungaltransreads",
              label = "Minimum gene read count",
              value = 5
            ),
            numericInput(
              inputId = "fungaltranscpm", 
              label = "Minimum gene CPM...", 
              value = 0
            ), 
            numericInput(
              inputId = "fungaltranscpmn", 
              label = "...in at least __ samples", 
              value = 0
            ), 
            actionButton(
              inputId = "fungaltranscalc", 
              label = "Apply filters", 
              icon = icon("circle-right"), 
              class = "btn btn-secondary btn-sm"
            ), 
            numericInput(
              inputId = "fungalvolcfc",
              label = "Minimum fold change",
              value = 2
            ),
            numericInput(
              inputId = "fungalvolcp",
              label = "FDR cut-off",
              value = 0.05
            ),
            selectizeInput(
              label = "Feature",
              inputId = "fungalfeature",
              choices = NULL)#,
            # h4("Proteomics settings"), 
            # radioButtons(
            #   inputId = "fungalprotnorm", 
            #   label = "Normalisation", 
            #   choices = c("None" = "none", 
            #               "Total fungal protein" = "total", 
            #               "Median" = "median", 
            #               "Fungal 18S/RDN25 qPCR" = "qPCR"), 
            #   selected = "total"
            # ), 
            # actionButton(
            #   inputId = "fungalprotnormgo", 
            #   label = "Apply normalisation", 
            #   icon = icon("calculator"),
            #   class = "btn btn-secondary btn-sm"
            # )
          ),
          layout_column_wrap(
            width = 1,  
            div(
              class = "navcard-custom", 
              navset_card_tab(
                id = "fungaltabs", 
                nav_panel(
                  title = "Volcanoes", 
                  layout_column_wrap(
                    width = "400px", 
                    gap = "10px", 
                    max_height = "400px", 
                    max_width = "600px", 
                    fixed_width = FALSE, 
                    heights_equal = "row", 
                    card(
                      class = "card-content", 
                      
                      card_body(
                        plotlyOutput("fungaltransvolcano", height = "340px")
                      )
                    ), 
                    card(
                      class = "card-content", 
                      card_body(
                        uiOutput("fungalint_content")
                      )
                    )#card 
                  )#layout_column_wrap
                ), 
                nav_panel(
                  title = "Transcriptomics stats", 
                  DTOutput(outputId = "fungaltransstattable")%>% withSpinner(type = 7), 
                  div(
                    width = "200px", 
                    downloadButton(
                      outputId = "fungaltransstatdown", 
                      label = "Download Results", 
                      width = NULL, 
                      class = "btn btn-secondary btn-sm"
                    )
                  )
                
                ),  #nav_panel
                nav_panel(
                  title = "Proteomics stats", 
                  DTOutput(outputId = "fungalprotstattable") %>% withSpinner(type = 7)
                )
           
          )#navset_card_tab
            )#div
          )#layout_column_wrap
          
        )#layout sidebar
      ),#nav_panel, 
  #mixOmics panel----
  nav_panel(
    title = "mixOmics",
    value = "mixomics",
    layout_sidebar(
      sidebar = sidebar(
        title = "Settings", 
        width = 250,
        position = "left", 
        h4("Settings to come") 
      ),
      h4("Content to come")
    )
  )
  
)

#SERVER --------
# Define server logic
server <- function(input, output, session) {
  
  #lazy loading data ---- 
  print("starting app")
  Sys.time()%>% print()
  trans.plotdata <- reactiveVal(NULL)
  all_genes <- reactiveVal(NULL)
  trans.statdata <- reactiveVal(NULL)
  prot.plotdata <- reactiveVal(NULL)
  all_prot_genes <- reactiveVal(NULL)
  all_proteins <- reactiveVal(NULL)
  prot.statdata <- reactiveVal(NULL)
  prot.lookup <- reactiveVal(NULL)
  metstd <- reactiveVal(NULL)
  # names <- reactiveVal(NULL)
  #met.lookup <- reactiveVal(NULL)
  met.data <- reactiveVal(NULL)
  # gene_aliases <- reactiveVAl(NULL)
  # gene_aliases_trans <- reactiveVal(NULL)
  # gene_aliases_prot <- reactiveVal(NULL)
  pathwaylist.all <- reactiveVal(NULL)
  # pathlist <- reactiveVal(NULL)
  # pathway.options <- reactiveVal(NULL)
  all_features <- reactiveVal(NULL)
  misc.data <- reactiveVal(NULL)
  misc.data.trans <- reactiveVal(NULL)
  fungal.trans.counts <- reactiveVal(NULL)
  
  observeEvent(TRUE, {
    waiter_show(
    spin_fading_circles(),
    html = tagList(
      h3("Loading app...")
    )
  )
    print("loading data")
    Sys.time()%>% print()
    print("trans.plotdata")
    trans.plotdata(qs_read(file = file.path("data", "HS17_transcriptomics_plotdata.qs2")))
    print("all_genes")
    all_genes(trans.plotdata()$gene %>% unique() %>% sort())
    print("trans.statdata")
    trans.statdata(qs_read(file = file.path("data", "HS17_transcriptomics_statdata.qs2")))
    prot.plotdata(qs_read(file = file.path("data", "HS17_proteomics_plotdata.qs2")) %>%
                    subset(PG.Organisms == "Mus musculus"))
    all_prot_genes(prot.plotdata()$gene %>% unique())
    all_proteins(prot.plotdata()$ProteinID %>% unique() %>% sort())
    prot.statdata(qs_read(file = file.path("data", "HS17_proteomics_statdata.qs2")) %>%
                    subset(PG.Organisms == "Mus musculus"))
    prot.lookup(setNames(prot.statdata()$ProteinID %>% unique(), prot.statdata()$gene %>% unique()))
    metstd(read.csv(file = file.path("data", "HS17_Metabolomics_Metabolite-standards.csv"), stringsAsFactors = FALSE))
    # names(read.csv("data/HS17_Metabolomics_Metabolite-names.csv", stringsAsFactors = FALSE) %>%
    #   merge.data.frame(metstd(), all.x = TRUE) %>%
    #   mutate(Standard = if_else(is.na(Standard), FALSE, TRUE)) %>%
    #   mutate(Label_Conf = paste0(Label, " (conf: ", Confidence, ")")))
    #met.lookup(setNames(names$Label, names$Metabolite))
    met.data(qs_read(file = file.path("data", "HS17_metabolomics_data.qs2")))
    # gene_aliases(qs_read(file = "data/HS17_omics_gene-aliases.qs2"))
    # gene_aliases_trans(subset(gene_aliases(), value %in% all_genes()))
    # gene_aliases_prot(subset(gene_aliases(), value %in% all_prot_genes()))
    pathwaylist.all(qs_read(file = file.path("data", "HS17_omics_pathwaylist.all.qs2")))      
    # pathlist(qs_read(file = "data/HS17_omics_pathlist.qs2"))
    # pathway.options(qs_read(file = "data/HS17_omics_pathway.options.qs2"))
    all_features(data.frame(data = "trans", Category = "Transcriptomics", 
                            Feature = all_genes(), label = all_genes()) %>%
                   rbind(data.frame(data = "prot", Category = "Proteomics", Feature = all_proteins(), label = all_proteins())) %>%
                   rbind(data.frame(data = "met", Category = "Metabolomics", Feature = names$Metabolite, label = names$Label)) %>% 
                   mutate(value = paste(data, Feature, sep = "_"), 
                          name = paste(Category, label, sep = ": ")))
    misc.data(qs_read(file = file.path("data", "HS17_omics_misc.data.qs2")))
    misc.data.trans(misc.data() %>% 
                      mutate(across(c(2:14, 16:28), log2)))
    fungal.trans.counts(qs_read(file = file.path("data", "HS17_transcriptomics_fungal-counts.qs2")))
    
    print("finished loading data")
    Sys.time() %>% print()
    waiter_hide()
    
  }, once = TRUE)
  
  shiny::observeEvent(input$main_nav, {
    print(input$main_nav)
    print(class(input$main_nav))
  })
  
  transprotupdated <- reactiveVal(FALSE)

  #updating the selectize input fields ----
  observeEvent(input$main_nav, {
    print("observer:update transprot selectize input start")
    print(transprotupdated())
    print(class(transprotupdated()))
    if(input$main_nav == "transprot" & transprotupdated() == FALSE){
      updateSelectizeInput(
        session, 
        inputId = "volcpathsel", 
        choices = setNames(c("", pathlist$value), c("", pathlist$Label)), 
        options = list(placeholder = 'Type to search...', create = FALSE), 
        server = FALSE
      )
      print('updating tpfeature selectize')
      updateSelectizeInput(
        session, 
        inputId = "tpfeature", 
        choices = setNames(c("", gene_aliases$value), c("Type to search...", gene_aliases$label)),
        selected = "",
        options = list(
          placeholder = "Type to search..."
        ), 
        server = FALSE)
      print("done updating tpfeature selectize")
      transprotupdated(TRUE)
    }
    print("observer:update transprot selectize input end")
  }, ignoreInit = TRUE, once = FALSE)
  
  metupdated <- reactiveVal(FALSE)
  observeEvent(input$main_nav, {
    print("observer:update met selectize input start")
    if(input$main_nav == "met" & metupdated() == FALSE){
      updateSelectizeInput(
        session, 
        inputId = "mfeature", 
        choices = setNames(c("", names$Metabolite), c("Type to search...", names$Label_Conf)),
        selected = "",
        options = list(
          placeholder = "Type to search..."
        ), 
        server = FALSE
      )
      metupdated(TRUE)
    }
    print("observer:update transprot selectize input end")
  }, ignoreInit = TRUE, once = FALSE)
  
  corrupdated <- reactiveVal(FALSE)
  observeEvent(input$main_nav, {
    print("observer:update corr selectize input start")
    if(input$main_nav == "corr" & corrupdated() == FALSE){
      print('updating corr selectizes')
      updateSelectizeInput(
        session, 
        inputId = "corrx_trans", 
        choices = setNames(c("", gene_aliases_trans$value), c("Type to search...", gene_aliases_trans$label)),
        selected = "",
        options = list(
          placeholder = "Type to search..."
        ), 
        server = FALSE 
      )
      updateSelectizeInput(
        session, 
        inputId = "corry_trans", 
        choices = setNames(c("", gene_aliases_trans$value), c("Type to search...", gene_aliases_trans$label)),
        selected = "",
        options = list(
          placeholder = "Type to search..."
        ), 
       server = FALSE 
       )
      updateSelectizeInput(
        session, 
        inputId = "corrx_prot", 
        choices = setNames(c("", gene_aliases_prot$value), c("Type to search...", gene_aliases_prot$label)),
        selected = "",
        options = list(
          placeholder = "Type to search..."
        ), 
        server = FALSE
      )
      updateSelectizeInput(
        session, 
        inputId = "corry_prot", 
        choices = setNames(c("", gene_aliases_prot$value), c("Type to search...", gene_aliases_prot$label)),
        selected = "",
        options = list(
          placeholder = "Type to search..."
        ), 
        server = FALSE
      )
      updateSelectizeInput(
        session, 
        inputId = "corrx_met",  
        choices = setNames(c("", names$Metabolite), c("Type to search...", names$Label_Conf)),
        selected = "",
        options = list(
          placeholder = "Type to search..."
        ), 
        server = FALSE
        
      )
      updateSelectizeInput(
        session, 
        inputId = "corry_met",  
        choices = setNames(c("", names$Metabolite), c("Type to search...", names$Label_Conf)),
        selected = "",
        options = list(
          placeholder = "Type to search..."
        ), 
        server = FALSE
        
      )
      updateSelectizeInput(
        session, 
        inputId = "corrx_misc",  
        choices = setNames(c("", misc.data.names$Measurement), c("Type to search...", misc.data.names$Name)),
        selected = "",
        options = list(
          placeholder = "Type to search..."
        ), 
        server = FALSE
        
      )
      updateSelectizeInput(
        session, 
        inputId = "corry_misc",  
        choices = setNames(c("", misc.data.names$Measurement), c("Type to search...", misc.data.names$Name)),
        selected = "",
        options = list(
          placeholder = "Type to search..."
        ), 
        server = FALSE
        
      )
      
      corrupdated(TRUE)
    }
    print("observer:update corr selectize input end")
  }, ignoreInit = TRUE, once = FALSE)
  

  
  #making sure the links on the landing page work----
  observeEvent(input$goabout, {
    updateNavbarPage(session, "main_nav", selected = "about")
  })
  observeEvent(input$gotransprot, {
    updateNavbarPage(session, "main_nav", selected = "transprot")
  })
  observeEvent(input$gomet, {
    updateNavbarPage(session, "main_nav", selected = "met")
  })
  observeEvent(input$gomet2, {
    updateNavbarPage(session, "main_nav", selected = "met")
  })
  observeEvent(input$godimred, {
    updateNavbarPage(session, "main_nav", selected = "dimred")
  })
  observeEvent(input$gopathway, {
    updateNavbarPage(session, "main_nav", selected = "pathway")
  })
  observeEvent(input$gocorr, {
    updateNavbarPage(session, "main_nav", selected = "corr")
  })
  observeEvent(input$goeuler, {
    updateNavbarPage(session, "main_nav", selected = "eulers")
  })
  observeEvent(input$gofungal, {
    updateNavbarPage(session, "main_nav", selected = "fungalomics")
  })  
  observeEvent(input$gomixomics, {
    updateNavbarPage(session, "main_nav", selected = "mixomics")
  })
  
  
  observeEvent(input$selectPCA, {
    updateRadioButtons(
      session, 
      inputId = "dimredtype", 
      selected = "pca"
    )
  })
  
  
  observeEvent(input$selectMDS, {
    updateRadioButtons(
      session, 
      inputId = "dimredtype", 
      selected = "mds"
    )
  })
  
  
  observeEvent(input$metvolcreset, {
    updateNumericInput(
      session, 
      inputId = "metvolcfc", 
      value = 1.5
    )
    updateNumericInput(
      session, 
      inputId = "metvolcp", 
      value = 0.05 
    )
  })
  
  observeEvent(input$volcreset, {
    updateNumericInput(
      session, 
      inputId = "volcfc", 
      value = 2
    )
    updateNumericInput(
      session, 
      inputId = "volcp", 
      value = 0.05 
    )
  })
  
  observeEvent(input$reseteulercolours, {
    updateColourInput(
      session, 
      inputId = "eulercol1", 
      value = "#106107AF"
    )
    updateColourInput(
      session, 
      inputId = "eulercol2", 
      value = "#1B7DBF9B"
    )
  })
  

  selected_tpfeature <- reactiveVal(NULL)
  selected_mfeature <- reactiveVal(NULL)
  clickupdate_tpfeature <- reactiveVal(FALSE)
  selected_fungalfeature <- reactiveVal("")

  observeEvent(input$tpfeature, {
    print("observer:input$tpfeature start")
    if(is.null(input$tpfeature) || input$tpfeature == "") return()
    if(clickupdate_tpfeature()) return()
    if(!identical(selected_tpfeature(), input$tpfeature)){
      
      selected_tpfeature(input$tpfeature)
    }
    print("observer:input$tpfeature end")

  })

  observeEvent(input$mfeature, {
    selected_mfeature(input$mfeature)
    
  })
  
  # observeEvent(input$tpfeature, {
  #  
  # })

  #observing plotly_click events---- 
  observeEvent(plotly::event_data("plotly_click", source = "transprott"), {
    click <- plotly::event_data("plotly_click", source = "transprott")
    req(click)
    req(click$key)
    
    clickupdate_tpfeature(TRUE)
    
    selected_tpfeature(as.character(click$key))

    updateSelectizeInput(
      session,
      inputId = "tpfeature",
      selected = click$key
    )
    
    clickupdate_tpfeature(FALSE)
  })
  observeEvent(plotly::event_data("plotly_click", source = "transprotp"), {
    click <- plotly::event_data("plotly_click", source = "transprotp")
    req(click)
    req(click$key)
    clickupdate_tpfeature(TRUE)

    gene <- prot.statdata()[prot.statdata()$ProteinID == click$key, "gene"] %>% unique() %>% .[1] %>% as.character()

    if(is.null(selected_tpfeature()) | !identical(selected_tpfeature(), gene)){

      selected_tpfeature(as.character(gene))
      
      updateSelectizeInput(
        session,
        inputId = "tpfeature",
        #choices = setNames(c("", gene_aliases$value), c("Type to search...", gene_aliases$label)),
        selected = gene
      )
    }
    
    clickupdate_tpfeature(FALSE)
  })
  
  observeEvent(plotly::event_data("plotly_click", source = "metvolc"), {
    click <- plotly::event_data("plotly_click", source = "metvolc")
    req(click)
    req(click$key)
    selected_mfeature(as.character(click$key))
    
    updateSelectizeInput(
      session,
      inputId = "mfeature",
      selected = click$key
    )
  })
  
  observeEvent(plotly::event_data("plotly_click", source = "fungaltransvolc"), {
    click <- plotly::event_data("plotly_click", source = "fungaltransvolc")
    req(click)
    req(click$key)
    selected_fungalfeature(as.character(click$key))
    
    # updateSelectizeInput(
    #   session,
    #   inputId = "mfeature",
    #   selected = click$key
    # )
  })
  
  ##### ABOUT ##### 
  
  #observe actionButton aboutmiscprev ----
  observeEvent(input$aboutmiscprev, {
    current <- match(input$selected_misc, misc.data.names$Measurement)
    if(nzchar(input$selected_misc) == FALSE){
      updateSelectizeInput(
        session = session, 
        inputId = "selected_misc", 
        selected = misc.data.names$Measurement[nrow(misc.data.names)]
      )
    }
    else if(current > 1){
      updateSelectizeInput(
        session = session, 
        inputId = "selected_misc", 
        selected = misc.data.names$Measurement[current  - 1]
      )
    }
    else if(current == 1){
      updateSelectizeInput(
        session = session, 
        inputId = "selected_misc", 
        selected = misc.data.names$Measurement[nrow(misc.data.names)]
      )
    }
  })
  #observe actionButton aboutmiscnext ----
  observeEvent(input$aboutmiscnext, {
    current <- match(input$selected_misc, misc.data.names$Measurement)
    if(nzchar(input$selected_misc) == FALSE){
      updateSelectizeInput(
        session = session, 
        inputId = "selected_misc", 
        selected = misc.data.names$Measurement[1]
      )
    }

    else if(current < nrow(misc.data.names)){
      updateSelectizeInput(
        session = session, 
        inputId = "selected_misc", 
        selected = misc.data.names$Measurement[current  + 1]
      )
    }
    else if(current == nrow(misc.data.names)){
      updateSelectizeInput(
        session = session, 
        inputId = "selected_misc", 
        selected = misc.data.names$Measurement[1]
      )
    }
  })

  #output$miscplot ----
  output$miscplot <- renderPlot({
    req(input$main_nav == "about")
    req(nzchar(input$selected_misc))
    
    ftr <- paste0("Misc_", input$selected_misc)
    
    data <- misc.data() %>%
      select(c("Mouse", ftr)) %>% 
      na.omit() %>% 
      dplyr::rename(Feature = ftr) %>% 
      merge.data.frame(metadata, by = "Mouse") 
      
    
   
    name <- misc.lookup[input$selected_misc]
    
    Max <- max(data[, "Feature"])
    Min <- min(data[, "Feature"])
    
    ytitle = parse(text = misc.data.names[misc.data.names$Measurement == input$selected_misc, "Label2"])
    
    if(nrow(data) < 24){
      data <- data %>%
        mutate(Group = factor(Group, levels = c("Con_Can", "T2D_Can")))

      stats <- data %>% 
        wilcox_test(Feature ~ Group, exact = TRUE, comparisons = list(c("Con_Can", "T2D_Can")))
      
      statplot <- stats %>%
        mutate(sig = unlist(lapply(X = p, FUN = asterisk)))%>% 
        mutate(ybarfact = 0, 
               size = if_else(sig == "ns", "small", "large"), 
               textnudge = if_else(size == "large", 0, 0.02))
      
      pbs <- prettybreaks_log(min = Min, max = Max)
      plotmin <- pbs[1]
      plotmax <- pbs[2]
      plotrange <- log10(plotmax)-log10(plotmin)
      
      p <- ggplot(data = data, mapping = aes(x = Group, y = Feature))+
        scale_x_discrete(expand=expansion(mult=0, add=0.67))+
        scale_y_log10(name = ytitle, limits = c(plotmin, plotmax),
                           expand = expansion(mult = 0, add = 0))+
        geom_boxplot(mapping = aes(fill = Phenotype), colour = "black", show.legend = F, fatten = 0.75, outliers = FALSE)+
        geom_jitter(mapping = aes(fill = Treatment, shape = Sex), stroke = 1.5/1.427, size = 2, width = 0.2, height = 0, show.legend = TRUE) +
        scale_fill_manual(values = c(coloursP, coloursT), labels = c(LabelsP, LabelsT)) +
        scale_shape_manual(values = c("Female" = 21, "Male" = 22))+
        scale_size_manual(values = c("large" = 6, "small" = 4), guide = "none")+
        theme(panel.background = element_rect(fill="transparent"),
              axis.line = element_line(colour="black", linewidth = 1/1.427),
              plot.title = element_text(hjust=0.5, size=13),
              plot.subtitle = element_text(hjust=0.5, size=11),
              axis.ticks.y = element_line(linewidth = 1/1.427, colour = "transparent"),
              axis.ticks.x = element_blank(),
              axis.text.y = element_text(colour="black", size = 11),
              axis.text.x = element_blank(),
              axis.title.y = element_text(size=12),
              axis.title.x = element_blank(),
              legend.key.height = unit(0.6, "cm"),
              legend.title = element_blank(),
              legend.key = element_blank(),
              legend.text = element_text(size=11, hjust = 0),
              legend.margin = ggplot2::margin(-5.5, 5.5, -5.5, 5.5),
              plot.caption = element_text (colour = "black", size = 10, hjust = 1),
              panel.grid = element_blank())+
        labs(title = name) +
        #coord_cartesian(clip = "off")+
        coord_fixed(clip = "off", ratio = (5.76/plotrange))+
        guides(fill = guide_legend(order = 1,
                                   override.aes = list(shape  = 22,
                                                       size = c(4, 4, 2))),
               shape = guide_legend(order = 2, override.aes = list(fill = "grey")))+
        annotation_logticks(sides = "l", size = 1/1.427, outside = TRUE, short = unit(0.05, "cm"),mid = unit(0.1, "cm"),long = unit(0.15, "cm"))+
        geom_segment(inherit.aes = TRUE, data = statplot, mapping = aes(x = group1, xend = group2, y = 10^(log10(plotmin)+plotrange*(0.98-ybarfact*0.07-0.02)), yend = 10^(log10(plotmin)+plotrange*(0.98-ybarfact*0.07-0.02))))+
        geom_text(inherit.aes = TRUE, data = statplot, mapping = aes(x = 1.5, label = sig, y = 10^(log10(plotmin)+plotrange*(0.98-ybarfact*0.07+textnudge)), size = size))
      
    }
    else {
      print("else")
      stats <- data %>% 
        wilcox_test(Feature ~ Group, exact = TRUE, 
                    comparisons = list(c("Con_PBS", "T2D_PBS"), 
                                       c("Con_Can", "T2D_Can"), 
                                       c("Con_PBS", "Con_Can"), 
                                       c("T2D_PBS", "T2D_Can")))
      
      statplot <- stats %>%
        mutate(Group1 = factor(group1, levels = Groups), 
               Group2 = factor(group2, levels = Groups)) %>% 
        mutate(barx1 = as.numeric(Group1),
               barx2 = as.numeric(Group2),
               starx = (barx1 + barx2)/2,
               ybarfact = mapply(FUN = function(grp1, grp2){
                 if ((grp1 == "Con_PBS" & grp2 == "T2D_PBS") | (grp2 == "Con_PBS" & grp1 == "T2D_PBS")){return(2)}
                 else if ((grp1 == "Con_Can" & grp2 == "T2D_Can") | (grp2 == "Con_Can" & grp1 == "T2D_Can")){return(2)}
                 else if ((grp1 == "Con_PBS" & grp2 == "Con_Can") | (grp2 == "Con_PBS" & grp1 == "Con_Can")){return(1)}
                 else if ((grp1 == "T2D_PBS" & grp2 == "T2D_Can") | (grp2 == "T2D_PBS" & grp1 == "T2D_Can")){return(0)}
               }, grp1 = Group1, grp2 = Group2), 
               size = if_else(p.adj.signif == "ns", "small", "large"), 
               textnudge = if_else(size == "large", 0, 0.02))
      
      

      if(grepl(pattern = "Perc.", x = ftr) == TRUE){
        plotrange <- 20
        plotmin <- 90
        plotmax <- 110
        breaks <- c(90, 95, 100, 105, 110)
      }
      else if(Min < 0){
        
        if(Max < 0){
          pbs <- prettybreaks_neg(min = Min, max = 1)
          plotrange <- pbs[2]-pbs[1]
          plotmin <- pbs[1]
          plotmax <- pbs[2]
          breaks <- seq(pbs[1], pbs[2], plotrange/pbs[3])
        }
        else{
          pbs <- prettybreaks_neg(min = Min, max = Max)
          plotrange <- pbs[2]-pbs[1]
          plotmin <- pbs[1]
          plotmax <- pbs[2]
          breaks <- seq(pbs[1], pbs[2], plotrange/pbs[3])
        }
        
      }
      else{
        plotmin <- 0
        pbs <- prettybreaks(Max)
        plotmax <- pbs[1]
        breaks <- seq(0, plotmax, plotmax/pbs[2])
        plotrange <- plotmax - plotmin
      }
      
      
      set.seed(seed = 793)
      
      
      
      p <- ggplot(data = data, mapping = aes(x = Group, y = Feature))+
        scale_x_discrete(expand=expansion(mult=0, add=0.67))+
        scale_y_continuous(name = ytitle, limits = c(plotmin, plotmax),
                           expand = expansion(mult = 0, add = 0), breaks = breaks)+
        geom_boxplot(mapping = aes(fill = Phenotype), colour = "black", show.legend = F, fatten = 0.75, outliers = FALSE)+
        geom_jitter(mapping = aes(fill = Treatment, shape = Sex), stroke = 1.5/1.427, size = 2, width = 0.2, height = 0, show.legend = TRUE) +
        scale_fill_manual(values = c(coloursP, coloursT), labels = c(LabelsP, LabelsT)) +
        scale_shape_manual(values = c("Female" = 21, "Male" = 22))+
        scale_size_manual(values = c("large" = 6, "small" = 4), guide = "none")+
        theme(panel.background = element_rect(fill="transparent"),
              axis.line = element_line(colour="black", linewidth = 1/1.427),
              plot.title = element_text(hjust=0.5, size=13),
              plot.subtitle = element_text(hjust=0.5, size=11),
              axis.ticks.y = element_line(linewidth = 1/1.427, colour = "black"),
              axis.ticks.x = element_blank(),
              axis.text.y = element_text(colour="black", size = 11),
              axis.text.x = element_blank(),
              axis.title.y = element_text(size=12),
              axis.title.x = element_blank(),
              legend.key.height = unit(0.6, "cm"),
              legend.title = element_blank(),
              legend.key = element_blank(),
              legend.text = element_text(size=11, hjust = 0),
              legend.margin = ggplot2::margin(-5.5, 5.5, -5.5, 5.5),
              plot.caption = element_text (colour = "black", size = 10, hjust = 1),
              panel.grid = element_blank())+
        labs(title = name) +
        #coord_cartesian(clip = "off")+
        coord_fixed(clip = "off", ratio = (5.76/plotrange))+
        guides(fill = guide_legend(order = 1,
                                   override.aes = list(shape  = 22,
                                                       size = c(4, 4, 2, 2))),
               shape = guide_legend(order = 2, override.aes = list(fill = "grey")))+
        geom_segment(inherit.aes = FALSE, data = statplot, mapping = aes(x = barx1, xend = barx2, y = plotrange*(0.98-ybarfact*0.07-0.02)+plotmin, yend = plotrange*(0.98-ybarfact*0.07-0.02)+plotmin))+
        geom_text(inherit.aes = FALSE, data = statplot, mapping = aes(x = starx, label = p.adj.signif, y = plotrange*(0.98-ybarfact*0.07)+plotmin+plotrange*textnudge, size = size))
    }
    
    
    p
  }, 
  res = 100 
  )

  ##### METABOLOMICS ##### 

  #calculating metabolomics data----
  met.calcdata <- eventReactive(input$recalcmet, {
    
    updateSelectizeInput(
      session = session, 
      inputId = "mfeature", 
      choices = setNames(c("", names$Metabolite), c("Type to search...", names$Label_Conf))
    )

    metnorm <- met.data() %>% 
      dplyr::select(-sample) %>%
      dplyr::left_join(met.data() %>%
                         dplyr::select(-sample) %>%
                         group_by(Mouse) %>%
                         summarise(median = median(Intensity)) %>%
                         as.data.frame(),
                       by = "Mouse") %>%
      mutate(Intensity = Intensity / median,
             Intlog2 = log2(Intensity)) %>%
      dplyr::select(-median)
    
    
    if(input$metfiltstd == TRUE){
      plot.data <- subset(metnorm, Standard == TRUE) 
    }
    else if(input$metfiltstd == FALSE){
      req(input$metfiltconf)
      plot.data <- subset(metnorm, Confidence >= input$metfiltconf)
    }
    else {
      plot.data <- metnorm}
    

    
    data.limma <- plot.data %>% 
      dplyr::select(Mouse, Metabolite, Intlog2) %>% 
      pivot_wider(names_from = "Mouse", values_from = "Intlog2") %>% 
      as.data.frame() %>% 
      `rownames<-`(.$Metabolite) %>% 
      dplyr::select(-Metabolite) 

    meta.limma <- data.frame(Mouse = c(colnames(data.limma))) %>%
      left_join(y= metadata %>% mutate(Mouse = as.character(Mouse))) %>%
      mutate(Group = factor(Group, levels = Groups))

    design <- model.matrix(~0 + meta.limma$Group) %>%
      `colnames<-`(levels(meta.limma$Group))
    
    fit <- lmFit (data.limma, design)#Fit linear model for each gene given a series of arrays
    
    #defining contrasts for pairwise comparisons 
    contrasts <- makeContrasts(
      PBS = T2D_PBS - Con_PBS,
      Can = T2D_Can - Con_Can,
      Con = Con_Can - Con_PBS,
      T2D = T2D_Can - T2D_PBS,
      levels = design
    )
    
    fit2 <- contrasts.fit(fit, contrasts) #Given a linear model fit to microarray data, compute estimated coefficients and standard errors for a given set of contrasts.
    fit2 <- eBayes(fit2) #Given a linear model fit from lmFit, compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression by empirical Bayes moderation of the standard errors towards a global value.
    
    stat.data <- cbind(Metabolite = rownames(data.limma), 
                     topTable(fit2, coef = "PBS", number = Inf, adjust.method = "BH", sort.by = "none")) %>%
      pivot_longer(cols = 2:ncol(.)) %>%
      mutate(Comparison = "PBS") %>%
      rbind(y = cbind(Metabolite = rownames(data.limma),
                      topTable(fit2, coef = "Can", number = Inf, adjust.method = "BH", sort.by = "none")) %>%
              pivot_longer(cols = 2:ncol(.)) %>%
              mutate(Comparison = "Can") ) %>%
      rbind(y = cbind(Metabolite = rownames(data.limma),
                      topTable(fit2, coef = "Con", number = Inf, adjust.method = "BH", sort.by = "none")) %>%
              pivot_longer(cols = 2:ncol(.)) %>%
              mutate(Comparison = "Con") )%>%
      rbind(y = cbind(Metabolite = rownames(data.limma),
                      topTable(fit2, coef = "T2D", number = Inf, adjust.method = "BH", sort.by = "none")) %>%
              pivot_longer(cols = 2:ncol(.)) %>%
              mutate(Comparison = "T2D") ) %>%
      pivot_wider(id_cols = c("Metabolite", "Comparison"),
                  names_from = c("name"),
                  values_from = "value") %>% 
      merge.data.frame(names)%>%
      merge.data.frame(y = data.frame(Comparison = c("PBS", "Can", "Con", "T2D"), 
                                      Group1 = c("Con_PBS", "Con_Can", "Con_PBS", "T2D_PBS"), 
                                      Group2 = c("T2D_PBS", "T2D_Can", "Con_Can", "T2D_Can"))) %>%
      mutate(Group1 = factor(Group1, levels = Groups, labels = Labs2), 
             Group2 = factor(Group2, levels = Groups, labels = Labs2))

      list(
        met.data = metnorm, 
        plot.data = plot.data, 
        stat.data = stat.data)
  }, ignoreNULL = FALSE, ignoreInit = FALSE )
  
  #output$mettable ----
  output$mettable <- renderDataTable({
    
  
    met.calcdata()$plot.data %>%
      DT::datatable(., rownames = FALSE,
                    options = list(
                      columnDefs = list(
                      list(
                        targets = c("Intensity", "Mass", "RT", "Intlog2"),
                        render = JS(
                          "function(data, type, row) {",
                          "  return type === 'display' ?",
                          "    Number(data).toPrecision(4) : data;",
                          "}"
                        )
                      ))
                    )
      )
  })
  
  #output$metstats ----
  output$metstats <- renderDataTable({
    met.calcdata()$stat.data %>% 
      DT::datatable(., rownames = FALSE,
                    options = list(
                      columnDefs = list(list(
                        targets = c("P.Value", "adj.P.Val"),
                        render = JS(
                          "function(data, type, row) {",
                          "  return type === 'display' ?",
                          "    Number(data).toExponential(3) : data;",
                          "}"
                        )
                      ),
                      
                      list(
                        targets = c("logFC", "AveExpr", "t", "B", "Mass", "RT"),
                        render = JS(
                          "function(data, type, row) {",
                          "  return type === 'display' ?",
                          "    Number(data).toPrecision(4) : data;",
                          "}"
                        )
                      ))
                    )
      )
  })

  #output$metvolcano ----
    output$metvolcano <- renderPlotly({
      req(input$main_nav == "met")
      req(input$metvolccomp)
      
      
      metvolcdf <- subset(met.calcdata()$stat.data, Comparison == input$metvolccomp)
      
      
      p <- ggplot(data = metvolcdf, mapping = aes(x = logFC, y = -log10(adj.P.Val), 
                                                  key = Metabolite, text = paste0(
                                                    "<b>", Label, "</b><br>",
                                                    "Confidence: ", Confidence, "<br>",
                                                    "log<sub>2</sub>(FC): ", signif(logFC, 3), "<br>", 
                                                    "FDR: ", signif(adj.P.Val, 3))))+
        geom_hline(yintercept = -log10(input$metvolcp), colour= "black", linetype = "dotted")+
        geom_vline(xintercept = c(-log2(input$metvolcfc), log2(input$metvolcfc)), colour= "black", linetype = "dotted")+
        geom_point(data = metvolcdf %>% subset(adj.P.Val >= input$metvolcp |
                                                  abs(logFC) < log2(input$metvolcfc)),
                   colour = "grey", alpha = 0.5)+
        geom_point(data = metvolcdf %>% subset(adj.P.Val < input$metvolcp &
                                                  logFC > log2(input$metvolcfc)),
                   colour = "#862164", alpha = 0.5)+
        geom_point(data = metvolcdf %>% subset(adj.P.Val < input$metvolcp &
                                                  logFC < -log2(input$metvolcfc)),
                   colour = "#006DAE", alpha = 0.5)+
        scale_x_continuous(name = expression("log"["2"]~"Fold Change (FC)"), limits = c(-10, 10), expand = expansion())+
        scale_y_continuous(name = expression("-log"["10"]~"False Discovery Rate (FDR)") , limits = c(0, 20), expand = expansion())+
        theme(panel.background = element_rect(fill="transparent"),
              axis.line = element_line(colour="black", linewidth = 1/1.427),
              plot.title = element_text(hjust=0.5, size=13),
              plot.subtitle = element_text(hjust=0.5, size=11),
              axis.ticks= element_line(linewidth = 1/1.427, colour = "black"),
              axis.text = element_text(colour="black", size = 11),
              axis.title = element_text(size=12),
              legend.key.height = unit(0.6, "cm"),
              legend.title = element_blank(),
              legend.key = element_blank(),
              legend.text = element_text(size=11, hjust = 0),
              legend.margin = ggplot2::margin(-5.5, 5.5, -5.5, 5.5),
              plot.caption = element_text (colour = "black", size = 10, hjust = 1),
              panel.grid = element_blank())+
        annotate(geom = "text", x = ((-10+log2(input$metvolcfc))/2)-log2(input$metvolcfc), y = 14.5/.75, label = metvolcdf[1, "Group1"], colour = "#006DAE", hjust = 0.5)+
        annotate(geom = "text", x = ((10-log2(input$metvolcfc))/2)+log2(input$metvolcfc), y = 14.5/.75, label =  metvolcdf[1, "Group2"], colour = "#862164", hjust = 0.5)+
        annotate(geom = "text", x = ((-10+log2(input$metvolcfc))/2)-log2(input$metvolcfc), y = .5/.75, label = metvolcdf %>% subset(adj.P.Val < input$metvolcp &
                                                                                          -logFC > log2(input$metvolcfc)) %>% nrow(),
                 colour = "#006DAE", hjust = 0.5, vjust = 1)+
        annotate(geom = "text", x = ((10-log2(input$metvolcfc))/2)+log2(input$metvolcfc), y = 0.5/.75, label = metvolcdf %>% subset(adj.P.Val < input$metvolcp &
                                                                                          logFC > log2(input$metvolcfc)) %>% nrow(),
                 colour = "#862164", hjust = 0.5, vjust = 1)+
        labs(title = "Metabolomics")#+
      #coord_fixed(clip = "off", ratio = 1.5)
      
      if (!is.null(selected_mfeature())){
        p <- p+
          geom_point (data = metvolcdf %>% subset(Metabolite == selected_mfeature()),
                      colour = "black", alpha = 0.8)
      } else {p <- p}
      
      
      
      ggplotly(p, tooltip = "text", source = "metvolc") %>%
        layout(
          xaxis = list(title = "log<sub>2</sub> Fold Change (FC)", titlefont = list(size = 14), tickfont = list(size = 10)),
          yaxis = list(title = "-log<sub>10</sub> False Discovery Rate (FDR)", titlefont = list(size = 14), tickfont = list(size = 10))
          ) %>%
        plotly::event_register("plotly_click")
    }
    ) #end of output$metvolcano
  
  #output$metintplot ----
  output$metintplot <- renderPlot({
    req(input$main_nav == "met")
    req(nzchar(selected_mfeature()))
    req(selected_mfeature() %in% met.calcdata()$stat.data$Metabolite)

    
    # boxplot of log2 intensity ---
    # with requested gene
    # 1. It is "reactive" and therefore should be automatically
    #    re-executed when inputs (input$gene) change
    # 2. Its output type is a plot
    metdf <- subset(met.calcdata()$plot.data, Metabolite == selected_mfeature()) %>%
      merge.data.frame(metadata, by = "Mouse") 
      
    
    name <- metdf[1, "Label"]
   
    Max <- max(metdf[, "Intlog2"])
    Min <- min(metdf[, "Intlog2"])
    # if(Min < 0){
    #   
    #   if(Max < 0){
    #     pbs <- prettybreaks_neg(min = Min, max = 1)
    #     plotrange <- pbs[2]-pbs[1]
    #     plotmin <- pbs[1]
    #     plotmax <- pbs[2]
    #     breaks <- seq(pbs[1], pbs[2], plotrange/pbs[3])
    #   }
    #   else{
    #     pbs <- prettybreaks_neg(min = Min, max = Max)
    #     plotrange <- pbs[2]-pbs[1]
    #     plotmin <- pbs[1]
    #     plotmax <- pbs[2]
    #     breaks <- seq(pbs[1], pbs[2], plotrange/pbs[3])
    #   }
    #   
    # }
    # else{
      pbs <- prettybreaks_log2(Min, Max)
      plotrange <- pbs[2]-pbs[1]
      plotmin <- pbs[1]
      plotmax <- pbs[2]
      breaks <- seq(pbs[1], pbs[2], plotrange/pbs[3])
      #plotrange <- plotmax - plotmin
    # }
    
    set.seed(seed = 793)
    
    
    p <- ggplot(data = metdf, mapping = aes(x = Group, y = Intlog2))+
      scale_x_discrete(expand=expansion(mult=0, add=0.67))+
      scale_y_continuous(name = expression("log"["2"]~"Intensity"), limits = c(plotmin, plotmax),
                         expand = expansion(mult = 0, add = 0), breaks = breaks)+
      geom_boxplot(mapping = aes(fill = Phenotype), colour = "black", show.legend = F, fatten = 0.75, outliers = FALSE)+
      geom_jitter(mapping = aes(fill = Treatment, shape = Sex), stroke = 1.5/1.427, size = 2, width = 0.2, height = 0, show.legend = TRUE) +
      scale_fill_manual(values = c(coloursP, coloursT), labels = c(LabelsP, LabelsT)) +
      scale_shape_manual(values = c("Female" = 21, "Male" = 22))+
      theme(panel.background = element_rect(fill="transparent"),
            axis.line = element_line(colour="black", linewidth = 1/1.427),
            plot.title = element_text(hjust=0.5, size=13),
            plot.subtitle = element_text(hjust=0.5, size=11),
            axis.ticks.y = element_line(linewidth = 1/1.427, colour = "black"),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(colour="black", size = 11),
            axis.text.x = element_blank(),
            axis.title.y = element_text(size=12),
            axis.title.x = element_blank(),
            legend.key.height = unit(0.6, "cm"),
            legend.title = element_blank(),
            legend.key = element_blank(),
            legend.text = element_text(size=11, hjust = 0),
            legend.margin = ggplot2::margin(-5.5, 5.5, -5.5, 5.5),
            plot.caption = element_text (colour = "black", size = 10, hjust = 1),
            panel.grid = element_blank())+
      labs(title = name) +
      #coord_cartesian(clip = "off")+
      coord_fixed(clip = "off", ratio = (5.76/plotrange))+
      guides(fill = guide_legend(order = 1,
                                 override.aes = list(shape  = c(22, 22, 21, 21),
                                                     size = c(4, 4, 2, 2))),
             shape = guide_legend(order = 2, override.aes = list(fill = "grey")))
    
    if(input$metshowstats){
      
      metsdf <- subset(met.calcdata()$stat.data, Metabolite == selected_mfeature()) %>%
        mutate(sig = unlist(lapply(X = adj.P.Val, FUN = function(p) {
          if (is.na(p) == TRUE){return(NA)}
          else if (p < 0.0001) {return("****")}
          else if (p < 0.001) {return ("***")}
          else if (p < 0.01) {return ("**")}
          else if (p < 0.05) {return("*")}
          else {return ("ns")}
        })), 
               barx1 = as.numeric(Group1),
               barx2 = as.numeric(Group2),
               starx = (barx1 + barx2)/2,
               ybarfact = mapply(FUN = function(grp1, grp2){
                 if ((grp1 == 1 & grp2 == 2) | (grp2 == 1 & grp1 == 2)){return(2)}
                 else if ((grp1 == 3 & grp2 == 4) | (grp2 == 3 & grp1 == 4)){return(2)}
                 else if ((grp1 == 1 & grp2 == 3) | (grp2 == 1 & grp1 == 3)){return(1)}
                 else if ((grp1 == 2 & grp2 == 4) | (grp2 == 2 & grp1 == 4)){return(0)}
               }, grp1 = as.numeric(Group1), grp2 = as.numeric(Group2)))
      
      p <- p+
        geom_segment(inherit.aes = FALSE, data = metsdf, mapping = aes(x = barx1, xend = barx2, y = plotrange*(0.98-ybarfact*0.07-0.025)+plotmin, yend = plotrange*(0.98-ybarfact*0.07-0.025)+plotmin))
      
      if(input$metstatrep == 1){
        p <- p +
          geom_text(inherit.aes = FALSE, data = metsdf, mapping = aes(x = starx, label = paste("p = ", signif(adj.P.Val, 3)), y = plotrange*(0.98-ybarfact*0.07)+plotmin), size = 3, parse = FALSE)
        
      }
      if(input$metstatrep == 2){
        p <- p +
          geom_text(inherit.aes = FALSE, data = metsdf, mapping = aes(x = starx, label = sig, y = plotrange*(0.98-ybarfact*0.07)+plotmin))
        
      }
      
    }
    
    
    
    p
    
  },
  res = 100
  ) #end of output metintplot
  
  #output$metintmessage ----
  output$metintmessage <- renderUI({
    if (!nzchar(selected_mfeature())) {
      tags$div(
        style = "
        height: 100%;
        display: flex;
        align-items: center;
        justify-content: center;
        text-align: center;
        padding: 50px;
      ",
        p(HTML("To display a boxplot of log<sub>2</sub>-transformed normalised intensity by group, <b>search for a 
          feature</b> in the box on the left or <b>select one from the volcano plot</b>"), 
          class = "msg")
      )
      
    } else if(!selected_mfeature() %in% met.calcdata()$stat.data$Metabolite) {
      tags$div(
        style = "
        height: 100%;
        display: flex;
        align-items: center;
        justify-content: center;
        text-align: center;
        padding: 50px;
      ",
        p("Metabolite doesn't meet confidence filtering", 
        class = "msg")
      )
    }  else{
      tags$div(
        style = "
        height: 100%;
        display: flex;
        align-items: center;
        justify-content: center;
        text-align: center;
        padding: 50px;
      ",
      ""
      )
    }
  }) #end of output$metintmessage
  
  #output$metint_content----
  output$metint_content <- renderUI({

    if (nzchar(selected_mfeature()) & selected_mfeature() %in% met.calcdata()$stat.data$Metabolite) {

      tags$div(
        #style = "min-height: 250px;",
        plotOutput("metintplot", height = "350px")
      )

    } else {

      tags$div(
        style = "
        height: 100%;
        display: flex;
        align-items: center;
        justify-content: center;
        text-align: center;
        padding: 50px;
      ",
        uiOutput("metintmessage")
      )

    }
  })
  
  ##### TRANCRIPTOMICS AND PROTEOMICS ##### 
  
  #defining data as basis for volcano plot ----
  tpvolc_base <- reactive({
    req(input$volccomp)
    tvolcdf <- subset(trans.statdata(), comp == input$volccomp)
    pvolcdf <- subset(prot.statdata(), comp == input$volccomp)
    tpfc = input$volcfc
    tpp = input$volcp
    
    transvolcplot <- ggplot(data = tvolcdf, mapping = aes(x = logFC, y = -log10(adj.P.Val)))+
      geom_hline(yintercept = -log10(tpp), colour= "black", linetype = "dotted")+
      geom_vline(xintercept = c(-log2(tpfc), log2(tpfc)), colour= "black", linetype = "dotted")+
      geom_point(data = tvolcdf %>% dplyr::slice_sample(n = 1000) %>% subset(adj.P.Val >= tpp |
                                                                                    abs(logFC) < log2(tpfc) ),
                 colour = "grey", alpha = 0.5)+
      geom_point(data = tvolcdf %>% subset(adj.P.Val < tpp &
                                                  logFC > log2(tpfc)),
                 mapping = aes(key = gene, text = paste0(
                   "<b><i>", gene, "</b></i><br>", 
                   "log<sub>2</sub>(FC): ", signif(logFC, 3), "<br>", 
                   "FDR: ", signif(adj.P.Val, 3))), 
                 colour = "#862164", alpha = 0.5)+
      geom_point(data = tvolcdf %>% subset(adj.P.Val < tpp &
                                                  logFC < -log2(tpfc)),
                 mapping = aes(key = gene, text = paste0(
                   "<b><i>", gene, "</b></i><br>", 
                   "log<sub>2</sub>(FC): ", signif(logFC, 3), "<br>", 
                   "FDR: ", signif(adj.P.Val, 3))), 
                 colour = "#006DAE", alpha = 0.5)+
      scale_x_continuous(name = expression("log"["2"]~"fold change"), limits = c(-15, 15), expand = expansion())+
      scale_y_continuous(name = expression("-log"["10"]~"p-value") , limits = c(0, 15), expand = expansion())+
      theme(panel.background = element_rect(fill="transparent"),
            axis.line = element_line(colour="black", linewidth = 1/1.427),
            plot.title = element_text(hjust=0.5, size=13),
            plot.subtitle = element_text(hjust=0.5, size=11),
            axis.ticks= element_line(linewidth = 1/1.427, colour = "black"),
            axis.text = element_text(colour="black", size = 11),
            axis.title = element_text(size=12),
            legend.key.height = unit(0.6, "cm"),
            legend.title = element_blank(),
            legend.key = element_blank(),
            legend.text = element_text(size=11, hjust = 0),
            legend.margin = ggplot2::margin(-5.5, 5.5, -5.5, 5.5),
            plot.caption = element_text (colour = "black", size = 10, hjust = 1),
            panel.grid = element_blank())+
      annotate(geom = "text", x =  ((-15+log2(input$volcfc))/2)-log2(input$volcfc), y = 14.5, label = tvolcdf[1, "Group1"], colour = "#006DAE", hjust = 0.5)+
      annotate(geom = "text", x = ((15-log2(input$volcfc))/2)+log2(input$volcfc), y = 14.5, label =  tvolcdf[1, "Group2"], colour = "#862164", hjust = 0.5)+
      annotate(geom = "text", x = ((-15+log2(input$volcfc))/2)-log2(input$volcfc), y = .5, label = tvolcdf %>% subset(adj.P.Val < tpp &
                                                                                  -logFC > log2(tpfc)) %>% nrow(),
               colour = "#006DAE", hjust = 0.5, vjust = 1)+
      annotate(geom = "text", x = ((15-log2(input$volcfc))/2)+log2(input$volcfc), y = 0.5, label = tvolcdf %>% subset(adj.P.Val < tpp &
                                                                                  logFC > log2(tpfc)) %>% nrow(),
               colour = "#862164", hjust = 0.5, vjust = 1)+
      labs(title = "Transcriptomics")
    
    
    message("sig can prot up")
    pvolcdf %>% 
      subset(adj.P.Val < tpp & logFC > log2(tpfc)) %>% 
      str() %>% 
      print() 
    protvolcplot <- ggplot(data = pvolcdf, mapping = aes(x = logFC, y = -log10(adj.P.Val), 
                                                   key = ProteinID, text = paste0(
                                                     "<b>", ProteinID, "</b><br>", 
                                                     "log<sub>2</sub>(FC): ", signif(logFC, 3), "<br>", 
                                                     "FDR: ", signif(adj.P.Val, 3))))+
      geom_hline(yintercept = -log10(tpp), colour= "black", linetype = "dotted")+
      geom_vline(xintercept = c(-log2(tpfc), log2(tpfc)), colour= "black", linetype = "dotted")+
      geom_point(data = pvolcdf %>% subset(adj.P.Val >= tpp |
                                                  abs(logFC) < log2(tpfc)),
                 colour = "grey", alpha = 0.5)+
      geom_point(data = pvolcdf %>% subset(adj.P.Val < tpp &
                                                  logFC > log2(tpfc)),
                 colour = "#862164", alpha = 0.5)+
      geom_point(data = pvolcdf %>% subset(adj.P.Val < tpp &
                                                  logFC < -log2(tpfc)),
                 colour = "#006DAE", alpha = 0.5)+
      scale_x_continuous(name = expression("log"["2"]~"fold change"), limits = c(-10, 10), expand = expansion())+
      scale_y_continuous(name = expression("-log"["10"]~"p-value") , limits = c(0, 20), expand = expansion())+
      theme(panel.background = element_rect(fill="transparent"),
            axis.line = element_line(colour="black", linewidth = 1/1.427),
            plot.title = element_text(hjust=0.5, size=13),
            plot.subtitle = element_text(hjust=0.5, size=11),
            axis.ticks= element_line(linewidth = 1/1.427, colour = "black"),
            axis.text = element_text(colour="black", size = 11),
            axis.title = element_text(size=12),
            legend.key.height = unit(0.6, "cm"),
            legend.title = element_blank(),
            legend.key = element_blank(),
            legend.text = element_text(size=11, hjust = 0),
            legend.margin = ggplot2::margin(-5.5, 5.5, -5.5, 5.5),
            plot.caption = element_text (colour = "black", size = 10, hjust = 1),
            panel.grid = element_blank())+
      annotate(geom = "text", x = ((-10+log2(input$volcfc))/2)-log2(input$volcfc), y = 14.5/.75, label = pvolcdf[1, "Group1"], colour = "#006DAE", hjust = 0.5)+
      annotate(geom = "text", x = ((10-log2(input$volcfc))/2)+log2(input$volcfc), y = 14.5/.75, label =  pvolcdf[1, "Group2"], colour = "#862164", hjust = .5)+
      annotate(geom = "text", x = ((-10+log2(input$volcfc))/2)-log2(input$volcfc), y = .5/.75, label = pvolcdf %>% subset(adj.P.Val < tpp &
                                                                                          -logFC > log2(tpfc)) %>% nrow(),
               colour = "#006DAE", hjust = 0.5, vjust = 1)+
      annotate(geom = "text", x = ((10-log2(input$volcfc))/2)+log2(input$volcfc), y = 0.5/.75, label = pvolcdf %>% subset(adj.P.Val < tpp &
                                                                                          logFC > log2(tpfc)) %>% nrow(),
               colour = "#862164", hjust = .5, vjust = 1)+
      labs(title = "Proteomics")#+
    #coord_fixed(clip = "off", ratio = 1.5)

    list(
      tvolcdf = tvolcdf, 
      pvolcdf = pvolcdf, 
      transvolcplot = transvolcplot, 
      protvolcplot = protvolcplot 
      
    )
    
  })
  #output$transvolcano ----
  output$transvolcano <- renderPlotly({
    req(
      input$main_nav == "transprot"
    )
    
    base <- tpvolc_base()

    p <- base$transvolcplot

    #+
    #coord_fixed(clip = "off", ratio = 1.5)

    if(input$volcpath & nzchar(input$volcpathsel)){
     

      tvolcgenes <- pathwaylist.all()[[input$volcpathsel]]

      p <- p+
        geom_point(data = base$tvolcdf %>% subset(gene %in% tvolcgenes),
                   mapping = aes(key = gene, text = paste0(
                     "<b><i>", gene, "</b></i><br>", 
                     "log<sub>2</sub>(FC): ", signif(logFC, 3), "<br>", 
                     "FDR: ", signif(adj.P.Val, 3))), 
                   colour = "black", alpha = 0.8)


    } else if (!is.null(selected_tpfeature)){
      p <- p+
        geom_point (data = base$tvolcdf %>% subset(gene == selected_tpfeature()),
                    mapping = aes(key = gene, text = paste0(
                      "<b><i>", gene, "</b></i><br>", 
                      "log<sub>2</sub>(FC): ", signif(logFC, 3), "<br>", 
                      "FDR: ", signif(adj.P.Val, 3))),
                    colour = "black", alpha = 0.8)
    } else {p <- p}



    ggplotly(p, tooltip = "text", source = "transprott") %>%
      layout(
        xaxis = list(title = "log<sub>2</sub> Fold Change (FC)", titlefont = list(size = 14), tickfont = list(size = 10)),
        yaxis = list(title = "-log<sub>10</sub> False Discovery Rate (FDR)", titlefont = list(size = 14), tickfont = list(size = 10))
      ) %>%
      # plotly::add_markers(
      #   x = c(NA_real_), 
      #   y = c(NA_real_), 
      #   type = "scatter", 
      #   mode = "markers", 
      #   marker = list(color = "black", size = 10, opacity = 0.8), 
      #   name = "selected", 
      #   hoveron = "points", 
      #   text = "selected",
      #   inherit = FALSE
      #   
      # ) %>%
      plotly::event_register("plotly_click")
    

  }
  ) #end of output$transvolcano

  
  #output$transcpmplot ----
  output$transcpmplot <- renderPlot({
    req(
      input$main_nav == "transprot",
      !is.null(selected_tpfeature()),
      length(selected_tpfeature()) == 1,
      nzchar(selected_tpfeature()), 
      selected_tpfeature() %in% trans.statdata()$gene
    )
  


    # boxplot of log2 cpm ---
    # with requested gene
    # 1. It is "reactive" and therefore should be automatically
    #    re-executed when inputs (input$gene) change
    # 2. Its output type is a plot
    transdf <- subset(trans.plotdata(), gene == selected_tpfeature())


    name <- paste(selected_tpfeature())
    Max <- max(transdf[, "CPMlog2"])
    Min <- min(transdf[, "CPMlog2"])
    # if(Min < 0){
    #   
    #   if(Max < 0){
    #     pbs <- prettybreaks_neg(min = Min, max = 1)
    #     plotrange <- pbs[2]-pbs[1]
    #     plotmin <- pbs[1]
    #     plotmax <- pbs[2]
    #     breaks <- seq(pbs[1], pbs[2], plotrange/pbs[3])
    #   }
    #   else{
    #     pbs <- prettybreaks_neg(min = Min, max = Max)
    #     plotrange <- pbs[2]-pbs[1]
    #     plotmin <- pbs[1]
    #     plotmax <- pbs[2]
    #     breaks <- seq(pbs[1], pbs[2], plotrange/pbs[3])
    #   }
    #   
    # }
    # else{
    pbs <- prettybreaks_log2(Min, Max)
    plotrange <- pbs[2]-pbs[1]
    plotmin <- pbs[1]
    plotmax <- pbs[2]
    breaks <- seq(pbs[1], pbs[2], plotrange/pbs[3])
    #plotrange <- plotmax - plotmin
    # }

    set.seed(seed = 793)


    p <- ggplot(data = transdf, mapping = aes(x = Group, y = CPMlog2))+
      scale_x_discrete(expand=expansion(mult=0, add=0.67))+
      scale_y_continuous(name = expression("log"["2"]~"CPM"), limits = c(plotmin, plotmax),
                         expand = expansion(mult = 0, add = 0), breaks = breaks)+
      geom_boxplot(mapping = aes(fill = Phenotype), colour = "black", show.legend = F, fatten = 0.75, outliers = FALSE)+
      geom_jitter(mapping = aes(fill = Treatment, shape = Sex), stroke = 1.5/1.427, size = 2, width = 0.2, height = 0, show.legend = TRUE) +
      scale_fill_manual(values = c(coloursP, coloursT), labels = c(LabelsP, LabelsT)) +
      scale_shape_manual(values = c("Female" = 21, "Male" = 22))+
      theme(panel.background = element_rect(fill="transparent"),
            axis.line = element_line(colour="black", linewidth = 1/1.427),
            plot.title = element_text(hjust=0.5, size=13),
            plot.subtitle = element_text(hjust=0.5, size=11),
            axis.ticks.y = element_line(linewidth = 1/1.427, colour = "black"),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(colour="black", size = 11),
            axis.text.x = element_blank(),
            axis.title.y = element_text(size=12),
            axis.title.x = element_blank(),
            legend.key.height = unit(0.6, "cm"),
            legend.title = element_blank(),
            legend.key = element_blank(),
            legend.text = element_text(size=11, hjust = 0),
            legend.margin = ggplot2::margin(-5.5, 5.5, -5.5, 5.5),
            plot.caption = element_text (colour = "black", size = 10, hjust = 1),
            panel.grid = element_blank())+
      labs(title = bquote(italic(.(selected_tpfeature())))) +
      #coord_cartesian(clip = "off")+
      coord_fixed(clip = "off", ratio = (5.76/plotrange))+
      guides(fill = guide_legend(order = 1,
                                 override.aes = list(shape  = c(22, 22, 21, 21),
                                                     size = c(4, 4, 2, 2))),
             shape = guide_legend(order = 2, override.aes = list(fill = "grey")))

    if(input$showstats){

      transsdf <- subset(trans.statdata(), gene == selected_tpfeature()) %>%
        mutate(barx1 = as.numeric(Group1),
               barx2 = as.numeric(Group2),
               starx = (barx1 + barx2)/2,
               ybarfact = mapply(FUN = function(grp1, grp2){
                 if ((grp1 == 1 & grp2 == 2) | (grp2 == 1 & grp1 == 2)){return(2)}
                 else if ((grp1 == 3 & grp2 == 4) | (grp2 == 3 & grp1 == 4)){return(2)}
                 else if ((grp1 == 1 & grp2 == 3) | (grp2 == 1 & grp1 == 3)){return(1)}
                 else if ((grp1 == 2 & grp2 == 4) | (grp2 == 2 & grp1 == 4)){return(0)}
               }, grp1 = as.numeric(Group1), grp2 = as.numeric(Group2)))

      p <- p+
        geom_segment(inherit.aes = FALSE, data = transsdf, mapping = aes(x = barx1, xend = barx2, y = plotrange*(0.98-ybarfact*0.07-0.025)+plotmin, yend = plotrange*(0.98-ybarfact*0.07-0.025)+plotmin))

      if(input$statrep == 1){
        p <- p +
          geom_text(inherit.aes = FALSE, data = transsdf, mapping = aes(x = starx, label = paste("p = ", signif(adj.P.Val, 3)), y = plotrange*(0.98-ybarfact*0.07)+plotmin), size = 3, parse = FALSE)

      }
      if(input$statrep == 2){
        p <- p +
          geom_text(inherit.aes = FALSE, data = transsdf, mapping = aes(x = starx, label = sig, y = plotrange*(0.98-ybarfact*0.07)+plotmin))

      }

    }



    p

  },
  res = 100
  ) #end of output transcpmplot
  #output$cpmmessage ----
  output$cpmmessage <- renderUI({
    req(input$main_nav == "transprot")
    if (is.null(selected_tpfeature()) || !nzchar(selected_tpfeature())) {
      tags$div(
        style = "
        height: 100%;
        display: flex;
        align-items: center;
        justify-content: center;
        text-align: center;
        padding: 50px;
      ",
        p(HTML("To display a boxplot of log<sub>2</sub>-transformed counts per million (cpm) by group, <b>search for a feature</b> 
          in the box on the left or <b>select one from the volcano plot</b>"), 
          class = "msg")
      )
     
    } else if(!selected_tpfeature() %in% trans.statdata()$gene) {
      tags$div(
        style = "
        height: 100%;
        display: flex;
        align-items: center;
        justify-content: center;
        text-align: center;
        padding: 50px;
      ",
        p("Feature not present in transcriptomics data", 
          class = "msg")
      )
    } else{
      tags$div(
        style = "
        height: 100%;
        display: flex;
        align-items: center;
        justify-content: center;
        text-align: center;
        padding: 50px;
      ", ""
      )
    }
  }) #end of output cpm message
  #output$transint_content----
  output$transint_content <- renderUI({
    req(input$main_nav == "transprot")
    if(is.null(selected_tpfeature()) || !nzchar(selected_tpfeature()) || !selected_tpfeature() %in% trans.statdata()$gene) {
      tags$div(
        style = "
        height: 100%;
        display: flex;
        align-items: center;
        justify-content: center;
        text-align: center;
        padding: 50px;
      ",
        uiOutput("cpmmessage")%>% withSpinner(type = 7)
      )
    }
    else {
      
      tags$div(
        #style = "min-height: 250px;",
        plotOutput("transcpmplot", height = "340px")%>% withSpinner(type = 7)
      )
      
    } 
  })

  

  #output$protvolcano ----
  output$protvolcano <- renderPlotly({
    req(input$main_nav == "transprot")
   
    base <- tpvolc_base()

    p <- base$protvolcplot

    

    if(input$volcpath & nzchar(input$volcpathsel)){
      req(input$volcpathsel)

      pvolcgenes <- pathwaylist.all()[[input$volcpathsel]]

      p <- p+
        geom_point(data = base$pvolcdf %>% subset(gene %in% pvolcgenes),
                   mapping = aes(x = logFC, y = -log10(adj.P.Val), 
                                 key = ProteinID, text = paste0(
                                   "<b>", ProteinID, "</b><br>", 
                                   "log<sub>2</sub>(FC): ", signif(logFC, 3), "<br>", 
                                   "FDR: ", signif(adj.P.Val, 3))), 
                   colour = "black", alpha = 0.8)


    } else if (!is.null(selected_tpfeature())){
      p <- p+
        geom_point (data = base$pvolcdf %>% subset(gene == selected_tpfeature()),
                    mapping = aes(x = logFC, y = -log10(adj.P.Val), 
                                  key = ProteinID, text = paste0(
                                    "<b>", ProteinID, "</b><br>", 
                                    "log<sub>2</sub>(FC): ", signif(logFC, 3), "<br>", 
                                    "FDR: ", signif(adj.P.Val, 3))), 
                    colour = "black", alpha = 0.8)
    } else {p <- p}



    ggplotly(p, tooltip = "text", source = "transprotp") %>%
      layout(
        xaxis = list(title = "log<sub>2</sub> Fold Change (FC)", titlefont = list(size = 14), tickfont = list(size = 10)), 
        yaxis = list(title = "-log<sub>10</sub> False Discovery Rate (FDR)", titlefont = list(size = 14), tickfont = list(size = 10))
      ) %>%
      plotly::event_register("plotly_click")
  }
  ) #end of output$protvolcano

  #output$protintplot ----
  output$protintplot <- renderPlot({
    req(
      input$main_nav == "transprot",
      !is.null(selected_tpfeature()),
      length(selected_tpfeature()) == 1,
      nzchar(selected_tpfeature()), 
      selected_tpfeature() %in% prot.statdata()$gene
    )

    # boxplot of log2 intensity ---
    # with requested gene
    # 1. It is "reactive" and therefore should be automatically
    #    re-executed when inputs (input$gene) change
    # 2. Its output type is a plot
    protdf <- subset(prot.plotdata(), gene == selected_tpfeature())


    name <- paste(selected_tpfeature())
    Max <- max(protdf[, "Intlog2"])
    Min <- min(protdf[, "Intlog2"])
    # if(Min < 0){
    #   
    #   if(Max < 0){
    #     pbs <- prettybreaks_neg(min = Min, max = 1)
    #     plotrange <- pbs[2]-pbs[1]
    #     plotmin <- pbs[1]
    #     plotmax <- pbs[2]
    #     breaks <- seq(pbs[1], pbs[2], plotrange/pbs[3])
    #   }
    #   else{
    #     pbs <- prettybreaks_neg(min = Min, max = Max)
    #     plotrange <- pbs[2]-pbs[1]
    #     plotmin <- pbs[1]
    #     plotmax <- pbs[2]
    #     breaks <- seq(pbs[1], pbs[2], plotrange/pbs[3])
    #   }
    #   
    # }
    # else{
    pbs <- prettybreaks_log2(Min, Max)
    plotrange <- pbs[2]-pbs[1]
    plotmin <- pbs[1]
    plotmax <- pbs[2]
    breaks <- seq(pbs[1], pbs[2], plotrange/pbs[3])
    #plotrange <- plotmax - plotmin
    # }

    set.seed(seed = 793)


    p <- ggplot(data = protdf, mapping = aes(x = Group, y = Intlog2))+
      scale_x_discrete(expand=expansion(mult=0, add=0.67))+
      scale_y_continuous(name = expression("log"["2"]~"Intensity"), limits = c(plotmin, plotmax),
                         expand = expansion(mult = 0, add = 0), breaks = breaks)+
      geom_boxplot(mapping = aes(fill = Phenotype), colour = "black", show.legend = F, fatten = 0.75, outliers = FALSE)+
      geom_jitter(mapping = aes(fill = Treatment, shape = Sex), stroke = 1.5/1.427, size = 2, width = 0.2, height = 0, show.legend = TRUE) +
      scale_fill_manual(values = c(coloursP, coloursT), labels = c(LabelsP, LabelsT)) +
      scale_shape_manual(values = c("Female" = 21, "Male" = 22))+
      theme(panel.background = element_rect(fill="transparent"),
            axis.line = element_line(colour="black", linewidth = 1/1.427),
            plot.title = element_text(hjust=0.5, size=13),
            plot.subtitle = element_text(hjust=0.5, size=11),
            axis.ticks.y = element_line(linewidth = 1/1.427, colour = "black"),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(colour="black", size = 11),
            axis.text.x = element_blank(),
            axis.title.y = element_text(size=12),
            axis.title.x = element_blank(),
            legend.key.height = unit(0.6, "cm"),
            legend.title = element_blank(),
            legend.key = element_blank(),
            legend.text = element_text(size=11, hjust = 0),
            legend.margin = ggplot2::margin(-5.5, 5.5, -5.5, 5.5),
            plot.caption = element_text (colour = "black", size = 10, hjust = 1),
            panel.grid = element_blank())+
      labs(title = str_to_upper(selected_tpfeature())) +
      #coord_cartesian(clip = "off")+
      coord_fixed(clip = "off", ratio = (5.76/plotrange))+
      guides(fill = guide_legend(order = 1,
                                 override.aes = list(shape  = c(22, 22, 21, 21),
                                                     size = c(4, 4, 2, 2))),
             shape = guide_legend(order = 2, override.aes = list(fill = "grey")))

    if(input$showstats){

      protsdf <- subset(prot.statdata(), gene == selected_tpfeature()) %>%
        mutate(barx1 = as.numeric(Group1),
               barx2 = as.numeric(Group2),
               starx = (barx1 + barx2)/2,
               ybarfact = mapply(FUN = function(grp1, grp2){
                 if ((grp1 == 1 & grp2 == 2) | (grp2 == 1 & grp1 == 2)){return(2)}
                 else if ((grp1 == 3 & grp2 == 4) | (grp2 == 3 & grp1 == 4)){return(2)}
                 else if ((grp1 == 1 & grp2 == 3) | (grp2 == 1 & grp1 == 3)){return(1)}
                 else if ((grp1 == 2 & grp2 == 4) | (grp2 == 2 & grp1 == 4)){return(0)}
               }, grp1 = as.numeric(Group1), grp2 = as.numeric(Group2)))

      p <- p+
        geom_segment(inherit.aes = FALSE, data = protsdf, mapping = aes(x = barx1, xend = barx2, y = plotrange*(0.98-ybarfact*0.07-0.025)+plotmin, yend = plotrange*(0.98-ybarfact*0.07-0.025)+plotmin))

      if(input$statrep == 1){
        p <- p +
          geom_text(inherit.aes = FALSE, data = protsdf, mapping = aes(x = starx, label = paste("p = ", signif(adj.P.Val, 3)), y = plotrange*(0.98-ybarfact*0.07)+plotmin), size = 3, parse = FALSE)

      }
      if(input$statrep == 2){
        p <- p +
          geom_text(inherit.aes = FALSE, data = protsdf, mapping = aes(x = starx, label = sig, y = plotrange*(0.98-ybarfact*0.07)+plotmin))

      }

    }



    p

  },
  res = 100
  ) #end of output protintplot

  #output$protintmessage ----
  output$protintmessage <- renderUI({
    req(input$main_nav == "transprot")
    if (is.null(selected_tpfeature()) || !nzchar(selected_tpfeature())) {
      tags$div(
        style = "
        height: 100%;
        display: flex;
        align-items: center;
        justify-content: center;
        text-align: center;
        padding: 50px;
      ",
      p(HTML("To display a boxplot of log<sub>2</sub>-transformed normalised intensity by group, <b>search for a feature</b> 
          in the box on the left or <b>select one from the volcano plot</b>"), 
        class = "msg")
      )
    } else if(!selected_tpfeature() %in% prot.statdata()$gene) {
      tags$div(
        style = "
        height: 100%;
        display: flex;
        align-items: center;
        justify-content: center;
        text-align: center;
        padding: 50px;
      ",
      p("Feature not present in proteomics data", 
        class = "msg")
      )
    } else{
      tags$div(
        style = "
        height: 100%;
        display: flex;
        align-items: center;
        justify-content: center;
        text-align: center;
        padding: 50px;
      ",
      ""
      )
    }
  }) #end of output cpm message
  #output$transint_content----
  output$protint_content <- renderUI({
    req(input$main_nav == "transprot")
    if(is.null(selected_tpfeature()) || !nzchar(selected_tpfeature()) || !selected_tpfeature() %in% prot.statdata()$gene) {
      tags$div(
        style = "
        height: 100%;
        display: flex;
        align-items: center;
        justify-content: center;
        text-align: center;
        padding: 50px;
      ",
        uiOutput("protintmessage")%>% withSpinner(type = 7)
      )
    }
    else {
      
      tags$div(
        #style = "min-height: 250px;",
        plotOutput("protintplot", height = "340px")%>% withSpinner(type = 7)
      )
      
    } 
  })
  
  
  ####DIMENSIONALITY REDUCTION ---- 
  #output$dimredplot ----
  output$dimredplot <- renderPlotly({
    req(input$main_nav == "dimred")
    req(!is.null(input$dimreddata))
    #print(input$dimreddata)
    trans_df <- trans.plotdata() %>%
      dplyr::select(Mouse, Geneid, CPMlog2) %>%
      pivot_wider(values_from = "CPMlog2", names_from = "Mouse") %>%
      as.data.frame() %>% 
      `rownames<-`(paste0("trans_", .$Geneid)) %>% 
      dplyr::select(-Geneid) %>%
      t() %>%
      scale() %>% 
      t()
    
    prot_df <- prot.plotdata() %>%
      dplyr::select(Mouse, gene, Intlog2) %>%
      pivot_wider(values_from = "Intlog2", names_from = "Mouse") %>%
      as.data.frame() %>% 
      `rownames<-`(paste0("prot_", .$gene)) %>% 
      dplyr::select(-gene) %>%
      t() %>%
      scale() %>% 
      t() %>%
      .[, colnames(trans_df)]
    
    met_df <- met.calcdata()$plot.data %>% 
      dplyr::select(Mouse, Metabolite, Intlog2) %>%
      pivot_wider(values_from = "Intlog2", names_from = "Mouse") %>%
      as.data.frame() %>% 
      `rownames<-`(paste0("met_", .$Metabolite)) %>% 
      dplyr::select(-Metabolite) %>%
      t() %>%
      scale() %>% 
      t()%>%
      .[, colnames(trans_df)]
    
    
    if(length(input$dimreddata) == 1){
      dimred_df <- paste0(input$dimreddata[1], "_df") %>% get()
    } else if (length(input$dimreddata) == 2){
      dimred_df <- paste0(input$dimreddata[1], "_df") %>% get() %>%
        rbind(paste0(input$dimreddata[2], "_df") %>% get())
    } else if (length(input$dimreddata) == 3){
      dimred_df <- paste0(input$dimreddata[1], "_df") %>% get() %>%
        rbind(paste0(input$dimreddata[2], "_df") %>% get())%>%
        rbind(paste0(input$dimreddata[3], "_df") %>% get())
    }
    
    #print(dim(dimred_df))


    
    if(input$dimredtype == "mds"){
        
      dist_mat <- as.dist(1-cor(dimred_df))

     
      #print(dist_mat)

      
      mds <- cmdscale(dist_mat, k = 23, eig = TRUE) 
      
      dimredplot <- mds$points %>%
        as.data.frame() %>% 
        `colnames<-`(paste0("Dim", as.character(1:ncol(.)))) %>%
        mutate(Mouse = rownames(.)) %>%
        merge.data.frame(metadata)
      
      var_exp <- mds$eig %>% 
        .[. > 0] %>%
        divide_by(sum(.))
      
      Dimname <- "Dim"
      Tit <- "Multi-dimensional scaling"
      #print(var_exp)
      
      #print(mds.plot)

    } else if(input$dimredtype == "pca"){
      
      pca <- prcomp(dimred_df %>% t(), center = FALSE, scale. = FALSE) %>%
        {dimredplot <<- .$x %>%
        as.data.frame() %>% 
        `colnames<-`(paste0("Dim", as.character(1:ncol(.)))) %>%
        mutate(Mouse = row.names(.)) %>%
        merge.data.frame(y = metadata); .} %>%
        {var_exp <<- summary(.)$importance[2,] %>% 
          as.data.frame() %>% 
          t(); .}
      
      Dimname <- "PC"
      Tit <- "Principal Component Analysis"
    }
    
    
    
    p <- dimredplot %>% 
      mutate(
        Dim1 = Dim1 * (if ("x" %in% input$dimredflip) -1 else 1),
        Dim2 = Dim2 * (if ("y" %in% input$dimredflip) -1 else 1)
      ) %>%
      ggplot(mapping = aes(x = Dim1, y = Dim2, 
                                                 text = paste0("Mouse ", Mouse, "<br><b>", Phenotype, "/", Treatment, "</b>"))) + 
      stat_ellipse(geom = "polygon", mapping = aes(group = Group, fill = Phenotype), colour = "black", linewidth = 1/1.427, alpha = 0.8)+
      geom_point(mapping = aes(fill = Treatment, shape = Phenotype), size = 2)+
      scale_shape_manual(values = c(22, 23), labels = LabelsP)+
      scale_fill_manual(values = c(coloursT, coloursP)) +
      scale_colour_manual(values = colours, labels = Labels) +
      scale_x_continuous(name = paste(Dimname, "1 (", (signif(var_exp[1]*100, digits = 3)), "% var explained)", sep = ""), 
                         expand = expansion(add = 0, mult = 0.2))+
      scale_y_continuous(name = paste(Dimname, "2 (", (signif(var_exp[2]*100, digits = 3)), "% var explained)", sep = ""), 
                         expand = expansion(add = 0, mult = 0.2))+
      theme(panel.background = element_rect(fill="transparent"),
            axis.line = element_line(colour="black", linewidth = 1/1.427),
            plot.title = element_text(hjust=0.5, size=13),
            plot.subtitle = element_text(hjust=0.5, size=11),
            axis.ticks = element_line(linewidth = 1/1.427, colour = "black"),
            axis.text = element_text(colour="black", size = 11),
            axis.title = element_text(size=12),
            legend.key.height = unit(0.6, "cm"),
            legend.title = element_blank(),
            legend.key = element_blank(),
            legend.text = element_text(size=11, hjust = 0),
            legend.margin = ggplot2::margin(-5.5, 5.5, -5.5, 5.5),
            plot.caption = element_text (colour = "black", size = 10, hjust = 1),
            panel.grid = element_blank())+
      labs(title = Tit) +
      coord_cartesian(clip = "off")+
      guides(fill = guide_legend(override.aes = list(shape = c(22, 22, 22, 22))))
  
 
    p<- ggplotly(p, tooltip = "text", source = "dimredplot") %>%
      layout(
        legend = list(title = list(text = "Legend (click to deselect)"))
      )
    
    legendlabels <- c("Control", "Diabetic", "Control/PBS", "Control/Candida", "Diabetic/PBS", "Diabetic/Candida")
    for(i in seq_along(p$x$data)){
      p$x$data[[i]]$name <- legendlabels[i]
    }
  
    
    p
  })
  # #output$dimredmessage ---- 
  # output$dimredmessage <- renderText({
  #   if (is.null(input$dimreddata)) {
  #     "\n\nSelect at least one data set to display a dimensionality reduction plot"
  #   } else{
  #     ""
  #   }
  # }) #end of output cpm message
  
  #output$dimred_content----
  output$dimred_content <- renderUI({
    req(input$main_nav == "dimred")
    if (is.null(input$dimreddata)) {
      tags$div(
        style = "
        height: 100%;
        display: flex;
        align-items: center;
        justify-content: center;
        text-align: center;
        padding: 50px;
      ",
        p("Select at least one data set to display a dimensionality reduction plot", 
          class = "msg")
      )
    }
    else {
      tags$div(
        style = "max-width: 600px;",
        plotlyOutput("dimredplot", height = "340px") %>% withSpinner(type = 7)
      )
    } 
  })
  
  ##### PATHWAY ANALYSIS ##### 
  #PAresults reactive ---- 
  PAresults <- eventReactive(input$pathwaycalc, {
    selection <- subset(pathway.options, id == input$pathwayid)
    
    if(nzchar(selection[1, "gs_subcollection"])){
      genesets <- msigdbr(db_species = selection[1, "dbs"], species = "Mus musculus", collection = selection[1,"gs_collection"], subcollection = selection[1,"gs_subcollection"])
    } else{
      genesets <- msigdbr(db_species = selection[1, "dbs"], species = "Mus musculus", collection = selection[1,"gs_collection"])
    }
    
    genelist <- split(genesets$gene_symbol, genesets$gs_name)
   
    if(input$pathwaymethod == "GSEA"){
      
      PAtransdata <- trans.statdata() %>% 
        subset(comp == input$pathwaycomp) 
      PAprotdata <- prot.statdata() %>% 
        subset(comp == input$pathwaycomp) 
      
      PAtranslist <- PAtransdata$logFC %>% 
        `names<-`(PAtransdata$Geneid) %>% 
        sort(decreasing = TRUE)
      PAprotlist <- PAprotdata$logFC %>% 
        `names<-`(PAprotdata$gene) %>% 
        sort(decreasing = TRUE)
      
      gseatrans <- GSEA(
        geneList = PAtranslist, 
        TERM2GENE = genesets[, c("gs_name", "gene_symbol")], 
        pvalueCutoff = 0.05,
        seed = TRUE,
        BPPARAM = BiocParallel::SerialParam()
      )
      gseaprot <- GSEA(
        geneList = PAprotlist, 
        TERM2GENE = genesets[, c("gs_name", "gene_symbol")], 
        pvalueCutoff = 0.05,
        seed = TRUE,
        BPPARAM = BiocParallel::SerialParam()
      )
      transresult <- gseatrans@result 
      protresult <- gseaprot@result
    }
    
    else if(input$pathwaymethod == "ORA"){
      PAtransdata <- trans.statdata() %>% 
        subset(comp == input$pathwaycomp & adj.P.Val < input$pathwayp) %>%
          filter(if(input$pathwaydir == "both"){abs(logFC) > log2(input$pathwayfc) }
                 else if (input$pathwaydir == "up"){ logFC > log2(input$pathwayfc)}
                 else if(input$pathwaydir == "down"){-logFC > log2(input$pathwayfc)})
      PAprotdata <- prot.statdata() %>% 
        subset(comp == input$pathwaycomp & adj.P.Val < input$pathwayp) %>%
        filter(if(input$pathwaydir == "both"){abs(logFC) > log2(input$pathwayfc) }
               else if (input$pathwaydir == "up"){ logFC > log2(input$pathwayfc)}
               else if(input$pathwaydir == "down"){-logFC > log2(input$pathwayfc)})
      
      bkgd.trans <- trans.statdata()[trans.statdata()$comp == input$pathwaycomp, "gene"] %>% unique()
      bkgd.prot <- prot.statdata()[prot.statdata()$comp == input$pathwaycomp, "gene"] %>% unique()

      oratrans <- enricher(
        gene = PAtransdata$gene, 
        TERM2GENE = genesets[, c("gs_name", "gene_symbol")], 
        universe = bkgd.trans, 
        pvalueCutoff = 0.05
      )
      oraprot <- enricher(
        gene = PAprotdata$gene, 
        TERM2GENE = genesets[, c("gs_name", "gene_symbol")], 
        universe = bkgd.prot, 
        pvalueCutoff = 0.05
      )
      
      transresult <- oratrans@result
      protresult <- oraprot@result
      
    }
    results <- rbind(transresult %>% mutate(Dataset = "Transcriptomics"), 
                     protresult %>% mutate(Dataset = "Proteomics")) %>%
      mutate(Geneset = mapply(FUN = str_extract, string = Description, MoreArgs = list(pattern = "^[^_]*") ), 
             Pathway = mapply(FUN = substr, x = Description, start = nchar(Geneset)+2, stop = nchar(Description)) 
             %>% gsub(pattern = "_", replacement = " ", x = .) %>% str_wrap(width = 40)) %>%
      mutate(Dataset = factor(Dataset, levels = c("Transcriptomics", "Proteomics")))
    
    
    if(input$pathwaymethod == "GSEA"){
      order <- results %>% 
        arrange(desc(NES)) %>% 
        .[, "Pathway"] %>% 
        unique() 
      
      sigorder <- results %>% 
        arrange(p.adjust) %>% 
        .[, "Pathway"] %>% 
        unique()
      
      results <- results %>%
        mutate(Count = str_count(core_enrichment, "/") + 1) 
      
      maxcount <- results %>% 
        group_by(Dataset) %>% 
        summarise(MaxCount = max(Count))

      
      p <- results %>%
        merge.data.frame(y = maxcount) %>%
        mutate(RelCount = Count/MaxCount, 
               RelCountSize = Count/setSize) %>% 
        subset(p.adjust < 0.05 & Pathway %in% sigorder[1:20]) %>% 
        mutate(Pathway = factor(Pathway, levels = rev(order))) %>% 
        ggplot(mapping = aes(x = NES, y = Pathway))+
        facet_wrap(facets = vars(Dataset), nrow = 1) + 
        geom_vline(xintercept = 0) +
        geom_segment(mapping = aes(y = Pathway, yend = Pathway, x = 0, xend = NES))+
        geom_point(mapping = aes(size = 100*RelCountSize, fill = -log10(p.adjust)),shape = 21, colour = "black")+
        scale_x_continuous(name = "Normalised enrichment score", expand = expansion(mult = c(0.2 ,0.2), add = c(0, 0)))+
        scale_size_continuous(range = c(2, 10), name = "Core genes (%set size)")+
        scale_fill_gradient(name = expression("-log"["10"]~"(adj. p-value)"), low = "white", high = "#862164")+
        theme(panel.background = element_rect(fill="transparent"),
              axis.line.x = element_line(colour="black", linewidth = 1/1.427),
              axis.line.y = element_blank(), 
              plot.title = element_text(hjust=0.5, size=14, face = "bold"),
              plot.subtitle = element_text(hjust=0.5, size=11),
              axis.ticks.x = element_line(linewidth = 1/1.427, colour = "black"),
              axis.ticks.y = element_blank(), 
              axis.text = element_text(colour="black", size = 11),
              axis.title.x = element_text(size=12),
              axis.title.y = element_blank(), 
              legend.key.height = unit(0.6, "cm"),
              legend.title = element_text(size=11, hjust = 0.5),
              legend.key = element_blank(),
              legend.text = element_text(size=11, hjust = 0),
              #legend.margin = margin(-5.5, 5.5, -5.5, 5.5),
              plot.caption = element_text (colour = "black", size = 10, hjust = 1),
              panel.grid = element_blank(), 
              panel.spacing = unit(2, "lines"), 
              strip.background = element_blank(), 
              strip.text = element_text(size = 13, colour = "black", face = "bold"))+
        coord_cartesian(clip = "off")+
        labs(title = paste0("Top ", results$Geneset[1], " gene sets"))
    }
    else if(input$pathwaymethod == "ORA"){
      
      order <- results %>% 
        arrange(desc(FoldEnrichment)) %>% 
        .[, "Pathway"] %>% 
        unique() 
      
      sigorder <- results %>% 
        arrange(p.adjust) %>% 
        .[, "Pathway"] %>% 
        unique()
      
      maxCount <- results %>% 
        dplyr::group_by(Dataset) %>% 
        dplyr::summarise(MaxCount = max(Count)) %>%
        as.data.frame()

       p <- results %>%
         separate(BgRatio, into = c("pathway_size", "background"),sep = "/")%>%
         mutate(pathway_size = as.numeric(pathway_size), 
                coverage = Count / pathway_size) %>% 
         subset(p.adjust < 0.05 & Pathway %in% sigorder[1:20]) %>% 
         merge.data.frame(y = maxCount) %>% 
         mutate(RelCount = Count/MaxCount) %>% 
         mutate(Pathway = factor(Pathway, levels = rev(order))) %>% 
         ggplot(mapping = aes(x = FoldEnrichment, y = Pathway))+
         facet_wrap(facets = vars(Dataset), nrow = 1) + 
         geom_vline(xintercept = 0) +
         geom_segment(mapping = aes(y = Pathway, yend = Pathway, x = 0, xend = FoldEnrichment))+
         geom_point(mapping = aes(size = coverage*100, fill = -log10(p.adjust)),shape = 21, colour = "black")+
         scale_x_continuous(name = "Fold Enrichment", expand = expansion(mult = c(0,0.2), add = c(0, 0)))+
         scale_size_continuous(range = c(2, 10), name = "Sig. genes (%set size)")+
         scale_fill_gradient(name = expression("-log"["10"]~"(adj. p-value)"), low = "white", high = "#862164")+
         theme(panel.background = element_rect(fill="transparent"),
               axis.line.x = element_line(colour="black", linewidth = 1/1.427),
               axis.line.y = element_blank(), 
               plot.title = element_text(hjust=0.5, size=13, face = "bold"),
               plot.subtitle = element_text(hjust=0.5, size=11),
               axis.ticks = element_line(linewidth = 1/1.427, colour = "black"),
               axis.text = element_text(colour="black", size = 11),
               axis.title.x = element_text(size=12),
               axis.title.y = element_blank(), 
               legend.key.height = unit(0.6, "cm"),
               legend.title = element_text(size=11, hjust = 0.5),
               legend.key = element_blank(),
               legend.text = element_text(size=11, hjust = 0),
               #legend.margin = margin(-5.5, 5.5, -5.5, 5.5),
               plot.caption = element_text (colour = "black", size = 10, hjust = 1),
               panel.grid = element_blank(), 
               panel.spacing = unit(2, "lines"), 
               strip.background = element_blank(), 
               strip.text = element_text(size = 12, colour = "black", face = "bold"))+
         coord_cartesian(clip = "off")+
         labs(title = paste0("Top ", results$Geneset[1], " gene sets"))
    }
    list(trans = transresult, 
         prot = protresult, 
         plot = p)
    }
  )
  
  #output$PAtable ----
  output$PAtable <- renderDataTable({
    req(input$main_nav == "pathway")
    
    if(input$pathwaymethod == "GSEA"){
      PAresults()$trans %>%
        mutate(Dataset = "Transcriptomics") %>% 
        rbind(PAresults()$prot %>% 
                mutate(Dataset = "Proteomics")) %>% 
        select(-ID) %>% 
        select(Dataset, everything()) %>% 
        dplyr::arrange(p.adjust) %>% 
        DT::datatable(., rownames = FALSE,
                      options = list(
                        columnDefs = list(list(
                          targets = 5:7,
                          render = JS(
                            "function(data, type, row) {",
                            "  return type === 'display' ?",
                            "    Number(data).toExponential(2) : data;",
                            "}"
                          )
                        ),
                        
                        list(
                          targets = 3:4,
                          render = JS(
                            "function(data, type, row) {",
                            "  return type === 'display' ?",
                            "    Number(data).toPrecision(4) : data;",
                            "}"
                          )
                        ))
                      )
        )
    } else if(input$pathwaymethod == "ORA"){
      PAresults()$trans %>%
        mutate(Dataset = "Transcriptomics") %>% 
        rbind(PAresults()$prot %>% 
                mutate(Dataset = "Proteomics")) %>% 
        select(-ID) %>% 
        select(Dataset, everything()) %>% 
        dplyr::arrange(p.adjust) %>% 
        DT::datatable(., rownames = FALSE,
                      options = list(
                        columnDefs = list(list(
                          targets = 7:9,
                          render = JS(
                            "function(data, type, row) {",
                            "  return type === 'display' ?",
                            "    Number(data).toExponential(2) : data;",
                            "}"
                          )
                        ),
                        
                        
                        list(
                          targets = 4:6,
                          render = JS(
                            "function(data, type, row) {",
                            "  return type === 'display' ?",
                            "    Number(data).toPrecision(4) : data;",
                            "}"
                          )
                        ))
                      )
        )
    }
    
    
  })
  
  #output$PAplot ---- 
  output$PAplot <- renderPlot({
    req(input$main_nav == "pathway", 
        input$pathwaytabs == "pathwayplots", 
        nrow(PAresults()$trans > 0),
        nrow(PAresults()$prot > 0))
    
    PAresults()$plot
  })
  
  #output$PAplotcont ---- 
  PAplotmsg <- reactiveVal(TRUE)
  observeEvent(input$pathwaycalc, {
    PAplotmsg(FALSE)
  })
  output$PAplotcont <- renderUI({
    if(PAplotmsg() == TRUE){
      tags$div(
        style = "
        height: 100%;
        display: flex;
        align-items: center;
        justify-content: center;
        text-align: center;
        padding: 50px;
      ",
       p(HTML('To perform pathway analysis, select a method, pathway collection and comparison of interest from the left and hit <i>Analyse</i>.<br>
          <b>Gene Set Enrichment Analysis (GSEA; slow)</b> ranks genes by their fold change and does not take p-values into account.<br>
           <b>Over-Representation Analysis (ORA; fast)</b> only considers signficant genes, and cut-offs for fold change and FDR as well as
              the direction of change can be specified.'), 
         class = "msg")
      )
      
    }
    else{
      tags$div(
        #style = "min-height: 250px;",
        plotOutput(outputId = "PAplot", height = "500px") %>% withSpinner(type = 7)
      )
      
    }
  })
  ##### CORRELATIONS ##### 
  #putting all three omics data tables plus misc data together once ----
  has_run <- reactiveVal(FALSE)
  
  all.data <- reactiveVal(NULL)
  observeEvent(input$main_nav, {
    if(input$main_nav == "corr" && !has_run()){
      has_run(TRUE)
      
      all.data <- trans.plotdata() %>% 
        dplyr::select(Mouse, Geneid, CPMlog2) %>% 
        pivot_wider(id_cols = "Mouse", names_from = "Geneid", values_from = "CPMlog2", names_prefix = "Trans_") %>% 
        merge.data.frame(prot.plotdata() %>% 
                dplyr::select(Mouse, gene, Intlog2) %>% 
                pivot_wider(id_cols = "Mouse", names_from = "gene", values_from = "Intlog2", names_prefix = "Prot_"), by = "Mouse") %>% 
        merge.data.frame(met.calcdata()$met.data %>% 
                dplyr::select(Mouse, Metabolite, Intlog2) %>% 
                pivot_wider(id_cols = "Mouse", names_from = "Metabolite", values_from = "Intlog2", names_prefix = "Met_"), by = "Mouse") %>% 
        merge.data.frame(misc.data.trans(), by = "Mouse") %>%
        as.data.frame()
      
      all.data(all.data)

    }
  })
  
  
  #creating reactives to store selected features 
  corrftr_x <- reactiveVal()
  corrtype_x <- reactiveVal()
  corrftr_y <- reactiveVal()
  corrtype_y <- reactiveVal()
  
  observeEvent(input$corrplotgo, {
    if(input$corrxtype == "Transcriptomics"){
      corrtype_x ("Transcriptomics")
      corrftr_x(input$corrx_trans)
    } else if(input$corrxtype == "Proteomics"){
      corrtype_x ("Proteomics")
      corrftr_x(input$corrx_prot)
    } else if (input$corrxtype == "Metabolomics"){
      corrtype_x ("Metabolomics")
      corrftr_x(input$corrx_met)
    } else if (input$corrxtype == "Miscellaneous"){
      corrtype_x ("Miscellaneous")
      corrftr_x(input$corrx_misc)
    }
    if(input$corrytype == "Transcriptomics"){
      corrtype_y ("Transcriptomics")
      corrftr_y(input$corry_trans)
    } else if(input$corrytype == "Proteomics"){
      corrtype_y ("Proteomics")
      corrftr_y(input$corry_prot)
    } else if (input$corrytype == "Metabolomics"){
      corrtype_y ("Metabolomics")
      corrftr_y(input$corry_met)
    } else if (input$corrytype == "Miscellaneous"){
      corrtype_y ("Miscellaneous")
      corrftr_y(input$corry_misc)
    }
  })
  
  
  #corrdata reactive----
  corr <- eventReactive(input$corrplotgo, {
    if (input$corrxtype == "Transcriptomics"){
      corrftrtypex <- "gene"
      corrdatax <- trans.plotdata() %>%
        subset(get(corrftrtypex) == corrftr_x())%>% 
        dplyr::rename("x.name" = corrftrtypex, "x.val" = "CPMlog2", "x.lab" = "Geneid") %>% 
        dplyr::select(c(x.name, x.val, x.lab, Mouse))
    }
    else if (input$corrxtype == "Proteomics"){
      corrftrtypex <- "gene"
      corrdatax <- prot.plotdata() %>%
        subset(get(corrftrtypex) == corrftr_x())%>% 
        dplyr::rename("x.name" = corrftrtypex, "x.val" = "Intlog2", "x.lab" = "ProteinID")%>% 
        dplyr::select(c(x.name, x.val, x.lab, Mouse))
    }
    else if (input$corrxtype == "Metabolomics"){
      corrftrtypex <- "Metabolite"
      corrdatax <- met.calcdata()$met.data %>%
        subset(get(corrftrtypex) == corrftr_x())%>% 
        dplyr::rename("x.name" = corrftrtypex, "x.val" = "Intlog2", "x.lab" = "Label")%>% 
        dplyr::select(c(x.name, x.val, x.lab, Mouse))
    }
    else if (input$corrxtype == "Miscellaneous"){
      corrftrtypex <- "Measurement"
      corrdatax <- misc.data.trans() %>% 
        pivot_longer(cols = Misc_KIF:Misc_Weightloss_abs, names_to = "Measurement", values_to = "Value", names_prefix = "Misc_") %>%
        merge.data.frame(misc.data.names) %>% 
        subset(get(corrftrtypex) == corrftr_x())%>%
        dplyr::rename("x.name" = corrftrtypex, "x.val" = "Value", "x.lab" = "Label")%>%
        dplyr::select(c(x.name, x.val, x.lab, Mouse))
    }
    
    
    if (input$corrytype == "Transcriptomics"){
      corrftrtypey <- "gene"
      corrdatay <- trans.plotdata() %>%
        subset(get(corrftrtypey) == corrftr_y())%>% 
        dplyr::rename("y.name" = corrftrtypey, "y.val" = "CPMlog2", "y.lab" = "Geneid") %>% 
        dplyr::select(c(y.name, y.val, y.lab, Mouse))
    }
    else if (input$corrytype == "Proteomics"){
      corrftrtypey <- "gene"
      corrdatay <- prot.plotdata() %>%
        subset(get(corrftrtypey) == corrftr_y())%>% 
        dplyr::rename("y.name" = corrftrtypey, "y.val" = "Intlog2", "y.lab" = "ProteinID")%>% 
        dplyr::select(c(y.name, y.val, y.lab, Mouse))
    }
    else if (input$corrytype == "Metabolomics"){
      corrftrtypey <- "Metabolite"
      corrdatay <- met.calcdata()$met.data %>%
        subset(get(corrftrtypey) == corrftr_y())%>% 
        dplyr::rename("y.name" = corrftrtypey, "y.val" = "Intlog2", "y.lab" = "Label")%>% 
        dplyr::select(c(y.name, y.val, y.lab, Mouse))
    }
    else if (input$corrytype == "Miscellaneous"){
      corrftrtypey <- "Measurement"
      corrdatay <- misc.data.trans() %>% 
        pivot_longer(cols = Misc_KIF:Misc_Weightloss_abs, names_to = "Measurement", values_to = "Value", names_prefix = "Misc_") %>%
        merge.data.frame(misc.data.names) %>% 
        subset(get(corrftrtypey) == corrftr_y())%>%
        dplyr::rename("y.name" = corrftrtypey, "y.val" = "Value", "y.lab" = "Label")%>% 
        dplyr::select(c(y.name, y.val, y.lab, Mouse))
    }
    
    corrdata <- merge.data.frame(x = metadata, y = corrdatax, by = "Mouse")%>% 
      merge.data.frame(corrdatay, by = "Mouse") %>%
      na.omit() 

    corrresult <- cor.test(x = corrdata$x.val, y = corrdata$y.val)
    
    list(corrdata = corrdata, corrresult = corrresult)
  })
  
  #output$corrplot---- 
  output$corrplot <- renderPlot({
    req(input$main_nav == "corr", 
        nrow(corr()$corrdata) > 1)
    
    corrdata <- corr()$corrdata
    corrresult <- corr()$corrresult
    
    if(corrtype_x() == "Transcriptomics"){xtitle = bquote(italic(.(corrdata[1, "x.lab"])) ~ "(log"[2] ~ "CPM)")}
    else if (corrtype_x() == "Proteomics" | corrtype_x() == "Metabolomics"){xtitle = bquote(.(corrdata[1, "x.lab"])~ "(log"[2]~"Intensity)")}
    else if(corrtype_x() == "Miscellaneous"){xtitle = parse(text = corrdata[1, "x.lab"])}
    else(xtitle = "xtitle")
    
    if(corrtype_y() == "Transcriptomics"){ytitle = bquote(italic(.(corrdata[1, "y.lab"])) ~ "(log"[2] ~ "CPM)")}
    else if (corrtype_y() == "Proteomics" | corrtype_y() == "Metabolomics"){ytitle = bquote(.(corrdata[1, "y.lab"])~ "(log"[2]~"Intensity)")}
    else if(corrtype_y() == "Miscellaneous"){ytitle = parse(text = corrdata[1, "y.lab"])}
    else(ytitle = "ytitle")
    
    
    
    p <- ggplot(data = corrdata, mapping = aes(x = x.val, y = y.val))+
      geom_smooth(formula = y ~ x, method = "lm", colour = "black", fill = "grey50", alpha = 0.4)+ 
      geom_point(mapping = aes(fill = Treatment, shape = Phenotype), size = 2)+
      scale_shape_manual(values = c(22, 23), labels = LabelsP)+
      scale_fill_manual(values = c(coloursT, coloursP)) +
      scale_colour_manual(values = colours, labels = Labels) +
      scale_x_continuous(name = xtitle, expand = expansion(add = 0, mult = 0.2))+
      scale_y_continuous(name = ytitle,  expand = expansion(add = 0, mult = 0.2))+
      theme(panel.background = element_rect(fill="transparent"),
            axis.line = element_line(colour="black", linewidth = 1/1.427),
            plot.title = element_text(hjust=0.5, size=13, face = "bold"),
            plot.subtitle = element_text(hjust=0.5, size=11),
            axis.ticks = element_line(linewidth = 1/1.427, colour = "black"),
            axis.text = element_text(colour="black", size = 11),
            axis.title = element_text(size=12),
            legend.key.height = unit(0.6, "cm"),
            legend.title = element_blank(),
            legend.key = element_blank(),
            legend.text = element_text(size=11, hjust = 0),
            legend.margin = ggplot2::margin(-5.5, 5.5, -5.5, 5.5),
            plot.caption = element_text (colour = "black", size = 10, hjust = 1),
            panel.grid = element_blank())+
      #coord_fixed(clip = "off", ratio = (max(corrdata$x.val) - min(corrdata$x.val))/(max(corrdata$y.val) - min(corrdata$y.val)))+
      labs(title = "Correlation plot", 
           subtitle = bquote(atop("R = "* .(signif(corrresult$estimate, 3)), italic("p")* "-value: " * .(signif(corrresult$p.value, 3))* " ("* .(asterisk(corrresult$p.value))* 
                                    ")")))+
      guides(fill = guide_legend(override.aes = list(shape = 22)))
    
    gb <- ggplot_build(p)
    
    min.y = min(gb$data[[1]]$ymin, gb$data[[2]]$y)
    max.y = max(gb$data[[1]]$ymax, gb$data[[2]]$y)
    
    min.x = min(gb$data[[1]]$x, gb$data[[2]]$x)
    max.x = max(gb$data[[1]]$x, gb$data[[2]]$x)
   
    p+
      coord_fixed(clip = "off", ratio = (max.x - min.x)/(max.y - min.y))
  
  },
  res = 100)
  


  #output$corrplotcont ----
  corrplotmsg <- reactiveVal(TRUE)

  observeEvent(input$corrplotgo, {
    corrplotmsg(FALSE)
  })
  output$corrplotcont <- renderUI({
    if(corrplotmsg() == TRUE){
      tags$div(
        style = "
        height: 100%;
        display: flex;
        align-items: center;
        justify-content: center;
        text-align: center;
        padding: 50px;
      ",
        p(HTML('Select <b>two parameters</b> to plot <b>against each other</b> from the options on the left (fast).<br>
        Once selected, you can calculate the <b>correlations</b> of either feature
        <b>against all other features</b> across data sets by selecting the corresponding <i>Corr. with Feature</i> tab (slow).'), 
          class = "msg")
      )
     
    }
    else{
      tags$div(
        #style = "min-height: 250px;",
        plotOutput(outputId = "corrplot", width = "500px")
      )
     
    }
  })
  
  #corrtable_x eventReactive ---- 
  corrtable_x <- eventReactive(corrftr_x(), {
      req(input$corrtabs == "xcorr", 
          nrow(corr()$corrdata) > 1, 
          nzchar(corrftr_x()))
     
      if(corrtype_x() == "Transcriptomics"){ftr = paste0("Trans_", corrftr_x())}
      else if(corrtype_x() == "Proteomics"){ftr = paste0("Prot_", corrftr_x())}
      else if(corrtype_x() == "Metabolomics"){ftr = paste0("Met_", corrftr_x())}
      else if(corrtype_x() == "Miscellaneous"){ftr = paste0("Misc_", corrftr_x())}
      

    
    mat <- all.data() %>% 
      `rownames<-`(.$Mouse) %>% 
      dplyr::select(-Mouse) %>%
      mutate(across(everything(), as.numeric)) %>% 
      as.matrix() 
    
    x <- mat[, ftr] 
    xrm <- x %>% na.omit() 
    if(length(xrm) < nrow(mat)){
      mat <- mat %>% 
        na.omit()
      x <- mat[, ftr]
    }
    else{mat <- mat %>%
      t() %>% 
      na.omit() %>% 
      t()}
      
      r <- cor(x, mat)
      n <- nrow(mat)
      p <- 2 * pt(-abs(r * sqrt((n - 2) / (1 - r^2))), df = n - 2)
      
      debug_df <- data.frame(
        Ftr = colnames(mat),
        R = as.numeric(r),
        P.Value = as.numeric(p)
      )
      
      data.frame(
        Ftr = colnames(mat), 
        R = as.numeric(r), 
        P.Value = as.numeric(p)
      ) %>% 
        subset(Ftr != ftr) %>% 
        mutate(
          Dataset = unlist(lapply(FUN = function(X){str_split_fixed(string = X, pattern = "_", n = 2) %>% .[,1]}, X = Ftr)), 
          Feature = unlist(lapply(FUN = function(X){str_split_fixed(string = X, pattern = "_", n = 2) %>% .[, 2]}, X = Ftr)),
          Rsquared = R*R, 
          adj.P.Value = p.adjust(P.Value, method = "BH"), 
          sig = unlist(lapply(FUN = asterisk, X = adj.P.Value)), 
     
        ) %>% 
        mutate(Dataset = unlist(lapply(FUN = function(X){if_else(X == "Trans", "Transcriptomics", if_else(X == "Prot", "Proteomics", if_else(X == "Met", "Metabolomics", "Miscellaneous")))}, X = Dataset)),
               Feature = case_when(
                 Dataset == "Transcriptomics" ~ Feature,
                 Dataset == "Proteomics" ~ prot.lookup()[Feature],
                 Dataset == "Metabolomics" ~ met.lookup[Feature],
                 Dataset == "Miscellaneous" ~ misc.lookup[Feature],
                 TRUE ~ Feature
               )
        ) %>%
        dplyr::select(Dataset, Feature, R, Rsquared, P.Value, adj.P.Value, sig)%>% 
        arrange(-Rsquared)
      
      
    }
  )
  
  #output$corrdatax ---- 
  output$corrdatax <- renderDataTable({
    req(input$main_nav == "corr", 
        nrow(corr()$corrdata) > 1, 
        input$corrtabs == "xcorr")
    
    
   results <- corrtable_x()
    
    DT::datatable(
      results,
      filter = "bottom", 
      rownames = FALSE,
      options = list(
        columnDefs = list(list(
          targets = 4:5,
          render = JS(
            "function(data, type, row) {",
            "  return type === 'display' ?",
            "    Number(data).toExponential(3) : data;",
            "}"
          )
        ),
        
        list(
          targets = 2:3,
          render = JS(
            "function(data, type, row) {",
            "  return type === 'display' ?",
            "    Number(data).toPrecision(4) : data;",
            "}"
          )
        ))
      ))
  })
  
  
  #corrtable_y eventReactive ---- 
  corrtable_y <- eventReactive(corrftr_y(), {
    req(input$corrtabs == "ycorr", 
        nrow(corr()$corrdata) > 1, 
        nzchar(corrftr_y()))
    
    if(corrtype_y() == "Transcriptomics"){ftr = paste0("Trans_", corrftr_y())}
    else if(corrtype_y() == "Proteomics"){ftr = paste0("Prot_", corrftr_y())}
    else if(corrtype_y() == "Metabolomics"){ftr = paste0("Met_", corrftr_y())}
    else if(corrtype_y() == "Miscellaneous"){ftr = paste0("Misc_", corrftr_y())}
    
    
    
    mat <- all.data() %>% 
      `rownames<-`(.$Mouse) %>% 
      dplyr::select(-Mouse) %>%
      mutate(across(everything(), as.numeric)) %>% 
      as.matrix() 
    
    x <- mat[, ftr] 
    xrm <- x %>% na.omit() 
    if(length(xrm) < nrow(mat)){
      mat <- mat %>% 
        na.omit()
      x <- mat[, ftr]
    }
    else{mat <- mat %>%
           t() %>% 
           na.omit() %>% 
           t()}
    
    r <- cor(x, mat)
    n <- nrow(mat)
    p <- 2 * pt(-abs(r * sqrt((n - 2) / (1 - r^2))), df = n - 2)
    
    debug_df <- data.frame(
      Ftr = colnames(mat),
      R = as.numeric(r),
      P.Value = as.numeric(p)
    )
    
    data.frame(
      Ftr = colnames(mat), 
      R = as.numeric(r), 
      P.Value = as.numeric(p)
    ) %>% 
      subset(Ftr != ftr) %>% 
      mutate(
        Dataset = unlist(lapply(FUN = function(X){str_split_fixed(string = X, pattern = "_", n = 2) %>% .[,1]}, X = Ftr)), 
        Feature = unlist(lapply(FUN = function(X){str_split_fixed(string = X, pattern = "_", n = 2) %>% .[, 2]}, X = Ftr)),
        Rsquared = R*R, 
        adj.P.Value = p.adjust(P.Value, method = "BH"), 
        sig = unlist(lapply(FUN = asterisk, X = adj.P.Value)), 
        
      ) %>% 
      mutate(Dataset = unlist(lapply(FUN = function(X){if_else(X == "Trans", "Transcriptomics", if_else(X == "Prot", "Proteomics", if_else(X == "Met", "Metabolomics", "Miscellaneous")))}, X = Dataset)),
             Feature = case_when(
               Dataset == "Transcriptomics" ~ Feature,
               Dataset == "Proteomics" ~ prot.lookup()[Feature],
               Dataset == "Metabolomics" ~ met.lookup[Feature],
               Dataset == "Miscellaneous" ~ misc.lookup[Feature],
               TRUE ~ Feature
             )
      ) %>% 
      dplyr::select(Dataset, Feature, R, Rsquared, P.Value, adj.P.Value, sig)%>% 
      arrange(-Rsquared)
    
    
  }
  )
  
  #output$corrdatay ---- 
  output$corrdatay <- renderDataTable({
    req(input$main_nav == "corr", 
        nrow(corr()$corrdata) > 1, 
        input$corrtabs == "ycorr")
    
    
    results <- corrtable_y()
    
    DT::datatable(
      results,
      filter = "bottom", 
      rownames = FALSE,
      options = list(
        columnDefs = list(list(
          targets = 4:5,
          render = JS(
            "function(data, type, row) {",
            "  return type === 'display' ?",
            "    Number(data).toExponential(3) : data;",
            "}"
          )
        ),
        
        list(
          targets = 2:3,
          render = JS(
            "function(data, type, row) {",
            "  return type === 'display' ?",
            "    Number(data).toPrecision(4) : data;",
            "}"
          )
        ))
      ))
  })
  
  
  output$corr_x_title <- renderUI({
    req(nzchar(corrftr_x()))
    tags$h3(paste("Correlation for", corrftr_x()))
  })
  output$corr_y_title <- renderUI({
    req(nzchar(corrftr_y()))
    tags$h3(paste("Correlation for", corrftr_y()))
  })
  ##### EULERS ##### 
  
  #reactive Eulertables ----

  edfs <- reactive({
    req(input$main_nav == "eulers", 
        nzchar(input$eulerdata1), 
        nzchar(input$eulerdata2), 
        nzchar(input$eulercomp1), 
        nzchar(input$eulercomp2), 
        nzchar(input$eulerfc1), 
        nzchar(input$eulerfc2), 
        nzchar(input$eulerp1),
        nzchar(input$eulerp2),
        nzchar(input$eulerdir1), 
        nzchar(input$eulerdir2))
    
    print("entered edfs reactive") 
    if(input$eulerdata1 == "prot") {ftr1 <- "gene"}
    else if(input$eulerdata1 == "trans") {ftr1 <- "Geneid"}
    else if (input$eulerdata1 == "met"){ftr1 <- "Metabolite"}
    if(input$eulerdata2 == "prot") {ftr2 <- "gene"}
    else if(input$eulerdata2 == "trans") {ftr2 <- "Geneid"}
    else if (input$eulerdata2 == "met"){ftr2 <- "Metabolite"}

    met.statdata <- met.calcdata()$stat.data %>%
      dplyr::rename(comp = Comparison) %>%
      mutate(sig = unlist(lapply(X = adj.P.Val, FUN = function(p) {
        if (is.na(p) == TRUE){return(NA)}
        else if (p < 0.0001) {return("****")}
        else if (p < 0.001) {return ("***")}
        else if (p < 0.01) {return ("**")}
        else if (p < 0.05) {return("*")}
        else {return ("ns")}
      })))

    statdata = list(
      trans = trans.statdata(),
      prot = prot.statdata(),
      met = met.statdata
    )



    edf1 <- statdata[[input$eulerdata1]]%>%
      subset(comp == input$eulercomp1) %>%
      subset(adj.P.Val < input$eulerp1) %>%
      filter(if(input$eulerdir1 == "both"){abs(logFC) > log2(input$eulerfc1) }
             else if (input$eulerdir1 == "up"){ logFC > log2(input$eulerfc1)}
             else if(input$eulerdir1 == "down"){-logFC > log2(input$eulerfc1)}
      ) %>%
      dplyr::select(c(ftr1, "logFC", "P.Value", "adj.P.Val", "sig", "Group1", "Group2")) %>%
      dplyr::rename("Feature" = ftr1) %>%
      arrange(adj.P.Val)

    edf2 <-  statdata[[input$eulerdata2]]%>%
      subset(comp == input$eulercomp2) %>%
      subset(adj.P.Val < input$eulerp2) %>%
      filter(if(input$eulerdir2 == "both"){abs(logFC) > log2(input$eulerfc2) }
             else if (input$eulerdir2 == "up"){ logFC > log2(input$eulerfc2)}
             else if(input$eulerdir2 == "down"){-logFC > log2(input$eulerfc2)}
      ) %>%
      dplyr::select(c(ftr2, "logFC", "P.Value", "adj.P.Val", "sig", "Group1", "Group2")) %>%
      dplyr::rename("Feature" = ftr2) %>%
      arrange(adj.P.Val)

    list(
      edf1 = edf1,
      edf2 = edf2
    )
  })

  #output$eulertable1 ----
  output$eulertable1 <- renderDataTable({
    req(input$main_nav == "eulers")

    edfs()$edf1 %>%
      DT::datatable(., rownames = FALSE,
                    options = list(
                      columnDefs = list(list(
                        targets = 2:3,
                        render = JS(
                          "function(data, type, row) {",
                          "  return type === 'display' ?",
                          "    Number(data).toExponential(3) : data;",
                          "}"
                        )
                      ),

                      list(
                        targets = 1,
                        render = JS(
                          "function(data, type, row) {",
                          "  return type === 'display' ?",
                          "    Number(data).toPrecision(4) : data;",
                          "}"
                        )
                      ))
                    )
      )
  })

  #output$eulertable2 ----
  output$eulertable2 <- renderDataTable({
    req(input$main_nav == "eulers")
    if(input$eulerdata2 == "trans" | input$eulerdata2 == "prot") {ftr <- "gene"}
    else if (input$eulerdata2 == "met"){ftr <- "Metabolite"}

    edfs()$edf2 %>%
      DT::datatable(., rownames = FALSE,
                    options = list(
                      columnDefs = list(list(
                        targets = 2:3,
                        render = JS(
                          "function(data, type, row) {",
                          "  return type === 'display' ?",
                          "    Number(data).toExponential(3) : data;",
                          "}"
                        )
                      ),

                      list(
                        targets = 1,
                        render = JS(
                          "function(data, type, row) {",
                          "  return type === 'display' ?",
                          "    Number(data).toPrecision(4) : data;",
                          "}"
                        )
                      ))
                    )
      )


  })

  #output$eulershared ----
  output$eulershared <- renderDataTable({
    req(input$main_nav == "eulers")


    edf <- merge.data.frame(edfs()$edf1, edfs()$edf2, by = "Feature", suffixes = c(".list1", ".list2")) %>%
      dplyr::select(c("Feature", "logFC.list1", "logFC.list2", "adj.P.Val.list1", "adj.P.Val.list2")) %>%
      arrange(adj.P.Val.list1)%>%
      DT::datatable(., rownames = FALSE,
                    options = list(
                      columnDefs = list(list(
                        targets = 3:4,
                        render = JS(
                          "function(data, type, row) {",
                          "  return type === 'display' ?",
                          "    Number(data).toExponential(3) : data;",
                          "}"
                        )
                      ),

                      list(
                        targets = 1:2,
                        render = JS(
                          "function(data, type, row) {",
                          "  return type === 'display' ?",
                          "    Number(data).toPrecision(4) : data;",
                          "}"
                        )
                      ))
                    )
      )

  })

  #output$eulerdiagram ----
  output$eulerdiagram <- renderPlot({
    p <- euler(list(list1 = edfs()$edf1$Feature,
                    list2 = edfs()$edf2$Feature))

    plot(p,
         fills = list(fill= c("list1" = input$eulercol1, "list2" = input$eulercol2), alpha = 0.8),
         labels = list(labels = c(paste("First set\n(n = ", nrow(edfs()$edf1), ")", sep = ""),
                                  paste("Second set\n(n = ", nrow(edfs()$edf2), ")", sep = "")),
                       fontsize = 12),
         quantities = list(cex = TRUE, fontsize = 12))
  })

  #output#eulerunique1
  output$eulerunique1 <- renderDataTable({
    req(input$main_nav == "eulers")

    edf1 <-  edfs()$edf1
    edf2 <- edfs()$edf2
    edf1[edf1$Feature %notin% edf2$Feature, ] %>%
      DT::datatable(., rownames = FALSE,
                    options = list(
                      columnDefs = list(list(
                        targets = 2:3,
                        render = JS(
                          "function(data, type, row) {",
                          "  return type === 'display' ?",
                          "    Number(data).toExponential(3) : data;",
                          "}"
                        )
                      ),

                      list(
                        targets = 1,
                        render = JS(
                          "function(data, type, row) {",
                          "  return type === 'display' ?",
                          "    Number(data).toPrecision(4) : data;",
                          "}"
                        )
                      ))
                    )
      )
  })

  #output#eulerunique2
  output$eulerunique2 <- renderDataTable({
    req(input$main_nav == "eulers")

    edf1 <-  edfs()$edf1
    edf2 <- edfs()$edf2
    edf2[edf2$Feature %notin% edf1$Feature, ] %>%
      DT::datatable(., rownames = FALSE,
                    options = list(
                      columnDefs = list(list(
                        targets = 2:3,
                        render = JS(
                          "function(data, type, row) {",
                          "  return type === 'display' ?",
                          "    Number(data).toExponential(3) : data;",
                          "}"
                        )
                      ),

                      list(
                        targets = 1,
                        render = JS(
                          "function(data, type, row) {",
                          "  return type === 'display' ?",
                          "    Number(data).toPrecision(4) : data;",
                          "}"
                        )
                      ))
                    )
      )
  })

  #output$eulercontent ---- 
 output$eulercontent <- renderUI({
   if(nzchar(input$eulerdata1) & 
      nzchar(input$eulerdata2) &
      nzchar(input$eulercomp1) &
      nzchar(input$eulercomp2) &
      nzchar(input$eulerdir1) &
      nzchar(input$eulerdir2)){
     div(
       #style = "display: grid; grid-template-columns: auto auto; grid-template-rows: auto auto; gap: 10px; width: 100%; ", 
       div(
         style = "width: 100%; ", 
         card(
           class = "card-content", 
           card_header(h3("Comparison results", class = "text-center", style = "margin: 0;")),
           # h4("First list"),
           navset_tab(
             selected = "Shared",
             nav_panel("First set", DTOutput(outputId = "eulertable1") %>% withSpinner(type = 7)),
             nav_panel("Second set", DTOutput(outputId = "eulertable2")%>% withSpinner(type = 7)),
             nav_panel("Shared", DTOutput(outputId = "eulershared")%>% withSpinner(type = 7)),
             nav_panel("Unique to first set", DTOutput(outputId = "eulerunique1")%>% withSpinner(type = 7)),
             nav_panel("Unique to second set", DTOutput(outputId = "eulerunique2")%>% withSpinner(type = 7))
           ) #end of navset_tab
         )  #end of card
       ), #table div 
       div(
         style = "width: 100%; display: grid; grid-template-columns: auto auto; ", 
         card(
           style = "width: fit-content; ", 
           max_height = "300px", 
           class = "card-content-sm", 
           card_header(h3("Euler diagram", class = "text-center", style = "margin: 0;")),
           div(
             style = "display: grid; grid-template-columns: auto auto; gap: 10px; align-items: center;",
             div(
               style = "width: 320px; height: 250px; display: flex; align-items: center; justify-content: center; ", 
               plotOutput(outputId = "eulerdiagram" , width = "300px", height = "230px")%>% withSpinner(type = 7)
             ),
             div(
               style = "max-width: 200px; display: flex; flex-direction: column; align-items: center;  ", 
               h4("Euler settings", style = "margin-bottom: 10px; " ), 
               colourInput(
                 inputId = "eulercol1",
                 label = "Colour for first set",
                 value = "#106107AF",
                 allowTransparent = TRUE,
                 closeOnClick = FALSE,
                 width = "150px"
               ),
               colourInput(
                 inputId = "eulercol2",
                 label = "Colour for second set",
                 value = "#1B7DBF9B",
                 allowTransparent = TRUE,
                 closeOnClick = FALSE,
                 width = "150px"
               ), 
               actionButton(
                 inputId = "reseteulercolours", 
                 label = "Reset colours",
                 class = "btn-outline-secondary btn-sm"
               )
             ) #colour settings div 
           )#whole card div 
         ), # end of card 
         card(
           tags$div(
             style = "
        height: 100%;
        display: flex;
        align-items: center;
        justify-content: center;
        text-align: center;
        padding: 50px;
      ",
             p(HTML('
              The tabs <b><i>First set</i></b> and <b><i>Second set</b></i> show <b>all significant features</b> for the chosen comparisons. 
              The tab <b><i>Shared</i></b> lists features that are <b>found in both sets</b>.
              The <b><i>Unique to...</i></b> tabs contain the features that are <b>not shared</b>.'), 
               class = "msg")
           ), #div
           class = "card-content-sm", 
         )#card
       )#second row div 
       
     )#end of largediv 
   } else{
     card(
       class = "card-content", 
       tags$div(
         style = "
        height: 100%;
        display: flex;
        align-items: center;
        justify-content: center;
        text-align: center;
        padding: 50px;
      ",
         p(HTML('A <b>Euler diagram</b>, similar to a Venn diagram, <b>shows relationships between two groups</b> (e.g., sets of significant features).
              The areas of the circles and their overlaps are relative to the sizes of each feature set and the number of overlapping features, respectively.
              <b>Choose two sets</b> using the setttings on the left and <b>see how they overlap</b>.'), 
           class = "msg")
       )
     )
     
   }
 })

  ##### Fungal OMICS ##### 
  fungal.trans <- reactiveVal(NULL)
  
  observeEvent(input$fungaltranscalc, {
    print("fungaltranscalc clicked")
    x <- fungal.trans.counts()
    id_cols <- c("Geneid", "gene", "product")
    counts <- x[, -c(1:length(id_cols))]
    design <- matrix(c(c(1,0,0,1,0,1,0,0,1,0,1,1),c(0,1,1,0,1,0,1,1,0,1,0,0)), ncol=2, dimnames=list(c('HS17_03','HS17_04','HS17_05','HS17_06','HS17_07','HS17_08','HS17_11','HS17_15','HS17_16','HS17_19','HS17_20','HS17_25'),c('Con_Can','T2D_Can')))
    grouping <- c("HS17_03" = "Con_Can", "HS17_04" = "T2D_Can", "HS17_05" = "T2D_Can", "HS17_06" = "Con_Can", 
                  "HS17_07" = "T2D_Can", "HS17_08" = "Con_Can", "HS17_11" = "T2D_Can", "HS17_15" = "T2D_Can", 
                  "HS17_16" = "Con_Can", "HS17_19" = "T2D_Can", "HS17_20" = "Con_Can", "HS17_25" = "Con_Can")
    canmice <- grouping[grouping %in% c("Con_Can", "T2D_Can")]
    meta <- data.frame(group = c("Con_PBS", "T2D_PBS", "Con_Can", "T2D_Can"),
                       Treatment = c("PBS", "PBS", "Candida", "Candida"), 
                       Phenotype = c("Control", "Diabetic", "Control", "Diabetic"))
    contrasts <- limma::makeContrasts(
      Can = T2D_Can - Con_Can, 
      levels = design
    )
    
    
    keepReads <- apply(counts, 1, max) >= input$fungaltransreads
    keepCpm <- rowSums(cpm(counts)> input$fungaltranscpm) >= input$fungaltranscpmn
    keep <- keepReads & keepCpm
    
    x <- x[keep, ]
    counts <- counts[keep, ]
    
    #normalising for library size 
    nf <- calcNormFactors(counts) #calculating the normalising factors for each sample 
    y <- voom(counts, design, plot=FALSE, lib.size=colSums(counts)*nf) #Transform count data to log2 counts-per-million (logCPM), using the normalisation factor
    
    fungal.trans.plot <- cbind(x[, 1:2], y$E) %>%
      pivot_longer(cols = 3:ncol(.), names_to = "Mouse", values_to = "CPMlog2") %>% 
      mutate(Mouse = gsub(pattern = "HS17_", replacement = "", x = Mouse) %>% as.numeric()) %>%
      merge.data.frame(y = metadata) 
    #doing the actual comparison stats 
    fit <- lmFit(y,design) #Fit linear model for each gene given a series of arrays
    fit2 <- contrasts.fit(fit, contrasts) #Given a linear model fit to microarray data, compute estimated coefficients and standard errors for a given set of contrasts.
    fit2 <- eBayes(fit2) #Given a linear model fit from lmFit, compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression by empirical Bayes moderation of the standard errors towards a global value.
    
    fungal.trans.stat <- cbind(x[, c("Geneid", "gene")], 
                               topTable(fit2, coef = "Can", number = Inf, adjust.method = "BH", sort.by = "none"))
    
  

    fungal.trans(
      list(
        plot = fungal.trans.plot, 
        stat = fungal.trans.stat
    )
    )
  })
  
  observeEvent(input$main_nav, {
    req(input$main_nav == "fungalomics")
      print("running fungal calculation")
      
      x <- fungal.trans.counts()
      id_cols <- c("Geneid", "gene", "product")
      counts <- x[, -c(1:length(id_cols))]
      design <- matrix(c(c(1,0,0,1,0,1,0,0,1,0,1,1),c(0,1,1,0,1,0,1,1,0,1,0,0)), ncol=2, dimnames=list(c('HS17_03','HS17_04','HS17_05','HS17_06','HS17_07','HS17_08','HS17_11','HS17_15','HS17_16','HS17_19','HS17_20','HS17_25'),c('Con_Can','T2D_Can')))
      grouping <- c("HS17_03" = "Con_Can", "HS17_04" = "T2D_Can", "HS17_05" = "T2D_Can", "HS17_06" = "Con_Can", 
                    "HS17_07" = "T2D_Can", "HS17_08" = "Con_Can", "HS17_11" = "T2D_Can", "HS17_15" = "T2D_Can", 
                    "HS17_16" = "Con_Can", "HS17_19" = "T2D_Can", "HS17_20" = "Con_Can", "HS17_25" = "Con_Can")
      canmice <- grouping[grouping %in% c("Con_Can", "T2D_Can")]
      meta <- data.frame(group = c("Con_PBS", "T2D_PBS", "Con_Can", "T2D_Can"),
                         Treatment = c("PBS", "PBS", "Candida", "Candida"), 
                         Phenotype = c("Control", "Diabetic", "Control", "Diabetic"))
      contrasts <- limma::makeContrasts(
        Can = T2D_Can - Con_Can, 
        levels = design
      )
      
      
      keepReads <- apply(counts, 1, max) >= input$fungaltransreads
      keepCpm <- rowSums(cpm(counts)> input$fungaltranscpm) >= input$fungaltranscpmn
      keep <- keepReads & keepCpm
      
      x <- x[keep, ]
      counts <- counts[keep, ]
      
      #normalising for library size 
      nf <- calcNormFactors(counts) #calculating the normalising factors for each sample 
      y <- voom(counts, design, plot=FALSE, lib.size=colSums(counts)*nf) #Transform count data to log2 counts-per-million (logCPM), using the normalisation factor
      
      fungal.trans.plot <- cbind(x[, 1:2], y$E) %>%
        pivot_longer(cols = 3:ncol(.), names_to = "Mouse", values_to = "CPMlog2") %>% 
        mutate(Mouse = gsub(pattern = "HS17_", replacement = "", x = Mouse) %>% as.numeric()) %>%
        merge.data.frame(y = metadata) 
      #doing the actual comparison stats 
      fit <- lmFit(y,design) #Fit linear model for each gene given a series of arrays
      fit2 <- contrasts.fit(fit, contrasts) #Given a linear model fit to microarray data, compute estimated coefficients and standard errors for a given set of contrasts.
      fit2 <- eBayes(fit2) #Given a linear model fit from lmFit, compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression by empirical Bayes moderation of the standard errors towards a global value.
      
      fungal.trans.stat <- cbind(x[, c("Geneid", "gene")], 
                       topTable(fit2, coef = "Can", number = Inf, adjust.method = "BH", sort.by = "none"))
      
      

      fungal.trans(
        list(
        plot = fungal.trans.plot, 
        stat = fungal.trans.stat
      )
      )

      
    })
  #output$fungaltransstattable ----
  output$fungaltransstattable <- renderDataTable({
    req(input$main_nav == "fungalomics")
    
   
    fungal.trans()$stat %>%
      arrange(adj.P.Val) %>%
      DT::datatable(., rownames = FALSE,
                    options = list(
                      columnDefs = list(list(
                        targets = 5:6,
                        render = JS(
                          "function(data, type, row) {",
                          "  return type === 'display' ?",
                          "    Number(data).toExponential(3) : data;",
                          "}"
                        )
                      ),
                      
                      list(
                        targets = c(2, 3, 4, 7),
                        render = JS(
                          "function(data, type, row) {",
                          "  return type === 'display' ?",
                          "    Number(data).toPrecision(4) : data;",
                          "}"
                        )
                      ))
                    )
      )
  })
  
  #output$fungaltransstatdown -----
  output$fungaltransstatdown <- downloadHandler(
    filename = function (){
      paste0(Sys.Date(), "_Stölting et al._Fungal transcriptomics_stats results.csv")
    }, 
    content = function(file){
      write.csv(
        fungal.trans()$stat %>%
          arrange(adj.P.Val) %>% 
          mutate(significant = if_else(adj.P.Val < 0.05, TRUE, FALSE), 
                 Dia_vs_Con = if_else(logFC > 0, "UP", "DOWN")), 
        file, 
        row.names = FALSE 
      )
    }
  )
  #output$fungaltransvolcano ----
  output$fungaltransvolcano <- renderPlotly({
    req(input$main_nav == "fungalomics")
    
    
    volcdf <- fungal.trans()$stat
    
    
    p <- ggplot(data = volcdf, mapping = aes(x = logFC, y = -log10(adj.P.Val), 
                                                key = Geneid, text = paste0(
                                                  "<b>", gene, "</b><br>",
                                                  "(", Geneid, ")<br>",
                                                  "log<sub>2</sub>(FC): ", signif(logFC, 3), "<br>", 
                                                  "FDR: ", signif(adj.P.Val, 3))))+
      geom_hline(yintercept = -log10(input$fungalvolcp), colour= "black", linetype = "dotted")+
      geom_vline(xintercept = c(-log2(input$fungalvolcfc), log2(input$fungalvolcfc)), colour= "black", linetype = "dotted")+
      geom_point(data = volcdf %>% subset(adj.P.Val >= input$fungalvolcp |
                                               abs(logFC) < log2(input$fungalvolcfc)),
                 colour = "grey", alpha = 0.5)+
      geom_point(data = volcdf %>% subset(adj.P.Val < input$fungalvolcp &
                                               logFC > log2(input$fungalvolcfc)),
                 colour = "#862164", alpha = 0.5)+
      geom_point(data = volcdf %>% subset(adj.P.Val < input$fungalvolcp &
                                               logFC < -log2(input$fungalvolcfc)),
                 colour = "#006DAE", alpha = 0.5)+
      scale_x_continuous(name = expression("log"["2"]~"Fold Change (FC)"), limits = c(-5, 5), expand = expansion())+
      scale_y_continuous(name = expression("-log"["10"]~"False Discovery Rate (FDR)") , limits = c(0, 5), expand = expansion())+
      theme(panel.background = element_rect(fill="transparent"),
            axis.line = element_line(colour="black", linewidth = 1/1.427),
            plot.title = element_text(hjust=0.5, size=13),
            plot.subtitle = element_text(hjust=0.5, size=11),
            axis.ticks= element_line(linewidth = 1/1.427, colour = "black"),
            axis.text = element_text(colour="black", size = 11),
            axis.title = element_text(size=12),
            legend.key.height = unit(0.6, "cm"),
            legend.title = element_blank(),
            legend.key = element_blank(),
            legend.text = element_text(size=11, hjust = 0),
            legend.margin = ggplot2::margin(-5.5, 5.5, -5.5, 5.5),
            plot.caption = element_text (colour = "black", size = 10, hjust = 1),
            panel.grid = element_blank())+
      annotate(geom = "text", x = ((-5+log2(input$fungalvolcfc))/2)-log2(input$fungalvolcfc), y = 14.5/3, label = "Control", colour = "#006DAE", hjust = 0.5)+
      annotate(geom = "text", x = ((5-log2(input$fungalvolcfc))/2)+log2(input$fungalvolcfc), y = 14.5/3, label =  "Diabetic", colour = "#862164", hjust = 0.5)+
      annotate(geom = "text", x = ((-5+log2(input$fungalvolcfc))/2)-log2(input$fungalvolcfc), y = .5/3, label = volcdf %>% subset(adj.P.Val < input$fungalvolcp &
                                                                                                                                    -logFC > log2(input$fungalvolcfc)) %>% nrow(),
               colour = "#006DAE", hjust = 0.5, vjust = 1)+
      annotate(geom = "text", x = ((5-log2(input$fungalvolcfc))/2)+log2(input$fungalvolcfc), y = 0.5/3, label = volcdf %>% subset(adj.P.Val < input$fungalvolcp &
                                                                                                                                    logFC > log2(input$fungalvolcfc)) %>% nrow(),
               colour = "#862164", hjust = 0.5, vjust = 1)+
      labs(title = "Fungal transcriptomics")#+
    #coord_fixed(clip = "off", ratio = 1.5)
    
    if (!is.null(selected_fungalfeature())){
      p <- p+
        geom_point (data = volcdf %>% subset(Geneid == selected_fungalfeature()),
                    colour = "black", alpha = 0.8)
    } else {p <- p}

    
    
    ggplotly(p, tooltip = "text", source = "fungaltransvolc") %>%
      layout(
        xaxis = list(title = "log<sub>2</sub> Fold Change (FC)", titlefont = list(size = 14), tickfont = list(size = 10)),
        yaxis = list(title = "-log<sub>10</sub> False Discovery Rate (FDR)", titlefont = list(size = 14), tickfont = list(size = 10))
      ) %>%
      plotly::event_register("plotly_click")
  }
  ) #end of output$metvolcano
  
  #output$fungalintplot ----
  output$fungalintplot <- renderPlot({
    req(input$main_nav == "fungalomics")
    req(nzchar(selected_fungalfeature()))
    # req(selected_mfeature() %in% met.calcdata()$stat.data$Metabolite)
    
    
    # boxplot of log2 intensity ---
    # with requested gene
    # 1. It is "reactive" and therefore should be automatically
    #    re-executed when inputs (input$gene) change
    # 2. Its output type is a plot
    df <- subset(fungal.trans()$plot, Geneid == selected_fungalfeature())
    
  
    name <- if_else(df[1, "gene"] == "NA", df[1, "Geneid"], df[1, "gene"]) 

    
    Max <- max(df[, "CPMlog2"])
    Min <- min(df[, "CPMlog2"])
    # if(Min < 0){
    #   
    #   if(Max < 0){
    #     pbs <- prettybreaks_neg(min = Min, max = 1)
    #     plotrange <- pbs[2]-pbs[1]
    #     plotmin <- pbs[1]
    #     plotmax <- pbs[2]
    #     breaks <- seq(pbs[1], pbs[2], plotrange/pbs[3])
    #   }
    #   else{
    #     pbs <- prettybreaks_neg(min = Min, max = Max)
    #     plotrange <- pbs[2]-pbs[1]
    #     plotmin <- pbs[1]
    #     plotmax <- pbs[2]
    #     breaks <- seq(pbs[1], pbs[2], plotrange/pbs[3])
    #   }
    #   
    # }
    # else{
    pbs <- prettybreaks_log2(Min, Max)
    plotrange <- pbs[2]-pbs[1]
    plotmin <- pbs[1]
    plotmax <- pbs[2]
    breaks <- seq(pbs[1], pbs[2], plotrange/pbs[3])
    #plotrange <- plotmax - plotmin
    # }
    
    set.seed(seed = 793)
    
    
    p <- ggplot(data = df, mapping = aes(x = Group, y = CPMlog2))+
      scale_x_discrete(expand=expansion(mult=0, add=0.67))+
      scale_y_continuous(name = expression("log"["2"]~"Intensity"), limits = c(plotmin, plotmax),
                         expand = expansion(mult = 0, add = 0), breaks = breaks)+
      geom_boxplot(mapping = aes(fill = Phenotype), colour = "black", show.legend = F, fatten = 0.75, outliers = FALSE)+
      geom_jitter(mapping = aes(fill = Treatment, shape = Sex), stroke = 1.5/1.427, size = 2, width = 0.2, height = 0, show.legend = TRUE) +
      scale_fill_manual(values = c(coloursP, coloursT), labels = c(LabelsP, LabelsT)) +
      scale_shape_manual(values = c("Female" = 21, "Male" = 22))+
      theme(panel.background = element_rect(fill="transparent"),
            axis.line = element_line(colour="black", linewidth = 1/1.427),
            plot.title = element_text(hjust=0.5, size=13),
            plot.subtitle = element_text(hjust=0.5, size=11),
            axis.ticks.y = element_line(linewidth = 1/1.427, colour = "black"),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(colour="black", size = 11),
            axis.text.x = element_blank(),
            axis.title.y = element_text(size=12),
            axis.title.x = element_blank(),
            legend.key.height = unit(0.6, "cm"),
            legend.title = element_blank(),
            legend.key = element_blank(),
            legend.text = element_text(size=11, hjust = 0),
            legend.margin = ggplot2::margin(-5.5, 5.5, -5.5, 5.5),
            plot.caption = element_text (colour = "black", size = 10, hjust = 1),
            panel.grid = element_blank())+
      labs(title = name) +
      #coord_cartesian(clip = "off")+
      coord_fixed(clip = "off", ratio = (5.76/plotrange))+
      guides(fill = guide_legend(order = 1,
                                 override.aes = list(shape  = 22,
                                                     size = c(4, 4, 2))),
             shape = guide_legend(order = 2, override.aes = list(fill = "grey")))
    
    # if(input$metshowstats){
    #   
    #   metsdf <- subset(met.calcdata()$stat.data, Metabolite == selected_mfeature()) %>%
    #     mutate(sig = unlist(lapply(X = adj.P.Val, FUN = function(p) {
    #       if (is.na(p) == TRUE){return(NA)}
    #       else if (p < 0.0001) {return("****")}
    #       else if (p < 0.001) {return ("***")}
    #       else if (p < 0.01) {return ("**")}
    #       else if (p < 0.05) {return("*")}
    #       else {return ("ns")}
    #     })), 
    #     barx1 = as.numeric(Group1),
    #     barx2 = as.numeric(Group2),
    #     starx = (barx1 + barx2)/2,
    #     ybarfact = mapply(FUN = function(grp1, grp2){
    #       if ((grp1 == 1 & grp2 == 2) | (grp2 == 1 & grp1 == 2)){return(2)}
    #       else if ((grp1 == 3 & grp2 == 4) | (grp2 == 3 & grp1 == 4)){return(2)}
    #       else if ((grp1 == 1 & grp2 == 3) | (grp2 == 1 & grp1 == 3)){return(1)}
    #       else if ((grp1 == 2 & grp2 == 4) | (grp2 == 2 & grp1 == 4)){return(0)}
    #     }, grp1 = as.numeric(Group1), grp2 = as.numeric(Group2)))
    #   
    #   p <- p+
    #     geom_segment(inherit.aes = FALSE, data = metsdf, mapping = aes(x = barx1, xend = barx2, y = plotrange*(0.98-ybarfact*0.07-0.025)+plotmin, yend = plotrange*(0.98-ybarfact*0.07-0.025)+plotmin))
    #   
    #   if(input$metstatrep == 1){
    #     p <- p +
    #       geom_text(inherit.aes = FALSE, data = metsdf, mapping = aes(x = starx, label = paste("p = ", signif(adj.P.Val, 3)), y = plotrange*(0.98-ybarfact*0.07)+plotmin), size = 3, parse = FALSE)
    #     
    #   }
    #   if(input$metstatrep == 2){
    #     p <- p +
    #       geom_text(inherit.aes = FALSE, data = metsdf, mapping = aes(x = starx, label = sig, y = plotrange*(0.98-ybarfact*0.07)+plotmin))
    #     
    #   }
    #   
    # }
    
    
    
    p
    
  },
  res = 100
  ) #end of output metintplot
  
  #output$fungalintmessage ----
  output$fungalintmessage <- renderUI({
    if (!nzchar(selected_fungalfeature())) {
      tags$div(
        style = "
        height: 100%;
        display: flex;
        align-items: center;
        justify-content: center;
        text-align: center;
        padding: 50px;
      ",
        p(HTML("To display a boxplot of log<sub>2</sub>-transformed normalised intensity by group, <b>search for a 
          feature</b> in the box on the left or <b>select one from the volcano plot</b>"), 
          class = "msg")
      )
      
    }  else{
      tags$div(
        style = "
        height: 100%;
        display: flex;
        align-items: center;
        justify-content: center;
        text-align: center;
        padding: 50px;
      ",
        ""
      )
    }
  }) #end of output$metintmessage
  
  #output$fungalint_content----
  output$fungalint_content <- renderUI({
    
    if (nzchar(selected_fungalfeature())) {
      
      tags$div(
        #style = "min-height: 250px;",
        plotOutput("fungalintplot", height = "350px")
      )
      
    } else {
      
      tags$div(
        style = "
        height: 100%;
        display: flex;
        align-items: center;
        justify-content: center;
        text-align: center;
        padding: 50px;
      ",
        uiOutput("fungalintmessage")
      )
      
    }
  })
  
}# end of server

shinyApp(ui = ui, server = server)


