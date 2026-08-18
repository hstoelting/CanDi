#### loading vectors and files #### 
print("loading vectors etc ")
print(Sys.time())
rm(list = ls())
# mainDir <- "C:/Users/hsto0009/OneDrive - Monash University/Experiments/HS17 - 1x10^5 C. albicans in BKS.db for omics/R_HS17"
# #mainDir <- "C:/Users/helen/OneDrive - Monash University/Experiments/HS17 - 1x10^5 C. albicans in BKS.db for omics/R_HS17"


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
