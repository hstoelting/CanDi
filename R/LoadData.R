library(dplyr)
library(qs2)
### loading vectors and files ####
print("loading vectors etc ")
print(Sys.time())

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

names <- read.csv(file = file.path("data", "HS17_Metabolomics_Metabolite-names.csv"), stringsAsFactors = FALSE)
met.lookup <- setNames(names$Label_Conf, names$Metabolite)

gene_aliases <- qs_read(file = file.path("data", "HS17_omics_gene-aliases.qs2"))
gene_aliases_trans <- qs_read(file = file.path("data", "HS17_omics_gene-aliases_trans.qs2"))
gene_aliases_prot <- qs_read(file = file.path("data", "HS17_omics_gene-aliases_prot.qs2"))
pathlist <- qs_read(file = file.path("data", "HS17_omics_pathlist.qs2"))
pathway.options <- qs_read(file = file.path("data", "HS17_omics_pathway.options.qs2"))
misc.data.names <- read.csv(file = file.path("data", "HS17_misc data_names.csv"), stringsAsFactors = FALSE)
misc.lookup <- setNames(misc.data.names$Name, misc.data.names$Measurement)
