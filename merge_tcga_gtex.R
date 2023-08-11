

############################### TCGA + GTEx ####################################

# 加载要用到的包，大家如果没有可以先自己安装一下

library(data.table)  # 用于高效处理大数据集的库
library(dplyr)       # 数据操作和转换的库
library(tidyverse)   # 数据处理和可视化的综合库



#################################### TCGA ######################################

# 读取表达矩阵
dlbc.fpkm <- fread("./download_data/TCGA-DLBC.htseq_fpkm.tsv", header = T, sep = '\t', data.table = F)

# 查看表达矩阵，我们发现第一列是ensembl id，后面的列是样本名
head(dlbc.fpkm)[1:5, 1:5]
#           Ensembl_ID TCGA-RQ-A6JB-01A TCGA-FF-8046-01A TCGA-FF-A7CW-01A TCGA-RQ-A68N-01A
# 1  ENSG00000242268.2       0.00000000         0.000000       0.00000000         0.000000
# 2  ENSG00000270112.3       0.04396508         0.000000       0.01294356         0.000000
# 3 ENSG00000167578.15       3.15192692         3.399737       3.43036482         3.332282
# 4  ENSG00000273842.1       0.00000000         0.000000       0.00000000         0.000000
# 5  ENSG00000078237.5       3.03661720         3.183184       3.04236324         2.787429


# 读取基因ID转换信息，目的是为了将ensembl id转换为gene symbol
dlbc.pro <- fread("./download_data/gencode.v22.annotation.gene.probeMap", header = T, sep = '\t', data.table = F)

# 查看基因ID转换信息，我们发现第一列是ensembl id，第二列是gene symbol
head(dlbc.pro)
#                  id         gene chrom chromStart chromEnd strand
# 1 ENSG00000223972.5      DDX11L1  chr1      11869    14409      +
# 2 ENSG00000227232.5       WASH7P  chr1      14404    29570      -
# 3 ENSG00000278267.1    MIR6859-3  chr1      17369    17436      -
# 4 ENSG00000243485.3 RP11-34P13.3  chr1      29554    31109      +
# 5 ENSG00000274890.1    MIR1302-9  chr1      30366    30503      +
# 6 ENSG00000237613.2      FAM138A  chr1      34554    36081      -


# 提取前两列用于进行转换
dlbc.pro <- dlbc.pro[ , c(1, 2)]
head(dlbc.pro)
#                  id         gene
# 1 ENSG00000223972.5      DDX11L1
# 2 ENSG00000227232.5       WASH7P
# 3 ENSG00000278267.1    MIR6859-3
# 4 ENSG00000243485.3 RP11-34P13.3
# 5 ENSG00000274890.1    MIR1302-9
# 6 ENSG00000237613.2      FAM138A

# 基因ID转换信息和表达矩阵合并
dlbc.fpkm.pro <- merge(dlbc.pro, dlbc.fpkm, by.y  = "Ensembl_ID", by.x = "id" )
head(dlbc.fpkm.pro)[1:5, 1:5]
#                   id     gene TCGA-RQ-A6JB-01A TCGA-FF-8046-01A TCGA-FF-A7CW-01A
# 1 ENSG00000000003.13   TSPAN6         1.026771       0.76890645       0.35573821
# 2  ENSG00000000005.5     TNMD         0.000000       0.07154948       0.02051499
# 3 ENSG00000000419.11     DPM1         5.134474       5.01494933       5.28879599
# 4 ENSG00000000457.12    SCYL3         1.534159       1.18431712       1.91343387
# 5 ENSG00000000460.15 C1orf112         1.894523       1.72454346       2.24504933


# 去重
dlbc.fpkm.pro <- distinct(dlbc.fpkm.pro,gene, .keep_all = T)

# 把基因名转换为行名
rownames(dlbc.fpkm.pro) <- dlbc.fpkm.pro$gene
dlbc.fpkm.pro <- dlbc.fpkm.pro[ , -c(1,2)]
dim(dlbc.fpkm.pro) # 58387  48
head(dlbc.fpkm.pro)[1:5, 1:5]
#          TCGA-RQ-A6JB-01A TCGA-FF-8046-01A TCGA-FF-A7CW-01A TCGA-RQ-A68N-01A TCGA-FM-8000-01A
# TSPAN6           1.026771       0.76890645       0.35573821        0.7864696        1.0575194
# TNMD             0.000000       0.07154948       0.02051499        0.0000000        0.1362624
# DPM1             5.134474       5.01494933       5.28879599        5.1281748        4.6654726
# SCYL3            1.534159       1.18431712       1.91343387        1.2004038        1.5156969
# C1orf112         1.894523       1.72454346       2.24504933        1.2145814        1.7478664

# 现在表达矩阵就构建完毕啦！
# 但是TCGA中不止有癌组织，还有癌旁组织，所以我们可以根据临床信息把它们分开

# 读取淋巴癌的临床信息
dlbc.phe <- fread("./download_data/TCGA-DLBC.GDC_phenotype.tsv", header = T, sep = '\t', data.table = F)
dlbc.phe$submitter_id.samples[1:5]
# [1] "TCGA-FA-A6HN-01A" "TCGA-GR-A4D4-01A" "TCGA-GS-A9TT-01A" "TCGA-FF-A7CQ-01A" "TCGA-FF-8046-01A"

# 把行名改成样本名
rownames(dlbc.phe) <- dlbc.phe$submitter_id.samples

# 查看临床信息的sample_type.samples列
table(dlbc.phe$sample_type.samples)
# Bone Marrow Normal      Primary Tumor 
#                  4                 48

# Primary Tumor为癌组织，我们将癌组织提取出来
dlbc.phe.t <- filter(dlbc.phe, sample_type.samples == "Primary Tumor")

# 临床信息与表达矩阵取交集
merge_phe_fpkm <- intersect(rownames(dlbc.phe.t), colnames(dlbc.fpkm.pro))

# 提取癌组织的表达矩阵
dlbc.exp <- dlbc.fpkm.pro[ , merge_phe_fpkm]
dim(dlbc.exp)
# [1] 58387    48

# TCGA DLBC 癌组织样本表达矩阵构建完成，保存一下！
saveRDS(dlbc.exp, file='./generated_data/dlbc.exp.rds')



#################################### GTEx ######################################

# 读取gtex的表达矩阵，解压和不解压其实都可以读取（如果超级慢或者出现error，不要慌，据说可能因为内存小）
gtex.exp <- fread("./download_data/gtex_RSEM_gene_fpkm", header = T, sep = '\t', data.table = F)
gtex.exp[1:5, 1:4]
#               sample GTEX-S4Q7-0003-SM-3NM8M GTEX-QV31-1626-SM-2S1QC GTEX-13QIC-0011-R1a-SM-5O9CJ
# 1  ENSG00000242268.2                 -3.8160                 -9.9658                      -9.9658
# 2  ENSG00000259041.1                 -9.9658                 -9.9658                      -9.9658
# 3  ENSG00000270112.3                 -4.0350                 -2.9324                      -2.1779
# 4 ENSG00000167578.16                  4.2958                  3.9204                       6.1528
# 5  ENSG00000278814.1                 -9.9658                 -9.9658                      -9.9658

dim(gtex.exp)
# [1] 60498  7863

# 注意这里包含GTEx中所有正常组织来源的样本，我们保存一下，后面继续提取自己需要的组织样本
# saveRDS(gtex.exp, file = './generated_data/gtex.exp.rds')


# 读取gtex的基因ID转换信息，将ensembl id转换为gene symbol
gtex.pro <- fread("./download_data/probeMap_gencode.v23.annotation.gene.probemap", header = T, sep = '\t', data.table = F)
head(gtex.pro)
#                  id         gene chrom chromStart chromEnd strand
# 1 ENSG00000223972.5      DDX11L1  chr1      11869    14409      +
# 2 ENSG00000227232.5       WASH7P  chr1      14404    29570      -
# 3 ENSG00000278267.1    MIR6859-1  chr1      17369    17436      -
# 4 ENSG00000243485.3 RP11-34P13.3  chr1      29554    31109      +
# 5 ENSG00000274890.1    MIR1302-2  chr1      30366    30503      +
# 6 ENSG00000237613.2      FAM138A  chr1      34554    36081      -

dim(dlbc.pro) # 60483     2
dim(gtex.pro) # 60498     6
# 我们顺便发现了v23和v22的差异有多少，60498 vs. 60483

# 提取前两列用于基因ID转换
gtex.pro <- gtex.pro[, c(1,2)]

# 基因ID转换信息和表达矩阵合并
gtex.fpkm.pro <- merge(gtex.pro, gtex.exp, by.y ="sample", by.x = "id" )

# 取交集，为了决定之后gtex和dlbc合并是按照gene symbol还是ensembl id进行
length(intersect(gtex.pro$id, dlbc.pro$id))  # 42566
length(intersect(rownames(dlbc.exp), gtex.fpkm.pro$gene)) # 57793

# 按照gene symbol合并会有57793个，按照ensembl id合并有42566个，所以这次我们选择按照gene symbol进行合并



# 接下来，提取正常淋巴组织的表达矩阵，根据gtex的临床信息来匹配淋巴组织的样本

# 读取gtex的临床信息
gtex.phe <- fread("./download_data/GTEX_phenotype", header = T, sep = '\t', data.table = F)
rownames(gtex.phe) <- gtex.phe$Sample

# 给它改一下列名，好看点
colnames(gtex.phe) <- c("Sample", "body_site_detail (SMTSD)", "primary_site", "gender", "patient", "cohort")
table(gtex.phe$primary_site)
# <not provided>  Adipose Tissue   Adrenal Gland         Bladder           Blood    Blood Vessel     Bone Marrow 
#              5             621             161              13             595             753             102 
#          Brain          Breast    Cervix Uteri           Colon       Esophagus  Fallopian Tube           Heart 
#           1426             221              11             384             805               7             493 
#         Kidney           Liver            Lung          Muscle           Nerve           Ovary        Pancreas 
#             38             141             381             478             335             112             203 
#      Pituitary        Prostate  Salivary Gland            Skin Small Intestine          Spleen         Stomach 
#            126             122              71             977             106             121             209 
#         Testis         Thyroid          Uterus          Vagina 
#            208             366              93              99 

# 上面展示的是GTEx中所有的组织来源，今天我们分析的是淋巴癌
# TCGA中数据集名为DLBC，在GTEx中对应的组织来源应该为blood，上面显示有595个

# 在本文最后，我附上了TCGA 与 GTEx 癌症类型与组织来源对应表，大家有需要可以去看哟！

# 筛选出blood组织来源的样本
gtex.phe.s <- filter(gtex.phe, primary_site == "Blood")

# 临床信息与表达矩阵取交集
merge_phe_fpkm_gtex <- intersect(rownames(gtex.phe.s), colnames(gtex.fpkm.pro)) # 444
gtex.s <- gtex.fpkm.pro[ , c("gene", merge_phe_fpkm_gtex)]

# 去重
gtex.s <- distinct(gtex.s, gene, .keep_all = T)
rownames(gtex.s) <- gtex.s$gene
gtex.s <- gtex.s[ , -1]
dim(gtex.s)
# [1] 58581   444

# 我们发现GTEx中有444个blood样本

# gtex的表达矩阵是按照log2(fpkm+0.001)处理的，dlbc是按照log2(fpkm+1)
# 所以合并之前，我们需要把它们的处理方式调整为相同的
gtex.s2 <- 2^gtex.s
gtex.s3 <- log2(gtex.s2-0.001+1)

# 现在数据的处理方式都相同，就有可比性啦！

# 正式合并开始！
# 把gtex的blood组织数据与tcga的dlbc数据合并
all.data <- merge(gtex.s3, dlbc.exp, by = 0)
all.data <- column_to_rownames(all.data, "Row.names")
dim(all.data)
# [1] 57793   492

head(all.data)[1:5, 1:4]
#           GTEX-111YS-0006-SM-5NQBE GTEX-1122O-0003-SM-5Q5DL GTEX-113IC-0006-SM-5NQ9C GTEX-113JC-0006-SM-5O997
# 5_8S_rRNA            -1.571525e-08            -1.571525e-08            -1.571525e-08            -1.571525e-08
# 5S_rRNA              -1.571525e-08            -1.571525e-08            -1.571525e-08            -1.571525e-08
# 7SK                  -1.571525e-08            -1.571525e-08            -1.571525e-08            -1.571525e-08
# A1BG                  2.641511e+00             3.407365e+00             3.005374e+00             1.545996e+00
# A1BG-AS1              1.580150e+00             2.618198e+00             2.657654e+00             6.690156e-01

# 前面444列是来自GTEx的正常blood组织样本，后面48列是来自TCGA的DLBC癌组织样本

# 浅浅保存一下
saveRDS(all.data, file = './generated_data/all_data_tcga_gtex.rds')

############################### 整合完毕 #######################################

# 去除批次效应
library(limma)
nromalized.data <- normalizeBetweenArrays(all.data)
nromalized.data <- as.data.frame(nromalized.data)

# 接下来，就可以进行我们的后续分析啦！



