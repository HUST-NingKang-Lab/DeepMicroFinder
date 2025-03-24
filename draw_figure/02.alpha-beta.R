library(tidyverse)
library(openxlsx)
library(ggplot2)
library(vegan)
library(ape)
library(ggpubr)
library(ggsci)

#读取文件
    #读取GGMP的糖尿病患者的metadata
        metadata_ggmp <- read.csv("/data4/wangnan/Shandong-DF/origin-data/ggmp/filter_metadata.csv")
    #读取GGMP的糖尿病患者的metadata
        data_ggmp <- read_tsv("/data4/wangnan/Shandong-DF/origin-data/ggmp/filter_data.tsv")
    #去掉Archaea的分类
        data_ggmp <- data_ggmp[4:nrow(data_ggmp),]

    #读取SGMP的糖尿病患者的metadata
        metadata_sgmp <- read.csv("/data4/wangnan/Shandong-DF/origin-data/sgmp/filter_metadata.csv")
    #读取SGMP的糖尿病患者的metadata
        data_sgmp <- read_tsv("/data4/wangnan/Shandong-DF/origin-data/sgmp/filter_data.tsv")

#计算α多样性
    #数据合并
        data <- full_join(data_ggmp,data_sgmp,by = "Sample ID")
        data[is.na(data)] <- 0
    #矩阵转置(行为样本，列为微生物分类)
        data <- as.data.frame(t(as.data.frame(data)))
        data <- data[-1,] #去掉转置后第一行的各种微生物分类名
        data[1:ncol(data)] <- lapply(data[1:ncol(data)],as.numeric) #将reads数变为数值型
    #Shannon指数
        Shannon <- diversity(data,index = "shannon",MARGIN = 1,base = exp(1))
    #Simpson 指数
        Simpson <- diversity(data, index = "simpson", MARGIN = 1, base = exp(1))
    #合并
        index <- as.data.frame(cbind(Shannon,Simpson))

    #添加分组信息
        #省份信息
            index$province <- "1"
            index$province[1:nrow(metadata_ggmp)] <- "Guangdong"
            index$province[(nrow(metadata_ggmp)+1):nrow(index)] <- "Shandong"
        #疾病状态信息
            index$group <- "1"
            index$group[1:nrow(metadata_ggmp)] <- metadata_ggmp$healthy
            index$group[(nrow(metadata_ggmp)+1):nrow(index)] <- metadata_sgmp$host_status
            index$group <- gsub("Health","Healthy",index$group)
            index$group <- gsub("TRUE","Healthy",index$group)
            index$group <- gsub("FALSE","T2D",index$group)
            index$group <- gsub("Type 2 diabetes","T2D",index$group)
    
    #画图
        comparision <- list(c("Healthy","T2D"))
        p <- ggplot(index,aes(x=group,y=Shannon,fill = group)) + 
                geom_boxplot()+
                theme_bw()+
                facet_grid(~ province,scales = "free")+ #箱线图分面
                geom_signif(comparisons = comparision,step_increase = 0.1,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                theme(
                    #panel.border=element_blank(),
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    axis.title = element_blank(),
                    axis.text.x = element_text(size = 12,colour = "black"),
                    axis.text.y = element_text(size = 12,colour = "black"),                   
                    axis.title.y = element_text(size = 15,colour = "black"),
                    axis.title.x = element_blank(),
                    legend.position = "none",
                    strip.text.x = element_text(size = 15,colour = "black"), #修改分面图标题的大小
                    #axis.line= element_line(colour = "black"),
                    plot.title = element_text(hjust = 0.5,size = 20,colour = "black"))+
                    labs(y = "Shannon diversity") +
                scale_fill_npg()
        ggsave("/home/wangnan/Shandong-DF/figure/alpha-diversity.pdf",p,height = 4.5,width = 6)


