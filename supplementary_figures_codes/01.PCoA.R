library(tidyverse)
library(openxlsx)
library(ggplot2)
library(vegan)
library(dplyr)
library(ggsci)
library(ggpubr)
library(aplot)

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


#GGMP糖尿病患者和健康人的PCoA
    #矩阵转置(行为样本，列为微生物分类)
        data <- as.data.frame(t(as.data.frame(data_ggmp)))
        data <- data[-1,] #去掉转置后第一行的各种微生物分类名
        data[1:ncol(data)] <- lapply(data[1:ncol(data)],as.numeric) #将reads数变为数值型
    #计算Bray-Curtis距离
        data <- vegdist(data, method = 'bray')
    #计算pcoa(计算的距离矩阵缺少一行，所以要减去样本名和额外一行)
        pcoa <- cmdscale(data, k = ncol(data_ggmp)-2, eig = TRUE)

    #提取样本点坐标（points记录了各样本在各排序轴中的坐标值）
        #前两轴
        plot_data <- data.frame({pcoa$points})[1:2]
        names(plot_data)[1:2] <- c('PCoA1', 'PCoA2') #命名为PCoA1和PCoA2
    #eig记录了PCoA排序结果中，主要排序轴的特征值（再除以特征值总和就是各轴的解释量）
        eig = pcoa$eig
    #为样本点坐标添加分组信息
        plot_data$group <- metadata_ggmp$healthy
        plot_data$group <- gsub("TRUE","Healthy",plot_data$group)
        plot_data$group <- gsub("FALSE","T2D",plot_data$group)

        p <- ggplot(plot_data,aes(x=PCoA1,y=PCoA2,color=group)) + 
                geom_point()+
                #添加置信椭圆
                stat_ellipse(level = 0.95, show.legend = F)+
                theme_bw()+
                theme(panel.border=element_blank(),
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    axis.title = element_blank(),
                    axis.text.x = element_text(size = 12,colour = "black"),
                    axis.text.y = element_text(size = 12,colour = "black"),                   
                    axis.title.y = element_text(size = 15,colour = "black"),
                    axis.title.x = element_text(size = 15,colour = "black"),
                    legend.text = element_text(size = 15,colour = "black"),
                    legend.title = element_text(size = 17,colour = "black"),
                    axis.line= element_line(colour = "black"),
                    plot.title = element_text(hjust = 0.5,size = 20,colour = "black"))+
                    labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
                        y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""))
        ggsave("/home/wangnan/Shandong-DF/figure/PCoA(Guangdong).pdf",p,height = 5,width = 6)

#SGMP糖尿病患者和健康人的PCoA
    #矩阵转置(行为样本，列为微生物分类)
        data <- as.data.frame(t(as.data.frame(data_sgmp)))
        data <- data[-1,] #去掉转置后第一行的各种微生物分类名
        data[1:ncol(data)] <- lapply(data[1:ncol(data)],as.numeric) #将reads数变为数值型
    #计算Bray-Curtis距离
        data <- vegdist(data, method = 'bray')
    #计算pcoa(计算的距离矩阵缺少一行，所以要减去样本名和额外一行)
        pcoa <- cmdscale(data, k = ncol(data_sgmp)-2, eig = TRUE)

    #提取样本点坐标（points记录了各样本在各排序轴中的坐标值）
        #前两轴
        plot_data <- data.frame({pcoa$points})[1:2]
        names(plot_data)[1:2] <- c('PCoA1', 'PCoA2') #命名为PCoA1和PCoA2
    #eig记录了PCoA排序结果中，主要排序轴的特征值（再除以特征值总和就是各轴的解释量）
        eig = pcoa$eig
    #为样本点坐标添加分组信息
        plot_data$group <- metadata_sgmp$host_status
        plot_data$group <- gsub("Health","Healthy",plot_data$group)
        plot_data$group <- gsub("Type 2 diabetes","T2D",plot_data$group)

        q <- ggplot(plot_data,aes(x=PCoA1,y=PCoA2,color=group)) + 
                geom_point()+
                #添加置信椭圆
                stat_ellipse(level = 0.95, show.legend = F)+
                theme_bw()+
                theme(panel.border=element_blank(),
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    axis.title = element_blank(),
                    axis.text.x = element_text(size = 12,colour = "black"),
                    axis.text.y = element_text(size = 12,colour = "black"),                   
                    axis.title.y = element_text(size = 15,colour = "black"),
                    axis.title.x = element_text(size = 15,colour = "black"),
                    legend.text = element_text(size = 15,colour = "black"),
                    legend.title = element_text(size = 17,colour = "black"),
                    axis.line= element_line(colour = "black"),
                    plot.title = element_text(hjust = 0.5,size = 20,colour = "black"))+
                    labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
                        y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""))
        ggsave("/home/wangnan/Shandong-DF/figure/PCoA(Shandong).pdf",q,height = 5,width = 6)

#GGMP和SGMP的样本的PCoA
    #数据合并
        data <- full_join(data_ggmp,data_sgmp,by = "Sample ID")
        data[is.na(data)] <- 0
    #矩阵转置(行为样本，列为微生物分类)
        data <- as.data.frame(t(as.data.frame(data)))
        data <- data[-1,] #去掉转置后第一行的各种微生物分类名
        data[1:ncol(data)] <- lapply(data[1:ncol(data)],as.numeric) #将reads数变为数值型
    #计算Bray-Curtis距离
        data <- vegdist(data, method = 'bray')
    #计算pcoa(计算的距离矩阵缺少一行，所以要减去样本名和额外一行)
        pcoa <- cmdscale(data, k = ncol(data_ggmp)+ncol(data_sgmp)-3, eig = TRUE)

    #提取样本点坐标（points记录了各样本在各排序轴中的坐标值）
        #前两轴
        plot_data <- data.frame({pcoa$points})[1:2]
        names(plot_data)[1:2] <- c('PCoA1', 'PCoA2') #命名为PCoA1和PCoA2
    #eig记录了PCoA排序结果中，主要排序轴的特征值（再除以特征值总和就是各轴的解释量）
        eig = pcoa$eig
    #为样本点坐标添加分组信息
        plot_data$group <- "1"
        plot_data$group[1:nrow(metadata_ggmp)] <- "Guangdong"
        plot_data$group[(nrow(metadata_ggmp)+1):nrow(plot_data)] <- "Shandong"
        plot_data$group <- factor(plot_data$group,levels = c("Guangdong","Shandong"))

        plot <- ggplot(plot_data,aes(x=PCoA1,y=PCoA2,color=group)) + 
                geom_point(size = 0.5)+
                #添加置信椭圆
                stat_ellipse(level = 0.95, show.legend = F)+
                theme_bw()+
                theme(
                    #panel.border=element_blank(),
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    axis.title = element_blank(),
                    axis.text.x = element_text(size = 12,colour = "black"),
                    axis.text.y = element_text(size = 12,colour = "black"),                   
                    axis.title.y = element_text(size = 15,colour = "black"),
                    axis.title.x = element_text(size = 15,colour = "black"),
                    legend.text = element_text(size = 15,colour = "black"),
                    legend.position = "none",
                    legend.title = element_blank(),
                    axis.line= element_line(colour = "black"),
                    plot.title = element_text(hjust = 0.5,size = 20,colour = "black"))+
                    labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
                        y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
                xlim(-0.6,0.6) +
                ylim(-0.45,0.6) +
                scale_color_npg()
        
        #PCoA 1的boxplot
            comparision <- list(c("Guangdong","Shandong")) # 这里一定要用list的数据形式，不然会报错
            b1 <- ggplot(plot_data,aes(x = group,y = PCoA1,fill=group)) +
                    geom_point(size = 0.3)+
                    geom_boxplot()+
                    geom_signif(comparisons = comparision,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                    scale_fill_npg()+
                    coord_flip()+ # 图形转90度
                    theme_bw()+
                    theme(
                        #panel.border=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        axis.title = element_blank(),
                        axis.text.x = element_blank(),
                        axis.text.y = element_text(size = 12,colour = "black"),                  
                        axis.title.y = element_blank(),
                        axis.title.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        legend.position = "none",
                        strip.text.x = element_text(size = 15,colour = "black"), #修改分面图标题的大小
                        #axis.line= element_line(colour = "black"),
                        plot.title = element_text(hjust = 0.5,size = 20,colour = "black")) +
                    ylim(-0.6,0.6)

        #PCoA 2的boxplot
            comparision <- list(c("Guangdong","Shandong")) # 这里一定要用list的数据形式，不然会报错
            b2 <- ggplot(plot_data,aes(x = group,y = PCoA2,fill=group)) +
                    geom_point(size = 0.3)+
                    geom_boxplot()+
                    geom_signif(comparisons = comparision,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                    scale_fill_npg()+
                    theme_bw()+
                    theme(
                        #panel.border=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        axis.title = element_blank(),
                        axis.text.x = element_blank(),
                        axis.text.y = element_blank(),                   
                        axis.title.y = element_blank(),
                        axis.title.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.ticks.y = element_blank(),
                        legend.position = "none",
                        strip.text.x = element_text(size = 15,colour = "black")) + #修改分面图标题的大小
                        ylim(-0.45,0.6)
                        
        
        # PCoA图和箱线图组合
            p <- plot %>% 
                insert_top(b1, height = 0.2) %>% 
                insert_right(b2, width = 0.4)

        ggsave("/home/wangnan/Shandong-DF/figure/PCoA_with_boxplot(Guangdong-Shandong).pdf",p,height = 5,width = 6)

