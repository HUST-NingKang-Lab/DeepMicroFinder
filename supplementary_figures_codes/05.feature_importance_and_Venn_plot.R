library(tidyverse)
library(ggplot2)
library(VennDiagram)

rm(list=ls())
gc()


#按照mean of absolute value对这四个模型筛选出来的feature进行排序
    #读取没有丰度的marker
        no_abundance_marker <- read.csv("/home/wangnan/Shandong-DF/table/no_abundance_marker.csv")
    #读取四个模型选择的marker
        sta_ind <- read.csv("/home/wangnan/Shandong-DF/table/markers selected by independent DNN model.csv")
            pos <- match(no_abundance_marker$no_abundance_marker,sta_ind$Marker)
            sta_ind <- sta_ind[setdiff(1:nrow(sta_ind),pos),]
            sta_ind <- sta_ind[order(sta_ind$Mean_of_absolute_value,decreasing = T),]
        sta_reg <- read.csv("/home/wangnan/Shandong-DF/table/markers selected by regional DNN model.csv")
            pos <- match(no_abundance_marker$no_abundance_marker,sta_reg$Marker)
            sta_reg <- sta_reg[setdiff(1:nrow(sta_reg),pos),]
            sta_reg <- sta_reg[order(sta_reg$Mean_of_absolute_value,decreasing = T),]
        sta_reg_extra <- read.csv("/home/wangnan/Shandong-DF/table/markers selected by regional+ DNN model.csv")
            pos <- match(no_abundance_marker$no_abundance_marker,sta_reg_extra$Marker)
            sta_reg_extra <- sta_reg_extra[setdiff(1:nrow(sta_reg_extra),pos),]
            sta_reg_extra <- sta_reg_extra[order(sta_reg_extra$Mean_of_absolute_value,decreasing = T),]
        sta_tra <- read.csv("/home/wangnan/Shandong-DF/table/markers selected by transfer DNN model.csv")
            pos <- match(no_abundance_marker$no_abundance_marker,sta_tra$Marker)
            sta_tra <- sta_tra[setdiff(1:nrow(sta_tra),pos),]
            sta_tra <- sta_tra[order(sta_tra$Mean_of_absolute_value,decreasing = T),]

    data_ggmp <- read_tsv("/data4/wangnan/Shandong-DF/origin-data/ggmp/table-filtered-feature-rarefied5k-L6.tsv")

    #替换feature
        for(i in 1:nrow(sta_ind))
        {
            strain_g <- strsplit(sta_ind$Marker[i],";",fixed = T)[[1]][6]
                strain_g <- strsplit(strain_g,"__",fixed =T)[[1]][2]
            sta_ind$Marker[i] <- strain_g
        }  
        for(i in 1:nrow(sta_reg))
        {
            strain_g <- strsplit(sta_reg$Marker[i],";",fixed = T)[[1]][6]
                strain_g <- strsplit(strain_g,"__",fixed =T)[[1]][2]
            sta_reg$Marker[i] <- strain_g
        }               
        for(i in 1:nrow(sta_reg_extra))
        {
            strain_g <- strsplit(sta_reg_extra$Marker[i],";",fixed = T)[[1]][6]
                strain_g <- strsplit(strain_g,"__",fixed =T)[[1]][2]
            sta_reg_extra$Marker[i] <- strain_g
        }  
        for(i in 1:nrow(sta_tra))
        {
            strain_g <- strsplit(sta_tra$Marker[i],";",fixed = T)[[1]][6]
                strain_g <- strsplit(strain_g,"__",fixed =T)[[1]][2]
            sta_tra$Marker[i] <- strain_g
        }  
    #画出来重要性的柱状图
        sta_ind <- sta_ind[1:20,]
        p <- ggplot(sta_ind, mapping = aes(x = reorder(Marker,Mean_of_absolute_value), y =Mean_of_absolute_value)) +
            geom_bar(stat = "identity",aes(fill = Mean_of_absolute_value))+
            scale_fill_gradient(low = "#66CCFF", high = "#003366",limits = c(0, 0.16),breaks = c(0.05,0.1,0.15)) + 
            theme_classic()+
            coord_flip()+
            theme(
                axis.text.x = element_text(size = 12,colour = "black"),
                axis.title.x = element_text(size = 12,colour = "black"),
                axis.title.y = element_text(size = 12,colour = "black"),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 12,colour = "black"),
                legend.text = element_blank(),
                legend.title = element_blank(),
                legend.position = "none",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                plot.title = element_text(hjust = 0.5,size = 15,colour = "black")) +
                #scale_fill_npg()+
                #geom_signif(comparisons = comparision,step_increase = 0.1,map_signif_level = F,test = wilcox.test)+
                #stat_compare_means(label = "p.format")+
                #facet_grid(~method,scales = "free")+
                xlab("Taxonomy in Genus level")+
                ylab("Mean absolute value of the change in AUROC")
         ggsave("/home/wangnan/Shandong-DF/figure/feature_importance/importance_independent.pdf",p,height =5,width = 8)

        sta_reg <- sta_reg[1:20,]
        p <- ggplot(sta_reg, mapping = aes(x = reorder(Marker,Mean_of_absolute_value), y =Mean_of_absolute_value)) +
            geom_bar(stat = "identity",aes(fill = Mean_of_absolute_value))+
            scale_fill_gradient(low = "#66CCFF", high = "#003366",limits = c(0, 0.20),breaks = c(0.05,0.1,0.15)) + 
            theme_classic()+
            coord_flip()+
            theme(
                axis.text.x = element_text(size = 12,colour = "black"),
                axis.title.x = element_text(size = 12,colour = "black"),
                axis.title.y = element_text(size = 12,colour = "black"),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 12,colour = "black"),
                legend.text = element_blank(),
                legend.title = element_blank(),
                legend.position = "none",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                plot.title = element_text(hjust = 0.5,size = 15,colour = "black")) +
                #scale_fill_npg()+
                #geom_signif(comparisons = comparision,step_increase = 0.1,map_signif_level = F,test = wilcox.test)+
                #stat_compare_means(label = "p.format")+
                #facet_grid(~method,scales = "free")+
                xlab("Taxonomy in Genus level")+
                ylab("Mean absolute value of the change in AUROC")
         ggsave("/home/wangnan/Shandong-DF/figure/feature_importance/importance_regional.pdf",p,height =5,width = 8)

        sta_reg_extra <- sta_reg_extra[1:20,]
        p <- ggplot(sta_reg_extra, mapping = aes(x = reorder(Marker,Mean_of_absolute_value), y =Mean_of_absolute_value)) +
            geom_bar(stat = "identity",aes(fill = Mean_of_absolute_value))+
            scale_fill_gradient(low = "#66CCFF", high = "#003366",limits = c(0, 0.25),breaks = c(0.05,0.1,0.15,0.20)) + 
            theme_classic()+
            coord_flip()+
            theme(
                axis.text.x = element_text(size = 12,colour = "black"),
                axis.title.x = element_text(size = 12,colour = "black"),
                axis.title.y = element_text(size = 12,colour = "black"),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 12,colour = "black"),
                legend.text = element_blank(),
                legend.title = element_blank(),
                legend.position = "none",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                plot.title = element_text(hjust = 0.5,size = 15,colour = "black")) +
                #scale_fill_npg()+
                #geom_signif(comparisons = comparision,step_increase = 0.1,map_signif_level = F,test = wilcox.test)+
                #stat_compare_means(label = "p.format")+
                #facet_grid(~method,scales = "free")+
                xlab("Taxonomy in Genus level")+
                ylab("Mean absolute value of the change in AUROC")
         ggsave("/home/wangnan/Shandong-DF/figure/feature_importance/importance_regional_extra.pdf",p,height =5,width = 8)

        sta_tra <- sta_tra[1:20,]
        p <- ggplot(sta_tra, mapping = aes(x = reorder(Marker,Mean_of_absolute_value), y =Mean_of_absolute_value)) +
            geom_bar(stat = "identity",aes(fill = Mean_of_absolute_value))+
            scale_fill_gradient(low = "#66CCFF", high = "#003366",limits = c(0, 0.2),breaks = c(0.05,0.1,0.15)) + 
            theme_classic()+
            coord_flip()+
            theme(
                axis.text.x = element_text(size = 12,colour = "black"),
                axis.title.x = element_text(size = 12,colour = "black"),
                axis.title.y = element_text(size = 12,colour = "black"),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 12,colour = "black"),
                legend.text = element_blank(),
                legend.title = element_blank(),
                legend.position = "none",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                plot.title = element_text(hjust = 0.5,size = 15,colour = "black")) +
                #scale_fill_npg()+
                #geom_signif(comparisons = comparision,step_increase = 0.1,map_signif_level = F,test = wilcox.test)+
                #stat_compare_means(label = "p.format")+
                #facet_grid(~method,scales = "free")+
                xlab("Taxonomy in Genus level")+
                ylab("Mean absolute value of the change in AUROC")
         ggsave("/home/wangnan/Shandong-DF/figure/feature_importance/importance_transfer.pdf",p,height =5,width = 8)


#画出来GGMP,SGMP,SD-DFI的韦恩图
#最终只筛选出来95个在三个队列中都存在，且有丰度的marker(之前是38个,现在已经更正)
    data_ggmp <- read_tsv("/data4/wangnan/Shandong-DF/origin-data/ggmp/table-filtered-feature-rarefied5k-L6.tsv")
    data_sgmp <- read_tsv("/data4/wangnan/Shandong-DF/origin-data/sgmp/table-filtered-feature-rarefied5k-L6.tsv")
    data_dfi <- read_tsv("/data4/wangnan/Shandong-DF/origin-data/sd-DFI/sd-dfi_abundance_processed.tsv")
        #metadata
            metadata <- read.csv("/data4/wangnan/Shandong-DF/origin-data/sd-DFI/DFI_metadata_250_used.csv")
            data_dfi <- data_dfi[match(metadata$SampleID,data_dfi$`Sample ID`),] #筛选出来参与膳食纤维干预的250个样本
            #筛选出来genus的行
                pos <- grep("g__",colnames(data_dfi))
                data_dfi <- data_dfi[,c(1,pos)] 
            #去掉genus中species的行    
                pos <- grep("s__",colnames(data_dfi))
                data_dfi <- data_dfi[,c(setdiff(1:ncol(data_dfi),pos))] 
            #筛选出来有丰度的marker
                no_abundance_marker <- read.csv("/home/wangnan/Shandong-DF/table/no_abundance_marker.csv")
                data_dfi <- data_dfi[,setdiff(1:ncol(data_dfi),match(no_abundance_marker[,2],colnames(data_dfi)))] 

    # 三个队列都存在，并且在Shandong-DFI中有丰度的marker
        tmp <- intersect(colnames(data_dfi)[2:ncol(data_dfi)],colnames(data_sgmp)[2:ncol(data_sgmp)])
        tmp <- intersect(tmp,colnames(data_ggmp)[2:ncol(data_ggmp)])

    #生成一个空列表
        x <- list()
        #导入各组的微生物分类
            x[["Guangdong"]] <- colnames(data_ggmp)[2:ncol(data_ggmp)]
            x[["Shandong"]] <- colnames(data_sgmp)[2:ncol(data_sgmp)]
            x[["Dietary fiber intervention"]] <- colnames(data_dfi)[2:ncol(data_dfi)]
        #设置颜色
            color <- c("#fb0007","#139177","#ed9e08")
            venn.plot <- venn.diagram(x,
                filename =NULL,  #保存路径
                height = 450, 
                width = 450,
                resolution =300, #imagetype="png", 
                col = "transparent",      #指定图形的圆周边缘颜色  transparent 透明           
                fill = color,  #填充颜色
                alpha = 0.50,  #透明度
                label.col = "black",
                cex = 1.2,    #每个区域label名称的大小
                fontfamily = "sans",  #字体
                fontface = "bold",     #字体格式
                cat.col = c("#fb0007","#139177","#ed9e08"),  #分类颜色 
                cat.cex = 1.2,      #每个分类名称大小
                cat.pos = c(100, 260, 0),        #
                cat.dist = c(0.14, 0.12, 0.02),    #
                cat.fontfamily = "sans",     #分类字体
                rotation.degree =180,        #旋转角度
                margin = 0.2               #在网格单元中给出图周围空白量的编号
                )
            ggsave("/home/wangnan/Shandong-DF/figure/Venn_plot.pdf",venn.plot,height = 5,width = 5)
