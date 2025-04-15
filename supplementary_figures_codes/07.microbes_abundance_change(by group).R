library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(gg.gap)
library(openxlsx)

rm(list=ls())
gc()
 
#筛选GGMP,SGMP和SD-DFI都存在的marker
    #GGMP的genus水平的丰度表
        data_ggmp <- read_tsv("/data4/wangnan/Shandong-DF/origin-data/ggmp/table-filtered-feature-rarefied5k-L6.tsv")
    #SGMP的genus水平的丰度表
        data_sgmp <- read_tsv("/data4/wangnan/Shandong-DF/origin-data/sgmp/table-filtered-feature-rarefied5k-L6.tsv")
    #SD-DFI的丰度表
        data_dfi <- read_tsv("/data4/wangnan/Shandong-DF/origin-data/sd-DFI/sd-dfi_abundance_processed.tsv")
        #筛选出来有丰度的marker
            no_abundance_marker <- read.csv("/home/wangnan/Shandong-DF/table/no_abundance_marker.csv")
            data_dfi <- data_dfi[,setdiff(1:ncol(data_dfi),match(no_abundance_marker[,2],colnames(data_dfi)))] 
    #筛选出来三个队列都存在的marker
        tmp <- intersect(colnames(data_dfi)[2:ncol(data_dfi)],colnames(data_sgmp)[2:ncol(data_sgmp)])
        tmp <- intersect(tmp,colnames(data_ggmp)[2:ncol(data_ggmp)])
        abundance <- data_dfi[,c(1,match(tmp,colnames(data_dfi)))]

#筛选所使用到的250个样本的丰度
    metadata <- read.csv("/data4/wangnan/Shandong-DF/origin-data/sd-DFI/DFI_metadata_250_used.csv")
    id <- intersect(abundance$`Sample ID`,metadata$SampleID)
    pos <- match(id,as.data.frame(abundance$`Sample ID`)[,1])
    abundance <- abundance[pos,]
    #丰度表行列转置
        abundance <- tibble::rownames_to_column(abundance)
        abundance <- as.data.frame(t(abundance))
        abundance <- tibble::rownames_to_column(abundance)
        abundance <- abundance[-1,]
        colnames(abundance) <- as.vector(abundance[1,])
        abundance <- as_tibble(abundance[-1,])
        abundance[,2:ncol(abundance)] <- lapply(abundance[,2:ncol(abundance)],as.numeric)
        #备份
            # abu <- abundance <- abu
            # meta <- metadata <- meta

#95个菌种在膳食纤维干预前后变化的箱线图
    #同一个病人同一个状态具有多次测序，因此需要去冗余
        # tmp <- abundance[,1]
        # pos <- c()
        # for(i in unique(metadata$Patient_id))
        # {
        #     for(j in unique(metadata$Intervention))
        #     {
        #         p <- which(metadata$Patient_id == i)
        #         p <- p[which(metadata[p,]$Intervention == j)]
        #         if(length(p) > 0)
        #         {
        #             tmp1 <- tibble(apply(abundance[,p+1],1,sum))
        #             tmp1 <- tibble(tmp1/length(p))
        #             colnames(tmp1) <- metadata$SampleID[p[1]]
        #             pos <- c(pos,p[1])
        #             tmp <- bind_cols(tmp,tmp1)
        #         }
        #     }  
        # }
        # abundance <- tmp
        # metadata <- metadata[pos,]
        #match(colnames(abundance)[2:ncol(abundance)],metadata$SampleID)
    #使丰度表转化为相对丰度表
        for(i in 2:ncol(abundance))
        {                
            sum <- sum(abundance[,i]) 
            for(j in 1:nrow(abundance))
            {     
                if(sum > 0)
                {                
                    abundance[j,i] <- abundance[j,i]/sum
                }
            } 
        }
    # write.csv(abundance,"/home/wangnan/Shandong-DF/abundance.csv")
    # write.csv(metadata,"/home/wangnan/Shandong-DF/metadata.csv")

# 按照transfer model选出的marker的重要性按由大到小进行排列
    marker_contribution_transfer <- read.csv("/home/wangnan/Shandong-DF/table/markers selected by transfer DNN model.csv")
        pos <- match(abundance$`Sample ID`,marker_contribution_transfer$Marker)
        marker_contribution_transfer <- marker_contribution_transfer[pos,]
    marker_contribution_transfer <- marker_contribution_transfer[order(marker_contribution_transfer$Mean_of_absolute_value,decreasing = T),]
    pos <- match(marker_contribution_transfer$Marker,abundance$`Sample ID`)
        abundance <- abundance[pos,] # 按照重要性由大到小进行排序 

# 按照GBTM模型的分类结果，将FBG和HbA1c响应者作为响应组，分别画出响应组和非响应组的marker丰度的差异
    group <- read.xlsx("/data4/wangnan/Shandong-DF/origin-data/sd-DFI/Combine_groups.xlsx")
    # 响应组
        id <- group$Patient_id[which(group$FBG == 1| group$HbA1c == 1)]
        metadata_response <- metadata[which(!is.na(match(metadata$Patient_id,id))),]
        abundance_response <- abundance[,c(1,match(metadata_response$SampleID,colnames(abundance)))]
    # 非响应组
        id_non <- setdiff(group$Patient_id,id)
        metadata_non_response <- metadata[which(!is.na(match(metadata$Patient_id,id_non))),]
        abundance_non_response <- abundance[,c(1,match(metadata_non_response$SampleID,colnames(abundance)))]

    # 响应组画图
        # 利用相对丰度表画图(pdf)
            # 画出前20个marker在膳食纤维干预的三个时间点的相对丰度变化
                sta  <- data.frame(sampleid = "1",strain = "1",abundance = 1,status = "1",stringsAsFactors = FALSE)
                for(i in 1:20)
                {
                    for(j in 2:ncol(abundance_response))
                    {
                        strain_g <- strsplit(abundance_response$`Sample ID`[i],";",fixed = T)[[1]][6]
                            strain_g <- strsplit(strain_g,"__",fixed =T)[[1]][2]
                        sta <- rbind(sta,c(metadata_response$Patient_id[j-1],strain_g,round(as.numeric(abundance_response[i,j])*100,4),metadata_response$Intervention[j-1]))
                    }
                    sta <- sta[-1,]
                    sta[,3] <- unlist(lapply(sta[,3],as.numeric))
                }
                sta$status <- factor(sta$status,levels = c("Early","Mid","Later")) # 添加因子水平
                sta$strain <- factor(sta$strain,levels = unique(sta$strain)) # 添加因子水平

                comparision <- list(c("Early","Mid"),c("Mid","Later"),c("Early","Later")) # 设置显著性检验的组合
                p <- ggplot(sta,aes(x=status,y=abundance)) + 
                    geom_point() +
                    geom_boxplot(aes(middle = mean(abundance)),position=position_dodge(0.7),width = 0.6, alpha = 0.1) +
                    # geom_line(aes(group=sampleid), color="black",linewidth = 0.1,alpha = 0.1)+ # alpha可以修改线的透明度
                    theme_bw() +
                    facet_wrap(~strain,scales = "free_y",ncol = 7,nrow = 3) +
                    #facet_grid(~ province,scales = "free")+ #箱线图分面
                    geom_signif(comparisons = comparision,step_increase = 0.08,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                    theme(strip.text.x = element_text(size = 16, face = "bold")) +  # 设置分面的字字体大小、颜色、背景、边框，)
                    theme(
                        #panel.border=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        axis.title = element_text(size = 10,colour = "black"),
                        axis.text.x = element_text(size = 15,colour = "black"),
                        axis.text.y = element_text(size = 15,colour = "black"),                   
                        axis.title.y = element_text(size = 20,colour = "black"),
                        axis.title.x = element_text(size = 20,colour = "black")) +
                        labs(x = "Intervention stage",y = "Relative abundance")
                ggsave(paste0("/home/wangnan/Shandong-DF/figure/GBTM_group/abundance_change/abundance_change_response(FBG+HbA1c)(前20个菌种).pdf"),p,height = 10.5,width = 21)

            # 画出后75个marker在膳食纤维干预的三个时间点的相对丰度变化
                sta  <- data.frame(sampleid = "1",strain = "1",abundance = 1,status = "1",stringsAsFactors = FALSE)
                for(i in 21:95)
                {
                    for(j in 2:ncol(abundance_response))
                    {
                        strain_g <- strsplit(abundance_response$`Sample ID`[i],";",fixed = T)[[1]][6]
                            strain_g <- strsplit(strain_g,"__",fixed =T)[[1]][2]
                        sta <- rbind(sta,c(metadata_response$Patient_id[j-1],strain_g,round(as.numeric(abundance_response[i,j])*100,4),metadata_response$Intervention[j-1]))
                    }
                    sta <- sta[-1,]
                    sta[,3] <- unlist(lapply(sta[,3],as.numeric))
                }
                sta$status <- factor(sta$status,levels = c("Early","Mid","Later")) # 添加因子水平
                sta$strain <- factor(sta$strain,levels = unique(sta$strain)) # 添加因子水平

                comparision <- list(c("Early","Mid"),c("Mid","Later"),c("Early","Later")) # 设置显著性检验的组合
                p <- ggplot(sta,aes(x=status,y=abundance)) + 
                    geom_point() +
                    geom_boxplot(position=position_dodge(0.7),width = 0.6, alpha = 0.1) +
                    # geom_line(aes(group=sampleid), color="black",linewidth = 0.1,alpha = 0.1)+ # alpha可以修改线的透明度
                    theme_bw() +
                    facet_wrap(~strain,scales = "free_y",ncol = 7,nrow = 11) +
                    #facet_grid(~ province,scales = "free")+ #箱线图分面
                    geom_signif(comparisons = comparision,step_increase = 0.08,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                    theme(strip.text.x = element_text(size = 16, face = "bold")) +  # 设置分面的字字体大小、颜色、背景、边框，)
                    theme(
                        #panel.border=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        axis.title = element_text(size = 10,colour = "black"),
                        axis.text.x = element_text(size = 15,colour = "black"),
                        axis.text.y = element_text(size = 15,colour = "black"),                   
                        axis.title.y = element_text(size = 20,colour = "black"),
                        axis.title.x = element_text(size = 20,colour = "black")) +
                        labs(x = "Intervention stage",y = "Relative abundance")
                ggsave(paste0("/home/wangnan/Shandong-DF/figure/GBTM_group/abundance_change/abundance_change_response(FBG+HbA1c)(后75个菌种).pdf"),p,height = 38.5,width = 21)

    # 非响应组画图
        # 利用相对丰度表画图(pdf)
            # 画出前20个marker在膳食纤维干预的三个时间点的相对丰度变化
                sta  <- data.frame(sampleid = "1",strain = "1",abundance = 1,status = "1",stringsAsFactors = FALSE)
                for(i in 1:20)
                {
                    for(j in 2:ncol(abundance_non_response))
                    {
                        strain_g <- strsplit(abundance_non_response$`Sample ID`[i],";",fixed = T)[[1]][6]
                            strain_g <- strsplit(strain_g,"__",fixed =T)[[1]][2]
                        sta <- rbind(sta,c(metadata_non_response$Patient_id[j-1],strain_g,round(as.numeric(abundance_non_response[i,j])*100,4),metadata_non_response$Intervention[j-1]))
                    }
                    sta <- sta[-1,]
                    sta[,3] <- unlist(lapply(sta[,3],as.numeric))
                }
                sta$status <- factor(sta$status,levels = c("Early","Mid","Later")) # 添加因子水平
                sta$strain <- factor(sta$strain,levels = unique(sta$strain)) # 添加因子水平

                comparision <- list(c("Early","Mid"),c("Mid","Later"),c("Early","Later")) # 设置显著性检验的组合
                p <- ggplot(sta,aes(x=status,y=abundance)) + 
                    geom_point() +
                    geom_boxplot(aes(middle = mean(abundance)),position=position_dodge(0.7),width = 0.6, alpha = 0.1) +
                    # geom_line(aes(group=sampleid), color="black",linewidth = 0.1,alpha = 0.1)+ # alpha可以修改线的透明度
                    theme_bw() +
                    facet_wrap(~strain,scales = "free_y",ncol = 7,nrow = 3) +
                    #facet_grid(~ province,scales = "free")+ #箱线图分面
                    geom_signif(comparisons = comparision,step_increase = 0.08,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                    theme(strip.text.x = element_text(size = 16, face = "bold")) +  # 设置分面的字字体大小、颜色、背景、边框，)
                    theme(
                        #panel.border=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        axis.title = element_text(size = 10,colour = "black"),
                        axis.text.x = element_text(size = 15,colour = "black"),
                        axis.text.y = element_text(size = 15,colour = "black"),                   
                        axis.title.y = element_text(size = 20,colour = "black"),
                        axis.title.x = element_text(size = 20,colour = "black")) +
                        labs(x = "Intervention stage",y = "Relative abundance")
                ggsave(paste0("/home/wangnan/Shandong-DF/figure/GBTM_group/abundance_change/abundance_change_non_response(FBG+HbA1c)(前20个菌种).pdf"),p,height = 10.5,width = 21)

            # 画出后75个marker在膳食纤维干预的三个时间点的相对丰度变化
                sta  <- data.frame(sampleid = "1",strain = "1",abundance = 1,status = "1",stringsAsFactors = FALSE)
                for(i in 21:95)
                {
                    for(j in 2:ncol(abundance_non_response))
                    {
                        strain_g <- strsplit(abundance_non_response$`Sample ID`[i],";",fixed = T)[[1]][6]
                            strain_g <- strsplit(strain_g,"__",fixed =T)[[1]][2]
                        sta <- rbind(sta,c(metadata_non_response$Patient_id[j-1],strain_g,round(as.numeric(abundance_non_response[i,j])*100,4),metadata_non_response$Intervention[j-1]))
                    }
                    sta <- sta[-1,]
                    sta[,3] <- unlist(lapply(sta[,3],as.numeric))
                }
                sta$status <- factor(sta$status,levels = c("Early","Mid","Later")) # 添加因子水平
                sta$strain <- factor(sta$strain,levels = unique(sta$strain)) # 添加因子水平

                comparision <- list(c("Early","Mid"),c("Mid","Later"),c("Early","Later")) # 设置显著性检验的组合
                p <- ggplot(sta,aes(x=status,y=abundance)) + 
                    geom_point() +
                    geom_boxplot(position=position_dodge(0.7),width = 0.6, alpha = 0.1) +
                    # geom_line(aes(group=sampleid), color="black",linewidth = 0.1,alpha = 0.1)+ # alpha可以修改线的透明度
                    theme_bw() +
                    facet_wrap(~strain,scales = "free_y",ncol = 7,nrow = 11) +
                    #facet_grid(~ province,scales = "free")+ #箱线图分面
                    geom_signif(comparisons = comparision,step_increase = 0.08,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                    theme(strip.text.x = element_text(size = 16, face = "bold")) +  # 设置分面的字字体大小、颜色、背景、边框，)
                    theme(
                        #panel.border=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        axis.title = element_text(size = 10,colour = "black"),
                        axis.text.x = element_text(size = 15,colour = "black"),
                        axis.text.y = element_text(size = 15,colour = "black"),                   
                        axis.title.y = element_text(size = 20,colour = "black"),
                        axis.title.x = element_text(size = 20,colour = "black")) +
                        labs(x = "Intervention stage",y = "Relative abundance")
                ggsave(paste0("/home/wangnan/Shandong-DF/figure/GBTM_group/abundance_change/abundance_change_non_response(FBG+HbA1c)(后75个菌种).pdf"),p,height = 38.5,width = 21)


# 按照GBTM模型的分类结果，按照union的分组结果，分别画出响应组和非响应组的marker丰度的差异
    group <- read.xlsx("/data4/wangnan/Shandong-DF/origin-data/sd-DFI/Combine_groups.xlsx")
    # 响应组
        id <- group$Patient_id[which(group$union == 1)]
        metadata_response <- metadata[which(!is.na(match(metadata$Patient_id,id))),]
        abundance_response <- abundance[,c(1,match(metadata_response$SampleID,colnames(abundance)))]
    # 非响应组
        id_non <- setdiff(group$Patient_id,id)
        metadata_non_response <- metadata[which(!is.na(match(metadata$Patient_id,id_non))),]
        abundance_non_response <- abundance[,c(1,match(metadata_non_response$SampleID,colnames(abundance)))]

    # 响应组画图
        # 利用相对丰度表画图(pdf)
            # 画出前20个marker在膳食纤维干预的三个时间点的相对丰度变化
                sta  <- data.frame(sampleid = "1",strain = "1",abundance = 1,status = "1",stringsAsFactors = FALSE)
                for(i in 1:20)
                {
                    for(j in 2:ncol(abundance_response))
                    {
                        strain_g <- strsplit(abundance_response$`Sample ID`[i],";",fixed = T)[[1]][6]
                            strain_g <- strsplit(strain_g,"__",fixed =T)[[1]][2]
                        sta <- rbind(sta,c(metadata_response$Patient_id[j-1],strain_g,round(as.numeric(abundance_response[i,j])*100,4),metadata_response$Intervention[j-1]))
                    }
                    sta <- sta[-1,]
                    sta[,3] <- unlist(lapply(sta[,3],as.numeric))
                }
                sta$status <- factor(sta$status,levels = c("Early","Mid","Later")) # 添加因子水平
                sta$strain <- factor(sta$strain,levels = unique(sta$strain)) # 添加因子水平

                comparision <- list(c("Early","Mid"),c("Mid","Later"),c("Early","Later")) # 设置显著性检验的组合
                p <- ggplot(sta,aes(x=status,y=abundance)) + 
                    geom_point() +
                    geom_boxplot(aes(middle = mean(abundance)),position=position_dodge(0.7),width = 0.6, alpha = 0.1) +
                    # geom_line(aes(group=sampleid), color="black",linewidth = 0.1,alpha = 0.1)+ # alpha可以修改线的透明度
                    theme_bw() +
                    facet_wrap(~strain,scales = "free_y",ncol = 7,nrow = 3) +
                    #facet_grid(~ province,scales = "free")+ #箱线图分面
                    geom_signif(comparisons = comparision,step_increase = 0.08,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                    theme(strip.text.x = element_text(size = 16, face = "bold")) +  # 设置分面的字字体大小、颜色、背景、边框，)
                    theme(
                        #panel.border=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        axis.title = element_text(size = 10,colour = "black"),
                        axis.text.x = element_text(size = 15,colour = "black"),
                        axis.text.y = element_text(size = 15,colour = "black"),                   
                        axis.title.y = element_text(size = 20,colour = "black"),
                        axis.title.x = element_text(size = 20,colour = "black")) +
                        labs(x = "Intervention stage",y = "Relative abundance")
                ggsave(paste0("/home/wangnan/Shandong-DF/figure/GBTM_group/abundance_change/abundance_change_response(union)(前20个菌种).pdf"),p,height = 10.5,width = 21)

            # 画出后75个marker在膳食纤维干预的三个时间点的相对丰度变化
                sta  <- data.frame(sampleid = "1",strain = "1",abundance = 1,status = "1",stringsAsFactors = FALSE)
                for(i in 21:95)
                {
                    for(j in 2:ncol(abundance_response))
                    {
                        strain_g <- strsplit(abundance_response$`Sample ID`[i],";",fixed = T)[[1]][6]
                            strain_g <- strsplit(strain_g,"__",fixed =T)[[1]][2]
                        sta <- rbind(sta,c(metadata_response$Patient_id[j-1],strain_g,round(as.numeric(abundance_response[i,j])*100,4),metadata_response$Intervention[j-1]))
                    }
                    sta <- sta[-1,]
                    sta[,3] <- unlist(lapply(sta[,3],as.numeric))
                }
                sta$status <- factor(sta$status,levels = c("Early","Mid","Later")) # 添加因子水平
                sta$strain <- factor(sta$strain,levels = unique(sta$strain)) # 添加因子水平

                comparision <- list(c("Early","Mid"),c("Mid","Later"),c("Early","Later")) # 设置显著性检验的组合
                p <- ggplot(sta,aes(x=status,y=abundance)) + 
                    geom_point() +
                    geom_boxplot(position=position_dodge(0.7),width = 0.6, alpha = 0.1) +
                    # geom_line(aes(group=sampleid), color="black",linewidth = 0.1,alpha = 0.1)+ # alpha可以修改线的透明度
                    theme_bw() +
                    facet_wrap(~strain,scales = "free_y",ncol = 7,nrow = 11) +
                    #facet_grid(~ province,scales = "free")+ #箱线图分面
                    geom_signif(comparisons = comparision,step_increase = 0.08,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                    theme(strip.text.x = element_text(size = 16, face = "bold")) +  # 设置分面的字字体大小、颜色、背景、边框，)
                    theme(
                        #panel.border=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        axis.title = element_text(size = 10,colour = "black"),
                        axis.text.x = element_text(size = 15,colour = "black"),
                        axis.text.y = element_text(size = 15,colour = "black"),                   
                        axis.title.y = element_text(size = 20,colour = "black"),
                        axis.title.x = element_text(size = 20,colour = "black")) +
                        labs(x = "Intervention stage",y = "Relative abundance")
                ggsave(paste0("/home/wangnan/Shandong-DF/figure/GBTM_group/abundance_change/abundance_change_response(union)(后75个菌种).pdf"),p,height = 38.5,width = 21)

    # 非响应组画图
        # 利用相对丰度表画图(pdf)
            # 画出前20个marker在膳食纤维干预的三个时间点的相对丰度变化
                sta  <- data.frame(sampleid = "1",strain = "1",abundance = 1,status = "1",stringsAsFactors = FALSE)
                for(i in 1:20)
                {
                    for(j in 2:ncol(abundance_non_response))
                    {
                        strain_g <- strsplit(abundance_non_response$`Sample ID`[i],";",fixed = T)[[1]][6]
                            strain_g <- strsplit(strain_g,"__",fixed =T)[[1]][2]
                        sta <- rbind(sta,c(metadata_non_response$Patient_id[j-1],strain_g,round(as.numeric(abundance_non_response[i,j])*100,4),metadata_non_response$Intervention[j-1]))
                    }
                    sta <- sta[-1,]
                    sta[,3] <- unlist(lapply(sta[,3],as.numeric))
                }
                sta$status <- factor(sta$status,levels = c("Early","Mid","Later")) # 添加因子水平
                sta$strain <- factor(sta$strain,levels = unique(sta$strain)) # 添加因子水平

                comparision <- list(c("Early","Mid"),c("Mid","Later"),c("Early","Later")) # 设置显著性检验的组合
                p <- ggplot(sta,aes(x=status,y=abundance)) + 
                    geom_point() +
                    geom_boxplot(aes(middle = mean(abundance)),position=position_dodge(0.7),width = 0.6, alpha = 0.1) +
                    # geom_line(aes(group=sampleid), color="black",linewidth = 0.1,alpha = 0.1)+ # alpha可以修改线的透明度
                    theme_bw() +
                    facet_wrap(~strain,scales = "free_y",ncol = 7,nrow = 3) +
                    #facet_grid(~ province,scales = "free")+ #箱线图分面
                    geom_signif(comparisons = comparision,step_increase = 0.08,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                    theme(strip.text.x = element_text(size = 16, face = "bold")) +  # 设置分面的字字体大小、颜色、背景、边框，)
                    theme(
                        #panel.border=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        axis.title = element_text(size = 10,colour = "black"),
                        axis.text.x = element_text(size = 15,colour = "black"),
                        axis.text.y = element_text(size = 15,colour = "black"),                   
                        axis.title.y = element_text(size = 20,colour = "black"),
                        axis.title.x = element_text(size = 20,colour = "black")) +
                        labs(x = "Intervention stage",y = "Relative abundance")
                ggsave(paste0("/home/wangnan/Shandong-DF/figure/GBTM_group/abundance_change/abundance_change_non_response(union)(前20个菌种).pdf"),p,height = 10.5,width = 21)

            # 画出后75个marker在膳食纤维干预的三个时间点的相对丰度变化
                sta  <- data.frame(sampleid = "1",strain = "1",abundance = 1,status = "1",stringsAsFactors = FALSE)
                for(i in 21:95)
                {
                    for(j in 2:ncol(abundance_non_response))
                    {
                        strain_g <- strsplit(abundance_non_response$`Sample ID`[i],";",fixed = T)[[1]][6]
                            strain_g <- strsplit(strain_g,"__",fixed =T)[[1]][2]
                        sta <- rbind(sta,c(metadata_non_response$Patient_id[j-1],strain_g,round(as.numeric(abundance_non_response[i,j])*100,4),metadata_non_response$Intervention[j-1]))
                    }
                    sta <- sta[-1,]
                    sta[,3] <- unlist(lapply(sta[,3],as.numeric))
                }
                sta$status <- factor(sta$status,levels = c("Early","Mid","Later")) # 添加因子水平
                sta$strain <- factor(sta$strain,levels = unique(sta$strain)) # 添加因子水平

                comparision <- list(c("Early","Mid"),c("Mid","Later"),c("Early","Later")) # 设置显著性检验的组合
                p <- ggplot(sta,aes(x=status,y=abundance)) + 
                    geom_point() +
                    geom_boxplot(position=position_dodge(0.7),width = 0.6, alpha = 0.1) +
                    # geom_line(aes(group=sampleid), color="black",linewidth = 0.1,alpha = 0.1)+ # alpha可以修改线的透明度
                    theme_bw() +
                    facet_wrap(~strain,scales = "free_y",ncol = 7,nrow = 11) +
                    #facet_grid(~ province,scales = "free")+ #箱线图分面
                    geom_signif(comparisons = comparision,step_increase = 0.08,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                    theme(strip.text.x = element_text(size = 16, face = "bold")) +  # 设置分面的字字体大小、颜色、背景、边框，)
                    theme(
                        #panel.border=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        axis.title = element_text(size = 10,colour = "black"),
                        axis.text.x = element_text(size = 15,colour = "black"),
                        axis.text.y = element_text(size = 15,colour = "black"),                   
                        axis.title.y = element_text(size = 20,colour = "black"),
                        axis.title.x = element_text(size = 20,colour = "black")) +
                        labs(x = "Intervention stage",y = "Relative abundance")
                ggsave(paste0("/home/wangnan/Shandong-DF/figure/GBTM_group/abundance_change/abundance_change_non_response(union)(后75个菌种).pdf"),p,height = 38.5,width = 21)

# 按照GBTM模型的分类结果，将FBG响应者作为响应组，分别画出响应组和非响应组的代谢物的差异
    group <- read.xlsx("/data4/wangnan/Shandong-DF/origin-data/sd-DFI/Combine_groups.xlsx")
    # 响应组
        id <- group$Patient_id[which(group$FBG == 1)]
        metadata_response <- metadata[which(!is.na(match(metadata$Patient_id,id))),]
        abundance_response <- abundance[,c(1,match(metadata_response$SampleID,colnames(abundance)))]
    # 非响应组
        id_non <- setdiff(group$Patient_id,id)
        metadata_non_response <- metadata[which(!is.na(match(metadata$Patient_id,id_non))),]
        abundance_non_response <- abundance[,c(1,match(metadata_non_response$SampleID,colnames(abundance)))]

    # 响应组画图
        # 利用相对丰度表画图(pdf)
            # 画出前20个marker在膳食纤维干预的三个时间点的相对丰度变化
                sta  <- data.frame(sampleid = "1",strain = "1",abundance = 1,status = "1",stringsAsFactors = FALSE)
                for(i in 1:20)
                {
                    for(j in 2:ncol(abundance_response))
                    {
                        strain_g <- strsplit(abundance_response$`Sample ID`[i],";",fixed = T)[[1]][6]
                            strain_g <- strsplit(strain_g,"__",fixed =T)[[1]][2]
                        sta <- rbind(sta,c(metadata_response$Patient_id[j-1],strain_g,round(as.numeric(abundance_response[i,j])*100,4),metadata_response$Intervention[j-1]))
                    }
                    sta <- sta[-1,]
                    sta[,3] <- unlist(lapply(sta[,3],as.numeric))
                }
                sta$status <- factor(sta$status,levels = c("Early","Mid","Later")) # 添加因子水平
                sta$strain <- factor(sta$strain,levels = unique(sta$strain)) # 添加因子水平

                comparision <- list(c("Early","Mid"),c("Mid","Later"),c("Early","Later")) # 设置显著性检验的组合
                p <- ggplot(sta,aes(x=status,y=abundance)) + 
                    geom_point() +
                    geom_boxplot(aes(middle = mean(abundance)),position=position_dodge(0.7),width = 0.6, alpha = 0.1) +
                    # geom_line(aes(group=sampleid), color="black",linewidth = 0.1,alpha = 0.1)+ # alpha可以修改线的透明度
                    theme_bw() +
                    facet_wrap(~strain,scales = "free_y",ncol = 7,nrow = 3) +
                    #facet_grid(~ province,scales = "free")+ #箱线图分面
                    geom_signif(comparisons = comparision,step_increase = 0.08,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                    theme(strip.text.x = element_text(size = 16, face = "bold")) +  # 设置分面的字字体大小、颜色、背景、边框，)
                    theme(
                        #panel.border=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        axis.title = element_text(size = 10,colour = "black"),
                        axis.text.x = element_text(size = 15,colour = "black"),
                        axis.text.y = element_text(size = 15,colour = "black"),                   
                        axis.title.y = element_text(size = 20,colour = "black"),
                        axis.title.x = element_text(size = 20,colour = "black")) +
                        labs(x = "Intervention stage",y = "Relative abundance")
                ggsave(paste0("/home/wangnan/Shandong-DF/figure/GBTM_group/abundance_change/abundance_change_response(FBG)(前20个菌种).pdf"),p,height = 10.5,width = 21)

            # 画出后75个marker在膳食纤维干预的三个时间点的相对丰度变化
                sta  <- data.frame(sampleid = "1",strain = "1",abundance = 1,status = "1",stringsAsFactors = FALSE)
                for(i in 21:95)
                {
                    for(j in 2:ncol(abundance_response))
                    {
                        strain_g <- strsplit(abundance_response$`Sample ID`[i],";",fixed = T)[[1]][6]
                            strain_g <- strsplit(strain_g,"__",fixed =T)[[1]][2]
                        sta <- rbind(sta,c(metadata_response$Patient_id[j-1],strain_g,round(as.numeric(abundance_response[i,j])*100,4),metadata_response$Intervention[j-1]))
                    }
                    sta <- sta[-1,]
                    sta[,3] <- unlist(lapply(sta[,3],as.numeric))
                }
                sta$status <- factor(sta$status,levels = c("Early","Mid","Later")) # 添加因子水平
                sta$strain <- factor(sta$strain,levels = unique(sta$strain)) # 添加因子水平

                comparision <- list(c("Early","Mid"),c("Mid","Later"),c("Early","Later")) # 设置显著性检验的组合
                p <- ggplot(sta,aes(x=status,y=abundance)) + 
                    geom_point() +
                    geom_boxplot(position=position_dodge(0.7),width = 0.6, alpha = 0.1) +
                    # geom_line(aes(group=sampleid), color="black",linewidth = 0.1,alpha = 0.1)+ # alpha可以修改线的透明度
                    theme_bw() +
                    facet_wrap(~strain,scales = "free_y",ncol = 7,nrow = 11) +
                    #facet_grid(~ province,scales = "free")+ #箱线图分面
                    geom_signif(comparisons = comparision,step_increase = 0.08,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                    theme(strip.text.x = element_text(size = 16, face = "bold")) +  # 设置分面的字字体大小、颜色、背景、边框，)
                    theme(
                        #panel.border=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        axis.title = element_text(size = 10,colour = "black"),
                        axis.text.x = element_text(size = 15,colour = "black"),
                        axis.text.y = element_text(size = 15,colour = "black"),                   
                        axis.title.y = element_text(size = 20,colour = "black"),
                        axis.title.x = element_text(size = 20,colour = "black")) +
                        labs(x = "Intervention stage",y = "Relative abundance")
                ggsave(paste0("/home/wangnan/Shandong-DF/figure/GBTM_group/abundance_change/abundance_change_response(FBG)(后75个菌种).pdf"),p,height = 38.5,width = 21)

    # 非响应组画图
        # 利用相对丰度表画图(pdf)
            # 画出前20个marker在膳食纤维干预的三个时间点的相对丰度变化
                sta  <- data.frame(sampleid = "1",strain = "1",abundance = 1,status = "1",stringsAsFactors = FALSE)
                for(i in 1:20)
                {
                    for(j in 2:ncol(abundance_non_response))
                    {
                        strain_g <- strsplit(abundance_non_response$`Sample ID`[i],";",fixed = T)[[1]][6]
                            strain_g <- strsplit(strain_g,"__",fixed =T)[[1]][2]
                        sta <- rbind(sta,c(metadata_non_response$Patient_id[j-1],strain_g,round(as.numeric(abundance_non_response[i,j])*100,4),metadata_non_response$Intervention[j-1]))
                    }
                    sta <- sta[-1,]
                    sta[,3] <- unlist(lapply(sta[,3],as.numeric))
                }
                sta$status <- factor(sta$status,levels = c("Early","Mid","Later")) # 添加因子水平
                sta$strain <- factor(sta$strain,levels = unique(sta$strain)) # 添加因子水平

                comparision <- list(c("Early","Mid"),c("Mid","Later"),c("Early","Later")) # 设置显著性检验的组合
                p <- ggplot(sta,aes(x=status,y=abundance)) + 
                    geom_point() +
                    geom_boxplot(aes(middle = mean(abundance)),position=position_dodge(0.7),width = 0.6, alpha = 0.1) +
                    # geom_line(aes(group=sampleid), color="black",linewidth = 0.1,alpha = 0.1)+ # alpha可以修改线的透明度
                    theme_bw() +
                    facet_wrap(~strain,scales = "free_y",ncol = 7,nrow = 3) +
                    #facet_grid(~ province,scales = "free")+ #箱线图分面
                    geom_signif(comparisons = comparision,step_increase = 0.08,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                    theme(strip.text.x = element_text(size = 16, face = "bold")) +  # 设置分面的字字体大小、颜色、背景、边框，)
                    theme(
                        #panel.border=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        axis.title = element_text(size = 10,colour = "black"),
                        axis.text.x = element_text(size = 15,colour = "black"),
                        axis.text.y = element_text(size = 15,colour = "black"),                   
                        axis.title.y = element_text(size = 20,colour = "black"),
                        axis.title.x = element_text(size = 20,colour = "black")) +
                        labs(x = "Intervention stage",y = "Relative abundance")
                ggsave(paste0("/home/wangnan/Shandong-DF/figure/GBTM_group/abundance_change/abundance_change_non_response(FBG)(前20个菌种).pdf"),p,height = 10.5,width = 21)

            # 画出后75个marker在膳食纤维干预的三个时间点的相对丰度变化
                sta  <- data.frame(sampleid = "1",strain = "1",abundance = 1,status = "1",stringsAsFactors = FALSE)
                for(i in 21:95)
                {
                    for(j in 2:ncol(abundance_non_response))
                    {
                        strain_g <- strsplit(abundance_non_response$`Sample ID`[i],";",fixed = T)[[1]][6]
                            strain_g <- strsplit(strain_g,"__",fixed =T)[[1]][2]
                        sta <- rbind(sta,c(metadata_non_response$Patient_id[j-1],strain_g,round(as.numeric(abundance_non_response[i,j])*100,4),metadata_non_response$Intervention[j-1]))
                    }
                    sta <- sta[-1,]
                    sta[,3] <- unlist(lapply(sta[,3],as.numeric))
                }
                sta$status <- factor(sta$status,levels = c("Early","Mid","Later")) # 添加因子水平
                sta$strain <- factor(sta$strain,levels = unique(sta$strain)) # 添加因子水平

                comparision <- list(c("Early","Mid"),c("Mid","Later"),c("Early","Later")) # 设置显著性检验的组合
                p <- ggplot(sta,aes(x=status,y=abundance)) + 
                    geom_point() +
                    geom_boxplot(position=position_dodge(0.7),width = 0.6, alpha = 0.1) +
                    # geom_line(aes(group=sampleid), color="black",linewidth = 0.1,alpha = 0.1)+ # alpha可以修改线的透明度
                    theme_bw() +
                    facet_wrap(~strain,scales = "free_y",ncol = 7,nrow = 11) +
                    #facet_grid(~ province,scales = "free")+ #箱线图分面
                    geom_signif(comparisons = comparision,step_increase = 0.08,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                    theme(strip.text.x = element_text(size = 16, face = "bold")) +  # 设置分面的字字体大小、颜色、背景、边框，)
                    theme(
                        #panel.border=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        axis.title = element_text(size = 10,colour = "black"),
                        axis.text.x = element_text(size = 15,colour = "black"),
                        axis.text.y = element_text(size = 15,colour = "black"),                   
                        axis.title.y = element_text(size = 20,colour = "black"),
                        axis.title.x = element_text(size = 20,colour = "black")) +
                        labs(x = "Intervention stage",y = "Relative abundance")
                ggsave(paste0("/home/wangnan/Shandong-DF/figure/GBTM_group/abundance_change/abundance_change_non_response(FBG)(后75个菌种).pdf"),p,height = 38.5,width = 21)


# 按照GBTM模型的分类结果，将HbA1c响应者作为响应组，分别画出响应组和非响应组的代谢物的差异
    group <- read.xlsx("/data4/wangnan/Shandong-DF/origin-data/sd-DFI/Combine_groups.xlsx")
    # 响应组
        id <- group$Patient_id[which(group$HbA1c == 1)]
        metadata_response <- metadata[which(!is.na(match(metadata$Patient_id,id))),]
        abundance_response <- abundance[,c(1,match(metadata_response$SampleID,colnames(abundance)))]
    # 非响应组
        id_non <- setdiff(group$Patient_id,id)
        metadata_non_response <- metadata[which(!is.na(match(metadata$Patient_id,id_non))),]
        abundance_non_response <- abundance[,c(1,match(metadata_non_response$SampleID,colnames(abundance)))]

    # 响应组画图
        # 利用相对丰度表画图(pdf)
            # 画出前20个marker在膳食纤维干预的三个时间点的相对丰度变化
                sta  <- data.frame(sampleid = "1",strain = "1",abundance = 1,status = "1",stringsAsFactors = FALSE)
                for(i in 1:20)
                {
                    for(j in 2:ncol(abundance_response))
                    {
                        strain_g <- strsplit(abundance_response$`Sample ID`[i],";",fixed = T)[[1]][6]
                            strain_g <- strsplit(strain_g,"__",fixed =T)[[1]][2]
                        sta <- rbind(sta,c(metadata_response$Patient_id[j-1],strain_g,round(as.numeric(abundance_response[i,j])*100,4),metadata_response$Intervention[j-1]))
                    }
                    sta <- sta[-1,]
                    sta[,3] <- unlist(lapply(sta[,3],as.numeric))
                }
                sta$status <- factor(sta$status,levels = c("Early","Mid","Later")) # 添加因子水平
                sta$strain <- factor(sta$strain,levels = unique(sta$strain)) # 添加因子水平

                comparision <- list(c("Early","Mid"),c("Mid","Later"),c("Early","Later")) # 设置显著性检验的组合
                p <- ggplot(sta,aes(x=status,y=abundance)) + 
                    geom_point() +
                    geom_boxplot(aes(middle = mean(abundance)),position=position_dodge(0.7),width = 0.6, alpha = 0.1) +
                    # geom_line(aes(group=sampleid), color="black",linewidth = 0.1,alpha = 0.1)+ # alpha可以修改线的透明度
                    theme_bw() +
                    facet_wrap(~strain,scales = "free_y",ncol = 7,nrow = 3) +
                    #facet_grid(~ province,scales = "free")+ #箱线图分面
                    geom_signif(comparisons = comparision,step_increase = 0.08,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                    theme(strip.text.x = element_text(size = 16, face = "bold")) +  # 设置分面的字字体大小、颜色、背景、边框，)
                    theme(
                        #panel.border=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        axis.title = element_text(size = 10,colour = "black"),
                        axis.text.x = element_text(size = 15,colour = "black"),
                        axis.text.y = element_text(size = 15,colour = "black"),                   
                        axis.title.y = element_text(size = 20,colour = "black"),
                        axis.title.x = element_text(size = 20,colour = "black")) +
                        labs(x = "Intervention stage",y = "Relative abundance")
                ggsave(paste0("/home/wangnan/Shandong-DF/figure/GBTM_group/abundance_change/abundance_change_response(HbA1c)(前20个菌种).pdf"),p,height = 10.5,width = 21)

            # 画出后75个marker在膳食纤维干预的三个时间点的相对丰度变化
                sta  <- data.frame(sampleid = "1",strain = "1",abundance = 1,status = "1",stringsAsFactors = FALSE)
                for(i in 21:95)
                {
                    for(j in 2:ncol(abundance_response))
                    {
                        strain_g <- strsplit(abundance_response$`Sample ID`[i],";",fixed = T)[[1]][6]
                            strain_g <- strsplit(strain_g,"__",fixed =T)[[1]][2]
                        sta <- rbind(sta,c(metadata_response$Patient_id[j-1],strain_g,round(as.numeric(abundance_response[i,j])*100,4),metadata_response$Intervention[j-1]))
                    }
                    sta <- sta[-1,]
                    sta[,3] <- unlist(lapply(sta[,3],as.numeric))
                }
                sta$status <- factor(sta$status,levels = c("Early","Mid","Later")) # 添加因子水平
                sta$strain <- factor(sta$strain,levels = unique(sta$strain)) # 添加因子水平

                comparision <- list(c("Early","Mid"),c("Mid","Later"),c("Early","Later")) # 设置显著性检验的组合
                p <- ggplot(sta,aes(x=status,y=abundance)) + 
                    geom_point() +
                    geom_boxplot(position=position_dodge(0.7),width = 0.6, alpha = 0.1) +
                    # geom_line(aes(group=sampleid), color="black",linewidth = 0.1,alpha = 0.1)+ # alpha可以修改线的透明度
                    theme_bw() +
                    facet_wrap(~strain,scales = "free_y",ncol = 7,nrow = 11) +
                    #facet_grid(~ province,scales = "free")+ #箱线图分面
                    geom_signif(comparisons = comparision,step_increase = 0.08,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                    theme(strip.text.x = element_text(size = 16, face = "bold")) +  # 设置分面的字字体大小、颜色、背景、边框，)
                    theme(
                        #panel.border=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        axis.title = element_text(size = 10,colour = "black"),
                        axis.text.x = element_text(size = 15,colour = "black"),
                        axis.text.y = element_text(size = 15,colour = "black"),                   
                        axis.title.y = element_text(size = 20,colour = "black"),
                        axis.title.x = element_text(size = 20,colour = "black")) +
                        labs(x = "Intervention stage",y = "Relative abundance")
                ggsave(paste0("/home/wangnan/Shandong-DF/figure/GBTM_group/abundance_change/abundance_change_response(HbA1c)(后75个菌种).pdf"),p,height = 38.5,width = 21)

    # 非响应组画图
        # 利用相对丰度表画图(pdf)
            # 画出前20个marker在膳食纤维干预的三个时间点的相对丰度变化
                sta  <- data.frame(sampleid = "1",strain = "1",abundance = 1,status = "1",stringsAsFactors = FALSE)
                for(i in 1:20)
                {
                    for(j in 2:ncol(abundance_non_response))
                    {
                        strain_g <- strsplit(abundance_non_response$`Sample ID`[i],";",fixed = T)[[1]][6]
                            strain_g <- strsplit(strain_g,"__",fixed =T)[[1]][2]
                        sta <- rbind(sta,c(metadata_non_response$Patient_id[j-1],strain_g,round(as.numeric(abundance_non_response[i,j])*100,4),metadata_non_response$Intervention[j-1]))
                    }
                    sta <- sta[-1,]
                    sta[,3] <- unlist(lapply(sta[,3],as.numeric))
                }
                sta$status <- factor(sta$status,levels = c("Early","Mid","Later")) # 添加因子水平
                sta$strain <- factor(sta$strain,levels = unique(sta$strain)) # 添加因子水平

                comparision <- list(c("Early","Mid"),c("Mid","Later"),c("Early","Later")) # 设置显著性检验的组合
                p <- ggplot(sta,aes(x=status,y=abundance)) + 
                    geom_point() +
                    geom_boxplot(aes(middle = mean(abundance)),position=position_dodge(0.7),width = 0.6, alpha = 0.1) +
                    # geom_line(aes(group=sampleid), color="black",linewidth = 0.1,alpha = 0.1)+ # alpha可以修改线的透明度
                    theme_bw() +
                    facet_wrap(~strain,scales = "free_y",ncol = 7,nrow = 3) +
                    #facet_grid(~ province,scales = "free")+ #箱线图分面
                    geom_signif(comparisons = comparision,step_increase = 0.08,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                    theme(strip.text.x = element_text(size = 16, face = "bold")) +  # 设置分面的字字体大小、颜色、背景、边框，)
                    theme(
                        #panel.border=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        axis.title = element_text(size = 10,colour = "black"),
                        axis.text.x = element_text(size = 15,colour = "black"),
                        axis.text.y = element_text(size = 15,colour = "black"),                   
                        axis.title.y = element_text(size = 20,colour = "black"),
                        axis.title.x = element_text(size = 20,colour = "black")) +
                        labs(x = "Intervention stage",y = "Relative abundance")
                ggsave(paste0("/home/wangnan/Shandong-DF/figure/GBTM_group/abundance_change/abundance_change_non_response(HbA1c)(前20个菌种).pdf"),p,height = 10.5,width = 21)

            # 画出后75个marker在膳食纤维干预的三个时间点的相对丰度变化
                sta  <- data.frame(sampleid = "1",strain = "1",abundance = 1,status = "1",stringsAsFactors = FALSE)
                for(i in 21:95)
                {
                    for(j in 2:ncol(abundance_non_response))
                    {
                        strain_g <- strsplit(abundance_non_response$`Sample ID`[i],";",fixed = T)[[1]][6]
                            strain_g <- strsplit(strain_g,"__",fixed =T)[[1]][2]
                        sta <- rbind(sta,c(metadata_non_response$Patient_id[j-1],strain_g,round(as.numeric(abundance_non_response[i,j])*100,4),metadata_non_response$Intervention[j-1]))
                    }
                    sta <- sta[-1,]
                    sta[,3] <- unlist(lapply(sta[,3],as.numeric))
                }
                sta$status <- factor(sta$status,levels = c("Early","Mid","Later")) # 添加因子水平
                sta$strain <- factor(sta$strain,levels = unique(sta$strain)) # 添加因子水平

                comparision <- list(c("Early","Mid"),c("Mid","Later"),c("Early","Later")) # 设置显著性检验的组合
                p <- ggplot(sta,aes(x=status,y=abundance)) + 
                    geom_point() +
                    geom_boxplot(position=position_dodge(0.7),width = 0.6, alpha = 0.1) +
                    # geom_line(aes(group=sampleid), color="black",linewidth = 0.1,alpha = 0.1)+ # alpha可以修改线的透明度
                    theme_bw() +
                    facet_wrap(~strain,scales = "free_y",ncol = 7,nrow = 11) +
                    #facet_grid(~ province,scales = "free")+ #箱线图分面
                    geom_signif(comparisons = comparision,step_increase = 0.08,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                    theme(strip.text.x = element_text(size = 16, face = "bold")) +  # 设置分面的字字体大小、颜色、背景、边框，)
                    theme(
                        #panel.border=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        axis.title = element_text(size = 10,colour = "black"),
                        axis.text.x = element_text(size = 15,colour = "black"),
                        axis.text.y = element_text(size = 15,colour = "black"),                   
                        axis.title.y = element_text(size = 20,colour = "black"),
                        axis.title.x = element_text(size = 20,colour = "black")) +
                        labs(x = "Intervention stage",y = "Relative abundance")
                ggsave(paste0("/home/wangnan/Shandong-DF/figure/GBTM_group/abundance_change/abundance_change_non_response(HbA1c)(后75个菌种).pdf"),p,height = 38.5,width = 21)


# 按照GBTM模型的分类结果，将FBG和HbA1c共同的响应者作为响应组，分别画出响应组和非响应组的代谢物的差异
    group <- read.xlsx("/data4/wangnan/Shandong-DF/origin-data/sd-DFI/Combine_groups.xlsx")
    # 响应组
        id <- group$Patient_id[which(group$HbA1c == 1 & group$FBG == 1)]
        metadata_response <- metadata[which(!is.na(match(metadata$Patient_id,id))),]
        abundance_response <- abundance[,c(1,match(metadata_response$SampleID,colnames(abundance)))]
    # 非响应组
        id_non <- setdiff(group$Patient_id,id)
        metadata_non_response <- metadata[which(!is.na(match(metadata$Patient_id,id_non))),]
        abundance_non_response <- abundance[,c(1,match(metadata_non_response$SampleID,colnames(abundance)))]

    # 响应组画图
        # 利用相对丰度表画图(pdf)
            # 画出前20个marker在膳食纤维干预的三个时间点的相对丰度变化
                sta  <- data.frame(sampleid = "1",strain = "1",abundance = 1,status = "1",stringsAsFactors = FALSE)
                for(i in 1:20)
                {
                    for(j in 2:ncol(abundance_response))
                    {
                        strain_g <- strsplit(abundance_response$`Sample ID`[i],";",fixed = T)[[1]][6]
                            strain_g <- strsplit(strain_g,"__",fixed =T)[[1]][2]
                        sta <- rbind(sta,c(metadata_response$Patient_id[j-1],strain_g,round(as.numeric(abundance_response[i,j])*100,4),metadata_response$Intervention[j-1]))
                    }
                    sta <- sta[-1,]
                    sta[,3] <- unlist(lapply(sta[,3],as.numeric))
                }
                sta$status <- factor(sta$status,levels = c("Early","Mid","Later")) # 添加因子水平
                sta$strain <- factor(sta$strain,levels = unique(sta$strain)) # 添加因子水平

                comparision <- list(c("Early","Mid"),c("Mid","Later"),c("Early","Later")) # 设置显著性检验的组合
                p <- ggplot(sta,aes(x=status,y=abundance)) + 
                    geom_point() +
                    geom_boxplot(aes(middle = mean(abundance)),position=position_dodge(0.7),width = 0.6, alpha = 0.1) +
                    # geom_line(aes(group=sampleid), color="black",linewidth = 0.1,alpha = 0.1)+ # alpha可以修改线的透明度
                    theme_bw() +
                    facet_wrap(~strain,scales = "free_y",ncol = 7,nrow = 3) +
                    #facet_grid(~ province,scales = "free")+ #箱线图分面
                    geom_signif(comparisons = comparision,step_increase = 0.08,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                    theme(strip.text.x = element_text(size = 16, face = "bold")) +  # 设置分面的字字体大小、颜色、背景、边框，)
                    theme(
                        #panel.border=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        axis.title = element_text(size = 10,colour = "black"),
                        axis.text.x = element_text(size = 15,colour = "black"),
                        axis.text.y = element_text(size = 15,colour = "black"),                   
                        axis.title.y = element_text(size = 20,colour = "black"),
                        axis.title.x = element_text(size = 20,colour = "black")) +
                        labs(x = "Intervention stage",y = "Relative abundance")
                ggsave(paste0("/home/wangnan/Shandong-DF/figure/GBTM_group/abundance_change/abundance_change_response(FBG&HbA1c)(前20个菌种).pdf"),p,height = 10.5,width = 21)

            # 画出后75个marker在膳食纤维干预的三个时间点的相对丰度变化
    # 非响应组画图
        # 利用相对丰度表画图(pdf)
            # 画出前20个marker在膳食纤维干预的三个时间点的相对丰度变化
                sta  <- data.frame(sampleid = "1",strain = "1",abundance = 1,status = "1",stringsAsFactors = FALSE)
                for(i in 1:20)
                {
                    for(j in 2:ncol(abundance_non_response))
                    {
                        strain_g <- strsplit(abundance_non_response$`Sample ID`[i],";",fixed = T)[[1]][6]
                            strain_g <- strsplit(strain_g,"__",fixed =T)[[1]][2]
                        sta <- rbind(sta,c(metadata_non_response$Patient_id[j-1],strain_g,round(as.numeric(abundance_non_response[i,j])*100,4),metadata_non_response$Intervention[j-1]))
                    }
                    sta <- sta[-1,]
                    sta[,3] <- unlist(lapply(sta[,3],as.numeric))
                }
                sta$status <- factor(sta$status,levels = c("Early","Mid","Later")) # 添加因子水平
                sta$strain <- factor(sta$strain,levels = unique(sta$strain)) # 添加因子水平

                comparision <- list(c("Early","Mid"),c("Mid","Later"),c("Early","Later")) # 设置显著性检验的组合
                p <- ggplot(sta,aes(x=status,y=abundance)) + 
                    geom_point() +
                    geom_boxplot(aes(middle = mean(abundance)),position=position_dodge(0.7),width = 0.6, alpha = 0.1) +
                    # geom_line(aes(group=sampleid), color="black",linewidth = 0.1,alpha = 0.1)+ # alpha可以修改线的透明度
                    theme_bw() +
                    facet_wrap(~strain,scales = "free_y",ncol = 7,nrow = 3) +
                    #facet_grid(~ province,scales = "free")+ #箱线图分面
                    geom_signif(comparisons = comparision,step_increase = 0.08,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                    theme(strip.text.x = element_text(size = 16, face = "bold")) +  # 设置分面的字字体大小、颜色、背景、边框，)
                    theme(
                        #panel.border=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        axis.title = element_text(size = 10,colour = "black"),
                        axis.text.x = element_text(size = 15,colour = "black"),
                        axis.text.y = element_text(size = 15,colour = "black"),                   
                        axis.title.y = element_text(size = 20,colour = "black"),
                        axis.title.x = element_text(size = 20,colour = "black")) +
                        labs(x = "Intervention stage",y = "Relative abundance")
                ggsave(paste0("/home/wangnan/Shandong-DF/figure/GBTM_group/abundance_change/abundance_change_non_response(FBG&HbA1c)(前20个菌种).pdf"),p,height = 10.5,width = 21)

            # 画出后75个marker在膳食纤维干预的三个时间点的相对丰度变化
                sta  <- data.frame(sampleid = "1",strain = "1",abundance = 1,status = "1",stringsAsFactors = FALSE)
                for(i in 21:95)
                {
                    for(j in 2:ncol(abundance_non_response))
                    {
                        strain_g <- strsplit(abundance_non_response$`Sample ID`[i],";",fixed = T)[[1]][6]
                            strain_g <- strsplit(strain_g,"__",fixed =T)[[1]][2]
                        sta <- rbind(sta,c(metadata_non_response$Patient_id[j-1],strain_g,round(as.numeric(abundance_non_response[i,j])*100,4),metadata_non_response$Intervention[j-1]))
                    }
                    sta <- sta[-1,]
                    sta[,3] <- unlist(lapply(sta[,3],as.numeric))
                }
                sta$status <- factor(sta$status,levels = c("Early","Mid","Later")) # 添加因子水平
                sta$strain <- factor(sta$strain,levels = unique(sta$strain)) # 添加因子水平

                comparision <- list(c("Early","Mid"),c("Mid","Later"),c("Early","Later")) # 设置显著性检验的组合
                p <- ggplot(sta,aes(x=status,y=abundance)) + 
                    geom_point() +
                    geom_boxplot(position=position_dodge(0.7),width = 0.6, alpha = 0.1) +
                    # geom_line(aes(group=sampleid), color="black",linewidth = 0.1,alpha = 0.1)+ # alpha可以修改线的透明度
                    theme_bw() +
                    facet_wrap(~strain,scales = "free_y",ncol = 7,nrow = 11) +
                    #facet_grid(~ province,scales = "free")+ #箱线图分面
                    geom_signif(comparisons = comparision,step_increase = 0.08,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                    theme(strip.text.x = element_text(size = 16, face = "bold")) +  # 设置分面的字字体大小、颜色、背景、边框，)
                    theme(
                        #panel.border=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        axis.title = element_text(size = 10,colour = "black"),
                        axis.text.x = element_text(size = 15,colour = "black"),
                        axis.text.y = element_text(size = 15,colour = "black"),                   
                        axis.title.y = element_text(size = 20,colour = "black"),
                        axis.title.x = element_text(size = 20,colour = "black")) +
                        labs(x = "Intervention stage",y = "Relative abundance")
                ggsave(paste0("/home/wangnan/Shandong-DF/figure/GBTM_group/abundance_change/abundance_change_non_response(FBG&HbA1c)(后75个菌种).pdf"),p,height = 38.5,width = 21)
