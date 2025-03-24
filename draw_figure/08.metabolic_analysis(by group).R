library(tidyverse)
library(ggplot2)
library(ggpubr)
library(psych) #psych包用于计算相关性、p值等信息
library(dplyr)
library(reshape2)
library(openxlsx)
#library(ggtree) #用于画进化树

rm(list=ls())
gc()

# 筛选GGMP,SGMP和SD-DFI都存在的marker
    #GGMP的genus水平的丰度表
        data_ggmp <- read_tsv("/data4/wangnan/Shandong-DF/origin-data/ggmp/table-filtered-feature-rarefied5k-L6.tsv")
    #SGMP的genus水平的丰度表
        data_sgmp <- read_tsv("/data4/wangnan/Shandong-DF/origin-data/sgmp/table-filtered-feature-rarefied5k-L6.tsv")
    #SD-DFI的丰度表
        data_dfi <- read_tsv("/data4/wangnan/Shandong-DF/origin-data/sd-DFI/sd-dfi_abundance_processed.tsv")
        #筛选出来有丰度的marker
            no_abundance_marker <- read.csv("/home/wangnan/Shandong-DF/table/no_abundance_marker.csv")
            data_dfi <- data_dfi[,setdiff(1:ncol(data_dfi),match(no_abundance_marker[,2],colnames(data_dfi)))] 
        #metadata
            metadata <- read.csv("/data4/wangnan/Shandong-DF/origin-data/sd-DFI/DFI_metadata_250_used.csv")
    #筛选出来三个队列都存在的marker
        tmp <- intersect(colnames(data_dfi)[2:ncol(data_dfi)],colnames(data_sgmp)[2:ncol(data_sgmp)])
        tmp <- intersect(tmp,colnames(data_ggmp)[2:ncol(data_ggmp)])
        abundance <- data_dfi[,c(1,match(tmp,colnames(data_dfi)))]

# 筛选所使用到的250个样本的丰度
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

# 按照transfer model选出的marker的重要性按由大到小进行排列
    marker_contribution_transfer <- read.csv("/home/wangnan/Shandong-DF/table/markers selected by transfer DNN model.csv")
        pos <- match(abundance$`Sample ID`,marker_contribution_transfer$Marker)
        marker_contribution_transfer <- marker_contribution_transfer[pos,]
    marker_contribution_transfer <- marker_contribution_transfer[order(marker_contribution_transfer$Mean_of_absolute_value,decreasing = T),]
    pos <- match(marker_contribution_transfer$Marker,abundance$`Sample ID`)
        abundance <- abundance[pos,] # 按照重要性由大到小进行排序    

# 使丰度表转化为相对丰度表
    no_abundance_marker <- c()
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

# 按照GBTM模型的分类结果，分别画出响应组和非响应组的代谢物的差异
    group <- read.xlsx("/data4/wangnan/Shandong-DF/origin-data/sd-DFI/Combine_groups.xlsx")
    # 响应组
        id <- group$Patient_id[which(group$union == 1)]
        metadata_response <- metadata[which(!is.na(match(metadata$Patient_id,id))),]
        abundance_response <- abundance[,c(1,match(metadata_response$SampleID,colnames(abundance)))]
    # 非响应组
        id_non <- setdiff(group$Patient_id,id)
        metadata_non_response <- metadata[which(!is.na(match(metadata$Patient_id,id_non))),]
        abundance_non_response <- abundance[,c(1,match(metadata_non_response$SampleID,colnames(abundance)))]

# 按照不同阶段画出来各种代谢物的变化情况
    # 响应组
        metabolism <- metadata_response[,c(2,23:50,52)]
        # 将na的值都取标准差，然后添加上一个随机的噪声，噪声的标准差为其余值的标准差
            for(i in 2:ncol(metabolism))
            {
                pos <- which(is.na(metabolism[,i]))
                if(length(pos) > 0)
                {
                    mean <- mean(metabolism[setdiff(1:nrow(metabolism),pos),i])
                    sd <- sd(metabolism[setdiff(1:nrow(metabolism),pos),i])
                    metabolism[pos,i] <- rnorm(length(pos),mean,sd)
                }
                metabolism[,i][which(metabolism[,i] < 0)] <- 0
            }
        # 将数据转换为长格式，"metabolism"指的是表头，"concentration"指的是值，"-SampleID"指的是不需要转换的列，"-Intervention"指的是不需要转换的列
            metabolism <- gather(metabolism,metabolism,concentration,-SampleID,-Intervention) 
            metabolism$concentration <- as.numeric(metabolism$concentration)
        # 设置因子
            metabolism$Intervention <- factor(metabolism$Intervention,levels = c("Early","Mid","Later"))
        # 设置显著性检验的组合   
            comparision <- list(c("Early","Later")) 
        # 计算显著性水平，选出来有差异的变化的代谢物和无差异变化的代谢物，以及按照增高和降低进行区分
            metabolism$p_value <- "1"
            for(m in unique(metabolism$metabolism))
            {
                tmp <- metabolism[which(metabolism$metabolism == m),] # 选出来此种代谢物的行
                    early <- tmp[which(tmp$Intervention == "Early"),] #选出来为阶段为"Early"的组
                    later <- tmp[which(tmp$Intervention == "Later"),] #选出来为阶段为"Later"的组
                pos <- which(tmp$Intervention == "Mid")
                    tmp <- tmp[setdiff(1:nrow(tmp),pos),] #去掉阶段为"Mid"的组
                wil_t <- wilcox.test(concentration ~ Intervention,tmp[,c(2,4)],var.equal = TRUE) #计算wilcoxon秩和检验的p值，判断是否显著
                p_value <- wil_t$p.value #提取p值
                #如果p值小于0.05，则判断为差异显著,否则为不显著
                    if(p_value <= 0.05)
                    {
                        if(mean(early$concentration) > mean(later$concentration))
                        {
                            metabolism$p_value[which(metabolism$metabolism == m)] <- "significant_lower"
                        }else {
                            metabolism$p_value[which(metabolism$metabolism == m)] <- "significant_higher"
                        }
                    }else 
                    {
                        metabolism$p_value[which(metabolism$metabolism == m)] <- "no_significant"
                    }
            }

        # 画图
            #显著降低的
                tmp <- metabolism[which(metabolism$p_value == "significant_lower"),]
                # (1)
                    # tmp <- tmp[which(tmp$metabolism == "AST_ALT" | tmp$metabolism == "Free_fatty_acids" | tmp$metabolism == "HDL"),]
                # (2)
                    #tmp <- tmp[-which(tmp$metabolism == "AST_ALT" | tmp$metabolism == "Free_fatty_acids" | tmp$metabolism == "HDL"),]
                
                # Alanine_aminotransferase很奇怪，平均数降低了但是中位数没有降低，
                    # t1 <- tmp[which(tmp$metabolism == "Alanine_aminotransferase" & tmp$Intervention=="Early"),]
                    # t2 <- tmp[which(tmp$metabolism == "Alanine_aminotransferase" & tmp$Intervention=="Later"),]
                    # median(t1$concentration)
                    # median(t2$concentration)
                    significant_lower <- 
                        ggplot(tmp,aes(x=Intervention,y=concentration)) + 
                        #geom_point() +
                        geom_boxplot(size = 1,aes(middle=mean(concentration),color = Intervention)) +
                        scale_color_manual(values=c("#606060","#99C0EA","#3757AF"),aesthetics = "color") + #手动添加颜色
                        geom_jitter(aes(colour=Intervention),width =0.1,size=0.2)+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
                        #geom_line(aes(group=sampleid), col="black",size = 0.1)+
                        theme_bw() +                        
                        facet_wrap(~metabolism,scales = "free_y",nrow = 1,ncol =6) +
                        geom_signif(comparisons = comparision,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                        theme(strip.text.x = element_text(size = 14, face = "bold")) +  # 设置分面的字字体大小、颜色、背景、边框，)
                        theme(
                            #panel.border=element_blank(),
                            panel.grid.major=element_blank(),
                            panel.grid.minor=element_blank(),
                            axis.title = element_text(size = 10,colour = "black"),
                            axis.text.x = element_blank(),
                            axis.ticks.x = element_blank(),
                            axis.text.y = element_text(size = 12,colour = "black"),                   
                            axis.title.y = element_text(size = 20,colour = "black"),
                            axis.title.x = element_text(size = 20,colour = "black"),
                            legend.position = "none") +
                            labs(x = "",y = "")
                    # ggsave(paste0("/home/wangnan/Shandong-DF/figure/GBTM_group/metabolism_change/metabolism_change_response(lower)(1).pdf"),significant_lower,height = 4.5,width = 12)
                    ggsave(paste0("/home/wangnan/Shandong-DF/figure/GBTM_group/metabolism_change/metabolism_change_response(lower)(2).pdf"),significant_lower,height = 4.5,width = 12)

            #显著增高的
                tmp <- metabolism[which(metabolism$p_value == "significant_higher"),]
                # (1)
                    # tmp <- tmp[which(tmp$metabolism == "Leukocyte" | tmp$metabolism == "Lymphocyte" | tmp$metabolism == "TG"),]
                # (2)
                    # tmp <- tmp[-which(tmp$metabolism == "Leukocyte" | tmp$metabolism == "Lymphocyte" | tmp$metabolism == "TG"),]
                # TG很奇怪，平均数升高了但是中位数没有升高，
                    t1 <- tmp[which(tmp$metabolism == "TG" & tmp$Intervention=="Early"),]
                    t2 <- tmp[which(tmp$metabolism == "TG" & tmp$Intervention=="Later"),]
                    mean(t1$concentration)
                    mean(t2$concentration)

                    significant_higher <- 
                        ggplot(tmp,aes(x=Intervention,y=concentration)) + 
                        #geom_point() +
                        geom_boxplot(size = 1,aes(middle=mean(concentration),colour = Intervention)) +
                        scale_color_manual(values=c("#606060","#D7A023","#F4722C"),aesthetics = "color") + #手动添加颜色
                        geom_jitter(aes(colour=Intervention),width =0.1,size=0.2)+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
                        #geom_line(aes(group=sampleid), col="black",size = 0.1)+
                        theme_bw() +
                        facet_wrap(~metabolism,scales = "free_y",ncol =6,nrow =1) +
                        #facet_grid(~ province,scales = "free")+ #箱线图分面
                        geom_signif(comparisons = comparision,step_increase = 0.08,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                        theme(strip.text.x = element_text(size = 14, face = "bold")) +  # 设置分面的字字体大小、颜色、背景、边框，)
                        theme(
                            #panel.border=element_blank(),
                            panel.grid.major=element_blank(),
                            panel.grid.minor=element_blank(),
                            axis.title = element_text(size = 10,colour = "black"),
                            axis.text.x = element_blank(),
                            axis.ticks.x = element_blank(),
                            axis.text.y = element_text(size = 12,colour = "black"),                   
                            axis.title.y = element_text(size = 20,colour = "black"),
                            axis.title.x = element_text(size = 20,colour = "black"),
                            legend.position = "none") +
                            labs(x = "",y = "")
                    # ggsave(paste0("/home/wangnan/Shandong-DF/figure/GBTM_group/metabolism_change/metabolism_change_response(higher)(1).pdf"),significant_higher,height = 4.5,width = 12)
                        ggsave(paste0("/home/wangnan/Shandong-DF/figure/GBTM_group/metabolism_change/metabolism_change_response(higher)(2).pdf"),significant_higher,height = 4.5,width = 8)

            #没有明显变化的
                tmp <- metabolism[which(metabolism$p_value == "no_significant"),]
                    no_significant <- 
                        ggplot(tmp,aes(x=Intervention,y=concentration)) + 
                        #geom_point() +
                        scale_color_manual(values=c("#606060","#606060","#606060"),aesthetics = "color") + #手动添加颜色
                        geom_boxplot(aes(middle=mean(concentration)),size = 1) +
                        geom_jitter(width =0.1,size=0.2)+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
                        #geom_line(aes(group=sampleid), col="black",size = 0.1)+
                        theme_bw() +
                        facet_wrap(~metabolism,scales = "free_y",nrow =3,ncol =6) +
                        #facet_grid(~ province,scales = "free")+ #箱线图分面
                        geom_signif(comparisons = comparision,step_increase = 0.08,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                        theme(strip.text.x = element_text(size = 14, face = "bold")) +  # 设置分面的字字体大小、颜色、背景、边框，)
                        theme(
                            #panel.border=element_blank(),
                            panel.grid.major=element_blank(),
                            panel.grid.minor=element_blank(),
                            axis.title = element_text(size = 10,colour = "black"),
                            axis.text.x = element_blank(),
                            axis.ticks.x = element_blank(),
                            axis.text.y = element_text(size = 12,colour = "black"),                   
                            axis.title.y = element_text(size = 20,colour = "black"),
                            axis.title.x = element_text(size = 20,colour = "black"),
                            legend.position = "top") +
                            labs(x = "",y = "")
                    ggsave(paste0("/home/wangnan/Shandong-DF/figure/GBTM_group/metabolism_change/metabolism_change_response(no_significant).pdf"),no_significant,height = 13.5,width = 24)

    # 非响应组
        metabolism <- metadata_non_response[,c(2,23:50,52)]
        # 将na的值都取标准差，然后添加上一个随机的噪声，噪声的标准差为其余值的标准差
            for(i in 2:ncol(metabolism))
            {
                pos <- which(is.na(metabolism[,i]))
                if(length(pos) > 0)
                {
                    mean <- mean(metabolism[setdiff(1:nrow(metabolism),pos),i])
                    sd <- sd(metabolism[setdiff(1:nrow(metabolism),pos),i])
                    metabolism[pos,i] <- rnorm(length(pos),mean,sd)
                }
                metabolism[,i][which(metabolism[,i] < 0)] <- 0
            }
        # 将数据转换为长格式，"metabolism"指的是表头，"concentration"指的是值，"-SampleID"指的是不需要转换的列，"-Intervention"指的是不需要转换的列
            metabolism <- gather(metabolism,metabolism,concentration,-SampleID,-Intervention) 
            metabolism$concentration <- as.numeric(metabolism$concentration)
        # 设置因子
            metabolism$Intervention <- factor(metabolism$Intervention,levels = c("Early","Mid","Later"))
        # 设置显著性检验的组合   
            comparision <- list(c("Early","Later")) 
        # 计算显著性水平，选出来有差异的变化的代谢物和无差异变化的代谢物，以及按照增高和降低进行区分
            metabolism$p_value <- "1"
            for(m in unique(metabolism$metabolism))
            {
                tmp <- metabolism[which(metabolism$metabolism == m),] # 选出来此种代谢物的行
                    early <- tmp[which(tmp$Intervention == "Early"),] #选出来为阶段为"Early"的组
                    later <- tmp[which(tmp$Intervention == "Later"),] #选出来为阶段为"Later"的组
                pos <- which(tmp$Intervention == "Mid")
                    tmp <- tmp[setdiff(1:nrow(tmp),pos),] #去掉阶段为"Mid"的组
                wil_t <- wilcox.test(concentration ~ Intervention,tmp[,c(2,4)],var.equal = TRUE) #计算wilcoxon秩和检验的p值，判断是否显著
                p_value <- wil_t$p.value #提取p值
                #如果p值小于0.05，则判断为差异显著,否则为不显著
                    if(p_value <= 0.05)
                    {
                        if(mean(early$concentration) > mean(later$concentration))
                        {
                            metabolism$p_value[which(metabolism$metabolism == m)] <- "significant_lower"
                        }else {
                            metabolism$p_value[which(metabolism$metabolism == m)] <- "significant_higher"
                        }
                    }else 
                    {
                        metabolism$p_value[which(metabolism$metabolism == m)] <- "no_significant"
                    }
            }

        # 画图
            #显著降低的
                tmp <- metabolism[which(metabolism$p_value == "significant_lower"),]
                    significant_lower <- 
                        ggplot(tmp,aes(x=Intervention,y=concentration)) + 
                        #geom_point() +
                        geom_boxplot(size = 1,aes(middle=mean(concentration),color = Intervention)) +
                        scale_color_manual(values=c("#606060","#99C0EA","#3757AF"),aesthetics = "color") + #手动添加颜色
                        geom_jitter(aes(colour=Intervention),width =0.1,size=0.2)+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
                        #geom_line(aes(group=sampleid), col="black",size = 0.1)+
                        theme_bw() +
                        facet_wrap(~metabolism,scales = "free_y",nrow = 1,ncol =6) +
                        #facet_grid(~ province,scales = "free")+ #箱线图分面
                        geom_signif(comparisons = comparision,step_increase = 0.08,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                        theme(strip.text.x = element_text(size = 14, face = "bold")) +  # 设置分面的字字体大小、颜色、背景、边框，)
                        theme(
                            #panel.border=element_blank(),
                            panel.grid.major=element_blank(),
                            panel.grid.minor=element_blank(),
                            axis.title = element_text(size = 10,colour = "black"),
                            axis.text.x = element_blank(),
                            axis.ticks.x = element_blank(),
                            axis.text.y = element_text(size = 12,colour = "black"),                   
                            axis.title.y = element_text(size = 20,colour = "black"),
                            axis.title.x = element_text(size = 20,colour = "black"),
                            legend.position = "none") +
                            labs(x = "",y = "")
                    ggsave(paste0("/home/wangnan/Shandong-DF/figure/GBTM_group/metabolism_change/metabolism_change_non_response(lower).pdf"),significant_lower,height = 4.5,width = 16)
            

            #显著增高的
                tmp <- metabolism[which(metabolism$p_value == "significant_higher"),]
                    significant_higher <- 
                        ggplot(tmp,aes(x=Intervention,y=concentration)) + 
                        #geom_point() +
                        geom_boxplot(size = 1,aes(middle=mean(concentration),colour = Intervention)) +
                        scale_color_manual(values=c("#606060","#D7A023","#F4722C"),aesthetics = "color") + #手动添加颜色
                        geom_jitter(aes(colour=Intervention),width =0.1,size=0.2)+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
                        #geom_line(aes(group=sampleid), col="black",size = 0.1)+
                        theme_bw() +
                        facet_wrap(~metabolism,scales = "free_y",ncol =6,nrow =1) +
                        #facet_grid(~ province,scales = "free")+ #箱线图分面
                        geom_signif(comparisons = comparision,step_increase = 0.08,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                        theme(strip.text.x = element_text(size = 14, face = "bold")) +  # 设置分面的字字体大小、颜色、背景、边框，)
                        theme(
                            #panel.border=element_blank(),
                            panel.grid.major=element_blank(),
                            panel.grid.minor=element_blank(),
                            axis.title = element_text(size = 10,colour = "black"),
                            axis.text.x = element_blank(),
                            axis.ticks.x = element_blank(),
                            axis.text.y = element_text(size = 12,colour = "black"),                   
                            axis.title.y = element_text(size = 20,colour = "black"),
                            axis.title.x = element_text(size = 20,colour = "black"),
                            legend.position = "none") +
                            labs(x = "",y = "")
                    ggsave(paste0("/home/wangnan/Shandong-DF/figure/GBTM_group/metabolism_change/metabolism_change_non_response(higher).pdf"),significant_higher,height = 4.5,width = 12)
            
            #没有明显变化的
                tmp <- metabolism[which(metabolism$p_value == "no_significant"),]
                    no_significant <- 
                        ggplot(tmp,aes(x=Intervention,y=concentration)) + 
                        #geom_point() +
                        scale_color_manual(values=c("#606060","#606060","#606060"),aesthetics = "color") + #手动添加颜色
                        geom_boxplot(aes(middle=mean(concentration)),size = 1) +
                        geom_jitter(width =0.1,size=0.2)+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
                        #geom_line(aes(group=sampleid), col="black",size = 0.1)+
                        theme_bw() +
                        facet_wrap(~metabolism,scales = "free_y",nrow =4,ncol =6) +
                        #facet_grid(~ province,scales = "free")+ #箱线图分面
                        geom_signif(comparisons = comparision,step_increase = 0.08,map_signif_level = T,test = wilcox.test) + #添加显著性检验
                        theme(strip.text.x = element_text(size = 14, face = "bold")) +  # 设置分面的字字体大小、颜色、背景、边框，)
                        theme(
                            #panel.border=element_blank(),
                            panel.grid.major=element_blank(),
                            panel.grid.minor=element_blank(),
                            axis.title = element_text(size = 10,colour = "black"),
                            axis.text.x = element_blank(),
                            axis.ticks.x = element_blank(),
                            axis.text.y = element_text(size = 12,colour = "black"),                   
                            axis.title.y = element_text(size = 20,colour = "black"),
                            axis.title.x = element_text(size = 20,colour = "black"),
                            legend.position = "top") +
                            labs(x = "",y = "")
                    ggsave(paste0("/home/wangnan/Shandong-DF/figure/GBTM_group/metabolism_change/metabolism_change_non_response(no_significant).pdf"),no_significant,height = 18,width = 24)

# 综合三个阶段进行微生物组(tansfer DNN model中排名前40个微生物)和代谢组的pearson相关性分析
    # 响应组
        #处理代谢物
            #  tmp1 <- abundance
            #  tmp2 <- metabolism    
            metabolism <- metadata_response[,c(2,23:50,52)]
            abundance <- abundance_response
    # 非响应组
        #处理代谢物
            #  tmp1 <- abundance
            #  tmp2 <- metabolism    
            metabolism <- metadata_non_response[,c(2,23:50,52)]
            abundance <- abundance_non_response

        # 将na的值都取标准差，然后添加上一个随机的噪声，噪声的标准差为其余值的标准差
            for(i in 2:ncol(metabolism))
            {
                pos <- which(is.na(metabolism[,i]))
                if(length(pos) > 0)
                {
                    mean <- mean(metabolism[setdiff(1:nrow(metabolism),pos),i])
                    sd <- sd(metabolism[setdiff(1:nrow(metabolism),pos),i])
                    metabolism[pos,i] <- rnorm(length(pos),mean,sd)
                }
                metabolism[,i][which(metabolism[,i] < 0)] <- 0
            }
        abundance <- abundance[1:40,]

    #提取列名(即各种菌种，并且对其处理，变成pylum+genus的形式)
        cnames <- c("SampleID",as.data.frame(abundance[,1])[,1]) 
        for(i in 2:length(cnames))
        {
            strain_g <- strsplit(cnames[i],";",fixed = T)[[1]][6]
                strain_g <- strsplit(strain_g,"__",fixed =T)[[1]][2]
            cnames[i] <- strain_g
        } 

    #丰度表转置
        abun <- as.data.frame(t(as.data.frame(abundance)))
        abun <- abun[-1,] #去掉转置后第一行的各种微生物分类名
        abun[1:ncol(abun)] <- lapply(abun[1:ncol(abun)],as.numeric) #将reads数变为数值型
        abun <- cbind(as.data.frame(rownames(abun)),abun)
        rownames(abun) <- NULL
        colnames(abun) <- cnames #将列名改为各种菌种
    
    #画出来菌种和代谢物的相关性热图
        #计算相关性矩阵（可选：”pearson”、”spearman”、”kendall”相关系数）、p值矩阵
            cor <- corr.test(abun[,2:ncol(abun)], metabolism[,2:(ncol(metabolism)-1)], method = "pearson",adjust="none")
            cmt <- cor$r #相关性矩阵
                cmt[is.na(cmt)] <- 0
            pmt <- cor$p #p值矩阵
                pmt[is.na(pmt)] <- 1
        #给微生物和代谢物分别进行聚类 
            row <- hclust(dist(cmt)) #对行进行聚类
            col <- hclust(dist(t(cmt))) #对列进行聚类
        #按照聚类结果进行排序
            cmt <- cmt[,col$order]
            cmt <- cmt[row$order,]
            pmt <- pmt[,col$order]
            pmt <- pmt[row$order,]

        #合并成长数据表
            df_cmt <- melt(cmt,value.name="cor")
            df_pmt <- melt(pmt,value.name="p")
            df <- bind_cols(df_cmt,df_pmt)
            df <- df[,c(1,2,3,6)]
            colnames(df) <- c("strain","metabolism","cor","p")

        #对p值进行校正(因为p<0.05被认为是显著相关的，但是对于很多个实验，所有都被认为是显著相关则错误率较高，因此需要进行FDR矫正)
            #df$p <- p.adjust(df$p,method = "fdr")

        #p值大于0.05的取"",p值小于0.05大于0.01的取"*",p值小于0.01的取"**",小于0.001的取"***"
            n_sig <- df$p > 0.05
            sig_1 <- df$p <= 0.05 & df$p > 0.01
            sig_2 <- df$p <= 0.01 & df$p > 0.001    
            sig_3 <- df$p <= 0.001
            df$p[which(n_sig == T)] <- ""
            df$p[which(sig_1 == T)] <- "*" 
            df$p[which(sig_2 == T)] <- "**"
            df$p[which(sig_3 == T)] <- "***"

        #把所有的代谢物名称中的"_"都替换成"
            levels(df$metabolism) <- gsub("_"," ",levels(df$metabolism))
    
        #画图
            heatmap <-
                ggplot(df,aes(x =strain,y = metabolism ,fill = cor)) +
                geom_tile() + #选择画热图
                geom_text(aes(label = p), color = "black", size = 5) + #添加显著性标记
                #coord_fixed()+
                scale_fill_gradient2(low = "#407EB1", high = "#D72E2A",limits = c(-0.7, 0.7),breaks = c(-1,-0.5,0,0.5,1)) + ###设置热图标尺的颜色、范围以及刻度###
                theme(axis.title.x = element_text(size = 17,colour = "black"),
                axis.text.x = element_text(size = 12,colour = "black",angle =90,vjust = 0.5),
                axis.title.y = element_text(size = 17,colour = "black"),
                axis.text.y = element_text(size = 12,colour = "black"),
                plot.title = element_text(hjust = 0.5,size = 20,colour = "black"),
                legend.title = element_blank(),
                legend.text = element_text(size = 12,colour = "black")) +
                xlab("")+
                ylab("")+
            labs(title = "")
            # ggsave("/home/wangnan/Shandong-DF/figure/GBTM_group/heatmap/heatmap_response(前40个微生物).pdf",heatmap,width = 15,height =7.5)
            ggsave("/home/wangnan/Shandong-DF/figure/GBTM_group/heatmap/heatmap_non_response(前40个微生物).pdf",heatmap,width = 15,height =7.5)

# 综合三个阶段进行微生物组(tansfer DNN model中排名后55个微生物)和代谢组的pearson相关性分析
    # 响应组
        #处理代谢物
            #  tmp1 <- abundance
            #  tmp2 <- metabolism    
            metabolism <- metadata_response[,c(2,23:50,52)]
            abundance <- abundance_response
    # 非响应组
        #处理代谢物
            #  tmp1 <- abundance
            #  tmp2 <- metabolism    
            metabolism <- metadata_non_response[,c(2,23:50,52)]
            abundance <- abundance_non_response

        # 将na的值都取标准差，然后添加上一个随机的噪声，噪声的标准差为其余值的标准差
            for(i in 2:ncol(metabolism))
            {
                pos <- which(is.na(metabolism[,i]))
                if(length(pos) > 0)
                {
                    mean <- mean(metabolism[setdiff(1:nrow(metabolism),pos),i])
                    sd <- sd(metabolism[setdiff(1:nrow(metabolism),pos),i])
                    metabolism[pos,i] <- rnorm(length(pos),mean,sd)
                }
                metabolism[,i][which(metabolism[,i] < 0)] <- 0
            }
        abundance <- abundance[41:95,]

    #提取列名(即各种菌种，并且对其处理，变成pylum+genus的形式)
        cnames <- c("SampleID",as.data.frame(abundance[,1])[,1]) 
        for(i in 2:length(cnames))
        {
            strain_g <- strsplit(cnames[i],";",fixed = T)[[1]][6]
                strain_g <- strsplit(strain_g,"__",fixed =T)[[1]][2]
            cnames[i] <- strain_g
        } 

    #丰度表转置
        abun <- as.data.frame(t(as.data.frame(abundance)))
        abun <- abun[-1,] #去掉转置后第一行的各种微生物分类名
        abun[1:ncol(abun)] <- lapply(abun[1:ncol(abun)],as.numeric) #将reads数变为数值型
        abun <- cbind(as.data.frame(rownames(abun)),abun)
        rownames(abun) <- NULL
        colnames(abun) <- cnames #将列名改为各种菌种
    
    #画出来菌种和代谢物的相关性热图
        #计算相关性矩阵（可选：”pearson”、”spearman”、”kendall”相关系数）、p值矩阵
            cor <- corr.test(abun[,2:ncol(abun)], metabolism[,2:(ncol(metabolism)-1)], method = "pearson",adjust="none")
            cmt <- cor$r #相关性矩阵
                cmt[is.na(cmt)] <- 0
            pmt <- cor$p #p值矩阵
                pmt[is.na(pmt)] <- 1
        #给微生物和代谢物分别进行聚类 
            row <- hclust(dist(cmt)) #对行进行聚类
            col <- hclust(dist(t(cmt))) #对列进行聚类
        #按照聚类结果进行排序
            cmt <- cmt[,col$order]
            cmt <- cmt[row$order,]
            pmt <- pmt[,col$order]
            pmt <- pmt[row$order,]

        #合并成长数据表
            df_cmt <- melt(cmt,value.name="cor")
            df_pmt <- melt(pmt,value.name="p")
            df <- bind_cols(df_cmt,df_pmt)
            df <- df[,c(1,2,3,6)]
            colnames(df) <- c("strain","metabolism","cor","p")

        #对p值进行校正(因为p<0.05被认为是显著相关的，但是对于很多个实验，所有都被认为是显著相关则错误率较高，因此需要进行FDR矫正)
            #df$p <- p.adjust(df$p,method = "fdr")

        #p值大于0.05的取"",p值小于0.05大于0.01的取"*",p值小于0.01的取"**",小于0.001的取"***"
            n_sig <- df$p > 0.05
            sig_1 <- df$p <= 0.05 & df$p > 0.01
            sig_2 <- df$p <= 0.01 & df$p > 0.001    
            sig_3 <- df$p <= 0.001
            df$p[which(n_sig == T)] <- ""
            df$p[which(sig_1 == T)] <- "*" 
            df$p[which(sig_2 == T)] <- "**"
            df$p[which(sig_3 == T)] <- "***"

        #把所有的代谢物名称中的"_"都替换成"
            levels(df$metabolism) <- gsub("_"," ",levels(df$metabolism))
    
        #画图
            heatmap <-
                ggplot(df,aes(x =strain,y = metabolism ,fill = cor)) +
                geom_tile() + #选择画热图
                geom_text(aes(label = p), color = "black", size = 5) + #添加显著性标记
                #coord_fixed()+
                scale_fill_gradient2(low = "#407EB1", high = "#D72E2A",limits = c(-0.9, 0.9),breaks = c(-1,-0.5,0,0.5,1)) + ###设置热图标尺的颜色、范围以及刻度###
                theme(axis.title.x = element_text(size = 17,colour = "black"),
                axis.text.x = element_text(size = 12,colour = "black",angle =90,vjust = 0.5),
                axis.title.y = element_text(size = 17,colour = "black"),
                axis.text.y = element_text(size = 12,colour = "black"),
                plot.title = element_text(hjust = 0.5,size = 20,colour = "black"),
                legend.title = element_blank(),
                legend.text = element_text(size = 12,colour = "black")) +
                xlab("")+
                ylab("")+
            labs(title = "")
            # ggsave("/home/wangnan/Shandong-DF/figure/GBTM_group/heatmap/heatmap_response(后55个微生物).pdf",heatmap,width = 15,height =7.5)
            ggsave("/home/wangnan/Shandong-DF/figure/GBTM_group/heatmap/heatmap_non_response(后55个微生物).pdf",heatmap,width = 15,height =7.5)

# 按照膳食干预的Early,Mid,Later阶段和前40个菌种分别进行pearson相关性分析
    # 响应组
        #处理代谢物
            metabolism <- metadata_response[,c(2,23:50,52)]   
                pos <- which(metadata_response$Intervention == "Later")
                metabolism <- metabolism[pos,]
            abundance <- abundance_response[,c(1,pos+1)]   
                abundance <- abundance[1:40,]

    # 非响应组
        #处理代谢物
            metabolism <- metadata_non_response[,c(2,23:50,52)]   
                pos <- which(metadata_non_response$Intervention == "Later")
                metabolism <- metabolism[pos,]
            abundance <- abundance_non_response[,c(1,pos+1)]   
                abundance <- abundance[1:40,]

    #提取列名(即各种菌种，并且对其处理，变成pylum+genus的形式)
        cnames <- c("SampleID",as.data.frame(abundance[,1])[,1]) 
        for(i in 2:length(cnames))
        {
            strain_g <- strsplit(cnames[i],";",fixed = T)[[1]][6]
                strain_g <- strsplit(strain_g,"__",fixed =T)[[1]][2]
            cnames[i] <- strain_g
        } 

    #丰度表转置
        abun <- as.data.frame(t(as.data.frame(abundance)))
        abun <- abun[-1,] #去掉转置后第一行的各种微生物分类名
        abun[1:ncol(abun)] <- lapply(abun[1:ncol(abun)],as.numeric) #将reads数变为数值型
        abun <- cbind(as.data.frame(rownames(abun)),abun)
        rownames(abun) <- NULL
        colnames(abun) <- cnames #将列名改为各种菌种
    
    #画出来菌种和代谢物的相关性热图
        #计算相关性矩阵（可选：”pearson”、”spearman”、”kendall”相关系数）、p值矩阵
            cor <- corr.test(abun[,2:ncol(abun)], metabolism[,2:(ncol(metabolism)-1)], method = "pearson",adjust="none")
            cmt <- cor$r #相关性矩阵
                cmt[is.na(cmt)] <- 0
            pmt <- cor$p #p值矩阵
                pmt[is.na(pmt)] <- 1
        #给微生物和代谢物分别进行聚类 
            row <- hclust(dist(cmt)) #对行进行聚类
            col <- hclust(dist(t(cmt))) #对列进行聚类
        #按照聚类结果进行排序
            cmt <- cmt[,col$order]
            cmt <- cmt[row$order,]
            pmt <- pmt[,col$order]
            pmt <- pmt[row$order,]

        #合并成长数据表
            df_cmt <- melt(cmt,value.name="cor")
            df_pmt <- melt(pmt,value.name="p")
            df <- bind_cols(df_cmt,df_pmt)
            df <- df[,c(1,2,3,6)]
            colnames(df) <- c("strain","metabolism","cor","p")

        #把所有的代谢物名称中的"_"都替换成"
            levels(df$metabolism) <- gsub("_"," ",levels(df$metabolism))

        #对p值进行校正(因为p<0.05被认为是显著相关的，但是对于很多个实验，所有都被认为是显著相关则错误率较高，因此需要进行FDR矫正)
            #df$p <- p.adjust(df$p,method = "fdr")

        #p值大于0.05的取"",p值小于0.05大于0.01的取"*",p值小于0.01的取"**",小于0.001的取"***"
            n_sig <- df$p > 0.05
            sig_1 <- df$p <= 0.05 & df$p > 0.01
            sig_2 <- df$p <= 0.01 & df$p > 0.001
            sig_3 <- df$p <= 0.001
            df$p[which(n_sig == T)] <- ""
            df$p[which(sig_1 == T)] <- "*" 
            df$p[which(sig_2 == T)] <- "**"
            df$p[which(sig_3 == T)] <- "***"
            df$cor[is.na(df$cor)] <- 0
            df$p[is.na(df$p)] <- ""

        #画图
            heatmap <-
                ggplot(df,aes(x =strain,y = metabolism ,fill = cor)) +
                geom_tile() + #选择画热图
                geom_text(aes(label = p), color = "black", size = 5) + #添加显著性标记
                #coord_fixed()+
                scale_fill_gradient2(low = "#407EB1", high = "#D72E2A",limits = c(-0.95, 0.95),breaks = c(-1,-0.5,0,0.5,1)) + ###设置热图标尺的颜色、范围以及刻度###
                theme(axis.title.x = element_text(size = 17,colour = "black"),
                axis.text.x = element_text(size = 12,colour = "black",angle =90,vjust = 0.5),
                axis.title.y = element_text(size = 17,colour = "black"),
                axis.text.y = element_text(size = 12,colour = "black"),
                plot.title = element_text(hjust = 0.5,size = 20,colour = "black"),
                legend.title = element_blank(),
                legend.text = element_text(size = 12,colour = "black")) +
                xlab("")+
                ylab("")+
            labs(title = "")
            ggsave("/home/wangnan/Shandong-DF/figure/GBTM_group/heatmap/heatmap_response_Later(前40个菌种).pdf",heatmap,width = 15,height =7.5)
            # ggsave("/home/wangnan/Shandong-DF/figure/GBTM_group/heatmap/heatmap_non_response_Later(前40个菌种).pdf",heatmap,width = 15,height =7.5)

# 按照膳食干预的Early,Mid,Later阶段和后55个菌种分别进行pearson相关性分析
    # 响应组
        #处理代谢物
            metabolism <- metadata_response[,c(2,23:50,52)]   
                pos <- which(metadata_response$Intervention == "Later")
                metabolism <- metabolism[pos,]
            abundance <- abundance_response[,c(1,pos+1)]   
                abundance <- abundance[41:95,]

    # 非响应组
        #处理代谢物
            metabolism <- metadata_non_response[,c(2,23:50,52)]   
                pos <- which(metadata_non_response$Intervention == "Later")
                metabolism <- metabolism[pos,]
            abundance <- abundance_non_response[,c(1,pos+1)]   
                abundance <- abundance[41:95,]

    #提取列名(即各种菌种，并且对其处理，变成pylum+genus的形式)
        cnames <- c("SampleID",as.data.frame(abundance[,1])[,1]) 
        for(i in 2:length(cnames))
        {
            strain_g <- strsplit(cnames[i],";",fixed = T)[[1]][6]
                strain_g <- strsplit(strain_g,"__",fixed =T)[[1]][2]
            cnames[i] <- strain_g
        } 

    #丰度表转置
        abun <- as.data.frame(t(as.data.frame(abundance)))
        abun <- abun[-1,] #去掉转置后第一行的各种微生物分类名
        abun[1:ncol(abun)] <- lapply(abun[1:ncol(abun)],as.numeric) #将reads数变为数值型
        abun <- cbind(as.data.frame(rownames(abun)),abun)
        rownames(abun) <- NULL
        colnames(abun) <- cnames #将列名改为各种菌种
    
    #画出来菌种和代谢物的相关性热图
        #计算相关性矩阵（可选：”pearson”、”spearman”、”kendall”相关系数）、p值矩阵
            cor <- corr.test(abun[,2:ncol(abun)], metabolism[,2:(ncol(metabolism)-1)], method = "pearson",adjust="none")
            cmt <- cor$r #相关性矩阵
                cmt[is.na(cmt)] <- 0
            pmt <- cor$p #p值矩阵
                pmt[is.na(pmt)] <- 1
        #给微生物和代谢物分别进行聚类 
            row <- hclust(dist(cmt)) #对行进行聚类
            col <- hclust(dist(t(cmt))) #对列进行聚类
        #按照聚类结果进行排序
            cmt <- cmt[,col$order]
            cmt <- cmt[row$order,]
            pmt <- pmt[,col$order]
            pmt <- pmt[row$order,]

        #合并成长数据表
            df_cmt <- melt(cmt,value.name="cor")
            df_pmt <- melt(pmt,value.name="p")
            df <- bind_cols(df_cmt,df_pmt)
            df <- df[,c(1,2,3,6)]
            colnames(df) <- c("strain","metabolism","cor","p")

        #把所有的代谢物名称中的"_"都替换成"
            levels(df$metabolism) <- gsub("_"," ",levels(df$metabolism))

        #对p值进行校正(因为p<0.05被认为是显著相关的，但是对于很多个实验，所有都被认为是显著相关则错误率较高，因此需要进行FDR矫正)
            #df$p <- p.adjust(df$p,method = "fdr")

        #p值大于0.05的取"",p值小于0.05大于0.01的取"*",p值小于0.01的取"**",小于0.001的取"***"
            n_sig <- df$p > 0.05
            sig_1 <- df$p <= 0.05 & df$p > 0.01
            sig_2 <- df$p <= 0.01 & df$p > 0.001
            sig_3 <- df$p <= 0.001
            df$p[which(n_sig == T)] <- ""
            df$p[which(sig_1 == T)] <- "*" 
            df$p[which(sig_2 == T)] <- "**"
            df$p[which(sig_3 == T)] <- "***"
            df$cor[is.na(df$cor)] <- 0
            df$p[is.na(df$p)] <- ""

        #画图
            heatmap <-
                ggplot(df,aes(x =strain,y = metabolism ,fill = cor)) +
                geom_tile() + #选择画热图
                geom_text(aes(label = p), color = "black", size = 5) + #添加显著性标记
                #coord_fixed()+
                scale_fill_gradient2(low = "#407EB1", high = "#D72E2A",limits = c(-0.9, 0.9),breaks = c(-1,-0.5,0,0.5,1)) + ###设置热图标尺的颜色、范围以及刻度###
                theme(axis.title.x = element_text(size = 17,colour = "black"),
                axis.text.x = element_text(size = 12,colour = "black",angle =90,vjust = 0.5),
                axis.title.y = element_text(size = 17,colour = "black"),
                axis.text.y = element_text(size = 12,colour = "black"),
                plot.title = element_text(hjust = 0.5,size = 20,colour = "black"),
                legend.title = element_blank(),
                legend.text = element_text(size = 12,colour = "black")) +
                xlab("")+
                ylab("")+
            labs(title = "")
            # ggsave("/home/wangnan/Shandong-DF/figure/GBTM_group/heatmap/heatmap_response_Later(后55个菌种).pdf",heatmap,width = 15,height =7.5)
            ggsave("/home/wangnan/Shandong-DF/figure/GBTM_group/heatmap/heatmap_non_response_Later(后55个菌种).pdf",heatmap,width = 15,height =7.5)
