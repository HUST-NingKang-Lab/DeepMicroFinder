library(tidyverse)
library(ggplot2)
library(ggsci)
library(vegan)
library(dplyr)
library(openxlsx)

# 清除所有变量并释放内存
    rm(list = ls())
    gc()

# 筛选GGMP,SGMP和SD-DFI都存在的marker
    # GGMP的genus水平的丰度表
        data_ggmp <- read_tsv("/data4/wangnan/Shandong-DF/origin-data/ggmp/table-filtered-feature-rarefied5k-L6.tsv")
    # SGMP的genus水平的丰度表
        data_sgmp <- read_tsv("/data4/wangnan/Shandong-DF/origin-data/sgmp/table-filtered-feature-rarefied5k-L6.tsv")
    # SD-DFI的丰度表
        # data_dfi <- read_tsv("/data4/wangnan/Shandong-DF/origin-data/sd-DFI/abundance_processed.tsv") 2023/07/12重新合并ASV ID之前的丰度表，现已更正
        data_dfi <- read_tsv("/data4/wangnan/Shandong-DF/origin-data/sd-DFI/sd-dfi_abundance_processed.tsv")
    # 干预之前的丰度表
        testdata <- read_tsv("/data4/wangnan/Shandong-DF/prediction/Guangdong-Shandong/exp1/exp1/testdata.tsv")
    # 筛选出来三个队列都存在的marker
        tmp <- intersect(colnames(data_dfi)[2:ncol(data_dfi)],colnames(data_sgmp)[2:ncol(data_sgmp)])
        tmp <- intersect(tmp,colnames(data_ggmp)[2:ncol(data_ggmp)])    
        data_dfi <- data_dfi[,c(1,match(tmp,colnames(data_dfi)))]
        marker <- tmp
        pos <- match(marker,as.data.frame(testdata[,1])[,1])
    # 读取没有丰度的marker
        no_abundance_marker <- read.csv("/home/wangnan/Shandong-DF/table/no_abundance_marker.csv")

    # 筛选四个模型分别找到的marker
        #     #膳食纤维干预之后的数据
        #         abundance <- read_tsv("/data4/wangnan/Shandong-DF/origin-data/sd-DFI/abundance_processed.tsv")
        #     #干预之前的
        #         testdata <- read_tsv("/data4/wangnan/Shandong-DF/prediction/Guangdong-Shandong/exp1/exp1/testdata.tsv")
        #     #筛选都存在的marker
        #         marker <- as.data.frame(testdata[,1])
        #         marker <- intersect(marker$`Sample ID`,colnames(abundance)[2:ncol(abundance)])
        #         pos <- match(marker,as.data.frame(testdata[,1])[,1])
    
    #Independent DNN model
        sta_ind <- data.frame(Number =1, Marker= "1",Exp_1 = 0,Exp_2 = 0,Exp_3 = 0,Exp_4 = 0,Exp_5 = 0,Mean_of_absolute_value = 0,Standard_Deviation = 0)
        flag = 0 
        for(i in pos)
        {   
            flag = flag + 1
            tmp <- c(flag,marker[flag])
            for(n in 1:5)
            {
                auc_pre <- read.csv(paste0("/data4/wangnan/Shandong-DF/prediction/Shandong/exp4/exp",n,"/Evaluation_Independent/overall.csv"))
                auc_now <- read.csv(paste0("/data4/wangnan/Shandong-DF/select_marker/exp",n,"/exp",i,"/Evaluation_Independent/overall.csv"))
                tmp <- c(tmp,round(c(mean(auc_pre$ROC.AUC)-mean(auc_now$ROC.AUC)),3))
            }
            tmp <- c(tmp, round(mean(abs(c(as.numeric(tmp[3:7])))),3))
            tmp <- c(tmp,round(sd(c(as.numeric(tmp[3:7]))),3))
            sta_ind[flag+1,] <- tmp
        }
        sta_ind <- sta_ind[-1,]
        sta_ind <- sta_ind[,-1]
        sta_ind <- sta_ind[setdiff(1:nrow(sta_ind),match(no_abundance_marker[,2],sta_ind$Marker)),]
        sta_ind <- sta_ind[order(sta_ind$Mean_of_absolute_value,decreasing = T),]
        write.csv(sta_ind,"/home/wangnan/Shandong-DF/table/markers selected by independent DNN model.csv",row.names = F)
        # tmp1 <- read_tsv(paste0("/data4/wangnan/Shandong-DF/prediction/Shandong/exp4/exp",n,"/trainmapper.tsv"))
        # tmp2 <- read_tsv(paste0("/data4/wangnan/Shandong-DF/select_marker/exp",n,"/exp",i,"/trainmapper_sgmp.tsv"))

    #Regional DNN model
        sta_reg <- data.frame(Number =1, Marker= "1",Exp_1 = 0,Exp_2 = 0,Exp_3 = 0,Exp_4 = 0,Exp_5 = 0,Mean_of_absolute_value = 0,Standard_Deviation = 0)
        flag = 0 
        for(i in pos)
        {   
            flag = flag + 1
            tmp <- c(flag,marker[flag])
            for(n in 1:5)
            {
                auc_pre <- read.csv(paste0("/data4/wangnan/Shandong-DF/prediction/Guangdong-Shandong/exp4/exp",n,"/Evaluation_Independent_Guangdong/overall.csv"))
                auc_now <- read.csv(paste0("/data4/wangnan/Shandong-DF/select_marker/exp",n,"/exp",i,"/Evaluation_Regional/overall.csv"))
                tmp <- c(tmp,round(c(mean(auc_pre$ROC.AUC)-mean(auc_now$ROC.AUC)),3))
            }
            tmp <- c(tmp, round(mean(abs(c(as.numeric(tmp[3:7])))),3))
            tmp <- c(tmp,round(sd(c(as.numeric(tmp[3:7]))),3))
            sta_reg[flag+1,] <- tmp
        }
        sta_reg <- sta_reg[-1,]
        sta_reg <- sta_reg[,-1]
        sta_reg <- sta_reg[setdiff(1:nrow(sta_reg),match(no_abundance_marker[,2],sta_reg$Marker)),]
        sta_reg <- sta_reg[order(sta_reg$Mean_of_absolute_value,decreasing = T),]
        write.csv(sta_reg,"/home/wangnan/Shandong-DF/table/markers selected by regional DNN model.csv",row.names = F)

    #Regional+ DNN model
        sta_reg_extra <- data.frame(Number =1, Marker= "1",Exp_1 = 0,Exp_2 = 0,Exp_3 = 0,Exp_4 = 0,Exp_5 = 0,Mean_of_absolute_value = 0,Standard_Deviation = 0)
        flag = 0 
        for(i in pos)
        {   
            flag = flag + 1
            tmp <- c(flag,marker[flag])
            for(n in 1:5)
            {
                auc_pre <- read.csv(paste0("/data4/wangnan/Shandong-DF/prediction/Guangdong-Shandong/exp4/exp",n,"/Evaluation_Independent/overall.csv"))
                auc_now <- read.csv(paste0("/data4/wangnan/Shandong-DF/select_marker/exp",n,"/exp",i,"/Evaluation_Regional_extra/overall.csv"))
                tmp <- c(tmp,round(c(mean(auc_pre$ROC.AUC)-mean(auc_now$ROC.AUC)),3))
            }
            tmp <- c(tmp, round(mean(abs(c(as.numeric(tmp[3:7])))),3))
            tmp <- c(tmp,round(sd(c(as.numeric(tmp[3:7]))),3))
            sta_reg_extra[flag+1,] <- tmp
        }
        sta_reg_extra <- sta_reg_extra[-1,]
        sta_reg_extra <- sta_reg_extra[,-1]
        sta_reg_extra <- sta_reg_extra[setdiff(1:nrow(sta_reg_extra),match(no_abundance_marker[,2],sta_reg_extra$Marker)),]
        sta_reg_extra <- sta_reg_extra[order(sta_reg_extra$Mean_of_absolute_value,decreasing = T),]
        write.csv(sta_reg_extra,"/home/wangnan/Shandong-DF/table/markers selected by regional+ DNN model.csv",row.names = F)

    #Transfer DNN model
        sta_tra <- data.frame(Number =1, Marker= "1",Exp_1 = 0,Exp_2 = 0,Exp_3 = 0,Exp_4 = 0,Exp_5 = 0,Mean_of_absolute_value = 0,Standard_Deviation = 0)
        flag = 0 
        for(i in pos)
        {   
            flag = flag + 1
            tmp <- c(flag,marker[flag])
            for(n in 1:5)
            {
                auc_pre <- read.csv(paste0("/data4/wangnan/Shandong-DF/prediction/Guangdong-Shandong/exp4/exp",n,"/Evaluation_Transfer/overall.csv"))
                auc_now <- read.csv(paste0("/data4/wangnan/Shandong-DF/select_marker/exp",n,"/exp",i,"/Evaluation_Transfer/overall.csv"))
                tmp <- c(tmp,round(c(mean(auc_pre$ROC.AUC)-mean(auc_now$ROC.AUC)),3))
            }
            tmp <- c(tmp, round(mean(abs(c(as.numeric(tmp[3:7])))),3))
            tmp <- c(tmp,round(sd(c(as.numeric(tmp[3:7]))),3))
            sta_tra[flag+1,] <- tmp
        }
        sta_tra <- sta_tra[-1,]
        sta_tra <- sta_tra[,-1]
        sta_tra <- sta_tra[setdiff(1:nrow(sta_tra),match(no_abundance_marker[,2],sta_tra$Marker)),]
        sta_tra <- sta_tra[order(sta_tra$Mean_of_absolute_value,decreasing = T),]
        write.csv(sta_tra,"/home/wangnan/Shandong-DF/table/markers selected by transfer DNN model.csv",row.names = F)


# 统计梯度选择marker进行建模预测的AUROC
    # 统计
      sta <- data.frame(markers_num = "0",exp = "1",model = "1",AUROC ="1",stringsAsFactors =FALSE)
        for(i in c("independent","regional","regional_extra","transfer"))
        {
            for(j in c(1:19)*5 )
            {
                for(n in 1:10)
                {
                    auc <-read.csv(paste0("/data4/wangnan/Shandong-DF/prediction_select_marker_95_features_GBTM-grouping/HbA1c/",i,"/",j,"/exp",n,"/Evaluation_Independent/overall.csv"))
                    if(j != 42)
                    {
                    sta <- rbind(sta,c(j,paste0("exp",n),i,mean(auc$ROC.AUC)))
                    }
                    # else
                    # {
                    #     sta <- rbind(sta,c(paste0("exp",n),"",mean(auc$ROC.AUC)))
                    # }
                }
            }
        }
        sta <- sta[-1,]
        sta$AUROC <- as.numeric(sta$AUROC)
        sta$model <- factor(sta$model,levels = c("independent","regional","regional_extra","transfer"))
        sta$markers_num <- factor(sta$markers_num,levels = c(1:19)*5)
        # comparison <- list(c("transfer","All 59 markers"),c("regional_extra","All 59 markers"),c("regional","All 59 markers"),c("independent","All 59 markers"))

    # 折线图
        q <- ggplot(sta, aes(x = markers_num,y = AUROC,fill = model)) +
            geom_boxplot(position=position_dodge(0.7),width = 0.3) +
            # geom_point() +
            scale_color_manual(values=c("#E64B35", "#4DBBD5", "#00A087","#3C5488"),aesthetics = "color") +  #手动调整线的颜色 
            stat_summary(fun=mean, geom="line", aes(colour = model,group = model)) + #添加均值线
            stat_summary(fun=mean, geom="point", aes(colour = model,group = model)) + #添加均值点
            # scale_color_manual(values=c("#E64B35", "#4DBBD5", "#00A087"),aesthetics = "fill")+ #箱线图的外框指定颜色
            # scale_color_manual(values=c("#000000", "#000000", "#000000"),aesthetics = "color")+ #箱线图的填充指定颜色
            theme_bw() +
            theme(
                axis.text.x = element_text(size = 15,colour = "black"),
                axis.text.y = element_text(size = 15,colour = "black"),
                axis.title.x = element_text(size = 15,colour = "black"),
                axis.title.y = element_text(size = 15,colour = "black"),
                legend.text = element_text(size = 15,colour = "black"),
                legend.title = element_text(size = 15,colour = "black"),
                # legend.position = "none",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                plot.title = element_text(hjust = 0.5,size = 15,colour = "black")) +
                scale_fill_npg()+
                # geom_signif(comparisons = comparison,step_increase = 0.2,map_signif_level = F,test = wilcox.test)+
                ylim(c(0.3,0.7))+
                # stat_compare_means(label = "p.format")+
                # facet_grid(~model,scales = "free")+
                xlab("Number of markers")+
                ylab("AUROC")

        ggsave("/home/wangnan/Shandong-DF/figure/AUROC_assessment_selected_markers_grouping(HbA1c)(line chart).pdf",q,height = 6,width = 12)

    # 拟合曲线
        q <- ggplot(sta, aes(x = markers_num,y = AUROC,fill = model)) +
            geom_boxplot(position=position_dodge(0.7),width = 0.3) +
            geom_smooth(aes(colour = model,group = model)) +
            # geom_point() +
            scale_color_manual(values=c("#E64B35", "#4DBBD5", "#00A087","#3C5488"),aesthetics = "color") +  #手动调整线的颜色 
            # stat_summary(fun=mean, geom="line", aes(colour = model,group = model)) + #添加均值线
            # stat_summary(fun=mean, geom="point", aes(colour = model,group = model)) + #添加均值点
            # scale_color_manual(values=c("#E64B35", "#4DBBD5", "#00A087"),aesthetics = "fill")+ #箱线图的外框指定颜色
            # scale_color_manual(values=c("#000000", "#000000", "#000000"),aesthetics = "color")+ #箱线图的填充指定颜色
            theme_classic() +
            theme(axis.text.x = element_text(size = 10,colour = "black"),
                axis.text.y = element_text(size = 15,colour = "black"),
                axis.title.x = element_text(size = 15,colour = "black"),
                axis.title.y = element_text(size = 15,colour = "black"),
                legend.text = element_text(size = 15,colour = "black"),
                legend.title = element_text(size = 15,colour = "black"),
                legend.position = "none",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                plot.title = element_text(hjust = 0.5,size = 15,colour = "black")) +
                scale_fill_npg()+
                # geom_signif(comparisons = comparison,step_increase = 0.2,map_signif_level = F,test = wilcox.test)+
                ylim(c(0.3,0.7))+
                # stat_compare_means(label = "p.format")+
                facet_grid(~model,scales = "free")+
            theme(strip.text.x = element_text(size = 15, face = "bold")) +
                xlab("Number of markers")+
                ylab("AUROC")

        ggsave("/home/wangnan/Shandong-DF/figure/AUROC_assessment_selected_markers_grouping(HbA1c)(fitted line).pdf",q,height = 6,width = 15)


# 因为grouping分组的结果都很差，怀疑干预早期患者的肠道微生物组成几乎没有差异，因此画一个PCoA图看一下
    #读取利用GBTM模型进行分类的结果(HbA1c/non-HbA1c)
        group <- read.xlsx("/data4/wangnan/Shandong-DF/origin-data/sd-DFI/Combine_groups.xlsx")    
    #metadata
        metadata <- read.csv("/data4/wangnan/Shandong-DF/origin-data/sd-DFI/DFI_metadata_250_used.csv")
    #只筛选出来膳食纤维干预早期的样本(因为膳食纤维干预后期，糖尿病症状已经不是很明显，所以可能没有很好的分类效果)
        metadata <- metadata[metadata$Intervention == "Early",]
    #丰度表
        abundance <- read_tsv("/data4/wangnan/Shandong-DF/origin-data/sd-DFI/sd-dfi_abundance_processed.tsv")
        #筛选所使用到的85个样本中的膳食纤维干预早期的样本的丰度'
            id <- intersect(abundance$`Sample ID`,metadata$SampleID)
            pos <- match(id,as.data.frame(abundance$`Sample ID`)[,1])
            abundance <- abundance[pos,]
            abundance[,2:ncol(abundance)] <- lapply(abundance[,2:ncol(abundance)],as.numeric)
        # 筛选出所使用到的95个marker
            sta_tra <- read.csv("/home/wangnan/Shandong-DF/table/markers selected by transfer DNN model.csv")
            marker <- sta_tra$Marker
            abundance <- abundance[,c(1,match(marker,colnames(abundance)))]
    # 存储样本的分组信息
        symbol <- "HbA1c"
        grouping <- c()
        for(i in 1:nrow(abundance)) 
        {
            patient_id <- metadata$Patient_id[which(metadata$SampleID == abundance$`Sample ID`[i])]
            pos <- which(colnames(group) == symbol)
            grouping <- c(grouping,group[,pos][which(group$Patient_id == patient_id)])
        }
        # 把标签进行替换
            grouping <- gsub(1,symbol,grouping)
            grouping <- gsub(0,paste0("non-",symbol),grouping)
        # 去除第一列的样本名
            abundance <- abundance[,-1]
        
    # 定义画PCoA的函数
        PCoA <- function(matrix,group){
            # 计算Bray-Curtis距离
                distance <- vegdist(matrix, method = 'bray')
            # 计算pcoa(计算的距离矩阵缺少一行，所以要减去样本名和额外一行)
                pcoa <- cmdscale(distance, k = nrow(matrix)-1, eig = TRUE)
            # 提取样本点坐标（points记录了各样本在各排序轴中的坐标值）
                #前两轴
                plot_data <- data.frame({pcoa$points})[1:2]
                names(plot_data)[1:2] <- c('PCoA1', 'PCoA2') #命名为PCoA1和PCoA2
            # eig记录了PCoA排序结果中，主要排序轴的特征值（再除以特征值总和就是各轴的解释量）
                eig = pcoa$eig
            # 为样本点坐标添加分组信息
                plot_data$group <- group
            
            # 画图
                p <- ggplot(plot_data,aes(x=PCoA1,y=PCoA2,color=group)) + 
                     geom_point(size = 0.5)+
                     stat_ellipse(level = 0.95, show.legend = F)+ #添加置信椭圆
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
                        legend.title = element_blank(),
                        axis.line= element_line(colour = "black"),
                        plot.title = element_text(hjust = 0.5,size = 20,colour = "black"))+
                        labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
                            y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
                    scale_color_npg()
            # 返回结果
                return(p)
        }

    # 画图
        p <- PCoA(abundance,grouping)
        ggsave(paste0("/home/wangnan/Shandong-DF/figure/PCoA(",symbol,").pdf"),p,height = 5,width = 7)


# 按照干预阶段画出来PCoA图   
    #metadata
        metadata <- read.csv("/data4/wangnan/Shandong-DF/origin-data/sd-DFI/DFI_metadata_250_used.csv")
    #丰度表
        abundance <- read_tsv("/data4/wangnan/Shandong-DF/origin-data/sd-DFI/sd-dfi_abundance_processed.tsv")
        #筛选所使用到的85个样本中的膳食纤维干预早期的样本的丰度'
            id <- intersect(abundance$`Sample ID`,metadata$SampleID)
            pos <- match(id,as.data.frame(abundance$`Sample ID`)[,1])
            abundance <- abundance[pos,]
            abundance[,2:ncol(abundance)] <- lapply(abundance[,2:ncol(abundance)],as.numeric)
        # 筛选出所使用到的95个marker
            sta_tra <- read.csv("/home/wangnan/Shandong-DF/table/markers selected by transfer DNN model.csv")
            marker <- sta_tra$Marker
            abundance <- abundance[,c(1,match(marker,colnames(abundance)))]
    # 存储样本的分组信息
        grouping <- c()
        for(i in 1:nrow(abundance)) 
        {
            grouping <- c(grouping,metadata$Intervention[which(metadata$SampleID == abundance$`Sample ID`[i])])
        }
        # 去除第一列的样本名
            abundance <- abundance[,-1]
        
    # 定义画PCoA的函数
        PCoA <- function(matrix,group){
            # 计算Bray-Curtis距离
                distance <- vegdist(matrix, method = 'bray')
            # 计算pcoa(计算的距离矩阵缺少一行，所以要减去样本名和额外一行)
                pcoa <- cmdscale(distance, k = nrow(matrix)-1, eig = TRUE)
            # 提取样本点坐标（points记录了各样本在各排序轴中的坐标值）
                #前两轴
                plot_data <- data.frame({pcoa$points})[1:2]
                names(plot_data)[1:2] <- c('PCoA1', 'PCoA2') #命名为PCoA1和PCoA2
            # eig记录了PCoA排序结果中，主要排序轴的特征值（再除以特征值总和就是各轴的解释量）
                eig = pcoa$eig
            # 为样本点坐标添加分组信息
                plot_data$group <- group
            
            # 画图
                p <- ggplot(plot_data,aes(x=PCoA1,y=PCoA2,color=group)) + 
                     geom_point(size = 0.5)+
                     stat_ellipse(level = 0.95, show.legend = F)+ #添加置信椭圆
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
                        legend.title = element_blank(),
                        axis.line= element_line(colour = "black"),
                        plot.title = element_text(hjust = 0.5,size = 20,colour = "black"))+
                        labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
                            y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
                    scale_color_npg()
            # 返回结果
                return(p)
        }

    # 画图
        p <- PCoA(abundance,grouping)
        ggsave("/home/wangnan/Shandong-DF/figure/PCoA(Intervention stage).pdf",p,height = 5,width = 7)
