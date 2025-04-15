library(tidyverse)
library(openxlsx)
library(psych) #psych包用于计算相关性、p值等信息
library(ggplot2)

# 清空之前的环境
    rm(list=ls())
    gc()

# 对SGMP的207个marker进行分类，包括和T2D无关的regional-shared，regional-specific，和疾病有关的regional-shared，regional-specific

    # 需要先将SGMP的marker和T2D做一个相关性分析，定量每个marker是否和T2D相关
        # SGMP的genus水平的丰度表
            # 这个是已经提前做好筛选的丰度表，所以不需要重新筛选有metadata信息的丰度表样本
            data_sgmp <- read_tsv("/data4/wangnan/Shandong-DF/origin-data/sgmp/filter_data.tsv")
        # 读取SGMP的metadata
            metadata_sgmp <- read.xlsx("/data4/wangnan/Shandong-DF/origin-data/sgmp/SGMP_metadata_700_selected.xlsx")
            abun <- data_sgmp
            meta <- metadata_sgmp

    # 将丰度表转为相对丰度表，用于计算后续的相关性
        for(i in 2:ncol(abun))
        {                
            sum <- sum(abun[,i]) 
            for(j in 1:nrow(abun))
            {     
                if(sum > 0)
                {                
                    abun[j,i] <- abun[j,i]/sum
                }
            } 
        }

    # 根据丰度信息筛选和疾病相关和不相关的菌种
        # 丰度表转置,改为行名为样本名，列名为各种菌属   
            cname <- c("SampleID",as.data.frame(abun[,1])[,1])
            abun <- as.data.frame(t(as.data.frame(abun)))
            abun <- abun[-1,] #去掉转置后第一行的各种微生物分类名
            abun[1:ncol(abun)] <- lapply(abun[1:ncol(abun)],as.numeric) # 将reads数变为数值型
            abun <- cbind(as.data.frame(rownames(abun)),abun)
            rownames(abun) <- NULL
            colnames(abun) <- cname # 将列名改为各种菌种
        
        # 计算相关性矩阵（可选：”pearson”、”spearman”、”kendall”相关系数）、p值矩阵
            # spearman相关性必须通过数值进行计算，所以需要把host_status转换成数值类型
                meta[,3] <- gsub("Type 2 diabetes","0",meta[,3])
                meta[,3] <- gsub("Health","1",meta[,3])
                meta[,3] <- as.numeric(meta[,3])

            # 计算相关性矩阵和p值矩阵
                cor <- corr.test(abun[,2:ncol(abun)], meta[,3], method = "spearman",adjust="none")
                cmt <- as.data.frame(cor$r) # 相关性矩阵
                    cmt[is.na(cmt)] <- 0
                pmt <- cor$p # p值矩阵
                    pmt[is.na(pmt)] <- 0
            
            #对p值进行校正(因为p<0.05被认为是显著相关的，但是对于很多个实验，所有都被认为是显著相关则错误率较高，因此需要进行FDR矫正)
                pmt[,1] <- p.adjust(pmt[,1],method = "fdr")


    # 通过图像判断应该划分的相关性的组别的阈值
        tmp <- cbind(cmt,1)
        colnames(tmp) <- c("cor","group")
        p <- ggplot(tmp,aes(x=group,y=cor)) +
            geom_point(color = "black")

        # 通过图像划分：cor>0.1 & p<0.05为正相关菌种，cor<-0.1 & p<0.05为负相关菌种，其余的为不明显相关菌种。
            related_pos <- c(intersect(which(cmt[,1] > 0.1),which(as.data.frame(pmt)[,1] < 0.05)),
                            intersect(which(cmt[,1] < -0.1),which(as.data.frame(pmt)[,1] < 0.05)))
            strain <- data.frame(strain = rownames(cmt), correlation = cmt[,1], p.value = as.data.frame(pmt)[,1],stringsAsFactors = FALSE)
            related_strain <- strain[related_pos,]
            unrelated_strain <- strain[-related_pos,]
            write.csv(related_strain,"/home/wangnan/Shandong-DF/table/feature_division/T2D related strain of SGMP.csv",row.names = FALSE)
            write.csv(unrelated_strain,"/home/wangnan/Shandong-DF/table/feature_division/T2D unrelated strain of SGMP.csv",row.names = FALSE)

# 和leave-one-out筛选的菌种重要性进行对比，分别筛选disease-related/regional-specific, disease-related/regional-shared
    sta_tra <- data.frame(Number =1, Marker= "1",Exp_1 = 0,Exp_2 = 0,Exp_3 = 0,Exp_4 = 0,Exp_5 = 0,Mean_of_absolute_value = 0,Standard_Deviation = 0)
    flag = 0 
    marker <- as.data.frame(data_sgmp[,1])
    for(i in 1:207)
    {               
        tmp <- c(i,marker[i,1])
        for(n in 1:5)
        {
            auc_pre <- read.csv(paste0("/data4/wangnan/Shandong-DF/prediction/Guangdong-Shandong/exp4/exp",n,"/Evaluation_Transfer/overall.csv"))
            auc_now <- read.csv(paste0("/data4/wangnan/Shandong-DF/select_marker/exp",n,"/exp",i,"/Evaluation_Transfer/overall.csv"))
            tmp <- c(tmp,round(c(mean(auc_pre$ROC.AUC)-mean(auc_now$ROC.AUC)),3))
        }
        tmp <- c(tmp, round(mean(abs(c(as.numeric(tmp[3:7])))),3))
        tmp <- c(tmp,round(sd(c(as.numeric(tmp[3:7]))),3))
        sta_tra[i+1,] <- tmp
    }
    sta_tra <- sta_tra[-1,]
    sta_tra <- sta_tra[order(sta_tra$Mean_of_absolute_value,decreasing = TRUE),]
    write.csv(sta_tra,"/home/wangnan/Shandong-DF/table/feature_division/selected feature by transfer model.csv",row.names = FALSE)

    # 重要性的筛选依据是Mean_of_absolute_value,值比较高的认为是在地域差异较大所以是regional-specific，反之则是regional-shared.
        #前20种菌被认为是regional_specific，后20种被认为是regioal_shared
        sta_tra <- read.csv("/home/wangnan/Shandong-DF/table/feature_division/selected feature by transfer model.csv")
        regional_specific_strain <- sta_tra[1:20,]        
        regional_shared_strain <- sta_tra[(nrow(sta_tra)-19):nrow(sta_tra),]

    # 和疾病相关的菌种进行匹配
        related_strain <- read.csv("/home/wangnan/Shandong-DF/table/feature_division/T2D related strain of SGMP.csv")
        unrelated_strain <- read.csv("/home/wangnan/Shandong-DF/table/feature_division/T2D unrelated strain of SGMP.csv")       
        # 和疾病有关的
            # regional-specific和disease-related
                pos1 <- match(intersect(regional_specific_strain$Marker,related_strain$strain),regional_specific_strain$Marker)
                pos2 <- match(intersect(regional_specific_strain$Marker,related_strain$strain),related_strain$strain)
                regional_specific_disease_related <- cbind(regional_specific_strain[pos1,],related_strain[pos2,2:3])
                write.csv(regional_specific_disease_related,"/home/wangnan/Shandong-DF/table/feature_division/regional-specific and disease-related strains.csv",row.names = FALSE)
            # regional-shared和disease-related
                pos1 <- match(intersect(regional_shared_strain$Marker,related_strain$strain),regional_shared_strain$Marker)
                pos2 <- match(intersect(regional_shared_strain$Marker,related_strain$strain),related_strain$strain)
                regional_shared_disease_related <- cbind(regional_shared_strain[pos1,],related_strain[pos2,2:3])
                write.csv(regional_shared_disease_related,"/home/wangnan/Shandong-DF/table/feature_division/regional-shared and disease-related strains.csv",row.names = FALSE)
        # 和疾病无关的
            # regional-specific和disease-unrelated
                pos1 <- match(intersect(regional_specific_strain$Marker,unrelated_strain$strain),regional_specific_strain$Marker)
                pos2 <- match(intersect(regional_specific_strain$Marker,unrelated_strain$strain),unrelated_strain$strain)
                regional_specific_disease_unrelated <- cbind(regional_specific_strain[pos1,],unrelated_strain[pos2,2:3])
                write.csv(regional_specific_disease_unrelated,"/home/wangnan/Shandong-DF/table/feature_division/regional-specific and disease-unrelated strains.csv",row.names = FALSE)
            # regional-shared和disease-unrelated
                pos1 <- match(intersect(regional_shared_strain$Marker,unrelated_strain$strain),regional_shared_strain$Marker)
                pos2 <- match(intersect(regional_shared_strain$Marker,unrelated_strain$strain),unrelated_strain$strain)
                regional_shared_disease_unrelated <- cbind(regional_shared_strain[pos1,],unrelated_strain[pos2,2:3])
                write.csv(regional_shared_disease_unrelated,"/home/wangnan/Shandong-DF/table/feature_division/regional-shared and disease-unrelated strains.csv",row.names = FALSE)
            



