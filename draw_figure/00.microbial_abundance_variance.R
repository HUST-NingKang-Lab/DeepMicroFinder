library(tidyverse)
library(openxlsx)
library(ggplot2)
library(reshape2)

# GGMP
    #取出GGMP的metadata中的T2D患者和健康人
        metadata_ggmp <- read.csv("/data4/wangnan/Shandong-DF/origin-data/ggmp/filter_metadata.csv")
    #读取phylum level的丰度表
        data_ggmp <- read_tsv("/data4/wangnan/Shandong-DF/origin-data/ggmp/table-filtered-feature-rarefied5k-L2.tsv")
        pos <- match(metadata_ggmp$SampleID,data_ggmp$`Sample ID`) #选出丰度表中的T2D患者和健康人
        data_ggmp <- data_ggmp[pos,]
        data_ggmp <- tibble::rownames_to_column(as.data.frame(t(data_ggmp),stringsAsFactors = FALSE),"V0") #行列转置
        colnames(data_ggmp) <- data_ggmp[1,]
        data_ggmp <- data_ggmp[-1,]
    #去掉k__Archaea;p__Euryarchaeota和k__Bacteria;__
        data_ggmp <- data_ggmp[-c(1,2),]

# SGMP
    #读取SGMP的metadata
        metadata_sgmp <- read.xlsx("/data4/wangnan/Shandong-DF/origin-data/sgmp/SGMP_metadata_700_selected.xlsx")
    #读取phylum level的丰度表
        data_sgmp <- read_tsv("/data4/wangnan/Shandong-DF/origin-data/sgmp/table-filtered-feature-rarefied5k-L2.tsv")
        pos <- match(metadata_sgmp$`SampleID`,data_sgmp$`Sample ID`)
        data_sgmp <- data_sgmp[pos,]
        data_sgmp <- tibble::rownames_to_column(as.data.frame(t(data_sgmp),stringsAsFactors = FALSE),"V0") #行列转置
        colnames(data_sgmp) <- data_sgmp[1,]
        data_sgmp <- data_sgmp[-1,]
    #去掉k__Bacteria;__
        data_sgmp <- data_sgmp[-1,]

# 将列表元素由character转换为numeric
    data_ggmp[,2:ncol(data_ggmp)] <- apply(data_ggmp[,2:ncol(data_ggmp)],2,as.numeric)
    data_sgmp[,2:ncol(data_sgmp)] <- apply(data_sgmp[,2:ncol(data_sgmp)],2,as.numeric)

# 更换phylum的名字
    data_ggmp[,1] <- gsub("k__Bacteria;p__","",data_ggmp[,1])
    data_sgmp[,1] <- gsub("k__Bacteria;p__","",data_sgmp[,1])

# 筛选出来GGMP和SGMP共有的phylum,GGMP剩余的phylum都被定义为Others
    library(dplyr)
    phylum_common <- intersect(data_ggmp[,1],data_sgmp[,1])
    pos_other <- setdiff(1:nrow(data_ggmp),match(phylum_common,data_ggmp[,1]))
    data_ggmp <- rbind(data_ggmp,c("Others",as.vector(colSums(data_ggmp[pos_other,2:ncol(data_ggmp)]))))
    data_ggmp <- data_ggmp[-pos_other,]
    data_ggmp[,2:ncol(data_ggmp)] <- apply(data_ggmp[,2:ncol(data_ggmp)],2,as.numeric)

# 计算出各phylum的相对丰度
    data_ggmp[,2:ncol(data_ggmp)] <- data_ggmp[,2:ncol(data_ggmp)]/colSums(data_ggmp[,2:ncol(data_ggmp)])
    data_sgmp[,2:ncol(data_sgmp)] <- data_sgmp[,2:ncol(data_sgmp)]/colSums(data_sgmp[,2:ncol(data_sgmp)])

# 不分T2D和control
    # 按照Bacteroidetes的丰度由高到低进行排序
        data_ggmp[,2:ncol(data_ggmp)] <- data_ggmp[,1+order(data_ggmp[2,2:ncol(data_ggmp)],decreasing = TRUE)]
        data_sgmp[,2:ncol(data_sgmp)] <- data_sgmp[,1+order(data_sgmp[2,2:ncol(data_sgmp)],decreasing = TRUE)]

    # 按照SGMP的phylum进行因子化，用于指定堆叠柱状图的phylum的顺序
        order <- data_sgmp[,1]
        data_ggmp[,1] <- factor(data_ggmp[,1],levels = c(order,"Others"))
        data_sgmp[,1] <- factor(data_sgmp[,1],levels = order)

    # 宽数据变为长数据
        data_ggmp <- data_ggmp %>% melt(id.vars = "Sample ID")
        colnames(data_ggmp) <- c("phylum","sampleid","value")

        data_sgmp <- data_sgmp %>% melt(id.vars = "Sample ID")
        colnames(data_sgmp) <- c("phylum","sampleid","value")
        
    # 绘制堆叠柱状图
        # 查看ggsci的帮助文档
            # vignette("ggsci")
        # 获取npg和nejm的二进制配色
            library(scales)
            library(ggsci)
            pal_npg(palette = c("nrc"), alpha = 1)(10)
            pal_nejm(palette = c("default"), alpha = 1)(8)
        # 展示选择的配色
            #show_col(c(pal_npg(palette = c("nrc"), alpha = 1)(10),pal_nejm(palette = c("default"), alpha = 1)(8)[7:8]))

        # 指定配色
            color <- c(pal_npg(palette = c("nrc"), alpha = 1)(10),pal_nejm(palette = c("default"), alpha = 1)(8)[7:8])[c(1,7,3,4,5,6,2,8,9,10,11,12)]

        # 画图
            # GGMP
                p_ggmp <- ggplot(data_ggmp,aes(sampleid,value,fill=phylum)) +
                    geom_bar(stat="identity", position = 'fill') +
                    xlab("") +
                    ylab("") +
                    #theme_classic(base_size = 7) +
                    scale_y_continuous(expand = c(0,0)) +
                    ggtitle('') +
                    guides(fill=guide_legend(title=NULL)) +
                    theme(axis.text.x=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.ticks.y = element_blank(),
                        legend.text = element_text(size = 12),
                        legend.position = "bottom") +
                    scale_fill_manual(values = c(color,pal_aaas(palette = c("default"), alpha = 1)(10)[9]))
                ggsave("/home/wangnan/Shandong-DF/figure/GGMP_stacked_bar.pdf",p_ggmp,height = 6,width = 15)

            # SGMP
                p_sgmp <- ggplot(data_sgmp,aes(sampleid,value,fill=phylum)) +
                    geom_bar(stat="identity", position = 'fill') +
                    xlab("") +
                    ylab("") +
                    #theme_classic(base_size = 7) +
                    scale_y_continuous(expand = c(0,0)) +
                    ggtitle('') +
                    guides(fill=guide_legend(title=NULL)) +
                    theme(axis.text.x=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.ticks.y = element_blank(),
                        legend.text = element_text(size = 12),
                        legend.position = "bottom") +
                    scale_fill_manual(values = color)
                ggsave("/home/wangnan/Shandong-DF/figure/SGMP_stacked_bar.pdf",p_sgmp,height = 6,width = 15)

# 分T2D和control
    # 筛选出来T2D和control 
        # GGMP
            pos_case <- which(metadata_ggmp$healthy == TRUE)
            pos_control <- which(metadata_ggmp$healthy == FALSE)
            data_ggmp_case <- data_ggmp[,c(1,pos_case+1)]
            data_ggmp_control <- data_ggmp[,c(1,pos_control+1)]
        # SGMP
            pos_case <- which(metadata_sgmp$host_status == "Type 2 diabetes")
            pos_control <- which(metadata_sgmp$host_status == "Health")
            data_sgmp_case <- data_sgmp[,c(1,pos_case+1)]
            data_sgmp_control <- data_sgmp[,c(1,pos_control+1)] 

    # 按照Bacteroidetes的丰度由高到低进行排序
        # GGMP
            data_ggmp_case[,2:ncol(data_ggmp_case)] <- data_ggmp_case[,1+order(data_ggmp_case[2,2:ncol(data_ggmp_case)],decreasing = TRUE)]
            data_ggmp_control[,2:ncol(data_ggmp_control)] <- data_ggmp_control[,1+order(data_ggmp_control[2,2:ncol(data_ggmp_control)],decreasing = TRUE)]
        # SGMP
            data_sgmp_case[,2:ncol(data_sgmp_case)] <- data_sgmp_case[,1+order(data_sgmp_case[2,2:ncol(data_sgmp_case)],decreasing = TRUE)]
            data_sgmp_control[,2:ncol(data_sgmp_control)] <- data_sgmp_control[,1+order(data_sgmp_control[2,2:ncol(data_sgmp_control)],decreasing = TRUE)]

    # 按照SGMP的phylum进行因子化，用于指定堆叠柱状图的phylum的顺序
        order <- data_sgmp_case[,1]
        # GGMP
            data_ggmp_case[,1] <- factor(data_ggmp_case[,1],levels = c(order,"Others"))
            data_ggmp_control[,1] <- factor(data_ggmp_control[,1],levels = c(order,"Others"))
        # SGMP
            data_sgmp_case[,1] <- factor(data_sgmp_case[,1],levels = order)
            data_sgmp_control[,1] <- factor(data_sgmp_control[,1],levels = order)

    # 宽数据变为长数据
        # GGMP
            data_ggmp_case <- data_ggmp_case %>% melt(id.vars = "Sample ID")
            colnames(data_ggmp_case) <- c("phylum","sampleid","value")
            data_ggmp_control <- data_ggmp_control %>% melt(id.vars = "Sample ID")
            colnames(data_ggmp_control) <- c("phylum","sampleid","value")            
        #SGMP
            data_sgmp_case <- data_sgmp_case %>% melt(id.vars = "Sample ID")
            colnames(data_sgmp_case) <- c("phylum","sampleid","value")
            data_sgmp_control <- data_sgmp_control %>% melt(id.vars = "Sample ID")
            colnames(data_sgmp_control) <- c("phylum","sampleid","value")            
        
    # 绘制堆叠柱状图
        # 查看ggsci的帮助文档
            # vignette("ggsci")
        # 获取npg和nejm的二进制配色
            library(scales)
            library(ggsci)
            pal_npg(palette = c("nrc"), alpha = 1)(10)
            pal_nejm(palette = c("default"), alpha = 1)(8)
        # 展示选择的配色
            #show_col(c(pal_npg(palette = c("nrc"), alpha = 1)(10),pal_nejm(palette = c("default"), alpha = 1)(8)[7:8]))

        # 指定配色
            color <- c(pal_npg(palette = c("nrc"), alpha = 1)(10),pal_nejm(palette = c("default"), alpha = 1)(8)[7:8])[c(1,7,3,4,5,6,2,8,9,10,11,12)]

        # 画图
            # GGMP
                p_ggmp_case <- ggplot(data_ggmp_case,aes(sampleid,value,fill=phylum)) +
                    geom_bar(stat="identity", position = 'fill') +
                    xlab("") +
                    ylab("") +
                    #theme_classic(base_size = 7) +
                    scale_y_continuous(expand = c(0,0)) +
                    ggtitle('') +
                    guides(fill=guide_legend(title=NULL)) +
                    theme(axis.text.x=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.ticks.y = element_blank(),
                        legend.text = element_text(size = 12),
                        legend.position = "bottom") +
                    scale_fill_manual(values = c(color,pal_aaas(palette = c("default"), alpha = 1)(10)[9]))
                ggsave("/home/wangnan/Shandong-DF/figure/GGMP_case_stacked_bar.pdf",p_ggmp_case,height = 6,width = 10)

                p_ggmp_control <- ggplot(data_ggmp_control,aes(sampleid,value,fill=phylum)) +
                    geom_bar(stat="identity", position = 'fill') +
                    xlab("") +
                    ylab("") +
                    #theme_classic(base_size = 7) +
                    scale_y_continuous(expand = c(0,0)) +
                    ggtitle('') +
                    guides(fill=guide_legend(title=NULL)) +
                    theme(axis.text.x=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.ticks.y = element_blank(),
                        legend.text = element_text(size = 12),
                        legend.position = "bottom") +
                    scale_fill_manual(values = c(color,pal_aaas(palette = c("default"), alpha = 1)(10)[9]))
                ggsave("/home/wangnan/Shandong-DF/figure/GGMP_control_stacked_bar.pdf",p_ggmp_control,height = 6,width = 6)

            # SGMP
                p_sgmp_case <- ggplot(data_sgmp_case,aes(sampleid,value,fill=phylum)) +
                    geom_bar(stat="identity", position = 'fill') +
                    xlab("") +
                    ylab("") +
                    #theme_classic(base_size = 7) +
                    scale_y_continuous(expand = c(0,0)) +
                    ggtitle('') +
                    guides(fill=guide_legend(title=NULL)) +
                    theme(axis.text.x=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.ticks.y = element_blank(),
                        legend.text = element_text(size = 12),
                        legend.position = "bottom") +
                    scale_fill_manual(values = color)
                ggsave("/home/wangnan/Shandong-DF/figure/SGMP_case_stacked_bar.pdf",p_sgmp_case,height = 6,width = 10)

                p_sgmp_control <- ggplot(data_sgmp_control,aes(sampleid,value,fill=phylum)) +
                    geom_bar(stat="identity", position = 'fill') +
                    xlab("") +
                    ylab("") +
                    #theme_classic(base_size = 7) +
                    scale_y_continuous(expand = c(0,0)) +
                    ggtitle('') +
                    guides(fill=guide_legend(title=NULL)) +
                    theme(axis.text.x=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.ticks.y = element_blank(),
                        legend.text = element_text(size = 12),
                        legend.position = "bottom") +
                    scale_fill_manual(values = color)
                ggsave("/home/wangnan/Shandong-DF/figure/SGMP_control_stacked_bar.pdf",p_sgmp_control,height = 6,width = 6)                
