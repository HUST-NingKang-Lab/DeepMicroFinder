library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)

#画图
#Guangdong
    sta <- data.frame(exp = "1" ,AUROC = "2",stringsAsFactors = FALSE)
    for(i in 1:10)
    {
        tmp <- read.csv(paste0("/data4/wangnan/Shandong-DF/prediction/Guangdong/exp",i,"/Evaluation_Independent/overall.csv"))
        sta <- rbind(sta,c(paste0("exp",i),mean(tmp$ROC.AUC)))
    }
    sta <- sta[2:nrow(sta),]
    sta$AUROC <- as.numeric(sta$AUROC)
    sta$exp <- factor(sta$exp,levels = c("exp1","exp2","exp3","exp4","exp5","exp6","exp7","exp8","exp9","exp10"))

    q <- ggplot(sta, mapping = aes(x = exp, y =AUROC)) +
        geom_bar(stat = "identity",aes(fill = AUROC))+
        theme_bw()+
        scale_fill_gradient(low = "#66CCFF", high = "#003366",limits = c(0.6, 0.85),breaks = c(0.6,0.7,0.8)) + 
        theme(axis.title.x = element_blank(),
            axis.text.x = element_text(size = 15,colour = "black"),
            axis.title.y = element_text(size = 15,colour = "black"),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 15,colour = "black"),
            legend.text = element_text(size = 15,colour = "black"),
            legend.title = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5,size = 15,colour = "black")) +
            #scale_fill_npg()+
            #geom_signif(comparisons = comparision,step_increase = 0.1,map_signif_level = F,test = wilcox.test)+
            #stat_compare_means(label = "p.format")+
            #facet_grid(~method,scales = "free")+
            ylab("AUROC")
    ggsave("/home/wangnan/Shandong-DF/figure/AUROC_GGMP.pdf",q,height = 4.5,width = 7)

#Shandong
    sta <- data.frame(method = "1" ,value = "2",stringsAsFactors = FALSE)
    for(n in 1:4)
    {
        for(i in 1:10)
        {
            tmp <- read.csv(paste0("/data4/wangnan/Shandong-DF/prediction/Shandong/exp",n,"/exp",i,"/Evaluation_Independent/overall.csv"))
            sta <- rbind(sta,c(paste0(n*20,"%"),tmp$ROC.AUC[1]))
            sta <- rbind(sta,c(paste0(n*20,"%"),tmp$ROC.AUC[2]))
        }
    }
    sta <- sta[2:nrow(sta),]
    sta$value <- as.numeric(sta$value)

    q <- ggplot(sta, aes(x = method,y = value,fill = method)) +
        geom_boxplot(aes(fill = method),position=position_dodge(0.7),width = 0.6) +
        #scale_color_manual(values=c("#E64B35", "#4DBBD5", "#00A087"),aesthetics = "fill")+ #箱线图的外框指定颜色
        #scale_color_manual(values=c("#000000", "#000000", "#000000"),aesthetics = "color")+ #箱线图的填充指定颜色
        theme_bw()+
        theme(axis.title.x = element_blank(),
            axis.text.x = element_text(size = 15,colour = "black"),
            axis.title.y = element_text(size = 15,colour = "black"),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 15,colour = "black"),
            legend.position = "none",
            legend.text = element_blank(),
            legend.title = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5,size = 15,colour = "black")) +
            scale_fill_npg()+
            #geom_signif(comparisons = comparision,step_increase = 0.1,map_signif_level = F,test = wilcox.test)+
            #stat_compare_means(label = "p.format")+
            #facet_grid(~method,scales = "free")+
            ylab("AUROC")


#Guangdong-Shandong
    sta <- data.frame(method = "1" ,proportion = "2",value = "2",stringsAsFactors = FALSE)
    for(n in 1:4)
    {
        for(i in 1:10)
        {
            tmp <- read.csv(paste0("/data4/wangnan/Shandong-DF/prediction/Shandong/exp",n,"/exp",i,"/Evaluation_Independent/overall.csv"))
            sta <- rbind(sta,c("Independent DNN model",paste0(n*20,"%"),tmp$ROC.AUC[1]))
            sta <- rbind(sta,c("Independent DNN model",paste0(n*20,"%"),tmp$ROC.AUC[1]))
            tmp <- read.csv(paste0("/data4/wangnan/Shandong-DF/prediction/Guangdong-Shandong/exp",n,"/exp",i,"/Evaluation_Independent/overall.csv"))
            sta <- rbind(sta,c("Regional+ DNN model",paste0(n*20,"%"),tmp$ROC.AUC[1]))
            sta <- rbind(sta,c("Regional+ DNN model",paste0(n*20,"%"),tmp$ROC.AUC[2]))
            tmp <- read.csv(paste0("/data4/wangnan/Shandong-DF/prediction/Guangdong-Shandong/exp",n,"/exp",i,"/Evaluation_Independent_Guangdong/overall.csv"))
            sta <- rbind(sta,c("Regional DNN model",paste0(n*20,"%"),tmp$ROC.AUC[1]))
            sta <- rbind(sta,c("Regional DNN model",paste0(n*20,"%"),tmp$ROC.AUC[2]))
            tmp <- read.csv(paste0("/data4/wangnan/Shandong-DF/prediction/Guangdong-Shandong/exp",n,"/exp",i,"/Evaluation_Transfer/overall.csv"))
            sta <- rbind(sta,c("Transfer DNN model",paste0(n*20,"%"),tmp$ROC.AUC[1]))
            sta <- rbind(sta,c("Transfer DNN model",paste0(n*20,"%"),tmp$ROC.AUC[2]))      
        }
    }
    sta <- sta[2:nrow(sta),]
    sta$value <- as.numeric(sta$value)

    # 计算各模型在各百分比数据下的平均准确性
        aggregate(sta$value,by = list(sta$method,sta$proportion),FUN =mean) 

    # 计算利用80%的数据进行迁移学习时模型的准确性
        pos <- intersect(which(sta$method == "Transfer DNN model"),which(sta$proportion == "80%"))
        print(paste0("Average AUROC of Transfer DNN model is: ", mean(sta$value[pos])))

    sta$method <- factor(sta$method,levels = c("Independent DNN model","Regional DNN model","Regional+ DNN model","Transfer DNN model"))
    comparison <- list(c("Independent DNN model","Transfer DNN model"))

    q <- ggplot(sta, aes(x = proportion,y = value,fill = method)) +
        geom_boxplot(aes(fill = method),position=position_dodge(0.55),width = 0.5) + # 箱线图取平均值线
        #scale_color_manual(values=c("#E64B35", "#4DBBD5", "#00A087"),aesthetics = "fill")+ #箱线图的外框指定颜色
        #scale_color_manual(values=c("#000000", "#000000", "#000000"),aesthetics = "color")+ #箱线图的填充指定颜色
        stat_summary(fun=mean, geom="line", aes(colour = method,group = method)) + #添加均值线
        stat_summary(fun=mean, geom="point", aes(colour = method,group = method)) + #添加均值点
        theme_bw()+
        theme(
            axis.title.x = element_text(size = 18,colour = "black"),
            axis.text.x = element_text(size = 18,colour = "black"),
            axis.title.y = element_text(size = 18,colour = "black"),
            axis.text.y = element_text(size = 18,colour = "black"),
            legend.text = element_text(size = 18,colour = "black"),
            legend.title = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5,size = 15,colour = "black")) +
            scale_fill_npg()+
            scale_color_npg()+
            ylim(c(0.25,1.2))+
            #stat_compare_means(label = "p.format")+
            #facet_grid(~method,scales = "free")+
            xlab("The proportion of training subset")+
            ylab("AUROC")
    
    ggsave("/home/wangnan/Shandong-DF/figure/AUROC_assessment.pdf",q,height = 6,width = 12)

# 流程图里面的一个小图
    sta <- data.frame(method = "1" ,proportion = "2",value = "2",stringsAsFactors = FALSE)
    for(n in 4)
    {
        for(i in 1:10)
        {
            tmp <- read.csv(paste0("/data4/wangnan/Shandong-DF/prediction/Shandong/exp",n,"/exp",i,"/Evaluation_Independent/overall.csv"))
            sta <- rbind(sta,c("Independent DNN model",paste0(n*20,"%"),tmp$ROC.AUC[1]))
            sta <- rbind(sta,c("Independent DNN model",paste0(n*20,"%"),tmp$ROC.AUC[1]))
            tmp <- read.csv(paste0("/data4/wangnan/Shandong-DF/prediction/Guangdong-Shandong/exp",n,"/exp",i,"/Evaluation_Independent/overall.csv"))
            sta <- rbind(sta,c("Regional+ DNN model",paste0(n*20,"%"),tmp$ROC.AUC[1]))
            sta <- rbind(sta,c("Regional+ DNN model",paste0(n*20,"%"),tmp$ROC.AUC[2]))
            tmp <- read.csv(paste0("/data4/wangnan/Shandong-DF/prediction/Guangdong-Shandong/exp",n,"/exp",i,"/Evaluation_Independent_Guangdong/overall.csv"))
            sta <- rbind(sta,c("Regional DNN model",paste0(n*20,"%"),tmp$ROC.AUC[1]))
            sta <- rbind(sta,c("Regional DNN model",paste0(n*20,"%"),tmp$ROC.AUC[2]))
            tmp <- read.csv(paste0("/data4/wangnan/Shandong-DF/prediction/Guangdong-Shandong/exp",n,"/exp",i,"/Evaluation_Transfer/overall.csv"))
            sta <- rbind(sta,c("Transfer DNN model",paste0(n*20,"%"),tmp$ROC.AUC[1]))
            sta <- rbind(sta,c("Transfer DNN model",paste0(n*20,"%"),tmp$ROC.AUC[2]))      
        }
    }
    sta <- sta[2:nrow(sta),]
    sta$value <- as.numeric(sta$value)

    sta$method <- factor(sta$method,levels = c("Independent DNN model","Regional DNN model","Regional+ DNN model","Transfer DNN model"))
    comparison <- list(c("Independent DNN model","Transfer DNN model"))

    q <- ggplot(sta, aes(x = method,y = value,fill = method)) +
        geom_boxplot(aes(fill = method),position=position_dodge(0.7),width = 0.6) +
        #scale_color_manual(values=c("#E64B35", "#4DBBD5", "#00A087"),aesthetics = "fill")+ #箱线图的外框指定颜色
        #scale_color_manual(values=c("#000000", "#000000", "#000000"),aesthetics = "color")+ #箱线图的填充指定颜色
        theme_classic()+
        theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y = element_text(size = 12,colour = "black"),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 12,colour = "black"),
            legend.text = element_blank(),
            legend.title = element_blank(),
            legend.position = "none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_blank()) +
            scale_fill_npg()+
            #stat_compare_means(label = "p.format")+
            #facet_grid(~method,scales = "free")+
            ylab("AUROC")
    
    ggsave("/home/wangnan/Shandong-DF/figure/figure 1B(1).pdf",q,height = 2,width = 3)
