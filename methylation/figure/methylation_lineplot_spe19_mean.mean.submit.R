library(dplyr)
library(data.table)
library (tidyr)
library(corrplot)
#library(ggsignif)
library(ggpubr)
library(svglite)
library(pheatmap)
library(ape)
library(ggpubr)

##########################Plotting different methylation types
#########CG
files =list.files("D:/OneDrive - zju.edu.cn/Quercus/adaptive introgression/methylation_analysis/summary/spe19/spe19_methlation_bin_mean/",pattern='*CG_methlation.genetype.bin.mean.gz$',full=T,rec=T)
#file1 <- fread("D:/OneDrive - zju.edu.cn/Quercus/adaptive introgression/methylation_analysis/summary/spe19/Qv.methylation.summary.bin.gz")  #Qa.methylation_union_type.anno.gz
species <- c("Q.acutissima","Q.aliena","Q.aliena var.acutiserrata","Q.chenii","Q.cocciferoides",
             "Q.dentata","Q.dolicholepis","Q.fabri","Q.franchetii","Q.glauca","Q.glaucoides",
             "Q.griffithii","Q.longispica","Q.mongolica","Q.phillyreoides","Q.robur",
             "Q.rehderiana","Q.serrata","Q.variabilis")
meth_spe19_CGtype=NULL
for (i in files){
  file = fread(i) 
  genetype <- file %>%  mutate(species=rep("19 species",nrow(file)))
  meth_spe19_CGtype <- rbind(meth_spe19_CGtype,genetype)
}
genetype <- meth_spe19_CGtype %>%     
  group_by(gene_type,Bintype) %>%
  summarise(
    mean = mean(mean_y),                 
    sd_y = sd(mean_y),                    
    n = n(),                          
    se = sd_y / sqrt(n),              
    ci_lower = mean - qt(0.975, df = n - 1) * se,    
    ci_upper = mean + qt(0.975, df = n - 1) * se     
    # max= max(mean_y),
    # min=min(mean_y),
  )
# ##Add Region
genetype$region <- genetype$Bintype
genetype$region <- substr(genetype$region,1,1) 
#data <- genetype %>% melt(id.vars=c("gene_type",'meth_type','Bintype'),variable.name="Bintype",value.name="methlation_ratio_mean")
genetype$Bintype <- factor(genetype$Bintype,levels = c("U1","U2","U3","U4","U5","U6","U7","U8","U9","U10","U11","U12","U13","U14","U15","U16","U17","U18","U19","U20","R1","R2","R3","R4","R5","R6","R7","R8","R9","R10",
                                                       "R11","R12","R13","R14","R15","R16","R17","R18","R19","R20","D1","D2","D3","D4","D5","D6","D7","D8","D9","D10","D11","D12","D13","D14","D15","D16","D17","D18","D19","D20"))
# breaks = seq(0,60,by=20)
# labels = c("-3K","TSS","TES","3K")
##Calculate significance
data <-compare_means(mean ~ gene_type,data = genetype,group.by="region",method = "t.test")
# region .y.   group1        group2                 p p.adj p.format p.signif method
# <chr>  <chr> <chr>         <chr>              <dbl> <dbl> <chr>    <chr>    <chr> 
#   1 D      mean  Introgression Non-introgression 0.0295 0.089 0.03     *        T-test
# 2 R      mean  Introgression Non-introgression 0.990  0.99  0.99     ns       T-test
# 3 U      mean  Introgression Non-introgression 0.250  0.5   0.25     ns       T-test
# Manually adjust the x and y positions of the p-value annotation
data$y.position <- rep(max(genetype$mean) + 0.05, nrow(data))  
p1 <- ggplot(genetype,aes(x=Bintype,y=mean,group =gene_type,color=gene_type))+
  geom_line(size=0.6)+
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper,fill = gene_type), 
              linetype = "blank" , alpha = 0.2) +  
  #facet_grid(meth_type~.)+
  #facet_grid(.~meth_type)+ 
  #facet_grid(rows = vars(region), cols = vars(meth_type), scales = "free_x", space = "free_x") +  
  #facet_grid(cols = vars(meth_type), scales = "free_x", space = "free_x") +
  scale_color_manual(values=c("#C75173","#358ACA")) +  
  scale_fill_manual(values=c("#C75173","#358ACA")) +  
  #scale_fill_manual(values=c("#34A471","#2169D3","#DE1A37")) +  
  #scale_fill_manual(values=c("#8ea1ca","#ec97cc")) +  
  #scale_fill_manual("",values=c("Internal branches"="#8ea1ca","External branches"="#ec97cc")) +  
  geom_vline(aes(xintercept="R1"),colour="#BB0000",linetype="dashed",size=0.6)+  
  geom_vline(aes(xintercept="R20"),colour="#BB0000",linetype="dashed",size=0.6)+  
  annotate("text", x = "U10", y = data$y.position, label = data[data$region == "U", ]$p.signif, size = 5, vjust =3.5) +  
  annotate("text", x = "R10", y = data$y.position, label = data[data$region == "R", ]$p.signif, size = 5, vjust =3.5) +  
  annotate("text", x = "D10", y = data$y.position, label = data[data$region == "D", ]$p.signif, size = 5, vjust =3.5) +   
  theme_bw() + 
  # theme(legend.position = c(0.94,0.95),legend.background = element_rect(fill = "NA", colour = "NA", size = 0.6),legend.key = element_rect(color = "NA", fill = "NA"),legend.key.size = unit(10, "pt"),legend.key.width=unit(20, "pt"),
  #       legend.text=element_text(size=9,family = "serif", color="black",face="bold"),
  #       legend.direction="vertical", legend.spacing.x = unit(0.5, 'cm'))+ 
  theme(panel.grid =element_blank()) +   
  guides(color=guide_legend(title=NULL))+ 
  theme(legend.position="none")+     
  scale_x_discrete(expand = c(0.01,0.01),breaks = c("U1","R1","R20","D20"),
                   labels = c("-3k","TSS","TES","+3k"))+
  #scale_x_discrete(expand = c(0.01,0.01))+
  scale_y_continuous(expand = c(0.01,0.00))+
  #theme(axis.text.x = element_blank()) +   
  theme(axis.title.x=element_blank())+ 
  theme(axis.text.x  = element_text(face="bold",family = "serif",color="black",size = 12,vjust = -0.4))+  #vjust负的越向外,face="bold"加粗
  theme(axis.text.y  = element_text(face="bold",family = "serif", color="black",size = 12))+ #angle = 90 #旋转字体
  theme(strip.background.y = element_rect(color="black", size=1, linetype="solid"))+    #theme(strip.background.x = element_rect(color="black", fill="#FC4E07", size=1.5, linetype="solid"))
  theme(strip.text.x = element_text(size = 12, family = "serif", color = "black", face = "bold"))+ 
  theme(panel.spacing.x = unit(0.5, "cm"))+ 
  theme(plot.margin=unit(rep(1,4),'lines'))+ 
  ylab("Methylation level")+theme(axis.title.y = element_text(size = 15, family = "serif", color = "black", face = "bold"))+
  labs(title = "CG")+theme(plot.title = element_text(size = 16, family = "serif", color = "black", face = "bold", hjust = 0.5, vjust = 0.5,angle = 0)) ##标题限制
#myfilename <- paste("D:/OneDrive - zju.edu.cn/Quercus/adaptive introgression/methylation_analysis/plot/lineplot/genetype_plot/", "Quercus variabilis",".png",sep ="") 
#ggsave(p, file = myfilename, dpi = 600,width=12, height=6)

#########CHG
files =list.files("D:/OneDrive - zju.edu.cn/Quercus/adaptive introgression/methylation_analysis/summary/spe19/spe19_methlation_bin_mean/",pattern='*CHG_methlation.genetype.bin.mean.gz$',full=T,rec=T)
#file1 <- fread("D:/OneDrive - zju.edu.cn/Quercus/adaptive introgression/methylation_analysis/summary/spe19/Qv.methylation.summary.bin.gz")  #Qa.methylation_union_type.anno.gz
species <- c("Q.acutissima","Q.aliena","Q.aliena var.acutiserrata","Q.chenii","Q.cocciferoides",
             "Q.dentata","Q.dolicholepis","Q.fabri","Q.franchetii","Q.glauca","Q.glaucoides",
             "Q.griffithii","Q.longispica","Q.mongolica","Q.phillyreoides","Q.robur",
             "Q.rehderiana","Q.serrata","Q.variabilis")
meth_spe19_CHGtype=NULL
for (i in files){
  file = fread(i) 
  genetype <- file %>%  mutate(species=rep("19 species",nrow(file)))
  meth_spe19_CHGtype <- rbind(meth_spe19_CHGtype,genetype)
}
genetype <- meth_spe19_CHGtype %>%     
  group_by(gene_type,Bintype) %>%
  summarise(
    mean = mean(mean_y),                 
    sd_y = sd(mean_y),                    
    n = n(),                          
    se = sd_y / sqrt(n),              
    ci_lower = mean - qt(0.975, df = n - 1) * se,    
    ci_upper = mean + qt(0.975, df = n - 1) * se     
    # max= max(mean_y),
    # min=min(mean_y),
  )
# ##Add Region
genetype$region <- genetype$Bintype
genetype$region <- substr(genetype$region,1,1) 
#data <- genetype %>% melt(id.vars=c("gene_type",'meth_type','Bintype'),variable.name="Bintype",value.name="methlation_ratio_mean")
genetype$Bintype <- factor(genetype$Bintype,levels = c("U1","U2","U3","U4","U5","U6","U7","U8","U9","U10","U11","U12","U13","U14","U15","U16","U17","U18","U19","U20","R1","R2","R3","R4","R5","R6","R7","R8","R9","R10",
                                                       "R11","R12","R13","R14","R15","R16","R17","R18","R19","R20","D1","D2","D3","D4","D5","D6","D7","D8","D9","D10","D11","D12","D13","D14","D15","D16","D17","D18","D19","D20"))
# breaks = seq(0,60,by=20)
# labels = c("-3K","TSS","TES","3K")
##Calculate significance
data <-compare_means(mean ~ gene_type,data = genetype,group.by="region",method = "t.test")
# region .y.   group1        group2                   p    p.adj p.format p.signif method
# <chr>  <chr> <chr>         <chr>                <dbl>    <dbl> <chr>    <chr>    <chr> 
#   1 D      mean  Introgression Non-introgression 2.32e- 3 4.6 e- 3 0.0023   **       T-test
# 2 R      mean  Introgression Non-introgression 8.02e-14 2.40e-13 8e-14    ****     T-test
# 3 U      mean  Introgression Non-introgression 2.11e- 1 2.1 e- 1 0.2111   ns       T-test
# Manually adjust the x and y positions of the p-value annotation
data$y.position <- rep(max(genetype$mean) + 0.05, nrow(data))  
p2 <- ggplot(genetype,aes(x=Bintype,y=mean,group =gene_type,color=gene_type))+
  geom_line(size=0.6)+
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper,fill = gene_type), 
              linetype = "blank" , alpha = 0.2) +  
  #facet_grid(meth_type~.)+
  #facet_grid(.~meth_type)+ 
  #facet_grid(rows = vars(region), cols = vars(meth_type), scales = "free_x", space = "free_x") +  
  #facet_grid(cols = vars(meth_type), scales = "free_x", space = "free_x") +
  scale_color_manual(values=c("#C75173","#358ACA")) +  
  scale_fill_manual(values=c("#C75173","#358ACA")) +  
  #scale_fill_manual(values=c("#34A471","#2169D3","#DE1A37")) +  
  #scale_fill_manual(values=c("#8ea1ca","#ec97cc")) +  
  #scale_fill_manual("",values=c("Internal branches"="#8ea1ca","External branches"="#ec97cc")) +  
  geom_vline(aes(xintercept="R1"),colour="#BB0000",linetype="dashed",size=0.6)+  
  geom_vline(aes(xintercept="R20"),colour="#BB0000",linetype="dashed",size=0.6)+  
  annotate("text", x = "U10", y = data$y.position, label = data[data$region == "U", ]$p.signif, size = 5, vjust =3.5) +  
  annotate("text", x = "R10", y = data$y.position, label = data[data$region == "R", ]$p.signif, size = 5, vjust =3.5) +  
  annotate("text", x = "D10", y = data$y.position, label = data[data$region == "D", ]$p.signif, size = 5, vjust =3.5) +   
  theme_bw() + 
  # theme(legend.position = c(0.94,0.95),legend.background = element_rect(fill = "NA", colour = "NA", size = 0.6),legend.key = element_rect(color = "NA", fill = "NA"),legend.key.size = unit(10, "pt"),legend.key.width=unit(20, "pt"),
  #       legend.text=element_text(size=9,family = "serif", color="black",face="bold"),
  #       legend.direction="vertical", legend.spacing.x = unit(0.5, 'cm'))+ 
  theme(panel.grid =element_blank()) +   
  guides(color=guide_legend(title=NULL))+ 
  theme(legend.position="none")+     
  # scale_x_discrete(expand = c(0.01,0.01),breaks = seq(0,53,by=15),
  #                  labels = c("-3K","TSS","TES","3K"))+
  scale_x_discrete(expand = c(0.01,0.01),breaks = c("U1","R1","R20","D20"),
                   labels = c("-3k","TSS","TES","+3k"))+
  #scale_x_discrete(expand = c(0.01,0.01))+
  scale_y_continuous(expand = c(0.01,0.00))+
  #theme(axis.text.x = element_blank()) +   
  theme(axis.title.x=element_blank())+ 
  theme(axis.text.x  = element_text(face="bold",family = "serif",color="black",size = 12,vjust = -0.4))+  #vjust负的越向外,face="bold"加粗
  theme(axis.text.y  = element_text(face="bold",family = "serif", color="black",size = 12))+ #angle = 90 #旋转字体
    theme(strip.background.y = element_rect(color="black", size=1, linetype="solid"))+    #theme(strip.background.x = element_rect(color="black", fill="#FC4E07", size=1.5, linetype="solid"))
  theme(strip.text.x = element_text(size = 12, family = "serif", color = "black", face = "bold"))+ 
  theme(panel.spacing.x = unit(0.5, "cm"))+ 
  theme(plot.margin=unit(rep(1,4),'lines'))+ 
  #ylab("Methylation level")+theme(axis.title.y = element_text(size = 15, family = "serif", color = "black", face = "bold"))+
  ylab("")+theme(axis.title.y = element_text(size = 15, family = "serif", color = "black", face = "bold"))+
  labs(title = "CHG")+theme(plot.title = element_text(size = 16, family = "serif", color = "black", face = "bold", hjust = 0.5, vjust = 0.5,angle = 0)) ##标题限制
#myfilename <- paste("D:/OneDrive - zju.edu.cn/Quercus/adaptive introgression/methylation_analysis/plot/lineplot/genetype_plot/", "Quercus variabilis",".png",sep ="") 
#ggsave(p, file = myfilename, dpi = 600,width=12, height=6)

#########CHH
files =list.files("D:/OneDrive - zju.edu.cn/Quercus/adaptive introgression/methylation_analysis/summary/spe19/spe19_methlation_bin_mean/",pattern='*CHH_methlation.genetype.bin.mean.gz$',full=T,rec=T)
#file1 <- fread("D:/OneDrive - zju.edu.cn/Quercus/adaptive introgression/methylation_analysis/summary/spe19/Qv.methylation.summary.bin.gz")  #Qa.methylation_union_type.anno.gz
species <- c("Q.acutissima","Q.aliena","Q.aliena var.acutiserrata","Q.chenii","Q.cocciferoides",
             "Q.dentata","Q.dolicholepis","Q.fabri","Q.franchetii","Q.glauca","Q.glaucoides",
             "Q.griffithii","Q.longispica","Q.mongolica","Q.phillyreoides","Q.robur",
             "Q.rehderiana","Q.serrata","Q.variabilis")

meth_spe19_CHHtype=NULL
for (i in files){
  file = fread(i) 
  genetype <- file %>%  mutate(species=rep("19 species",nrow(file)))
  meth_spe19_CHHtype <- rbind(meth_spe19_CHHtype,genetype)
}
genetype <- meth_spe19_CHHtype %>%     
  group_by(gene_type,Bintype) %>%
  summarise(
    mean = mean(mean_y),                 
    sd_y = sd(mean_y),                    
    n = n(),                          
    se = sd_y / sqrt(n),              
    ci_lower = mean - qt(0.975, df = n - 1) * se,    
    ci_upper = mean + qt(0.975, df = n - 1) * se     
    # max= max(mean_y),
    # min=min(mean_y),
  )
# ##Add Region
genetype$region <- genetype$Bintype
genetype$region <- substr(genetype$region,1,1) 
#data <- genetype %>% melt(id.vars=c("gene_type",'meth_type','Bintype'),variable.name="Bintype",value.name="methlation_ratio_mean")
genetype$Bintype <- factor(genetype$Bintype,levels = c("U1","U2","U3","U4","U5","U6","U7","U8","U9","U10","U11","U12","U13","U14","U15","U16","U17","U18","U19","U20","R1","R2","R3","R4","R5","R6","R7","R8","R9","R10",
                                                       "R11","R12","R13","R14","R15","R16","R17","R18","R19","R20","D1","D2","D3","D4","D5","D6","D7","D8","D9","D10","D11","D12","D13","D14","D15","D16","D17","D18","D19","D20"))
# breaks = seq(0,60,by=20)
# labels = c("-3K","TSS","TES","3K")
##Calculate significance
data <-compare_means(mean ~ gene_type,data = genetype,group.by="region",method = "t.test")
# region .y.   group1        group2                   p         p.adj p.format p.signif method
# <chr>  <chr> <chr>         <chr>                <dbl>         <dbl> <chr>    <chr>    <chr> 
#   1 D      mean  Introgression Non-introgression 6.81e- 1 0.68          0.68     ns       T-test
# 2 R      mean  Introgression Non-introgression 3.14e-10 0.00000000094 3.1e-10  ****     T-test
# 3 U      mean  Introgression Non-introgression 1.29e- 1 0.26          0.13     ns       T-test
R_intro<- genetype %>% filter(region=='R'& gene_type=='Introgression')
R_nointro<- genetype %>% filter(region=='R'& gene_type=='Non-introgression')
t.test(R_intro$mean,R_nointro$mean)#p-value = 3.136e-10
# Manually adjust the x and y positions of the p-value annotation
data$y.position <- rep(max(genetype$mean) + 0.05, nrow(data))  
p3 <- ggplot(genetype,aes(x=Bintype,y=mean,group =gene_type,color=gene_type))+
  geom_line(size=0.6)+
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper,fill = gene_type), 
              linetype = "blank" , alpha = 0.2) +  
  #facet_grid(meth_type~.)+
  #facet_grid(.~meth_type)+ 
  #facet_grid(rows = vars(region), cols = vars(meth_type), scales = "free_x", space = "free_x") +  
  #facet_grid(cols = vars(meth_type), scales = "free_x", space = "free_x") +
  scale_color_manual(values=c("#C75173","#358ACA")) +  
  scale_fill_manual(values=c("#C75173","#358ACA")) +  
  #scale_fill_manual(values=c("#34A471","#2169D3","#DE1A37")) +  
  #scale_fill_manual(values=c("#8ea1ca","#ec97cc")) +  
  #scale_fill_manual("",values=c("Internal branches"="#8ea1ca","External branches"="#ec97cc")) +  
  geom_vline(aes(xintercept="R1"),colour="#BB0000",linetype="dashed",size=0.6)+  
  geom_vline(aes(xintercept="R20"),colour="#BB0000",linetype="dashed",size=0.6)+  
  annotate("text", x = "U10", y = data$y.position, label = data[data$region == "U", ]$p.signif, size = 5, vjust =3.5) +  
  annotate("text", x = "R10", y = data$y.position, label = data[data$region == "R", ]$p.signif, size = 5, vjust =3.5) +  
  annotate("text", x = "D10", y = data$y.position, label = data[data$region == "D", ]$p.signif, size = 5, vjust =3.5) +   
  theme_bw() + 
  theme(legend.position = c(0.8,0.95),legend.background = element_rect(fill = "NA", colour = "NA", size = 0.6),legend.key = element_rect(color = "NA", fill = "NA"),legend.key.size = unit(10, "pt"),legend.key.width=unit(20, "pt"),
        legend.text=element_text(size=9,family = "serif", color="black",face="bold"),
        legend.direction="vertical", legend.spacing.x = unit(0.5, 'cm'))+ 
  theme(panel.grid =element_blank()) +   
  #guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))  # 合并图例
  guides(color=guide_legend(title=NULL))+ 
  guides(fill=guide_legend(title=NULL))+
  scale_x_discrete(expand = c(0.01,0.01),breaks = c("U1","R1","R20","D20"),
                 labels = c("-3k","TSS","TES","+3k"))+
  #scale_x_discrete(expand = c(0.01,0.01))+
  scale_y_continuous(expand = c(0.01,0.00))+
  #theme(axis.text.x = element_blank()) +   
  theme(axis.title.x=element_blank())+ 
  theme(axis.text.x  = element_text(face="bold",family = "serif",color="black",size = 12,vjust = -0.4))+  #vjust负的越向外,face="bold"加粗
  theme(axis.text.y  = element_text(face="bold",family = "serif", color="black",size = 12))+ #angle = 90 #旋转字体
    theme(strip.background.y = element_rect(color="black", size=1, linetype="solid"))+    #theme(strip.background.x = element_rect(color="black", fill="#FC4E07", size=1.5, linetype="solid"))
  theme(strip.text.x = element_text(size = 12, family = "serif", color = "black", face = "bold"))+ 
  theme(panel.spacing.x = unit(0.5, "cm"))+ 
  theme(plot.margin=unit(rep(1,4),'lines'))+ 
  #ylab("Methylation level")+theme(axis.title.y = element_text(size = 15, family = "serif", color = "black", face = "bold"))+
  ylab("")+theme(axis.title.y = element_text(size = 15, family = "serif", color = "black", face = "bold"))+
  labs(title = "CHH")+theme(plot.title = element_text(size = 16, family = "serif", color = "black", face = "bold", hjust = 0.5, vjust =0.5,angle = 0)) ##标题限制,hjust = 0.5, vjust = 0.5居中对齐，朝外
#myfilename <- paste("D:/OneDrive - zju.edu.cn/Quercus/adaptive introgression/methylation_analysis/plot/lineplot/genetype_plot/", "Quercus variabilis",".png",sep ="") 
#ggsave(p, file = myfilename, dpi = 600,width=12, height=6)
library(ggpubr)
##Set the margins of the plot
p1 <- p1 + theme(plot.margin = unit(c(0.8, 0, 0.5, 0.3), "cm")) #分别表示上、右、下、左四侧边距
p2 <- p2 + theme(plot.margin = unit(c(0.8, 0, 0.5, 0), "cm"))
p3 <- p3 + theme(plot.margin = unit(c(0.8, 0.4, 0.5, 0), "cm"))
p4 <- ggarrange(p1,p2,p3,nrow=1,ncol=3)
#p4 <- annotate_figure(p3,top = text_grob(title_text, color = "black", family = "serif",face = "bold", size = 16,vjust =3))
p5 <- annotate_figure(p4,top = text_grob("The average methylation levels for 19 species", color = "black", family = "serif",face = "bold", size = 16,vjust =2))
myfilename <- paste("D:/OneDrive - zju.edu.cn/Quercus/adaptive introgression/methylation_analysis/plot/lineplot/genetype_plot/update/", "spe19_permean_mean_final",".png",sep ="") 
myfilename2 <- paste("D:/OneDrive - zju.edu.cn/Quercus/adaptive introgression/methylation_analysis/plot/lineplot/genetype_plot/update/", "spe19_permean_mean_final",".svg",sep ="") 
ggsave(p5, file = myfilename, dpi = 600,width=16, height=6)
ggsave(file=myfilename2,plot=p4, width=16,height=6) 





