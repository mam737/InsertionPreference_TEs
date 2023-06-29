.libPaths('/home/brandvai/mmunasin/Rlibs4')

library(reshape2)
library(ggplot2)
library(gridExtra)
library(stringr)
library(tidyr)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)
library(svglite)
library(nord)
library(ggpubr)
library(ggnewscale)
library(viridis)

source('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/plotting_script_functions.R')

# The following versions reference the final design/script used for each model
# Model 1 - v3.2 - One Chromosome - No Excision - No Non-Autonomous
# Model 2 - v4.1 - One Chromosome - With Excision - No Non-Autonomous
# Model 3 - v7.1 - One Chromosome - No Excision - With Non-Autonomous
# Model 4 - v8 - One Chromosome - With Excision - With Non-Autonomous
# Model 5 - v3.3 - Five Chromosome - No Excision - No Non-Autonomous
# Model 6 - v10 - Five Chromosome - No Excision - No Non-Autonomous - With HR

## Load summarized data of # of TE Loss, Host Pop Extinction, and Dual Survival Outcomes
v3.2_transformed_log.df <- format_summary_df(read.csv('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/v3.2_results_updated/summary_outputs_scaled.csv',header=T,stringsAsFactors=F))
v4.1_transformed_log.df <- format_summary_df(read.csv('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/v4.1_results_updated/summary_outputs_scaled.csv',header=T,stringsAsFactors=F))
v7.1_transformed_log.df <- format_summary_df(read.csv('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/v7.1_results_updated/summary_outputs_scaled.csv',header=T,stringsAsFactors=F))
v8_transformed_log.df <- format_summary_df(read.csv('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/v8_results_updated/summary_outputs_scaled.csv',header=T,stringsAsFactors=F))
v3.3_transformed_log.df <- format_summary_df(read.csv('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/v3.3_results_updated/summary_outputs_scaled.csv',header=T,stringsAsFactors=F))
v10_transformed_log.df <- format_summary_df(read.csv('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/v10_results_updated/summary_outputs_scaled.csv',header=T,stringsAsFactors=F))

## Load info on TE Loss outcomes across all models
v3.2_TE_Loss.df <- read.csv('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/v3.2_results_updated/TE_Loss.csv',header=T,stringsAsFactors=F)
v3.2_TE_Loss.df$Model <- 'Model 1'

v4.1_TE_Loss.df <- read.csv('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/v4.1_results_updated/TE_Loss.csv',header=T,stringsAsFactors=F)
v4.1_TE_Loss.df$Model <- 'Model 2'

v7.1_TE_Loss.df <- read.csv('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/v7.1_results_updated/TE_Loss.csv',header=T,stringsAsFactors=F)
v7.1_TE_Loss.df <- v7.1_TE_Loss.df[,c(1,2,3,4)]
v7.1_TE_Loss.df$Model <- 'Model 3'

v8_TE_Loss.df <- read.csv('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/v8_results_updated/TE_Loss.csv',header=T,stringsAsFactors=F)
v8_TE_Loss.df <- v8_TE_Loss.df[,c(1,2,3,4)]
v8_TE_Loss.df$Model <- 'Model 4'

v3.3_TE_Loss.df <- read.csv('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/v3.3_results_updated/TE_Loss.csv',header=T,stringsAsFactors=F)
v3.3_TE_Loss.df$Model <- 'Model 5'

v10_TE_Loss.df <- read.csv('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/v10_results_updated/TE_Loss.csv',header=T,stringsAsFactors=F)
v10_TE_Loss.df$Model <- 'Model 6'

TE_Loss.df <- rbind(v3.2_TE_Loss.df,v4.1_TE_Loss.df,v7.1_TE_Loss.df,v8_TE_Loss.df,v3.3_TE_Loss.df,v10_TE_Loss.df)

## For Model 1 (v3.2) + Model 5 (v3.3) - read stored data on full trajectories across param space
v3.3_sim_data.df <- format_summary_df(read.csv('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/v3.3_results_updated/sim_data.csv',header=T,stringsAsFactors=F))
v10_sim_data.df <- format_summary_df(read.csv('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/v10_results_updated/sim_data.csv',header=T,stringsAsFactors=F))

## For Model 5 (v3.3) - Extract Final Allele Frequencies For Dual Survival Events
v3.3_DS_FAF.df <- read.csv('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/v3.3_results_updated/DS_FAF.csv',header=T,stringsAsFactors=F)
v10_DS_FAF.df <- read.csv('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/v10_results_updated/DS_FAF.csv',header=T,stringsAsFactors=F)


## For all models, extract mean final copy # and TE allele freqs in final gen
v3.2_final.df <- format_summary_df(read.csv('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/v3.2_results_updated/v3.2_final_gen.csv',header=T,stringsAsFactors=F))
v4.1_final.df <- format_summary_df(read.csv('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/v4.1_results_updated/v4.1_final_gen.csv',header=T,stringsAsFactors=F))
v7.1_final.df <- format_summary_df(read.csv('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/v7.1_results_updated/v7.1_final_gen.csv',header=T,stringsAsFactors=F))
v8_final.df <- format_summary_df(read.csv('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/v8_results_updated/v8_final_gen.csv',header=T,stringsAsFactors=F))
v3.3_final.df <- format_summary_df(read.csv('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/v3.3_results_updated/v3.3_final_gen.csv',header=T,stringsAsFactors=F))
v10_final.df <- format_summary_df(read.csv('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/v10_results_updated/v10_final_gen.csv',header=T,stringsAsFactors=F))

###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################

# Fig S1
p1 <- TE_loss_heatmaps(v3.3_transformed_log.df,'Model 1',TRUE)
p2 <- TE_loss_heatmaps(v10_transformed_log.df,'Model 2',TRUE)
p3 <- TE_loss_heatmaps(v3.2_transformed_log.df,'Model 3',TRUE)
p4 <- TE_loss_heatmaps(v4.1_transformed_log.df,'Model 4',TRUE)
p5 <- TE_loss_heatmaps(v7.1_transformed_log.df,'Model 5',TRUE)
p6 <- TE_loss_heatmaps(v8_transformed_log.df,'Model 6',TRUE)

FigS1 <- ggarrange(p1,p2,p3,p4,p5,p6,ncol=2,nrow=3,common.legend=T,legend='bottom')
ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/FigS1.pdf',width=7,height=8)

###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################

# Fig S2
p1 <- ext_ds_heatmaps(v3.3_transformed_log.df,'Model 1',TRUE)
p2 <- ext_ds_heatmaps(v10_transformed_log.df,'Model 2',TRUE)
p3 <- ext_ds_heatmaps(v3.2_transformed_log.df,'Model 3',TRUE)
p4 <- ext_ds_heatmaps(v4.1_transformed_log.df,'Model 4',TRUE)
p5 <- ext_ds_heatmaps(v7.1_transformed_log.df,'Model 5',TRUE)
p6 <- ext_ds_heatmaps(v8_transformed_log.df,'Model 6',TRUE)

FigS2 <- ggarrange(p1,p2,p3,p4,p5,p6,ncol=2,nrow=3,common.legend=T,legend='bottom')
ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/FigS2.png',width=8,height=10.5,bg='white')

###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################

# Fig S3
p1 <- establishment_prop_plot(v3.3_transformed_log.df,FALSE)
p2 <- establishment_prop_plot(v10_transformed_log.df,FALSE)
p3 <- establishment_prop_plot(v3.2_transformed_log.df,FALSE)
p4 <- establishment_prop_plot(v4.1_transformed_log.df,FALSE)
p5 <- establishment_prop_plot(v7.1_transformed_log.df,FALSE)
p6 <- establishment_prop_plot(v8_transformed_log.df,FALSE)

FigS3 <- ggarrange(p1,p2,p3,p4,p5,p6,ncol=6,nrow=1,common.legend=T,legend='bottom')
ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/FigS3.pdf',width=7,height=8.5)

###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################

# Fig S4
FigS4 <- fitness_trajectories(v3.3_sim_data.df,'Model 1')
ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/FigS4.pdf',plot=FigS4,width=18,height=18,bg='white')

FigS5 <- extinction_time_plot(v3.3_final.df,FALSE)
ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/FigS5.pdf',plot=FigS5,width=7,height=5)

FigS6 <- extinction_neutral_count_plot(v3.3_final.df, FALSE)
ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/FigS6.pdf',plot=FigS6,width=7,height=5)

FigS7 <- extinction_del_count_plot(v3.3_final.df, FALSE)
ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/FigS7.pdf',plot=FigS7,width=7,height=5)

FigS8 <- final_copy_number_ranges(v3.3_final.df, v10_final.df,v3.2_final.df,v4.1_final.df,v7.1_final.df,v8_final.df)
ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/FigS8.pdf',plot=FigS8,width=7,height=8)

FigS9 <- final_copy_number_bins(v3.3_final.df, v10_final.df,v3.2_final.df,v4.1_final.df,v7.1_final.df,v8_final.df,'by_cn')
ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/FigS9.pdf',plot=FigS9,width=18,height=18,bg='white')

FigS10 <- final_copy_number_bins(v3.3_final.df, v10_final.df,v3.2_final.df,v4.1_final.df,v7.1_final.df,v8_final.df,'by_outcome')
ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/FigS10.pdf',plot=FigS10,width=18,height=18,bg='white')

FigS11 <- TE_cn_trajectories(v10_sim_data.df)
ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/FigS11.png',plot=FigS11,width=18,height=18,bg='white')

FigS12 <- cn_delta_trajectories('/scratch.global/mmunasin/HR_results/yb_eq.csv')
ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/FigS12.pdf',plot=FigS12,width=18,height=18,bg='white')

del_traject <- v10_sim_data.df %>% subset(class=='Maint')

del_traject$TotalDel  <- del_traject$m2_mcn + del_traject$m3_mcn + del_traject$m4_mcn


FigS13 <- ggplot(del_traject,aes(x=Gen,y=TotalDel,group=tag,color=m1_mcn))+
  geom_line()+facet_grid(fct_rev(teJumpP)~neutP,scales='free_x') + 
  scale_color_viridis(option="viridis",end=0.75,direction=-1,name='Neutral TE \nCopy Number')+
  xlab('Generation') + ylab(expression(paste("Deleterious Mean \nCopy Number"))) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                          axis.title.x = element_text(size=16,margin=margin(t=10)),
                          axis.title.y = element_text(size=16,margin=margin(r=15,l=10)),
                          plot.title = element_blank(),
                          strip.text.x = element_text(size=12),
                          panel.grid=element_blank())+ theme(plot.margin=margin(t=25,l=10))


ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/FigS13.png',plot=FigS13,width=18,height=18,bg='white')

FigS14 <- ggplot(del_traject,aes(x=Gen,y=P1_MeanFit,group=tag,color=m1_mcn))+
  geom_line()+facet_grid(fct_rev(teJumpP)~neutP,scales='free_x') + 
  scale_color_viridis(option="viridis",end=0.75,direction=-1,name='Neutral TE \nCopy Number')+
  ggtitle(expression(paste("Host Population Mean Fitness Across Replicates for a Subset of Population Extinction Outcomes \n"))) +
  xlab('Generation') + ylab(expression(paste("Population \nMean Fitness"))) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                          axis.title.x = element_text(size=16,margin=margin(t=10)),
                          axis.title.y = element_text(size=16,margin=margin(r=15,l=10)),
                          plot.title = element_text(size = 22,face='bold'),
                          strip.text.x = element_text(size=12),
                          panel.grid=element_blank())+ theme(plot.margin=margin(t=25,l=10))

ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/FigS14.png',plot=FigS14,width=18,height=18,bg='white')


###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################

#Fig 2

Fig2_p1 <- TE_loss_heatmaps(v3.3_transformed_log.df,'Model 1',FALSE)
Fig2_p2 <- TE_loss_heatmaps(v10_transformed_log.df,'Model 2',FALSE)

Fig2_Row1 <- ggarrange(Fig2_p1,Fig2_p2,ncol=2,nrow=1,common.legend=T,legend='right')
ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/Fig2_Row1.png',plot=Fig2_Row1,width=6,height=2.5,bg='white')

Fig2_p3 <- ext_ds_heatmaps(v3.3_transformed_log.df,'Model 1',FALSE)
Fig2_p4 <- ext_ds_heatmaps(v10_transformed_log.df,'Model 2',FALSE)

Fig2_Row2 <- ggarrange(Fig2_p3,Fig2_p4,ncol=2,nrow=1,common.legend=T,legend='right')
ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/Fig2_Row2.png',plot=Fig2_Row2,width=10,height=3.5,bg='white')

Fig2_p5 <- establishment_prop_plot(v3.3_transformed_log.df,TRUE)

ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/Fig2_Row3.png',plot=Fig2_p5,width=6,height=2,bg='white')

###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################

# Fig 3

Fig3_p1 <- extinction_time_plot(v3.3_final.df,TRUE)
ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/Fig3_Row1.png',plot=Fig3_p1,width=6.5,height=1.5,bg='white')

Fig3_p2 <- extinction_neutral_count_plot(v3.3_final.df, TRUE)
ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/Fig3_Row2.png',plot=Fig3_p2,width=6.5,height=1.5,bg='white')

Fig3_p3 <- extinction_del_count_plot(v3.3_final.df, TRUE)
ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/Fig3_Row3.png',plot=Fig3_p3,width=6.5,height=1.5,bg='white')


Fig3_Row4 <- del_fit_trajects(v3.3_sim_data.df,'del','Ext')

ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/Fig3_Row4.png',plot=Fig3_Row4,width=16,height=3,bg='white')

Fig3_Row5 <- del_fit_trajects(v3.3_sim_data.df,'fit','Ext')

ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/Fig3_Row5.png',plot=Fig3_Row5,width=16,height=3,bg='white')

###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################

test_1 <- v3.3_final.df %>%
subset(class=='Maint') %>%
subset(!(teJumpP %in% c('1.0e-01','7.5e-02','5.0e-02'))) %>%
dplyr::summarise()



test_1 <- final_TEcn_df(v3.3_final.df)
test_2 <- final_TEcn_df(v10_final.df)


plot <- ggplot() +
geom_point(data=test_1,aes(x=neutP,y=teJumpP,color=MeanCN),shape="\u25E4",size=7)+
scale_color_gradient(low='grey100',high='#16697A',name='',limits=c(1,150000),trans='log10')+
scale_y_discrete(limits=levels(test_1$teJumpP)[1:10],drop=FALSE)+
new_scale_color()+
geom_point(data=test_2,aes(x=neutP,y=teJumpP,color=MeanCN),shape="\u25E2",size=7)+
scale_color_gradient(low='grey100',high='#16697A',name='',limits=c(1,150000),trans='log10')+
theme(panel.grid.major = element_line(linetype = "blank"),
    panel.grid.minor = element_line(linetype = "blank")) 
ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/FigX_triangle.png',plot=plot,width=4.8,height=3,bg='white')

plot2 <- ggplot() +
geom_point(data=test_1,aes(x=neutP,y=teJumpP,color=MeanFit),shape="\u25E4",size=7)+
scale_color_gradient(low='grey100',high='#7DBA84',name='',limits=c(0.6,1))+
scale_y_discrete(limits=levels(test_1$teJumpP)[1:10],drop=FALSE)+
new_scale_color()+
geom_point(data=test_2,aes(x=neutP,y=teJumpP,color=MeanFit),shape="\u25E2",size=7)+
scale_color_gradient(low='grey100',high='#7DBA84',name='',limits=c(0.6,1))+
theme(panel.grid.major = element_line(linetype = "blank"),
    panel.grid.minor = element_line(linetype = "blank")) 
ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/FigX_triangle2.png',plot=plot2,width=4.8,height=3,bg='white')

cn_bin_v3.3 <- v3.3_final.df %>% 
dplyr::mutate(total_mcn=m1_mcn+m2_mcn+m3_mcn+m4_mcn) %>%
dplyr::mutate(final_cn_bin=case_when(
  total_mcn < 10 ~ '< 10',
  total_mcn >= 10 & total_mcn < 100 ~ '10 - 100',
  total_mcn >= 100 & total_mcn < 1000 ~ '100 - 1,000',
  total_mcn >= 1000 & total_mcn < 10000 ~ '1,000 - 10,000',
  total_mcn >= 10000 ~ '> 10,000')) %>%
group_by(final_cn_bin,class) %>%
dplyr::summarise(n=n()) %>%
dplyr::mutate(Model='1')

cn_bin_v10 <- v10_final.df %>% 
dplyr::mutate(total_mcn=m1_mcn+m2_mcn+m3_mcn+m4_mcn) %>%
dplyr::mutate(final_cn_bin=case_when(
  total_mcn < 10 ~ '< 10',
  total_mcn >= 10 & total_mcn < 100 ~ '10 - 100',
  total_mcn >= 100 & total_mcn < 1000 ~ '100 - 1,000',
  total_mcn >= 1000 & total_mcn < 10000 ~ '1,000 - 10,000',
  total_mcn >= 10000 ~ '> 10,000')) %>%
group_by(final_cn_bin,class) %>%
dplyr::summarise(n=n()) %>%
dplyr::mutate(Model='2')


cn_bin <- rbind(cn_bin_v3.3,cn_bin_v10) %>%
dplyr::mutate(final_cn_bin=factor(final_cn_bin,levels=c('< 10','10 - 100','100 - 1,000','1,000 - 10,000','> 10,000')))

p <- ggplot(cn_bin, aes(fill=class, y=n, x=Model)) + 
    geom_bar(position="stack", stat="identity",width=0.9)+
    facet_wrap(~final_cn_bin,ncol=5) +
    scale_fill_manual(values=c('#16697A','#7DBA84'))+
    labs(x='',y='')+
    scale_y_continuous(breaks=c(0,2500,5000))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8),
      legend.position = "none",
      strip.text.x = element_text(size=6),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=8))
ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/FigX_stackedbar.png',plot=p,width=2.75,height=1.5,bg='white')

a <- v3.3_final.df %>% dplyr::mutate(total_mcn=m1_mcn+m2_mcn+m3_mcn+m4_mcn) %>% subset(class=='Ext') %>% pull(total_mcn)
b <- v3.3_final.df %>% dplyr::mutate(total_mcn=m1_mcn+m2_mcn+m3_mcn+m4_mcn) %>% subset(class=='Maint') %>% pull(total_mcn)
c <- v10_final.df %>% dplyr::mutate(total_mcn=m1_mcn+m2_mcn+m3_mcn+m4_mcn) %>% subset(class=='Ext') %>% pull(total_mcn)
d <- v10_final.df %>% dplyr::mutate(total_mcn=m1_mcn+m2_mcn+m3_mcn+m4_mcn) %>% subset(class=='Maint') %>% pull(total_mcn)

v3.3_DS_FAF.df$teJumpP <- formatC(as.numeric(v3.3_DS_FAF.df$teJumpP),format='e',digits=1)

FAF_teJumpP_factor_levels <- formatC(c(1.0e-04,5.0e-04,1.0e-03),format='e',digits=1)

v3.3_DS_FAF.df$teJumpP <- factor(v3.3_DS_FAF.df$teJumpP,levels=FAF_teJumpP_factor_levels) 

v3.3_DS_FAF.df$neutP <- as.character(v3.3_DS_FAF.df$neutP)
  
v3.3_DS_FAF.df$neutP <- factor(v3.3_DS_FAF.df$neutP,levels=c(0.010,0.025,0.050,0.075,0.10,0.25,0.50,0.75,0.90,0.925,0.950,0.975,0.99))

v3.3_DS_FAF.df$Rebin_neutP <- unlist(lapply(v3.3_DS_FAF.df$neutP,rebin_neutP))
v3.3_DS_FAF.df$Rebin_neutP <- factor(v3.3_DS_FAF.df$Rebin_neutP,
  levels=c("Low Neutral Insertion Preference","Moderate Neutral Insertion Preference","High Neutral Insertion Preference"))
v3.3_DS_FAF.df$Model <- '1'

v10_DS_FAF.df$teJumpP <- formatC(as.numeric(v10_DS_FAF.df$teJumpP),format='e',digits=1)

v10_DS_FAF.df$teJumpP <- factor(v10_DS_FAF.df$teJumpP,levels=FAF_teJumpP_factor_levels) 

v10_DS_FAF.df$neutP <- as.character(v10_DS_FAF.df$neutP)
  
v10_DS_FAF.df$neutP <- factor(v10_DS_FAF.df$neutP,levels=c(0.010,0.025,0.050,0.075,0.10,0.25,0.50,0.75,0.90,0.925,0.950,0.975,0.99))

v10_DS_FAF.df$Rebin_neutP <- unlist(lapply(v10_DS_FAF.df$neutP,rebin_neutP))
v10_DS_FAF.df$Rebin_neutP <- factor(v10_DS_FAF.df$Rebin_neutP,
  levels=c("Low Neutral Insertion Preference","Moderate Neutral Insertion Preference","High Neutral Insertion Preference"))
v10_DS_FAF.df$Model <- '2'

DS_FAF.df <- rbind(v3.3_DS_FAF.df,v10_DS_FAF.df)

Fig4_Row2_C1 <- DS_FAF.df %>% 
  subset(Rebin_neutP=='Low Neutral Insertion Preference')  %>%
  ggplot(aes(x=MutFreq)) + geom_histogram(binwidth=0.05,fill='#7DBA84') +
  facet_wrap(~fct_rev(teJumpP),nrow=3,scales='free_y')+
  scale_x_continuous(breaks=c(0,0.5,1.0))+
  labs(title='',x='',y='') + 
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    legend.position = "none",
    axis.title.x = element_text(size=14,margin=margin(t=10)),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid=element_blank())

ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/Fig4_Row2_C1.png',plot=Fig4_Row2_C1,width=1.5,height=2.5,bg='white')

Fig4_Row2_C2 <- DS_FAF.df %>% 
  subset(Rebin_neutP=='Moderate Neutral Insertion Preference') %>%
  ggplot(aes(x=MutFreq)) + geom_histogram(binwidth=0.05,fill='#7DBA84') +
  facet_wrap(~fct_rev(teJumpP),nrow=3,scales='free_y')+
  scale_x_continuous(breaks=c(0,0.5,1.0))+
  labs(title='',x='',y='') + 
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    legend.position = "none",
    axis.title.x = element_text(size=14,margin=margin(t=10)),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid=element_blank())

ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/Fig4_Row2_C2.png',plot=Fig4_Row2_C2,width=1.5,height=2.5,bg='white')


Fig4_Row2_C3 <- DS_FAF.df %>% 
  subset(Rebin_neutP=='High Neutral Insertion Preference')  %>%
  ggplot(aes(x=MutFreq)) + geom_histogram(binwidth=0.05,fill='#7DBA84') +
  facet_wrap(~fct_rev(teJumpP),nrow=3,scales='free_y',drop=F)+
  scale_x_continuous(breaks=c(0,0.5,1.0))+
  labs(title='',x='',y='') + 
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    legend.position = "none",
    axis.title.x = element_text(size=14,margin=margin(t=10)),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid=element_blank())

ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/Fig4_Row2_C3.png',plot=Fig4_Row2_C3,width=1.5,height=2.5,bg='white')




###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################

#Fig 4

v3.3_final.df$TE_Bin <- lapply(v3.3_final.df$TotalTE,bin_final_TE_vec)
v3.3_final.df$TE_Bin <- factor(v3.3_final.df$TE_Bin,
  levels=c('< 10','10 - 100','100 - 1,000','1,000 - 10,000','> 10,000'))


maint_v3.3_Total <- v3.3_final.df %>% subset(class=='Maint') %>% droplevels()

binned_maint <- maint_v3.3_Total %>% 
  subset(teJumpP %in% c('1.0e-03','5.0e-04','1.0e-04')) %>%
  group_by(teJumpP,neutP,TE_Bin) %>% dplyr::summarise(BinCount=n()) %>% droplevels()

label_df <- binned_maint %>% group_by(teJumpP,neutP) %>% dplyr::summarise(TotalObs=paste0('N = ',sum(BinCount)),PercObs=1.0)

Fig4_Row1 <- ggplot(binned_maint,aes(x=neutP,y=BinCount)) + facet_wrap(~fct_rev(teJumpP),nrow=3) +
  geom_bar(aes(fill=TE_Bin),position='fill',stat='identity') +
  geom_text(data=label_df,aes(x=neutP,y=1.25,label=TotalObs),size=3)+
  scale_y_continuous(breaks=seq(0,1.0,by=0.2),expand = expansion(mult = c(0, 0.3)))+
  scale_fill_manual(values=c("#6D2F20","#BA5649","#DF8162","#E3A15A"),name='Final TE Copy Number')+
  ggtitle(expression(paste("Mean Final TE Family Copy Number For Dual Survival Outcomes \n"))) +
  xlab('Neutral Insertion Preference') + ylab(expression(paste("\tProportion of \nDual Survival Outcomes"))) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                          axis.title.x = element_text(size=14,margin=margin(t=10)),
                          axis.title.y = element_text(size=14,margin=margin(r=15,l=10)),
                          plot.title = element_text(size = 18,face='bold'),
                          strip.text.x = element_text(size=12),
                          panel.grid=element_blank())+ theme(plot.margin=margin(t=25,l=10))

ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/03.03.22/Fig4_Row1.png',plot=Fig4_Row1,width=10,height=5,bg='white')


v3.3_DS_FAF.df$teJumpP <- formatC(as.numeric(v3.3_DS_FAF.df$teJumpP),format='e',digits=1)

FAF_teJumpP_factor_levels <- formatC(c(1.0e-04,5.0e-04,1.0e-03),format='e',digits=1)

v3.3_DS_FAF.df$teJumpP <- factor(v3.3_DS_FAF.df$teJumpP,levels=FAF_teJumpP_factor_levels) 

v3.3_DS_FAF.df$neutP <- as.character(v3.3_DS_FAF.df$neutP)
  
v3.3_DS_FAF.df$neutP <- factor(v3.3_DS_FAF.df$neutP,levels=c(0.010,0.025,0.050,0.075,0.10,0.25,0.50,0.75,0.90,0.925,0.950,0.975,0.99))

v3.3_DS_FAF.df$Rebin_neutP <- unlist(lapply(v3.3_DS_FAF.df$neutP,rebin_neutP))
v3.3_DS_FAF.df$Rebin_neutP <- factor(v3.3_DS_FAF.df$Rebin_neutP,
  levels=c("Low Neutral Insertion Preference","Moderate Neutral Insertion Preference","High Neutral Insertion Preference"))

Fig4_Row2_C1 <- v3.3_DS_FAF.df %>% 
  subset(Rebin_neutP=='Low Neutral Insertion Preference')  %>%
  ggplot(aes(x=MutFreq)) + geom_histogram(binwidth=0.01,fill='#7DBA84') +
  facet_wrap(~fct_rev(teJumpP),nrow=3,scales='free_y')+
  labs(title='',x='',y='') + 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                          axis.title.x = element_text(size=14,margin=margin(t=10)),
                          strip.background = element_blank(),
                          strip.text.x = element_blank(),
                          panel.grid=element_blank())

ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/03.03.22/Fig4_Row2_C1.png',plot=Fig4_Row2_C1,width=2.5,height=4,bg='white')

Fig4_Row2_C2 <- v3.3_DS_FAF.df %>% 
  subset(Rebin_neutP=='Moderate Neutral Insertion Preference') %>%
  ggplot(aes(x=MutFreq)) + geom_histogram(binwidth=0.01,fill='#7DBA84') +
  facet_wrap(~fct_rev(teJumpP),nrow=3,scales='free_y')+
  labs(title='',x='',y='') + 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                          axis.title.x = element_text(size=14,margin=margin(t=10)),
                          strip.background = element_blank(),
                          strip.text.x = element_blank(),
                          panel.grid=element_blank())

ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/03.03.22/Fig4_Row2_C2.png',plot=Fig4_Row2_C2,width=2.5,height=4,bg='white')


Fig4_Row2_C3 <- v3.3_DS_FAF.df %>% 
  subset(Rebin_neutP=='High Neutral Insertion Preference')  %>%
  ggplot(aes(x=MutFreq)) + geom_histogram(binwidth=0.01,fill='#7DBA84') +
  facet_wrap(~fct_rev(teJumpP),nrow=3,scales='free_y',drop=F)+
  labs(title='',x='',y='') + 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                          axis.title.x = element_text(size=14,margin=margin(t=10)),
                          strip.background = element_blank(),
                          strip.text.x = element_blank(),
                          panel.grid=element_blank())

ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/03.03.22/Fig4_Row2_C3.png',plot=Fig4_Row2_C3,width=2.5,height=4,bg='white')



Fig4_Row2 <- ggplot(v3.3_DS_FAF.df,aes(x=MutFreq)) + geom_histogram(binwidth=0.01,fill='#7DBA84') +
  facet_grid(fct_rev(teJumpP)~Rebin_neutP) +
  ggtitle(expression(paste("Final Allele Frequency of Individual TEs For Dual Survival Outcomes \n"))) +
  xlab('Final Allele Frequency') + ylab(expression(paste("Count"))) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                          axis.title.x = element_text(size=14,margin=margin(t=10)),
                          axis.title.y = element_text(size=14,margin=margin(r=15,l=10)),
                          plot.title = element_text(size = 18,face='bold'),
                          strip.text.x = element_text(size=12),
                          panel.grid=element_blank())+ theme(plot.margin=margin(t=25,l=10))


ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/03.03.22/Fig4_Row2.png',plot=Fig4_Row2,width=10,height=4,bg='white')


del_traject <- v10_sim_data.df %>% subset(neutP %in% c('0.01','0.075','0.5','0.925','0.99')) %>% 
  subset(teJumpP=='2.5e-03') %>% droplevels() %>% subset(class=='Maint')

del_traject$TotalDel  <- del_traject$m2_mcn + del_traject$m3_mcn + del_traject$m4_mcn


Fig4_Row2 <- ggplot(del_traject,aes(x=Gen,y=TotalDel,group=tag,color=m1_mcn))+
  geom_line()+facet_grid(teJumpP~neutP,scales='free_x') + 
  scale_color_viridis(option="viridis",end=0.75,direction=-1,name='Neutral TE \nCopy Number')+
  xlab('Generation') + ylab(expression(paste("Deleterious Mean \nCopy Number"))) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                          axis.title.x = element_text(size=16,margin=margin(t=10)),
                          axis.title.y = element_text(size=16,margin=margin(r=15,l=10)),
                          plot.title = element_blank(),
                          strip.text.x = element_text(size=12),
                          panel.grid=element_blank())+ theme(plot.margin=margin(t=25,l=10))



ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/Fig4_Row2.png',plot=Fig4_Row2,width=16,height=3,bg='white')

Fig4_Row3 <- ggplot(del_traject,aes(x=Gen,y=P1_MeanFit,group=tag,color=m1_mcn))+
  geom_line()+facet_grid(teJumpP~neutP,scales='free_x') + 
  scale_color_viridis(option="viridis",end=0.75,direction=-1,name='Neutral TE \nCopy Number')+
  ggtitle(expression(paste("Host Population Mean Fitness Across Replicates for a Subset of Population Extinction Outcomes \n"))) +
  xlab('Generation') + ylab(expression(paste("Population \nMean Fitness"))) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                          axis.title.x = element_text(size=16,margin=margin(t=10)),
                          axis.title.y = element_text(size=16,margin=margin(r=15,l=10)),
                          plot.title = element_text(size = 22,face='bold'),
                          strip.text.x = element_text(size=12),
                          panel.grid=element_blank())+ theme(plot.margin=margin(t=25,l=10))



ggsave('/home/brandvai/mmunasin/nonWF_SLiM_TE_NoRep/Rplots/Figures/FinalFigures_04.18.23/Fig4_Row3.png',plot=Fig4_Row3,width=16,height=3,bg='white')

