### Functions Used to Plot Results of Simulations

format_summary_df <- function(summary.df) {
	summary.df$teJumpP <- formatC(summary.df$teJumpP,format='e',digits=1)

	teJumpP_factor_levels <- formatC(c(1.0e-04,2.5e-04,5.0e-04,7.5e-04,1.0e-03,2.5e-03,5.0e-03,7.5e-03,1.0e-02,2.5e-02,5.0e-02,7.5e-02,1.0e-01),format='e',digits=1)

	summary.df$teJumpP <- factor(summary.df$teJumpP,levels=teJumpP_factor_levels)	
	
	summary.df$neutP <- as.character(summary.df$neutP)
	
	summary.df$neutP <- factor(summary.df$neutP,levels=c(0.010,0.025,0.050,0.075,0.10,0.25,0.50,0.75,0.90,0.925,0.950,0.975,0.99))

	return(summary.df)

}

ext_df <- function(transformed_log.df) {
	melted_transformed_log.df <- melt(transformed_log.df,id.vars=c("teJumpP",'neutP'),measure.vars=c("Ext_Events","Maint_Events"))

	teJumpP_levels <- levels(melted_transformed_log.df$teJumpP)
	
	reordered_levels <- formatC(sort(as.numeric(as.character(teJumpP_levels)),decreasing=T),format='e',digits=1)
	
	melted_transformed_log.df$teJumpP <- factor(melted_transformed_log.df$teJumpP,levels=reordered_levels)

	ext_melted_transformed_log.df <- melted_transformed_log.df[melted_transformed_log.df$variable=='Ext_Events',]

	ext_teJumpP_levels <- levels(ext_melted_transformed_log.df$teJumpP)
	
	ext_reordered_levels <- formatC(sort(as.numeric(as.character(ext_teJumpP_levels)),decreasing=F),format='e',digits=1)
	
	ext_melted_transformed_log.df$teJumpP <- factor(ext_melted_transformed_log.df$teJumpP,levels=ext_reordered_levels)

	return(ext_melted_transformed_log.df)

}

maint_df <- function(transformed_log.df) {
	melted_transformed_log.df <- melt(transformed_log.df,id.vars=c("teJumpP",'neutP'),measure.vars=c("Ext_Events","Maint_Events"))

	teJumpP_levels <- levels(melted_transformed_log.df$teJumpP)
	
	reordered_levels <- formatC(sort(as.numeric(as.character(teJumpP_levels)),decreasing=T),format='e',digits=1)
	
	melted_transformed_log.df$teJumpP <- factor(melted_transformed_log.df$teJumpP,levels=reordered_levels)

	maint_melted_transformed_log.df <- melted_transformed_log.df[melted_transformed_log.df$variable!='Ext_Events',]

	maint_teJumpP_levels <- levels(maint_melted_transformed_log.df$teJumpP)
	
	maint_reordered_levels <- formatC(sort(as.numeric(as.character(maint_teJumpP_levels)),decreasing=F),format='e',digits=1)
	
	maint_melted_transformed_log.df$teJumpP <- factor(maint_melted_transformed_log.df$teJumpP,levels=maint_reordered_levels)

	return(maint_melted_transformed_log.df)

}

bin_final_TE_df <- function(total_te.vec) {
  final_total_TE <- tail(total_te.vec,n=1)

  if (final_total_TE > 10000) {
    final_TE_bin <- '> 10,000'
  } else if (final_total_TE > 1000 & final_total_TE <= 10000) {
    final_TE_bin <- '1,000 - 10,000'
  } else if (final_total_TE > 100 & final_total_TE <= 1000) {
    final_TE_bin <- '100 - 1,000'
  } else if (final_total_TE > 10 & final_total_TE <= 100) {
    final_TE_bin <- '10 - 100'
  } else {
    final_TE_bin <- '< 10'
  }

  return(final_TE_bin)
}

bin_final_TE_vec <- function(final_total_TE) {

  if (final_total_TE > 10000) {
    final_TE_bin <- '> 10,000'
  } else if (final_total_TE > 1000 & final_total_TE <= 10000) {
    final_TE_bin <- '1,000 - 10,000'
  } else if (final_total_TE > 100 & final_total_TE <= 1000) {
    final_TE_bin <- '100 - 1,000'
  } else if (final_total_TE > 10 & final_total_TE <= 100) {
    final_TE_bin <- '10 - 100'
  } else {
    final_TE_bin <- '< 10'
  }

  return(final_TE_bin)
}

rebin_neutP <- function(neutP_value) {
	if (neutP_value %in% c(0.010,0.025,0.050,0.075)) {
		rebinned_neutP <- 'Low Neutral Insertion Preference'
	} else if (neutP_value %in% c(0.10,0.25,0.50,0.75,0.90)) {
		rebinned_neutP <- 'Moderate Neutral Insertion Preference'
	} else {
		rebinned_neutP <- 'High Neutral Insertion Preference'
	}
}

TE_loss_heatmaps <- function(transformed_log.df,model,full) {
  if (full==TRUE) {
    plot <- ggplot(transformed_log.df) + geom_tile(aes(x=neutP,y=teJumpP,fill=TE_Losses)) +
    labs(title=model,x='Neutral Insertion Preference',
       y='Transposition Probability',fill='TE Loss Outcome')+
    scale_fill_distiller(limits=c(0,1e6),palette='Reds',direction=1,na.value='grey88') +
    theme_minimal() + 
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_blank(),
      strip.text.x = element_text(size=8),
      legend.key.size=unit(0.5,'cm'),
      legend.key.height=unit(0.5,'cm'),
      legend.key.width=unit(1.0,'cm'),
      legend.title=element_text(size=8),
      legend.text=element_text(size=8),
      panel.grid=element_blank())
  } else {

    plot <- transformed_log.df %>% 
    subset(!(teJumpP %in% c('1.0e-01','7.5e-02','5.0e-02'))) %>%
    ggplot() + 
    geom_tile(aes(x=neutP,y=teJumpP,fill=TE_Losses)) +
    labs(title=model,x='Neutral Insertion Preference',
       y='Transposition Probability',fill='TE Loss Outcome')+
    scale_fill_distiller(limits=c(0,1e6),palette='Reds',direction=1,na.value='grey88') +
    theme_minimal() + 
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_blank(),
      strip.text.x = element_text(size=8),
      legend.key.size=unit(0.5,'cm'),
      legend.key.height=unit(0.5,'cm'),
      legend.key.width=unit(0.35,'cm'),
      legend.title=element_text(size=8),
      legend.text=element_text(size=8),
      panel.grid=element_blank())   
  }
  return(plot)
}

ext_ds_heatmaps <- function(transformed_log.df,model,full) {
  ext.df <- ext_df(transformed_log.df)
  maint.df <- maint_df(transformed_log.df)

  sub_ext.df <- ext.df %>% 
  subset(!(teJumpP %in% c('1.0e-01','7.5e-02','5.0e-02')))
  
  sub_maint.df <- maint.df %>% 
  subset(!(teJumpP %in% c('1.0e-01','7.5e-02','5.0e-02')))

  if (full==TRUE) {
    plot <- ggplot() +
    geom_point(data=ext.df,aes(x=neutP,y=teJumpP,color=value),shape="\u25E4",size=7)+
    scale_color_gradient(low='grey100',high='#16697A',name='')+
    new_scale_color()+
    geom_point(data=maint.df,aes(x=neutP,y=teJumpP,color=value),shape="\u25E2",size=7) +
    scale_color_gradient(low='grey100',high='#7DBA84',name='')+
    labs(title='',x='',y='') + 
    theme_minimal() + 
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_blank(),
      strip.text.x = element_text(size=8),
      legend.key.size=unit(0.5,'cm'),
      legend.key.height=unit(0.5,'cm'),
      legend.key.width=unit(1.0,'cm'),
      legend.title=element_text(size=8),
      legend.text=element_text(size=8),
      panel.grid=element_blank())
  } else {
    plot <- ggplot() +
    geom_point(data=sub_ext.df,aes(x=neutP,y=teJumpP,color=value),shape="\u25E4",size=9)+
    scale_color_gradient(low='grey100',high='#16697A',name='')+
    new_scale_color()+
    geom_point(data=sub_maint.df,aes(x=neutP,y=teJumpP,color=value),shape="\u25E2",size=9) +
    scale_color_gradient(low='grey100',high='#7DBA84',name='')+
    labs(title='',x='',y='') + 
    theme_minimal() + 
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_blank(),
      strip.text.x = element_text(size=8),
      legend.key.size=unit(0.5,'cm'),
      legend.key.height=unit(0.5,'cm'),
      legend.key.width=unit(0.5,'cm'),
      legend.title=element_text(size=8),
      legend.text=element_text(size=8),
      panel.grid=element_blank())
  }
  return(plot)
}

establishment_prop_plot <- function(transformed_log.df,main) {
	est_prop.df <- transformed_log.df[,c(1,2,8,9)] %>%
	pivot_longer(cols=c(ExtEvent_Prop,MaintEvent_Prop),names_to='Outcome',values_to='Proportion')
	
	if (main==FALSE) {
		plot <- ggplot(est_prop.df,aes(x=neutP,y=Proportion,fill=Outcome))+
		geom_bar(position='stack',stat='identity')+
		facet_wrap(~rev(teJumpP),nrow=13,ncol=1,scales='free_y')+
		geom_hline(yintercept=1/2000,linetype='dashed', color = "#fc7857") +
		geom_hline(data = . %>% 
			group_by(teJumpP)%>% 
			summarise(tejumpp = mean(as.numeric(as.character(teJumpP)))),
			aes(yintercept = 2*tejumpp), color = "#fc7857")+
		labs(x='',y='',title='')+
		scale_fill_manual(values=c("#16697A","#7DBA84"),labels=c("Population Extinction",'Dual Survival'))+ 
  	theme_bw()+
  	theme(
  		axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8),
    	axis.text.y = element_text(size=8),
    	axis.title.x = element_blank(),
    	axis.title.y = element_blank(),
    	plot.title = element_blank(),
    	strip.background = element_blank(),
    	strip.text.x = element_blank(),
    	legend.title=element_text(size=8),
    	legend.text=element_text(size=8),
    	panel.grid=element_blank())
  	} else {
  		est_prop.df <- est_prop.df %>% subset(teJumpP %in% c('1.0e-04','2.5e-04','5.0e-04','7.5e-04','1.0e-03','2.5e-03'))
		plot <- ggplot(est_prop.df,aes(x=neutP,y=Proportion,fill=Outcome))+
		geom_bar(position='stack',stat='identity')+
		facet_wrap(~rev(teJumpP),nrow=2,scales='free_y')+
		geom_hline(yintercept=1/2000,linetype='dashed', color = "#fc7857") +
		geom_hline(data = . %>% 
			group_by(teJumpP)%>% 
			summarise(tejumpp = mean(as.numeric(as.character(teJumpP)))),
			aes(yintercept = 2*tejumpp), color = "#fc7857")+
		labs(x='',y='',title='')+
		scale_fill_manual(values=c("#16697A","#7DBA84"),labels=c("Population Extinction",'Dual Survival'))+ 
		theme_bw()+
		theme(
			axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8),
    		axis.text.y = element_text(size=6),
    		axis.title.x = element_blank(),
    		axis.title.y = element_blank(),
    		plot.title = element_blank(),
    		legend.position='none',
    		panel.grid=element_blank())  		
  	}
  	return(plot)
}

extinction_time_plot <- function(final.df,main) {
  ext_time.df <- final.df %>% 
  subset(class=='Ext') %>% 
  group_by(teJumpP,neutP) %>% 
  dplyr::summarise(MeanGen=mean(Gen),
    MeanNeut=mean(m1_mcn),
    MeanDel=mean(m2_mcn+m3_mcn+m4_mcn))

  if (main==TRUE) {   
    plot <- ext_time.df %>%
    subset(teJumpP %in% c('2.5e-04','2.5e-03','2.5e-02')) %>%
    ggplot(aes(x=neutP,y=MeanGen)) + 
    geom_bar(stat='identity',fill='#16697A')+
    facet_wrap(~teJumpP,nrow=1)+
    labs(x='',y='')+
    scale_y_continuous(breaks=c(0,25000,50000))+
    theme_bw()+
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_blank(),
      strip.text.x = element_text(size=8),
      panel.grid=element_blank())+ theme(plot.margin=margin(t=25))
  } else {
    plot <- ggplot(ext_time.df,aes(x=neutP,y=MeanGen)) + 
    geom_bar(stat='identity',fill='#16697A')+
    facet_wrap(~teJumpP,scales='free_y',nrow=4)+
    labs(x='',y='')+
    theme_bw()+
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_blank(),
      strip.text.x = element_text(size=12),
      panel.grid=element_blank())+ theme(plot.margin=margin(t=25))
  }

  return(plot)
}

extinction_neutral_count_plot <- function(final.df,main) {
  ext_time.df <- final.df %>% 
  subset(class=='Ext') %>% 
  group_by(teJumpP,neutP) %>% 
  dplyr::summarise(MeanGen=mean(Gen),
    MeanNeut=mean(m1_mcn),
    MeanDel=mean(m2_mcn+m3_mcn+m4_mcn))


  if (main==TRUE) {   
    plot <- ext_time.df %>%
    subset(teJumpP %in% c('2.5e-04','2.5e-03','2.5e-02')) %>%
    ggplot(aes(x=neutP,y=MeanNeut)) + 
    geom_bar(stat='identity',fill='#16697A')+
    facet_wrap(~teJumpP,nrow=1)+
    labs(x='',y='')+
    theme_bw()+
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_blank(),
      strip.text.x = element_text(size=8),
      panel.grid=element_blank())+ theme(plot.margin=margin(t=25))
  } else {
    plot <- ggplot(ext_time.df,aes(x=neutP,y=MeanNeut)) + 
    geom_bar(stat='identity',fill='#16697A')+
    facet_wrap(~teJumpP,scales='free_y',nrow=4)+
    labs(x='',y='')+
    theme_bw()+
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_blank(),
      strip.text.x = element_text(size=12),
      panel.grid=element_blank())+ theme(plot.margin=margin(t=25))
  }

  return(plot)
}

extinction_del_count_plot <- function(final.df,main) {
  ext_time.df <- final.df %>% 
  subset(class=='Ext') %>% 
  group_by(teJumpP,neutP) %>% 
  dplyr::summarise(MeanGen=mean(Gen),
    MeanNeut=mean(m1_mcn),
    MeanDel=mean(m2_mcn+m3_mcn+m4_mcn))


  if (main==TRUE) {   
    plot <- ext_time.df %>%
    subset(teJumpP %in% c('2.5e-04','2.5e-03','2.5e-02')) %>%
    ggplot(aes(x=neutP,y=MeanDel)) + 
    geom_bar(stat='identity',fill='#16697A')+
    facet_wrap(~teJumpP,nrow=1)+
    labs(x='',y='')+
    scale_y_continuous(breaks=c(0,75,150))+
    theme_bw()+
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_blank(),
      strip.text.x = element_text(size=8),
      panel.grid=element_blank())+ theme(plot.margin=margin(t=25))
  } else {
    plot <- ggplot(ext_time.df,aes(x=neutP,y=MeanDel)) + 
    geom_bar(stat='identity',fill='#16697A')+
    facet_wrap(~teJumpP,scales='free_y',nrow=4)+
    labs(x='',y='')+
    theme_bw()+
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_blank(),
      strip.text.x = element_text(size=12),
      panel.grid=element_blank())+ theme(plot.margin=margin(t=25))
  }

  return(plot)
}

del_fit_trajects <- function(sim_data.df,plot_type,class_type) {
  trajects.df <- sim_data.df %>%
  subset(neutP %in% c('0.01','0.075','0.5','0.925','0.99')) %>% 
  subset(teJumpP=='5.0e-03') %>% droplevels() %>% subset(class==class_type) %>%
  dplyr::mutate(TotalDel=m2_mcn+m3_mcn+m4_mcn)

  if (plot_type=='del') {
    plot <- ggplot(trajects.df,aes(x=Gen,y=TotalDel,group=tag,color=m1_mcn))+
    geom_line()+f
    acet_grid(teJumpP~neutP,scales='free_x') +
    scale_color_viridis(option="viridis",end=0.75,direction=-1,name='Neutral TE \nCopy Number')+
    ggtitle(expression(paste("Mean Deleterious TE Copy Number Across Replicates for a Subset of Population Extinction Outcomes \n"))) +
    xlab('Generation') +
    ylab(expression(paste("Deleterious Mean \nCopy Number"))) +
    theme_bw()+
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      axis.title.x = element_text(size=16,margin=margin(t=10)),
      axis.title.y = element_text(size=16,margin=margin(r=15,l=10)),
      plot.title = element_text(size = 22,face='bold'),
      strip.text.x = element_text(size=12),
      panel.grid=element_blank())+ 
    theme(plot.margin=margin(t=25,l=10))
  }

  if (plot_type=='fit') {
    plot <- ggplot(trajects.df,aes(x=Gen,y=P1_MeanFit,group=tag,color=m1_mcn))+
    geom_line()+
    facet_grid(teJumpP~neutP,scales='free_x') +
    scale_color_viridis(option="viridis",end=0.75,direction=-1,name='Neutral TE \nCopy Number')+
    ggtitle(expression(paste("Host Population Mean Fitness Across Replicates for a Subset of Population Extinction Outcomes \n"))) +
    xlab('Generation') + 
    ylab(expression(paste("Population \nMean Fitness"))) +
    theme_bw()+
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      axis.title.x = element_text(size=16,margin=margin(t=10)),
      axis.title.y = element_text(size=16,margin=margin(r=15,l=10)),
      plot.title = element_text(size = 22,face='bold'),
      strip.text.x = element_text(size=12),
      panel.grid=element_blank()) + 
    theme(plot.margin=margin(t=25,l=10))
  }
  return(plot)
}

fitness_trajectories <- function(sim_data.df,model) {
  plot <- ggplot(sim_data.df,aes(x=Gen,y=P1_MeanFit,color=class)) + 
  geom_line() + 
  facet_grid(fct_rev(teJumpP)~neutP) +
  scale_color_manual(values=c('#16697A','#7DBA84'),labels=c("Population Extinction","Dual Survival")) + 
  labs(title=model,
    x='Generation',y='Mean Population Fitness',colour='Outcome') +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.title.x = element_text(size=14,margin=margin(t=10)),
    axis.title.y = element_text(size=14,margin=margin(r=15,l=10)),
    plot.title = element_text(size = 16),
    panel.grid=element_blank(),
    legend.text = element_text(size =  12),
    legend.title = element_text(size = 13, face = "bold"))

  return(plot)
}

final_copy_number_ranges <- function(final_1.df,final_2.df,final_3.df,final_4.df,final_5.df,final_6.df) {
  final_1.df <- final_1.df %>%
  dplyr::mutate(Tot_mcn=m1_mcn+m2_mcn+m3_mcn+m4_mcn) %>%
  select(class,tag,Tot_mcn) %>%
  dplyr::mutate(Model='Model 1')

  final_2.df <- final_2.df %>%
  dplyr::mutate(Tot_mcn=m1_mcn+m2_mcn+m3_mcn+m4_mcn) %>%
  select(class,tag,Tot_mcn) %>%
  dplyr::mutate(Model='Model 2')

  final_3.df <- final_3.df %>%
  dplyr::mutate(Tot_mcn=m1_mcn+m2_mcn+m3_mcn+m4_mcn) %>%
  select(class,tag,Tot_mcn) %>%
  dplyr::mutate(Model='Model 3')

  final_4.df <- final_4.df %>%
  dplyr::mutate(Tot_mcn=m1_mcn+m2_mcn+m3_mcn+m4_mcn) %>%
  select(class,tag,Tot_mcn) %>%
  dplyr::mutate(Model='Model 4')

  final_5.df <- final_5.df %>%
  dplyr::mutate(Tot_mcn=m1_mcn+m2_mcn+m3_mcn+m4_mcn+m5_mcn+m6_mcn+m7_mcn+m8_mcn) %>%
  select(class,tag,Tot_mcn) %>%
  dplyr::mutate(Model='Model 5')

  final_6.df <- final_6.df %>%
  dplyr::mutate(Tot_mcn=m1_mcn+m2_mcn+m3_mcn+m4_mcn+m5_mcn+m6_mcn+m7_mcn+m8_mcn) %>%
  select(class,tag,Tot_mcn) %>%
  dplyr::mutate(Model='Model 6')

  final_all.df <- rbind(final_1.df,final_2.df,final_3.df,final_4.df,final_5.df,final_6.df) %>%
  dplyr::mutate(class=case_when(
    class=='Ext' ~ 'Population Extinction',
    TRUE ~ 'Dual Survival'))

  plot <- ggplot(final_all.df,aes(x=Tot_mcn,fill=class))+
  geom_histogram() +
  scale_fill_manual(values=c('#7DBA84','#16697A')) +
  scale_x_continuous(trans='log10')+
  facet_grid(Model~class)+
  labs(x='Final TE Copy Number',y='Count',fill='Outcome')
  return(plot)

}

final_copy_number_bins <- function(final_1.df,final_2.df,final_3.df,final_4.df,final_5.df,final_6.df,tag) {
  final_1.df <- final_1.df %>%
  dplyr::mutate(Tot_mcn=m1_mcn+m2_mcn+m3_mcn+m4_mcn) %>%
  select(class,tag,Tot_mcn,neutP,teJumpP) %>%
  dplyr::mutate(Model='1')

  final_2.df <- final_2.df %>%
  dplyr::mutate(Tot_mcn=m1_mcn+m2_mcn+m3_mcn+m4_mcn) %>%
  select(class,tag,Tot_mcn,neutP,teJumpP) %>%
  dplyr::mutate(Model='2')

  final_3.df <- final_3.df %>%
  dplyr::mutate(Tot_mcn=m1_mcn+m2_mcn+m3_mcn+m4_mcn) %>%
  select(class,tag,Tot_mcn,neutP,teJumpP) %>%
  dplyr::mutate(Model='3')

  final_4.df <- final_4.df %>%
  dplyr::mutate(Tot_mcn=m1_mcn+m2_mcn+m3_mcn+m4_mcn) %>%
  select(class,tag,Tot_mcn,neutP,teJumpP) %>%
  dplyr::mutate(Model='4')

  final_5.df <- final_5.df %>%
  dplyr::mutate(Tot_mcn=m1_mcn+m2_mcn+m3_mcn+m4_mcn+m5_mcn+m6_mcn+m7_mcn+m8_mcn) %>%
  select(class,tag,Tot_mcn,neutP,teJumpP) %>%
  dplyr::mutate(Model='5')

  final_6.df <- final_6.df %>%
  dplyr::mutate(Tot_mcn=m1_mcn+m2_mcn+m3_mcn+m4_mcn+m5_mcn+m6_mcn+m7_mcn+m8_mcn) %>%
  select(class,tag,Tot_mcn,neutP,teJumpP) %>%
  dplyr::mutate(Model='6')

  if (tag=='by_cn') {
  	final_all.df <- rbind(final_1.df,final_2.df,final_3.df,final_4.df,final_5.df,final_6.df) %>%
  	dplyr::mutate(class=case_when(
  		class=='Ext' ~ 'Population Extinction',
  		TRUE ~ 'Dual Survival')) %>%
  	 dplyr::mutate(final_cn_bin=case_when(
  	 	Tot_mcn < 10 ~ '< 10',
  	 	Tot_mcn >= 10 & Tot_mcn < 100 ~ '10 - 100',
  		Tot_mcn >= 100 & Tot_mcn < 1000 ~ '100 - 1,000',
  		Tot_mcn >= 1000 & Tot_mcn < 10000 ~ '1,000 - 10,000',
  		Tot_mcn >= 10000 ~ '> 10,000')) %>%
  	 dplyr::mutate(final_cn_bin=factor(final_cn_bin,levels=c('< 10','10 - 100','100 - 1,000','1,000 - 10,000','> 10,000'))) %>%
  	 group_by(Model, neutP, teJumpP,final_cn_bin,class) %>%
  	 dplyr::summarise(n=n())

 	plot <- ggplot(final_all.df,aes(x=Model,y=n,fill=final_cn_bin))+
 	geom_bar(position="stack", stat="identity") +
 	scale_fill_manual(values=c("#6D2F20","#BA5649","#DF8162","#E3A15A","#CCC080"),name='Final TE Copy Number')+
 	facet_grid(fct_rev(teJumpP)~neutP) +
 	labs(x='',y='')
  }

  if (tag=='by_outcome') {
  	 final_all.df <- rbind(final_1.df,final_2.df,final_3.df,final_4.df,final_5.df,final_6.df) %>%
  	 dplyr::mutate(class=case_when(
  	 	class=='Ext' ~ 'Population Extinction',
  	 	TRUE ~ 'Dual Survival')) %>%
  	 group_by(Model, neutP, teJumpP,class) %>%
  	 dplyr::summarise(n=n())
  	  
  	 plot <- ggplot(final_all.df,aes(x=Model,y=n,fill=class))+
  	 geom_bar(position="stack", stat="identity") +
  	 scale_fill_manual(values=c('#7DBA84','#16697A')) +
  	 facet_grid(fct_rev(teJumpP)~neutP) +
  	 labs(x='',y='')  	
  }

  return(plot)

}

TE_cn_trajectories <- function(sim_data.df) {
  TE_cn.df <- sim_data.df %>% subset(class=='Maint') %>%
  group_by(tag) %>% 
  dplyr::mutate(Tot_mcn=last(mean_cn)) %>%
  dplyr::mutate(final_cn_bin=case_when(
      Tot_mcn < 10 ~ '< 10',
      Tot_mcn >= 10 & Tot_mcn < 100 ~ '10 - 100',
      Tot_mcn >= 100 & Tot_mcn < 1000 ~ '100 - 1,000',
      Tot_mcn >= 1000 & Tot_mcn < 10000 ~ '1,000 - 10,000',
      Tot_mcn >= 10000 ~ '> 10,000')) %>%
  dplyr::mutate(final_cn_bin=factor(final_cn_bin,levels=c('< 10','10 - 100','100 - 1,000','1,000 - 10,000','> 10,000')))

  p1 <- ggplot(TE_cn.df,aes(x=Gen,y=mean_cn,group=tag,colour=final_cn_bin)) +
  geom_line() + 
  facet_grid(fct_rev(teJumpP) ~ neutP,scales='free_y')+
  scale_colour_manual(values=c("#6D2F20","#BA5649","#DF8162","#E3A15A","#CCC080"),name='Final TE Copy Number')+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(p1)
}

cn_delta_trajectories <- function(eq.csv) {
	teJumpP_factor_levels <- formatC(c(1.0e-04,2.5e-04,5.0e-04,7.5e-04,1.0e-03,2.5e-03,5.0e-03,7.5e-03,1.0e-02,2.5e-02,5.0e-02,7.5e-02,1.0e-01),format='e',digits=1)


	eq <- read_csv(eq.csv)%>%
	separate(col = file,sep="_", into = c("Rep","teJumpP","neutP"))%>% 
	mutate(
		teJumpP = str_remove(teJumpP,"teJumpP="),
		neutP = str_remove(neutP,"neutP="),
		neutP = str_remove(neutP,".txt")) %>% 
	dplyr::mutate(teJumpP2 = factor(formatC(as.numeric(teJumpP),format='e',digits=1),levels=teJumpP_factor_levels))

	plot <- eq %>% 
	group_by(generation, teJumpP ,neutP ) %>%
	summarise(delta = mean(delta)) %>%
    ungroup()%>%
    mutate(teJumpP= fct_rev(teJumpP))%>%
    ggplot(aes(x = generation, y = delta))+
    geom_line()+
  	facet_grid(teJumpP ~ neutP)+
  	geom_smooth(method = "lm",se = FALSE)+
  	geom_hline(yintercept = 0, lty = 2, color = "red")+
  	scale_y_continuous(limits = c(-1.5,1.5))+
  	scale_x_continuous(breaks = c(45000,46000,47000,48000,49000, 50000))+
  	theme(axis.text.x = element_text(angle = 90))+
  	labs(x = "generation (thousands)", y = "change in mean TE copy number")

  	return(plot)
}

final_TEcn_df <- function(final.df) {
  final_TEcn.df <- final.df %>%
  subset(class=='Maint') %>%
  group_by(teJumpP,neutP) %>%
  subset(!(teJumpP %in% c('1.0e-01','7.5e-02','5.0e-02'))) %>%
  dplyr::summarise(MeanCN=mean(m1_mcn+m2_mcn+m3_mcn+m4_mcn),
    MeanFit=mean(P1_MeanFit))

  teJumpP_levels <- levels(final_TEcn.df$teJumpP)
  reordered_levels <- formatC(sort(as.numeric(as.character(teJumpP_levels)),decreasing=F),format='e',digits=1)
  final_TEcn.df$teJumpP <- factor(final_TEcn.df$teJumpP,levels=reordered_levels,ordered=T)
  return(final_TEcn.df)
}