#this shows a graphical illustration of model simplification.
#I made this myself for this specific application, but a generalizable version of it might be useful
my_model_simp_figure<-function(model_type="brms", 
                               models=NULL, return=NULL, title=NULL){
  
  
  models_df<-data.frame()
  for(i in 1:length(models)){
    this_df<-plot_my_coefs(models[[i]], option = "df")
    this_df$model<-paste0("model_", i)
    this_df$model_order<-i
    this_df$N<-nrow(models[[i]]$data)
    models_df<-rbind(models_df, this_df)
  }
  models_df$model_order
  models_df$par<-gsub("b_|_scale", "", models_df$par)
  models_df$par<-gsub("ee_", "STG_", models_df$par)
  models_df$par<-gsub("ee", "STG", models_df$par)
  models_df$par<-gsub(":", " x ", models_df$par)
  models_df$par<-gsub("_", " ", models_df$par)
  models_df$par<-gsub("Torpor2 levelyes", "Torpor", models_df$par)
  models_df$par<-gsub("tbd" ,"bout duration", models_df$par)
  models_df$sig_pol<-"no"
  models_df$sig_pol[models_df$sig=="significant" & models_df$Q50<0]<-"neg_sig"
  models_df$sig_pol[models_df$sig=="significant" & models_df$Q50>0]<-"pos_sig"
  models_df$sig_pol[models_df$sig=="marginal" & models_df$Q50<0]<-"neg_mar"
  models_df$sig_pol[models_df$sig=="marginal" & models_df$Q50>0]<-"pos_mar"
  models_df<-models_df %>% arrange(model, abs(Q50))
  models_df$par<-factor(models_df$par, levels = unique(models_df$par))  
  models_df$par_num<-as.numeric(models_df$par)  
  models_df$sig_pol<-factor(models_df$sig_pol, levels = c("pos_sig","pos_mar" , "no", "neg_mar", "neg_sig"))
  cols=c("pos_sig"="blue",  "pos_mar"="#8599ff", "no"="grey50", "neg_mar"="#ff8d85", "neg_sig"="red")
  
  
  if(return=="plot"){
  return(ggplot(models_df)+
           geom_label(aes(y=model_order, x=par_num, label=round(Q50, 2), color=sig_pol))+
           geom_label(data = models_df %>% group_by(model_order) %>% summarise(N=first(N)), 
                      aes(y=model_order, label=N),  x=1, color="black")+
    scale_color_manual(values=cols, name="", labels=c("significant\n+", "marginal\n+", "not significant", "marginal\n-", "significant\n-"), 
                       guide=guide_legend(override.aes = list(label="*", size=10)))+
      geom_label(x=1, y=2, label="N", color="black")+
    theme_tufte()+
    theme(axis.line.x.top = element_line(color="black"),
          axis.line.y.left = element_line(color="black"), 
          axis.text.x.top = element_text(angle=90),
          legend.text = element_text(size=10, hjust = 0),
          axis.text.x.bottom = element_blank(),
          axis.ticks.x.bottom =  element_blank(), legend.position = "top", plot.margin = margin(0, 1, 0, 0, unit = "in"))+
      coord_cartesian(clip = "off")+
    xlab("")+ylab("")+
      scale_x_reverse(name = "", breaks = length(unique(models_df$par)):1, labels=rev(unique(models_df$par)), sec.axis = dup_axis())+
      #scale_x_reverse(name = "", sec.axis = dup_axis())+
      scale_y_reverse(name = "", breaks = max(models_df$model_order):1, labels=paste0("M", max(models_df$model_order):1))+
      ggtitle(title)
      
  )
  }
  if(return=="df"){return(models_df)}
}