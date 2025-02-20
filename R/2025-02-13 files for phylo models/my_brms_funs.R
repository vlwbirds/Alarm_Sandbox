#set some conditions
iter <- 25000
warmup <- 15000
thin <- 10
adapt <- 0.99
chains <- 4
cores <- 4
treedepth <- 11
#and make a useful summary function
plot_my_coefs<-function(x, option){
  df<-posterior_samples(x) %>% 
    select(names(.)[grepl("b_", names(.))]) %>% 
    apply(., 2, FUN=function(x){quantile(x, c(0.025, 0.05, 0.5, 0.9, 0.95))}) %>% 
    t %>% 
    round(., 3)
  colnames(df)<-c("Q2.5", "Q5", "Q50", "Q95", "Q97.5")
  
  df<-df %>% as_tibble() %>% mutate(par=rownames(df)) %>%
    #filter(grepl(par, pattern = "b_"))%>%
    #filter(!grepl(par, pattern = "b_Intercept"))%>%
    mutate(f_sig=ifelse(Q2.5>0 &Q97.5>0|Q2.5<0 & Q97.5<0, 1, 0),
           m_sig=ifelse(Q5>0 &Q95>0|Q5<0 & Q95<0, 1, 0),
           delta_zero=abs(Q50)) %>%
    mutate(sig="not") %>%
    arrange(delta_zero)
  df$sig[df$m_sig==1]<-"marginal"
  df$sig[df$f_sig==1]<-"significant"
  
  
  
  p=ggplot(df, aes(y=Q50, x=factor(par, levels=par)))+
    geom_pointrange(aes(ymax=Q97.5, ymin=Q2.5, color=as.factor(sig)))+
    geom_pointrange(aes(ymax=Q95, ymin=Q5, color=as.factor(sig)), size=1.5)+
    scale_color_manual(name=NULL, values = c("not"="grey80", "significant"="firebrick", "marginal"="#ff9991"))+
    coord_flip()+theme_tufte()+
    theme(legend.position = "none")+
    geom_hline(yintercept = 0, lty=3)+xlab('')+ylab("Posterior Median")
  #+ggtitle(paste0("WAIC = ", round(as.numeric(waic(x)$estimates[3,]),3)))
  if(option=="plot"){return(p)}
  if(option=="df"){return(df[,c(6,3,1,5,2,4, 10)])}
}





my_coef_comp_plot<-function(mod1, mod2, options=c("1v1")){
  if(options=="1vMany"){
  #mod1 and mod2 are brms models
  #mod1 is normal, mod2 is multiphylo
  df=rbind(
    posterior_summary(mod1) %>% 
      as.data.frame() %>% 
      mutate(par=row.names(.),
             Tree=paste0("Consensus Tree")) %>% 
      dplyr::filter(grepl(par, pattern = "b_")),
    
    posterior_summary(mod2) %>% 
      as.data.frame() %>% 
      mutate(par=row.names(.),
             Tree=paste0("Multiphylo")) %>% 
      dplyr::filter(grepl(par, pattern = "b_"))
  ) %>% mutate(par=gsub(par, pattern = "b_|_scaled|_scale", replacement = ""))
  }
  if(options=="1v1"){
    #mod1 and mod2 are brms models
    #mod1 is normal, mod2 is normal
    df=rbind(
      posterior_summary(mod1) %>% 
        as.data.frame() %>% 
        mutate(par=row.names(.),
               Tree=paste0("Tree 1")) %>% 
        dplyr::filter(grepl(par, pattern = "b_")),
      
      posterior_summary(mod2) %>% 
        as.data.frame() %>% 
        mutate(par=row.names(.),
               Tree=paste0("Tree 2")) %>% 
        dplyr::filter(grepl(par, pattern = "b_"))
    ) %>% mutate(par=gsub(par, pattern = "b_|_scaled|_scale", replacement = ""))
  }
  
  
  p=ggplot(df, aes(x=par, y=Estimate, color=as.factor(Tree)))+
    geom_pointrange(aes(ymin=Q2.5, ymax=Q97.5), position = position_dodge(width=0.1))+
    theme_tufte()+
    theme(axis.text.x.bottom = element_text(angle = 90))+
    coord_flip()+theme(legend.title = element_blank())
  return(p)
}





#from https://discourse.mc-stan.org/t/sampling-over-many-phylogenetic-trees-in-brms/5927/24

#this function takes a brms model and a bunch of trees
#it refits the brms model on all of the trees and then combines the results

#this is useful when you've done model selection in brms on a consensus tree
#and then want to refit that final model on a bunch of trees to appreciate
#phylogenetic uncertainty
#choose the trees before fitting the model

my_tree_brms<-function(trees = NA,
                       brms_model=NA, 
                       ncores=NA){
  #loop model
  m.fits <- vector("list", length(trees)) 
  for (i in seq_along(m.fits)) {
    m.fits[[i]] <- update(brms_model,
                          data2 = list(A=ape::vcv.phylo(trees[[i]]),
                                       cores = ncores)
    )
  }
  
  ### combine using combine_models()
  m.fits_comb <- combine_models(m.fits[[i]])
  return(m.fits_comb)
}




#forinvestagting 2 predictors of a brms model at once
my_heatmap_plot<-function(x, vars){
  #x is my brms model
  #vars are my two vars
  coefs<-posterior_summary(x) %>% 
    as.data.frame() %>% mutate(par=row.names(.)) %>% 
    filter(grepl(par, pattern = "b_"))
  intercept<-coefs$Estimate[which(coefs$par=="b_Intercept")]
  coef1<-coefs$Estimate[which(coefs$par==paste0("b_", vars[1]))]
  coef2<-coefs$Estimate[which(coefs$par==paste0("b_", vars[2]))]
  coef3<-coefs$Estimate[grepl(coefs$par, 
                              pattern = paste0(paste(vars, collapse = ":"),
                                               "|", 
                                               paste(rev(vars), collapse = ":")),
  )]
  df=expand.grid(seq(from=min(x$data[,vars[1]], na.rm = T),
                     to=max(x$data[,vars[1]], na.rm = T),
                     by=0.01),
                 seq(from=min(x$data[,vars[2]], na.rm = T),
                     to=max(x$data[,vars[2]], na.rm = T),
                     by=0.01))
  names(df)[1:2]<-c("xvar", "yvar")
  df$z<-intercept + df$xvar*coef1 + df$yvar*coef2 + df$xvar*df$yvar*coef3
  
  p=ggplot(data=df,aes(x=xvar,y=yvar)) +
    geom_tile(aes(fill=z)) + 
    scale_fill_viridis_c(option="magma")+
    xlab(vars[1])+ylab(vars[2])
  return(p)
}
