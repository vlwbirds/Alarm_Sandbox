my_brms_acr_plot_phylosig_fun<-function(model=NULL, tree, data, trait=NULL, trait_se=NULL, my_ylab="", motmot_bayesian=F){
  
  if(is.null(trait) & !is.null(model)){
  ri<-posterior_summary(model) %>% as.data.frame() %>% mutate(par=rownames(.)) %>% filter(par=="b_Intercept") %$% Estimate+
    posterior_summary(model) %>% as.data.frame() %>% mutate(par=rownames(.)) %>% filter(grepl("r_taxon", par)) %$% Estimate
  names(ri)<- sort(unique(model$data$taxon))
  }
  if(!is.null(trait) & is.null(model)){ri<-trait}
  
  tree<-keep.tip(tree, names(ri))
  ri_fit <- phytools::fastAnc(tree,ri,vars=TRUE,CI=TRUE)
  td <- data.frame(node = nodeid(tree, names(ri)),ri = ri)
  nd <- data.frame(node = names(ri_fit$ace), ri = ri_fit$ace)
  d <- rbind(td, nd)
  d$node <- as.numeric(d$node)
  tree_ri <- full_join(tree, d, by = 'node')
  tree_ri@phylo$tip.label<-gsub(tree_ri@phylo$tip.label,pattern = "_",replacement = " ")
  ri_tree_plot<-ggtree(tree_ri, aes(color=ri),ladderize = T, continuous=T, size=3) +
    geom_tiplab(color="black", size=3, align=TRUE, linesize=.5, hjust=-0.1) +
    scale_color_viridis_c(option="inferno", name=my_ylab, direction = 1)+
    theme(legend.position = "left", legend.title = element_text(family = "serif"))+
    xlim(0,50)
  
  
  #phylosig? BLOM K
  
  #now with se
  BlomSE<-NA
  if(is.null(trait) & !is.null(model)){
  ri_se<-posterior_summary(model) %>% as.data.frame() %>% mutate(par=rownames(.)) %>% filter(grepl("r_taxon", par)) %$%  Est.Error 
  names(ri_se)<- sort(unique(model$data$taxon))
  BlomSE<-phylosig(tree=tree, x=ri ,method = "K", test=T, nsim=1000, se=ri_se)
  }
  if(!is.null(trait_se) & is.null(model)){
    BlomSE<-phylosig(tree=tree, x=ri ,method = "K", test=T, nsim=1000, se=trait_se)
  }
 
  
  
  motmot_bayesian_result<-NA
  if(motmot_bayesian){
    set.seed(666) 
    ri_m<-as.matrix(c(ri))
    rownames(ri_m)<-names(ri)
    motmot_bayesian_result <- transformPhylo.MCMC(y = ri_m, phy = tree, 
                                       model = "lambda", mcmc.iteration = 20000, burn.in = 0.5, 
                                       random.start = FALSE, sample.every = 10)
    mcmc.plot(motmot_bayesian_result)
  }
  
  return(list(plot=ri_tree_plot, 
              Blom=phylosig(tree=tree, x=ri ,method = "K", test=T, nsim=1000), 
              BlomSE=BlomSE, 
              motmot_bayesian_result=motmot_bayesian_result[1:4]
))
}