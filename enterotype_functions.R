library(tidyverse)
library(ggtree)
library(compositions)
library(umap)
library(MASS)

### functions 

### enterotype
JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
KLD <- function(x,y) sum(x * log(x/y))
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}


### bulk analysis

ordinate_and_cluster <- function(dat, minPts=20, input_type, n_neighbors = 25){
  require(Rtsne)
  require(umap)
  require(dbscan)
  require(cluster)
  if(input_type == "dist"){
    input = "dist"
    is_distance = TRUE
  }else if(input_type == "matrix"){
    input = "data"
    is_distance = FALSE
  } else {print("input_type must be 'dist' or 'matrix'"); break()}
  dat = as.matrix(dat)
  outlist = list()
  ### UMAP
  print("Calculating UMAP...")
  outlist[["umap"]] = list()
  outlist[["umap"]][["raw"]] = umap(dat, input = input, n_neighbors = n_neighbors)
  outlist[["umap"]][["ordination"]] = outlist[["umap"]][["raw"]]$layout %>% data.frame(., row.names = rownames(dat))
  colnames(outlist[["umap"]][["ordination"]]) = c("Dim1","Dim2")
  outlist[["umap"]][["hdbscan"]] = hdbscan(outlist[["umap"]][["ordination"]], minPts=minPts)
  
  ### tSNE
  print("Calculating tSNE...")
  outlist[["tsne"]] = list()
  outlist[["tsne"]][["raw"]] = Rtsne(dat, is_distance=is_distance)
  outlist[["tsne"]][["ordination"]] = outlist[["tsne"]][["raw"]]$Y %>% data.frame(., row.names = rownames(dat))
  colnames(outlist[["tsne"]][["ordination"]]) = c("Dim1","Dim2")
  outlist[["tsne"]][["hdbscan"]] = hdbscan(outlist[["tsne"]][["ordination"]], minPts=minPts)
  
  ### PCoA - 2 Dim
  print("Calculating PCoA...")
  outlist[["pcoa"]] = list()
  if(is_distance == T){
    outlist[["pcoa"]][["raw"]] = cmdscale(dat)
  }else{
    ### euclidean distance
    outlist[["pcoa"]][["raw"]] = cmdscale(dist(dat))
  }
  outlist[["pcoa"]][["ordination"]] = outlist[["pcoa"]][["raw"]] %>% data.frame(., row.names = rownames(dat))
  colnames(outlist[["pcoa"]][["ordination"]]) = c("Dim1","Dim2")
  outlist[["pcoa"]][["hdbscan"]] = hdbscan(outlist[["pcoa"]][["ordination"]], minPts=minPts)
  
  ### NMDS - 2 Dim
  print("Calculating NMDS...")
  outlist[["nmds"]] = list()
  if(is_distance == T){
    outlist[["nmds"]][["raw"]] = isoMDS(dat, trace=F)
  }else{
    ### euclidean distance
    outlist[["nmds"]][["raw"]] = cmdscale(dist(dat))
  }
  outlist[["nmds"]][["ordination"]] = outlist[["nmds"]][["raw"]]$points %>% data.frame(., row.names = rownames(dat))
  colnames(outlist[["nmds"]][["ordination"]]) = c("Dim1","Dim2")
  outlist[["nmds"]][["hdbscan"]] = hdbscan(outlist[["nmds"]][["ordination"]], minPts=minPts)
  
  ### PAM
  print("Calculating PAM...")
  outlist[["pam"]] = list()
  if(is_distance == T){
    outlist[["pam"]][["all"]] = lapply(2:12, function(k) cluster::pam(dat, k=k))
  }else{
    ### euclidean distance
    outlist[["pam"]][["all"]] = lapply(2:12, function(k) cluster::pam(dist(dat), k=k))
  }
  pam_frame = data.frame(n_cluster=2:12, id=(2:12)-1, avg_width = unlist(lapply(outlist[["pam"]][["all"]], function(x) x$silinfo$avg.width)))
  pam_frame$cand = sapply(1:11, function(x) pam_frame$avg_width[x] < pam_frame$avg_width[x+1])
  pam_best = ifelse(any(pam_frame$cand, na.rm = T), min(which(pam_frame$cand))+1, max(pam_frame$n_cluster)) 
  if(is_distance == T){
    outlist[["pam"]][["best"]] = cluster::pam(dat, k=pam_best)
  }else{
    ### euclidean distance
    outlist[["pam"]][["best"]] = cluster::pam(dist(dat), k=pam_best)
  }
  outlist[["params"]] = list(dat = dat, minPts = minPts, input_type = input_type)
  
  return(outlist)
}

plot_ord <- function(ord_dat, abu_dat, ordination = "umap", abu_thresh = 1, prev_thresh=0.2, maxplot = 20, clustering="pam", k_pam = "best", focus="all", format = "heatmap"){
  colset = data.frame(cluster = paste0("Cluster", 0:12), colors = c("lightgrey", RColorBrewer::brewer.pal(12,"Paired")[c(2,4,6,8,10,12,1,3,5,7,9,11)]))
  if( length(clustering) > 1){
    cluster = clustering
  }else if(clustering=="hdbscan") {cluster = ord_dat[[ordination]][["hdbscan"]]$cluster
  }else if (clustering=="pam") {
    if(k_pam == "best"){ cluster = ord_dat[["pam"]][["best"]]$clustering
    }else{print(paste0("k_pam = ",k_pam)); cluster = ord_dat[["pam"]][["all"]][[k_pam - 1]]$clustering}
  }else {print("ERROR: 'clustering' must be 'pam' or 'hbdscan' or external clustering object"); return(NULL)}
  ord_dat2 = ord_dat[[ordination]] %>% (function(this) bind_cols(this[["ordination"]], cluster=paste0("Cluster",cluster))) %>% rownames_to_column("sample") 
  ord_plot = ord_dat2 %>%
    ggplot(aes(x=Dim1, y=Dim2, col=cluster))  + 
    theme_bw() + theme(panel.grid = element_blank()) + scale_color_manual(values=colset$color, breaks=colset$cluster) + 
    geom_hline(yintercept = 0, lty=2, color="lightgrey") + geom_vline(xintercept = 0, lty=2, color="lightgrey") +
    geom_point() 
  abu_dat_long = abu_dat %>% data.frame() %>% rownames_to_column("sample") %>% reshape2::melt(.) %>% left_join(ord_dat2) %>% filter(cluster!="Cluster0")
  prev_dat = abu_dat_long %>% group_by(cluster, variable) %>% summarize(prev = mean(value>0), abu = mean(value))
  test_vars = prev_dat %>% filter(prev>=!!prev_thresh, abu>=!!abu_thresh) %>% pull(variable) %>% unique %>% as.character()
  print("Calculating cluster-defining differences")
  print(focus)
  if(focus == "all"){
    test_res = lapply(test_vars, function(x) {
      abu_dat_long_sub = abu_dat_long %>% filter(variable == x)
      lapply(unique(abu_dat_long_sub$cluster), function(cl) (abu_dat_long_sub %>% mutate(cluster2 = cluster==cl) %>% lm(value ~ cluster2, data=.) %>% summary() %>% coefficients())[2,] %>% t %>% data.frame(tax=x, cluster=cl) %>% return())
    }) %>% do.call(bind_rows,.)
    if(format == "boxplot") {
      plot_vars = test_res %>% mutate(dir=ifelse(sign(Estimate) < 0, "Depleted","Enrichted")) %>% group_by(cluster, dir) %>% top_n(n=3,wt = -Pr...t..)
      plot_dat_long = abu_dat_long %>% left_join(plot_vars, by=c("variable" = "tax")) %>% filter(!is.na(cluster.y))
      tax_plot_list = lapply(unique(plot_vars$cluster) %>% sort, function(cl) plot_dat_long %>% filter(cluster.y==cl) %>% 
                               ggplot(aes(x=cluster.x, y=value, fill=cluster.x)) + geom_boxplot(outlier.colour = NA) + geom_jitter(width=.25,size=.75, shape=21, stroke=0.1) + 
                               facet_grid(cluster.y ~ dir + variable, scales="free") +
                               theme_bw()+ theme(panel.grid = element_blank(), axis.text.x = element_blank(), legend.position = "none", axis.ticks.x = element_blank()) + 
                               scale_fill_manual(values=colset$color, breaks=colset$cluster) + ylab("Abundance") + xlab(""))
    } else {
      plot_vars = test_res %>% mutate(dir=ifelse(sign(Estimate) < 0, "Depleted","Enrichted")) %>% filter(Estimate > 0) %>% group_by(cluster, dir) %>% top_n(n=10,wt = -Pr...t..)
      plot_dat_long = abu_dat_long %>% filter(variable %in% plot_vars$tax) 
      plot_dat2 = plot_dat_long %>% group_by(cluster, variable) %>% summarise(mean=mean(value)) %>% (function(df) df %>% left_join(df %>% group_by(variable) %>% 
                                                                                                                                     summarise(mean_all = mean(mean), sd_all=sd(mean)))) %>% mutate(Z = (mean - mean_all) / sd_all) %>% 
        left_join(test_res %>% mutate(q = p.adjust(Pr...t.., "bonferroni")), by=c("variable"="tax", "cluster")) %>% 
        mutate(q_stars = case_when(q < 10^-3 ~ "***", q<10^-2 ~ "**", q < 5*10^-2 ~ "*", TRUE ~ ""))
      hc =  plot_dat2 %>% reshape2::dcast(cluster ~ variable, value.var="Z") %>% data.frame(row.names=1) %>% t %>% dist() %>% hclust()
      tax_plot_list = list(plot_dat2 %>% left_join(data.frame(variable = hc$labels, order = hc$order)) %>% mutate(variable = factor(variable, levels=hc$labels[hc$order])) %>%
                             ggplot(aes(y=variable, x=cluster, fill=Z)) + geom_tile() + scale_fill_gradient2() + theme_bw() + geom_text(aes(label=q_stars)))
    }
  }else{
    test_res = lapply(test_vars, function(x) {
      abu_dat_long_sub = abu_dat_long %>% filter(variable == x)
      lapply(focus, function(cl) (abu_dat_long_sub %>% mutate(cluster2 = cluster==cl) %>% lm(value ~ cluster2, data=.) %>% summary() %>% coefficients())[2,] %>% t %>% data.frame(tax=x, cluster=cl) %>% return())}) %>% do.call(bind_rows,.) 
    plot_vars = test_res %>% mutate(dir=ifelse(sign(Estimate) < 0, "Depleted","Enrichted")) %>% group_by(cluster,dir) %>% top_n(n=12,wt = -Pr...t..)
    plot_dat_long = abu_dat_long %>% left_join(plot_vars, by=c("variable" = "tax")) %>% filter(!is.na(cluster.y))
    tax_plot_list = lapply(unique(plot_dat_long$dir) %>% sort, function(tdir) plot_dat_long %>% filter(dir==tdir) %>% 
                             ggplot(aes(x=cluster.x, y=value, fill=cluster.x)) + geom_boxplot(outlier.colour = NA) + geom_jitter(width=.25,size=.75, shape=21, stroke=0.1) + 
                             facet_wrap( ~ variable, scales="free", ncol=6) +
                             theme_bw()+ theme(panel.grid = element_blank(), axis.text.x = element_blank(), legend.position = "none", axis.ticks.x = element_blank()) + 
                             scale_fill_manual(values=colset$color, breaks=colset$cluster) + ylab("Abundance") + xlab(""))
    
  }
  return(cowplot::plot_grid(ord_plot, cowplot::plot_grid(plotlist=tax_plot_list, ncol=1)))
}

###

