jrcatcell_idx=which(l1_test$V1==1)
Tsample_idx=setdiff(1:length(l1_test$V1),jrcatcell_idx)

  
  
  
  
  
  
  
  gethetamap_sample<-function(test_cell){
    
    #test_cell0=`test.(1)`
    
    
    
    jrcatsample=list()
    Tsample=list()
    for(i in 1:length(jrcatcell_idx)){
      jrcatsample[[i]]=test_cell[(jrcatcell_idx[i]*2-1):(jrcatcell_idx[i]*2),]
    }
    for(i in 1:length(Tsample_idx)){
      Tsample[[i]]=test_cell[(Tsample_idx[i]*2-1):(Tsample_idx[i]*2),]
    }
    
    s1=0
    s2=0
    for(i in length(jrcatcell_idx)){
      s1=s1+jrcatsample[[i]][1,]
      s2=s2+jrcatsample[[i]][2,]
    }
    s_jrcat=rbind(s1/length(jrcatcell_idx),s2/length(jrcatcell_idx))
    
    s1=0
    s2=0
    for(i in length(Tsample_idx)){
      s1=s1+Tsample[[i]][1,]
      s2=s2+Tsample[[i]][2,]
    }
    s_Tsample=rbind(s1/length(Tsample_idx),s2/length(Tsample_idx))
    
    p=s_jrcat +s_Tsample
    #p[1,]=as.numeric(log(p[1,]/(1-p[1,])))
    #p[2,]=as.numeric(log(p[2,]/(1-p[2,])))
    
    rownames(p)<-c( '293T','Jurkat cell')
    Heatmap(as.matrix(p), name="coupling_coeff",col = structure(brewer.pal(9, "YlOrRd")),rect_gp = gpar(col = "White", lwd = 0.2),cluster_rows = FALSE,cluster_columns = FALSE, show_column_names = TRUE,show_row_names = TRUE, column_names_gp = gpar(fontsize = 6),  row_title = NULL, column_names_side="bottom",show_heatmap_legend = "TRUE")
    
    
    
    # a=array()
    # for(i in 1:10){
    #   a=append(a,as.matrix(rep(uniq_cell[i],13)))
    # }
    # a=a[2:length(a)]
    # a=as.data.frame(a)
    # colnames(a)<-c('sample')
    
    
    
    h1=Heatmap(as.matrix(mat_sample), name="coupling_coeff",col = structure(brewer.pal(9, "YlOrRd")),rect_gp = gpar(col = "White", lwd = 0.2),cluster_rows = FALSE,cluster_columns = FALSE, show_column_names = TRUE,show_row_names = FALSE, column_names_gp = gpar(fontsize = 6),  row_title = NULL, column_names_side="bottom",show_heatmap_legend = "FALSE",heatmap_legend_param = list(title_position = "topcenter",color_bar = "continuous", legend_direction = "horizontal",title_gp = gpar(fontsize = 9)),split=rep(1:10, each =13 ))#h2=Heatmap(as.matrix(a),name="cell_type",col = c("Eryth"=colo[1],"NK"=colo[2],"CD14+ Mono"=colo[3],"Mk"=colo[4], "CD34+" =colo[5], "DC"=colo[6], "Memory CD4 T"=colo[7],"CD8 T"=colo[8],"CD16+ Mono"=colo[9], "B"=colo[10], "T/Mono doublets"=colo[11], "pDCs"=colo[12], "Naive CD4 T"=colo[13]),show_heatmap_legend = "TRUE",column_names_side = "bottom",column_names_rot = 45,column_names_gp = gpar(fontsize = 8),heatmap_legend_param = list(title_position = "topcenter",nrow = 4,title_gp = gpar(fontsize = 9), labels_gp = gpar(col = "black", fontsize = 8)))
    
    #h1=Heatmap(as.matrix(mat_sample[1:13,]), name="coupling_coeff",col = structure(brewer.pal(9, "YlOrRd")),rect_gp = gpar(col = "White", lwd = 0.2),cluster_rows = FALSE,cluster_columns = FALSE, show_column_names = TRUE,show_row_names = FALSE, column_names_gp = gpar(fontsize = 6),  row_title = NULL, column_names_side="bottom",show_heatmap_legend = "FALSE",heatmap_legend_param = list(title_position = "topcenter",color_bar = "continuous", legend_direction = "horizontal",title_gp = gpar(fontsize = 9)))#h2=Heatmap(as.matrix(a),name="cell_type",col = c("Eryth"=colo[1],"NK"=colo[2],"CD14+ Mono"=colo[3],"Mk"=colo[4], "CD34+" =colo[5], "DC"=colo[6], "Memory CD4 T"=colo[7],"CD8 T"=colo[8],"CD16+ Mono"=colo[9], "B"=colo[10], "T/Mono doublets"=colo[11], "pDCs"=colo[12], "Naive CD4 T"=colo[13]),show_heatmap_legend = "TRUE",column_names_side = "bottom",column_names_rot = 45,column_names_gp = gpar(fontsize = 8),heatmap_legend_param = list(title_position = "topcenter",nrow = 4,title_gp = gpar(fontsize = 9), labels_gp = gpar(col = "black", fontsize = 8)))
    
    h2=Heatmap(as.matrix(a1),col = c("Eryth"=colo[1],"NK"=colo[2],"CD14+ Mono"=colo[3],"Mk"=colo[4], "CD34+" =colo[5], "DC"=colo[6], "Memory CD4 T"=colo[7],"CD8 T"=colo[8],"CD16+ Mono"=colo[9], "B"=colo[10], "T/Mono doublets"=colo[11], "pDCs"=colo[12], "Naive CD4 T"=colo[13],"NA"="white"),rect_gp = gpar(col = "Black", lwd = 0.3),show_heatmap_legend = FALSE,column_names_side = "bottom",column_names_rot = 45,column_names_gp = gpar(fontsize = 8))
    
    ht=h1+h2
    
    draw(ht, heatmap_legend_side = "bottom",ht_gap = unit(2, "mm"))
    
    return(mat_sample)
    
  }
  
  
  
  get_matsample<-function(test_cell){
    
    # new_test1=matrix(0,nrow(test_cell),ncol(test_cell))
    # for(i in 1:nrow(test_cell)){#nrow(`test.(1)`)){
    #   print(i)
    #   p=test_cell[i,]
    #   new_test1[i,]=p#as.numeric(log(p/(1-p)))
    # }
    
    
    #new_test1_mat=new_test1[1:nrow(test_cell),]  
    
    
    mat_sample=matrix(0,nrow(test_cell),32)
    for(i in 1:nrow(test_cell)){
      k=1
      for(j in 1:32){
        mat_sample[i,j]=max(test_cell[i,k:(j*73)])
        k=j*73+1
      }
    }
    colnames(mat_sample)<-1:32
    rownames(mat_sample)<-1:nrow(test_cell)
    
    return(mat_sample)
  }
  
  
  
  Heatmap(as.matrix(s_jrcat), name="coupling_coeff",col = structure(brewer.pal(9, "YlOrRd")),rect_gp = gpar(col = "White", lwd = 0.2),cluster_rows = FALSE,cluster_columns = FALSE, show_column_names = TRUE,show_row_names = TRUE, column_names_gp = gpar(fontsize = 6),  row_title = NULL, column_names_side="bottom",show_heatmap_legend = "TRUE")
  
  
  coup_gene_all_jrcat= coup_gene_all_jrcat_upto100#rbind(rbind(coup_gene_all_jrcat_upto100,coup_gene_all_jrcat100to500),coup_gene_all_jrcat500to1000)
  
  k=1
  pricap_15val1=vector()
 for(i in 1:(nrow(coup_gene_all_jrcat)/2)){
   pricap_15val1[i]=coup_gene_all_jrcat[k,8]
   k=k+2
 }
  
  k=2
  pricap_15val2=vector()
  for(i in 1:(nrow(coup_gene_all_jrcat)/2)){
    pricap_15val2[i]=coup_gene_all_jrcat[k,8]
    k=k+2
  }
  pricap_15val=as.data.frame(append(pricap_15val1,pricap_22val2))
  pricap_22val$grp=c(rep('293T',1000),rep('jurkat',1000))
  colnames(pricap_22val)<-c('coupling_coeff','cell_type')
  ggplot(pricap_22val,aes(cell_type,coupling_coeff))+geom_boxplot()
  
  pricap_15val1=as.data.frame(pricap_15val1)
  colnames(pricap_15val1)<-c('var')
  ggplot(pricap_15val1, aes(1:1000,var))+geom_point(col="steelblue")+theme_bw()
  
  
  n_class=3
  s=list()
  mat_final=matrix(0,1,32)
  for(i in 1:n_class){
 s[[i]]=matrix(0,1,32)
    for(j in i:seq(1,(nrow(test)/n_class),3)){
    s[[i]]=s[[i]]+mat_sample[j,]
    }
    s[[i]]=s[[i]]/nrow(test)
  mat_final=rbind(mat_final,s[[i]])
  }
  mat_final=mat_final[2:nrow(mat_final),]
    
  mat_final=log(mat_final/(1-mat_final))
  