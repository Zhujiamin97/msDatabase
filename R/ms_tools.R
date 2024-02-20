plot_match_spec <- function(search.result = NULL,
                            row_num = c(1),
                            ggplotly = FALSE,
                            shiny = FALSE){
  options(warn = -1)
  library(RColorBrewer)
  library(ggrepel)
  if(is.null(search.result$match_matrix)){
    message("No Match!!")
    return(NA)
  }
  if(row_num > length(search.result$match_matrix)){
    message(paste0("The max number of ms/ms is ",length(search.result$match_matrix),"\nreturn row_num = 1"))
    row_num = 1
  }
  match_matrix <- search.result$match_matrix[row_num][[1]]
  if(nrow(match_matrix)==0){
    message("No Match!!")
    return(NA)
  }
  lib.spec = match_matrix %>% filter(int.lib!=0)
  lib.spec = data.frame(mz=lib.spec$mz.lib,
                        int=lib.spec$int.lib,
                        lab=rep("lib.spec",nrow(lib.spec)))
  exp.spec = match_matrix %>% filter(int.exp!=0)
  exp.spec = data.frame(mz=exp.spec$mz.exp,
                        int=-exp.spec$int.exp, # 负值 放到对面
                        lab=rep("exp.spec",nrow(exp.spec)))
  spec.table = rbind(lib.spec, exp.spec)
  match.num = sum(match_matrix$int.exp!=0 & match_matrix$int.lib!=0)
  MFR = paste0(match.num,"/",nrow(lib.spec))
  if(nrow(spec.table)>0){
    p <- ggplot(data=spec.table,mapping=aes(x=mz,y=int))+
      # ylim(-105,105)+
      geom_hline(yintercept = 0, color = "black", linewidth = 0.5) + # add y=0
      # geom_point(aes(color = lab), size = 1) +
      geom_segment(aes(xend = mz, yend = 0,colour = lab), linewidth = 1.2, lineend = "butt")+
      # geom_bar(aes(fill = lab), stat="identity", width=0.15)+
      labs(x = "mz",y = "Relative intensity",title = "Mass spectra of matched fragment ions")+
      theme_bw(base_size = 20,)+
      # geom_ribbon(aes(ymin=0, ymax=int), fill="gray")+ 填充
      theme(legend.position = "none")+
      # scale_y_continuous(limits = c(-105,105),breaks = seq(-105,105,50),labels = c(100,50,0,50,100))+
      scale_color_manual(values = c("#E41A1C","#377EB8"))+
      geom_text_repel(aes(label = round(mz,4), vjust = -0.5, hjust = 0.5), show.legend = TRUE)+
      annotate("text",x= max(spec.table$mz)-5,y= -104,label="Experimental",colour="#E41A1C",size = 8)+
      annotate("text",x= max(spec.table$mz)-5,y= 104,label="Reference",colour="#377EB8",size = 8)+
      annotate("text",x= min(spec.table$mz)+5,y= 104,label=paste0("MFR:",MFR),colour="#377EB8",size = 8)

    p
    
  }
  if(ggplotly){
    return(plotly::ggplotly(p))
  }
  return(p)
}

plot_spec <- function(search.result = NULL,
                      row_num = c(1)){
    if(all(is.na(search.result))){
      message("No search result!")
      return(NA)
    }
    ms2_spec = search.result$ms2
    if(length(ms2_spec)==0){
      message("No MS/MS")
      return(NA)
    }
    if(!is.numeric(row_num)){
      stop("row_num must be number!")
    }
    if(row_num > length(ms2_spec)){
      message(paste0("The max number of ms/ms is ",length(ms2_spec),"\nreturn row_num = 1"))
      row_num = 1
    }
    ms2_spec = ms2_spec[row_num][[1]]
    ms2_spec$intensity = ms2_spec$intensity/max(ms2_spec$intensity)*100
    p <- ggplot(data=ms2_spec,mapping=aes(x=mz,y=intensity))+
         ylim(0,105)+
         geom_segment(aes(xend = mz, yend = 0,colour = "red"), linewidth = 1, lineend = "butt")+
         labs(x = "mz",y = "Relative intensity",title = "Mass spectra of matched fragment ions")+
         theme_bw(base_size = 20,)+
         theme(legend.position = "none")
    return(plotly::ggplotly(p))
}


ms2_matrix_demo <- function(){
  
  ms2_matrix <- matrix(c(
    273.096,22,
    289.086,107,
    290.118,14,
    291.2316,999,
    292.2342,162,
    293.054,34,
    579.169,37,
    580.179,15), ncol=2, byrow=TRUE)
  colnames(ms2_matrix) = c("mz","intensity")
  return(ms2_matrix)
}


matchLibScore <- function(exp.spec,
                          lib.spec,
                          ms2.mz.tol = 0.02,
                          mz.tol.type.ms2 = "Da",
                          sn.thresh = 3,
                          noise = 100,
                          match.type = "DP",
                          min.frag.int = 100, 
                          fraction.weight = 0.2,
                          dp.forward.weight = 0.7,
                          dp.reverse.weight = 0.1){
  library(RaMS)
  # # data.frame
  exp.spec[,1] <- as.numeric(exp.spec[,1])
  exp.spec[,2] <- as.numeric(exp.spec[,2])
  exp.spec <- data.table(exp.spec) %>% filter(intensity >= min.frag.int)
  
  # exp.spec = tools.centroid_spec(exp.spec,ms2_ppm = 10)
  lib.spec[,1] <- as.numeric(lib.spec[,1])
  lib.spec[,2] <- as.numeric(lib.spec[,2])
  # normalize
  exp.spec[,2] <- exp.spec$intensity/max(exp.spec$intensity)*100
  lib.spec[,2] <- lib.spec$intensity/max(lib.spec$intensity)*100
  # clean and centroid
  # filter relative int lower than 1
  lib.spec <- as.data.frame(lib.spec) %>% filter(intensity > 1) # 去除相对强度低于 1%的离子
  if(nrow(exp.spec)==0|nrow(lib.spec)==0){
    score = 0
    match.matrix = data.frame()
  }else{
    # get match matrix
    match.matrix <- matchMS2(exp.spec = exp.spec,
                             lib.spec = lib.spec,
                             mz.tol.type.ms2 = mz.tol.type.ms2,
                             ms2.mz.tol = ms2.mz.tol)
    if(nrow(match.matrix)>0){
      # match fragment number
      fraction <- sum(match.matrix$int.exp!=0 & match.matrix$int.lib!=0)/nrow(lib.spec)
      # Calculate weight score
      dp.forward <- getDP(exp.int = match.matrix$int.exp,
                          lib.int = match.matrix$int.lib)
      dp.reverse <- getDP(exp.int = match.matrix$int.exp[which(match.matrix$int.lib > 0)],
                          lib.int = match.matrix$int.lib[which(match.matrix$int.lib > 0)])
      dp.forward[is.na(dp.forward)] <- 0
      dp.reverse[is.na(dp.reverse)] <- 0
      score <- round(as.numeric(dp.forward * dp.forward.weight + dp.reverse * dp.reverse.weight + fraction * fraction.weight),4)   
    }else{
      score = 0
    }
  }
  return(list(score,match.matrix))
}


# tools for match msms
# 找出数据库和实验数据相匹配的离子
matchMS2 <- function(exp.spec,
                     lib.spec,
                     ms2.mz.tol = 0.02,
                     mz.tol.type.ms2 = "Da"){
  lib_num = 1:nrow(lib.spec)
  match_idx <- lapply(lib_num, function(x){
    lib.mz = lib.spec[x,1]
    if(mz.tol.type.ms2=="Da"){
      mz_range <-  lib.mz + c(-ms2.mz.tol,+ms2.mz.tol)
    }else{
      mz_diff <- ms2.mz.tol*10^(-6)*lib.mz
      mz_range <- lib.mz+c(-mz_diff, +mz_diff)
    }
    # 数据库中单个离子在实验数据中匹配的离子数
    match_exp_spec = exp.spec[exp.spec$mz%between%mz_range,]
    
    if(nrow(match_exp_spec)==1){
      macthResult = cbind(lib.spec[x,],lib_idx = x,
                          match_exp_spec,exp_idx = which(exp.spec$mz==match_exp_spec$mz),
                          mz.error = round(lib.spec[x,1]-match_exp_spec[,1],digits = 4))
      colnames(macthResult) = c("mz.lib","int.lib","idx.lib",
                                "mz.exp","int.exp","idx.exp",
                                "mz.error")
      macthResult
    }else if(nrow(match_exp_spec)>1){
      match_exp_spec = match_exp_spec[which(match_exp_spec$intensity==max(match_exp_spec$intensity)),]
      macthResult = cbind(lib.spec[x,],lib_idx = x,
                          match_exp_spec,exp_idx = which(exp.spec$mz==match_exp_spec$mz),
                          mz.error = round(lib.spec[x,1]-match_exp_spec[,1],digits = 5))
      colnames(macthResult) = c("mz.lib","int.lib","idx.lib",
                                "mz.exp","int.exp","idx.exp",
                                "mz.error")
      macthResult
    }
  }) %>% do.call(rbind,.)
  # 这里的match_idx会出现没有匹配上离子的情况
  # 根据匹配上的序号，就没有匹配上离子的提取出,并合并起来
  # 使用nrow>0判断不行
  if(length(match_idx)>0){
    match.matrix.lib = data.frame(lib.spec[-match_idx$idx.lib,],lib.spec[-match_idx$idx.lib,])
    # 这里去掉了匹配上的离子后，会存在没有离子剩下的情况
    if(nrow(match.matrix.lib)>0){
      match.matrix.lib[,4] = 0
    }
    match.matrix.exp = data.frame(exp.spec[-match_idx$idx.exp,],exp.spec[-match_idx$idx.exp,])
    if(nrow(match.matrix.exp)>0){
      match.matrix.exp[,2] = 0
    }
    match.matrix <- rbind(match.matrix.lib,match.matrix.exp)
    colnames(match.matrix) <- c("mz.lib","int.lib","mz.exp","int.exp")
    macthResult = rbind(match_idx %>% select("mz.lib","int.lib","mz.exp","int.exp"),match.matrix)
  }else{
    macthResult <- data.frame()
  }
  return(macthResult)
}


getDP <- function (exp.int, lib.int){
  exp.weight <- lapply(exp.int, function(x) {
    1/(1 + x/(sum(exp.int) - 0.5))
  }) %>% unlist()
  lib.weight <- lapply(lib.int, function(x) {
    1/(1 + x/(sum(lib.int) - 0.5))
  }) %>% unlist()
  x <- exp.weight * exp.int
  y <- lib.weight * lib.int
  return(sum(x * y)^2/(sum(x^2) * sum(y^2)))
}
