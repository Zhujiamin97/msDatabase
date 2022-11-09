
#file <- "F:/HMDB.msp"
GNPS_ms2 <- function (file){
  library(tidyverse)
  message("正在读取msp数据")
  msp.data <- readr::read_lines(file, progress = FALSE)
  #检查一下文件中名称的形式
  test_index<- grep("NAME",msp.data[1])
  if(length(test_index)>0){
    colname <- c("NAME","PRECURSORMZ","PRECURSORTYPE","FORMULA","IONMODE","COLLISIONENERGY")
    ms2_type <- "2"
  }else{
    colname <- c("Name","PrecursorMZ","Precursor_type","Formula","Ion_mode","Collision_energy")
    ms2_type <- "1"
  }
  message("读取完成,转换中")
  msp.data <- unlist(list("", msp.data))
  blank_index <- which(msp.data == "")
  if (length(blank_index)%%2 == !0) {
    stop("信息缺失")
  }
  num <- 1:(length(blank_index)/2)
  pb <- txtProgressBar(0, length(num), style = 3)
  gnps_ms2 <- lapply(num, function(x) {
    setTxtProgressBar(pb, x)
    compound_data <- msp.data[blank_index[-1 + 2 * x]:blank_index[0 +                                                                2 * x]]
    info <- compound_data[grep(": ", compound_data)]
    info_index <- grep(": ", compound_data)
    #ms2 <- compound_data[grep("\t", compound_data)]
   
    info_all <- lapply(1:length(info), FUN = function(x) {
      j <- strsplit(info[x], ": ")
      info <- sapply(j, "[", 1)
      value <- sapply(j, "[", 2)
      info_all <- cbind(info, value)
    }) %>% do.call(rbind, .) %>% as.data.frame()
    compound_info <- lapply(colname, function(x){
      a <- filter(info_all, info_all$info == x)[[2]]
      if(length(a)==0){
        a="NA"
      }
      a
    }) %>% do.call(cbind,.) %>% as.data.frame()
    colnames(compound_info) <- colname
    if(ms2_type==1){
      ms2 <- compound_data[c(-1,-length(compound_data),-info_index)]
      msms <- lapply(1:length(ms2), FUN = function(x) {
        j <- strsplit(ms2[x], " ")
        mz <- sapply(j, "[", 1)
        intensity <- sapply(j, "[", 2)
        cbind(mz, intensity)
      }) %>% do.call(rbind, .) %>% as.data.frame()
    }else{
      ms2 <- compound_data[grep("\t", compound_data)]
      msms <- lapply(1:length(ms2), FUN = function(x) {
        j <- strsplit(ms2[x], "\t")
        mz <- sapply(j, "[", 1)
        intensity <- sapply(j, "[", 2)
        cbind(mz, intensity)
      }) %>% do.call(rbind, .) %>% as.data.frame()
    }
    
    gnps_ms2 <- compound_info %>% mutate(ms2 = list(data.frame(mz = msms$mz, 
                                                               intensity = msms$intensity)))
  }) %>% do.call(rbind, .)
  if (dir.exists("./GNPS") == FALSE) {
    dir.create("./GNPS")
  }
  save(gnps_ms2, file = "./GNPS/gnps_ms2.Rdata")
  message("数据库构建完成,已保存至默认路径")
  gnps_ms2
}
