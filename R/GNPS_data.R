
#file <- "F:/桌面/R开发/数据库操作/HMDB.msp"
gnps_ms2 <- function(file){
  library(tidyverse)
  message("正在将msp数据转换Rdata格式")
  msp.data <- readr::read_lines(file,progress = FALSE)
  #show tail data
  #tail(msp.data,100)
  #head(msp.data,100)
  #blank_index <- which(msp.data=="")
  #header with whitespace
  msp.data <- unlist(list("",msp.data))
  
  #检查一下blank是否为偶数
 if(length(blank_index)%%2==!0 ){
   stop("信息缺失")
 }
  num <- 1:(length(blank_index)/2)##代表一共有多少个化合物
  pb <- txtProgressBar(0, length(num), style = 3)
  gnps_ms2 <- lapply(num, function(x){
    setTxtProgressBar(pb, x)
    compound_data <-msp.data[blank_index[-1+2*x]:blank_index[0+2*x]]
    info<- compound_data[grep(": ",compound_data)]
    ms2 <- compound_data[grep("\t",compound_data)]
    #字符串提取info
    
    info_all <- lapply(1:length(info), FUN = function(x){
      
      j <- strsplit(info[x], ": ")
      info <- sapply(j, "[",1)
      value <- sapply(j, "[",2)
      info_all <- cbind(info,value)
    }) %>% do.call(rbind,.) %>% as.data.frame() 
    NAME <- filter(info_all,info_all$info==c("NAME"))[[2]]
    PRECURSORMZ <- filter(info_all,info_all$info==c("PRECURSORMZ"))[[2]]
    PRECURSORTYPE <- filter(info_all,info_all$info==c("PRECURSORTYPE"))[[2]]
    FORMULA <- filter(info_all,info_all$info==c("FORMULA"))[[2]]
    IONMODE <- filter(info_all,info_all$info==c("IONMODE"))[[2]]
    #COLLISIONENERGY <- filter(info_all,info_all$info==c("COLLISIONENERGY"))[[2]]
    compound_info <- data.frame(nameE=NAME ,precursor=PRECURSORMZ,
                                type=PRECURSORTYPE,formula=FORMULA,
                                ionmode=IONMODE)
    ##字符串提取 ms2
    msms <- lapply(1:length(ms2), FUN = function(x){
      j <- strsplit(ms2[x], "\t")
      mz <- sapply(j, "[",1) 
      intensity <- sapply(j, "[",2) 
      #list(data.frame(mz=mz,intensity=intensity))
      cbind(mz,intensity)
    }) %>% do.call(rbind,.) %>% as.data.frame() 
    gnps_ms2 <- compound_info %>% mutate(ms2=list(data.frame(mz=msms$mz,intensity=msms$intensity)))
  }) %>% do.call(rbind,.)
  
  if(dir.exists("./GNPS")==FALSE){
    dir.create("./GNPS")
  }
  save(gnps_ms2,file = "./GNPS/gnps_ms2.Rdata")
  message("数据库构建完成,已保存至默认路径")
  gnps_ms2
}

