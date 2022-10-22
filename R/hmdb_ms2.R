#提高来自hmdb数据库的xml文件
#example: id <- get_hmdb_id(file)
#example: MSMS <- hmdb_ms2(id,2)
file <- "D:/桌面/R/HMDB/hmdb.xml"
file <- "D:/桌面/R/HMDB/csf.xml"
get_hmdb_id <- function(file){
  library(xml2)
  library(tidyverse)
  message("正在读取xml文件")
  hmdb_metabolites <- read_xml(file)
  met.nodes <- xml_find_all(hmdb_metabolites, './/d1:metabolite')
  message("正在获取xml文件中代谢物的HMDBID号")
 # pb <- txtProgressBar(0, length(met.nodes), style = 3)
  HMDB_ID <- lapply(1:length(met.nodes), FUN = function(x){
    #setTxtProgressBar(pb, x)
    HMDB_ID <- xml_find_all(met.nodes[x],'./d1:accession')%>% xml_text()
    message("获取",HMDB_ID,"成功")
    HMDB_ID
  }) %>% do.call(rbind,.)
  HMDB_ID
}

#提取hmdb二级数据

hmdb_ms2 <- function(id_file,sleep_time){
  library(RCurl)
  library(rvest)
  library(curl)
  #先测试ID对应的HMDB数据库网页是否可以
  hmdbid <- id_file
  #获取二级数据网页对应的序列号
  hmdb_info <- lapply(1:length(hmdbid), function(x){
    #setTxtProgressBar(pb, x)
    Sys.sleep(sleep_time)
    message(paste0("正在获取HMDB_ID为：",hmdbid[x]," 的二级数据"))
    url_1 <- "https://hmdb.ca/metabolites/"
    info_1 <- paste0(hmdbid[x],"#spectra")
    new_url_1 <- paste0(url_1,info_1)
    ##获取html网页
    web <- try(read_html(curl(new_url_1,handle = curl::new_handle("useragent" = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/105.0.0.0 Safari/537.36 Edg/105.0.1343.33"))),silent=TRUE)
    if('try-error' %in% class(web)==FALSE){
    ##提取文档中指定元素,提取二级谱图对应的序列号
    news <- web %>% html_nodes("tr td tbody tr td")
    index <- grep("/spectra/ms_ms",news,value = TRUE)
##存在二级谱图
    if(length(index)>0){
      ms2_num <- lapply(1:length(index), FUN = function(i){
        str_num <- strsplit(index[i],"ms/")[[1]][2] %>% strsplit(.,"\"")
        id_num <- str_num[[1]][1]
        #先测试T3BD000001的二级序列号URL是否可用，去除不可用的url
        web <-try(read_html(paste0("http://www.hmdb.ca/spectra/ms_ms/",id_num)),silent=TRUE) 
        if('try-error' %in% class(web)==TRUE){
          id_num <- NULL
        }else{
          id_num <- id_num
        }
      }) %>% unlist()
    
#获取ms2谱图对应的网页url
  ms2_url <- lapply(1:length(ms2_num), FUN = function(j){
        #setTxtProgressBar(pb, x)
        #Sys.sleep(3)
        url_2 <- "http://www.hmdb.ca/spectra/ms_ms/"
        num <- ms2_num[j]#改
        new_url_2 <- paste0(url_2,num)
        #获取HTML网页
        web <- read_html(curl(new_url_2,handle = curl::new_handle("useragent" = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/105.0.0.0 Safari/537.36 Edg/105.0.1343.33")))
        #提取文档中指定元素
        news <- web %>% html_nodes("tbody tr td")
        index <- grep(".txt",news,value = TRUE)[1]  #对应的下载二级数据地址 这边选择第一个
        str_ms <- strsplit(index,"\">Download")[[1]][1]
        ms_url <- strsplit(str_ms,"href=\"")[[1]][2]
      })
      
      ###################获取Experimental Conditions#############################
      ms2_info <- lapply(1:length(ms2_num), FUN = function(k){
        
        url_2 <- "http://www.hmdb.ca/spectra/ms_ms/"
        num <- ms2_num[k]#改
        new_url_2 <- paste0(url_2,num)
        #获取HTML网页
        #web <- read_html(new_url_2)
        
        web <- read_html(curl(new_url_2,handle = curl::new_handle("useragent" = "Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/86.0.4240.198 Safari/537.36")))
        #提取文档中指定元素
        news <- web %>% html_nodes("div table tr")
        titles <- news %>% html_text()
        
        HMDB_ID <- hmdbid[x]
        
        str_name <- grep("Compound name",titles,value = TRUE)
        if(length(str_name)>0){
          Compound_name <- strsplit(str_name,":")[[1]][2]
        }else{
          Compound_name <- "NA"
        }
        
        str_mode <- grep("Ionization Mode",titles,value = TRUE)
        if(length(str_mode)>0){
          Ionization_Mode <- strsplit(str_mode,":")[[1]][2]
        }else{
          Ionization_Mode <- "NA"
        }
        
        str_Energy <- grep("Collision Energy Voltage",titles,value = TRUE)
        if(length(str_Energy)>0){
          Collision_Energy <- strsplit(str_Energy,":")[[1]][2]
        }else{
          Collision_Energy <- "NA"
        }
        
        str_type <- grep("Instrument Type",titles,value = TRUE)
        if(length(str_type)>0){
          Instrument_Type <- strsplit(str_type,":")[[1]][2]
        }else{
          Instrument_Type <- "NA"
        }
        data.frame(HMDB_ID=HMDB_ID,Compound_name=Compound_name,Ionization_Mode=Ionization_Mode,
                   Collision_Energy=Collision_Energy,Instrument_Type=Instrument_Type)
      }) %>% do.call(rbind,.)
      ##获取信息成功
      #获取二级数据text
      ################下载二级碎片离子信息(mz，intensity)
      ##################最终结果为一个csv列表，并保存到默认路径
      msms <- lapply(1:length(ms2_url), function(s){
        num <- ms2_url[[s]]
        if(is.na(num)==TRUE){    ######出现NA是因为网页中文档不可用
          new_ms2 <- "Web page documents are not available"
        }else{
          web_ms2 <- read_html(num)
          new_ms2 <- web_ms2 %>% html_nodes("body") %>% html_text()
        }
        data.frame(ms2=new_ms2)
      }) %>% do.call(rbind,.)
      info_final<- cbind(ms2_info,msms) 
    }else{
      #如果不存在二级谱图
      info_finla <- NULL
      message(paste0(hmdbid[x]),"没有二级谱图,跳过")
    }
  }#第一个try函数 判断url是否存在
    
  }) %>% do.call(rbind,.)
  if(dir.exists("./HMDB")==FALSE){
    dir.create("./HMDB")
  }
  write.csv(hmdb_info,"./HMDB/MS2.csv")
  save(hmdb_info,file = "./HMDB/MS2.Rdata")
 message("获取全部数据完成") 
}
