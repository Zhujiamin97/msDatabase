# start_page_num:The page number you want to start crawling
# end_page_num:The page number you want to end the crawl
# name:The class of substance, such as "pollutant"

get_t3db <- function(start_page_num=c("1"),end_page_num=c("2"),
                     name=c("pollutant"),sleep_time=c("2")){
  library(tidyverse)
  library(rvest)
  library(curl)
  library(RCurl)
  message("In order to avoid website denials, it may take some time!")
  name_list <- c("pollutant","Fungal Toxin","Synthetic Toxin","Uremic Toxin")
  if(name %in% name_list==FALSE){
    stop("This compound category does not exist with the T3DB database")
  }
  url <- "http://www.t3db.ca/categories?c=title&d=up&filter=true"
  new_url <- paste0(url,"&",name,"=1")
  
  num <- c(start_page_num:end_page_num)
  
  pb <- txtProgressBar(0, length(num), style = 3)
  t3db_data <- lapply(num,FUN = function(x){
    setTxtProgressBar(pb, x)
    Sys.sleep(sleep_time)
    #url <- "http://www.t3db.ca/categories?c=title&d=up&filter=true&pollutant=1"
    page_num <- paste0("&page=",x)
    new_url <- paste0(new_url,page_num)
    #读取html网页
    #web <- read_html(new_url)
    web <- read_html(curl(new_url,handle = curl::new_handle("useragent" = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/105.0.0.0 Safari/537.36 Edg/105.0.1343.33")))
    #提取文档中指定元素。节点的部分
    news <- web %>% html_nodes("tr td")
    #提取标签内的文本
    titles <- news %>% html_text()
    #length(title)
    t3dbid_index <- grep("T3D",titles)
    name_index <- t3dbid_index+1
    formula_index <- t3dbid_index+2
    data.frame(t3dbid=titles[t3dbid_index],name=titles[name_index],formula=titles[formula_index])
  }) %>% do.call(rbind,.)

  close(pb)
  if(dir.exist("./T3DB")==FALSE){
dir.create("./T3DB")
  }
  file.name <- paste0(name,"_t3dID")
  write.csv(t3db_data,"./T3DB/file.name.csv")
  message(paste0("获取",name,"的T3D_id成功，已保存到默认路径"))
#出问题 是T3D0119这个网页无法打开 待修复  在输出id时先测试对应网页是否可用
#测试可以用的T3DB_ID网页
id<- t3db_data$t3dbid
message("正在测试可查询的T3DB_ID网页")
t3dbid <- lapply(1:length(id), function(x){
    url_0 <- "http://www.t3db.ca/toxins/"
    info_0 <- paste0(id[x],"#spectra")
    new_url_0 <- paste0(url_0,info_0)
    if(url.exists(new_url_0)==TRUE){
      url <- id[x]
    }
  }) %>% unlist()
message("测试可查询的T3DB_ID网页结束")
#获取二级数据的url
  #t3dbid <- t3db_data$t3dbid
  #pb <- txtProgressBar(0, length(t3dbid), style = 3)
t3d_info <- lapply(1:length(t3dbid), function(x){
    #setTxtProgressBar(pb, x)
    Sys.sleep(sleep_time)
    message(paste0("正在获取T3DB_ID为：",t3dbid[x]," 的二级数据"))
    url_1 <- "http://www.t3db.ca/toxins/"
    info_1 <- paste0(t3dbid[x],"#spectra")
    new_url_1 <- paste0(url_1,info_1)
    ##获取html网页
    #System.time()  #中间放进去的时间单位为秒

    web <- read_html(curl(new_url_1,handle = curl::new_handle("useragent" = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/105.0.0.0 Safari/537.36 Edg/105.0.1343.33")))

    # web <- read_html(new_url_1)
    ##提取文档中指定元素,提取二级谱图对应的序列号
    news <- web %>% html_nodes("tr td tbody tr td")
    index <- grep("/spectra/ms_ms",news,value = TRUE)
    ##存在二级谱图
    if(length(index)>0){
ms2_num <- lapply(1:length(index), FUN = function(i){
        str_num <- strsplit(index[i],"ms/")[[1]][2] %>% strsplit(.,"\"")
        id_num <- str_num[[1]][1]
        #先测试T3BD000001的二级序列号URL是否可用，去除不可用的url
        if(url.exists(paste0("http://www.t3db.ca/spectra/ms_ms/",id_num))==FALSE){
          id_num <- NULL
        }else{
          id_num <- id_num
        }
      }) %>% unlist()
##去除ms2_num中为NA的项
#   if(length(which(ms2_num=="NULL"))>0){
#            ms2_num <- ms2_num[-which(ms2_num=="NULL")]  ##500问题所在 因为网页确实无法打开
#            }else{
#             ms2_num <- ms2_num
#        } 
#获取T3BD000001的二级序列号成功，为ms2_num
#获取二级碎片的URL
ms2_url <- lapply(1:length(ms2_num), FUN = function(j){
  #setTxtProgressBar(pb, x)
  #Sys.sleep(3)
  url_2 <- "http://www.t3db.ca/spectra/ms_ms/"
  num <- ms2_num[j]#改
  new_url_2 <- paste0(url_2,num)
  #获取HTML网页
  web <- read_html(curl(new_url_2,handle = curl::new_handle("useragent" = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/105.0.0.0 Safari/537.36 Edg/105.0.1343.33")))
  #提取文档中指定元素
  news <- web %>% html_nodes("tbody tr td")
  index <- grep("http",news,value = TRUE)[1]
  str_ms <- strsplit(index,"\">Download")[[1]][1]
  ms_url <- strsplit(str_ms,"href=\"")[[1]][2]
})
  #获取T3BD000001的二级url成功，为ms2_url
  ###################获取Experimental Conditions#############################
  #把所有有二级数据的T3DB号找出来，最终id号作为媒介与物质信息表进行合并
  ##数量过大 网站容易拒绝访问
  ms2_info <- lapply(1:length(ms2_num), FUN = function(k){
    
    url_2 <- "http://www.t3db.ca/spectra/ms_ms/"
    num <- ms2_num[k]#改
    new_url_2 <- paste0(url_2,num)
    #获取HTML网页
    #web <- read_html(new_url_2)
    
    web <- read_html(curl(new_url_2,handle = curl::new_handle("useragent" = "Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/86.0.4240.198 Safari/537.36")))
    #提取文档中指定元素
    news <- web %>% html_nodes("div table tr")
    titles <- news %>% html_text()
    
    T3DB_ID <- t3dbid[x]
    
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
    
    str_Energy <- grep("Collision",titles,value = TRUE)
    if(length(str_Energy)>0){
      Collision_Energy <- strsplit(str_Energy,":")[[1]][2]
    }else{
      Collision_Energy <- "NA"
    }
    
    str_Formula <- grep("Molecular Formula",titles,value = TRUE)
    if(length(str_Formula)>0){
      Molecular_Formula <- strsplit(str_Formula,":")[[1]][2]
    }else{
      Molecular_Formula <- "NA"
    }
    
    str_mass <- grep("Monoisotopic Mass",titles,value = TRUE)
    if(length(str_mass)>0){
      Monoisotopic_Mass <- strsplit(str_mass,":")[[1]][2]
    }else{
      Monoisotopic_Mass <- "NA"
    }
    
    str_type <- grep("Instrument Type",titles,value = TRUE)
    if(length(str_type)>0){
      Instrument_Type <- strsplit(str_type,":")[[1]][2]
    }else{
      Instrument_Type <- "NA"
    }
    data.frame(T3DB_ID=T3DB_ID,Compound_name=Compound_name,Ionization_Mode=Ionization_Mode,
               Collision_Energy=Collision_Energy,Molecular_Formula=Molecular_Formula,
               Monoisotopic_Mass=Monoisotopic_Mass,Instrument_Type=Instrument_Type)
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
message(paste0(t3dbid[x]),"没有二级谱图,跳过")
      
   }
}) %>% do.call(rbind,.)
message("全部数据获取完成")
write.csv(t3d_info,"./file.name.csv")
save(t3d_info,file = "./file.name.Rdata")
t3d_info
}





