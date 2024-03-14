#'@title Search HMDB Database
#'
#'@author JiaMin ZHU
#'@description
#'
#'@param HMDB_ID   "HMDB0000001"
#'@param Exact_Mass
#'@param mz_error
#'@param Instrument_Type
#'@param plot_ms2_match
#'@param ms2_matrix
#'
#'@example 
#' use name seacrh : search.hmdb(HMDB_ID = "HMDB0000001")

search.hmdb <-function(HMDB_ID = "HMDB0000001",
                       Exact_Mass = NULL,
                       mz_ppm = 5,
                       ms2.mz.tol = 0.02,
                       mz.tol.type.ms2 = "Da",
                       min.frag.int = 5,
                       plot_ms2_match = FALSE,
                       ms2_matrix = NULL,
                       author = "ZJM"){
  
  url = "https://hmdb.ca/metabolites/"
  # search by exact_mass
  if(is.numeric(Exact_Mass)&is.numeric(mz_ppm)){
    
    mz_diff <- mz_ppm*10^(-6)*Exact_Mass
    hmdb_id <- c()
    for (page in 1:10000) {
      
      mass_url <- paste0("https://hmdb.ca/structures/search/metabolites/mass?",
                         "&page=",page,
                         "&query_from=",Exact_Mass-mz_diff,
                         "&query_to=",Exact_Mass+mz_diff,
                         "&search_type=monoisotopic")
      web <- try(rvest::read_html(mass_url))
      if(any(class(web) %in% "try-error")){
        message("please check your network!!")
        return(NA)
      }
      # If the connection time is exceeded ,return NA
      # get compounds HMDB ID
      news <- web %>% rvest::html_nodes("a") %>% rvest::html_attr("href")
      str_ids = grep("/metabolites/HMDB",news,value = TRUE,fixed = TRUE)
      if(length(str_ids) == 0){
        break
      }
      ids <- lapply(str_ids, function(x){
                 strsplit(x,"/metabolites/")[[1]][2]
             }) %>% unlist()
      hmdb_id <- c(hmdb_id,ids)
    }
 
  }else if(!is.null(HMDB_ID)){
    hmdb_id <- HMDB_ID
    # hmdb_id <- "HMDB0000001"  # test
  }
  message("The number of results: ",length(hmdb_id))
  # after getting hmdb ids
  id_url <- paste0("https://hmdb.ca/metabolites/",hmdb_id,"#spectra")
  spec_ids <- 
    lapply(id_url, function(x){
    web <- try(rvest::read_html(x))
    if(any(class(web) %in% "try-error")){
      message("please check your network or input information!!")
      return(NULL)
    }
    # If the connection time is exceeded ,return NULL
    # get compounds HMDB ID
    news <- web %>% rvest::html_nodes("a") %>% rvest::html_attr("href")
    ## get ms2
    ms2_url <- grep("/spectra/ms_ms/",news,value = TRUE)[-1]
  }) %>% unlist()
  if(length(spec_ids) == 0){
    message("No MS/MS !")
    return(NA)
  }
  #-------get ms2---------#
  ms2_urls <- paste0("https://hmdb.ca",spec_ids)
  message("Thu number of ms2_url : ", length(ms2_urls))
  ms2_spectra <- lapply(ms2_urls, function(temp_id) {
    print(paste0("ms2_url: ",temp_id))
    html_document <- rvest::read_html(temp_id)
    link <- html_document %>% rvest::html_element("tr:nth-child(1) a") %>%
      rvest::html_attr("href")
    # 存在谱图但没有对应的文件
    ms2 <- tryCatch(read.table(link, header = FALSE), error = function(e) {
      #get msms
      web_text <- html_document %>% html_nodes("body")%>% html_text()
      msms <- strsplit(web_text,"x\":")[[1]][-1]
      #number of fragments
      num_f <- length(msms)
      if(num_f>0){
        MS2 <- lapply(seq(num_f), function(x){
          f <- msms[x]
          #split mz and intensity
          str1 <- strsplit(f,",\"y\":")
          mz <- as.numeric(str1[[1]][1]) %>% round(.,digits = 4)
          intensity <- as.numeric(strsplit(str1[[1]][2],",\"")[[1]][1])%>% round(.,digits = 4)
          ms2 <- data.frame(mz,intensity)
        }) %>% do.call(rbind,.) %>% unique()  #还需要去重复的
      }else{
        NULL
      }
    })
    if (!is.null(ms2)) {
      colnames(ms2) <- c("mz", "intensity")
      # return(NULL)
    }
    componud_info <- html_document %>% rvest::html_table()
    componud_info <- rbind(componud_info[[1]], componud_info[[2]])
    ms1_info <- matrix(data = componud_info$X2 ,nrow = 1,ncol = nrow(componud_info)) %>% as.data.frame()
    colnames(ms1_info) <- lapply(componud_info$X1,function(x){
      name_str <- strsplit(x,":")[[1]][1]
      gsub(" ","_",name_str)
    }) %>% unlist()
    #
    if(!is.null(ms1_info$Collision_Energy)){
      Collision_Energy = ms1_info$Collision_Energy
    }else if(!is.null(ms1_info$Collision_Energy_Voltage)){
      Collision_Energy = ms1_info$Collision_Energy_Voltage
    }else{
      Collision_Energy = NA
    }

 
    ms1_info = data.frame(HMDB_ID = ms1_info$HMDB_ID,
                          Compound_name = ms1_info$Compound_name,
                          Spectrum_type = ifelse(is.null(ms1_info$Spectrum_type),
                                                 NA,ms1_info$Spectrum_type),
                          Ionization_Mode = ifelse(is.null(ms1_info$Ionization_Mode),
                                                   NA,ms1_info$Ionization_Mode),
                          Instrument_type = ifelse(is.null(ms1_info$Instrument_Type),
                                                   NA,ms1_info$Instrument_Type),
                          Collision_Energy = Collision_Energy
                          )
    as.data.frame(ms1_info) %>% mutate(ms2=list(ms2))
  }) %>% do.call(rbind,.) 
  
  if(is.null(ms2_matrix)){
    message("No MS/MS matrix provided!")
    return(ms2_spectra)
  }
  #match ms2
  message("Matching ms/ms...")
  ms2_spectra_match <- data.table::data.table()
  for (lib.spec in ms2_spectra$ms2) {
    
    scoreAndmatch.matrix <- matchLibScore(exp.spec = ms2_matrix,
                                          lib.spec = lib.spec,
                                          ms2.mz.tol = ms2.mz.tol,
                                          mz.tol.type.ms2 = mz.tol.type.ms2,
                                          min.frag.int = min.frag.int)
    match_result <- data.frame(dp.score = scoreAndmatch.matrix[[1]]) %>% 
      mutate(match_matrix = list(scoreAndmatch.matrix[[2]]))
    ms2_spectra_match <- rbind(ms2_spectra_match,match_result)
  }
  return(cbind(ms2_spectra,ms2_spectra_match))
}
