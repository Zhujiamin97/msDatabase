#'@title Search massbank_eu Database
#'
#'@author JiaMin ZHU
#'@description
#' MASSBANK_EU: https://massbank.eu/MassBank/
#'@param compound_name 
#'@param formula
#'@param InChIKey
#'@param Exact_Mass
#'@param mz_error
#'@param Instrument_Type
#'@param plot_ms2_match
#'@param ms2_matrix for example use: ms2_matrix_demo()
#'
#'@example 
#' use name seacrh : search.massbank_eu(compound_name = "Dihydrotestosterone",
#'                                      formula = NULL,
#'                                      Instrument_Type = "LC-ESI-QTOF",
#'                                      ms2_matrix = ms2_matrix_demo())

search.massbank_eu <- function(compound_name = "Dihydrotestosterone",
                               formula = NULL,
                               InChIKey = NULL,
                               Exact_Mass = NULL,
                               mz_error = 0.3,
                               Ion_Mode = c("all","P","N"),
                               MS_Type = "MS2",
                               Instrument_Type = "ESI",
                               ms2.mz.tol = 0.02,
                               mz.tol.type.ms2 = "Da",
                               min.frag.int = 5,
                               ms2_matrix = NULL,
                               author = "ZJM"){
  # Ion_Mode
  if(all(Ion_Mode %in% c("all","P","N"))){
    Ion_Mode ="all"
  }
  if(Ion_Mode=="all"){
    ion = 0
  }else if(Ion_Mode=="P"){
    ion = 1
  }else if(Ion_Mode=="N"){
    ion = -1
  }
  
  if(is.null(InChIKey)){
    # 将空白替换为+号以满足url搜索需求
    compound_name = gsub("\\s", "+", compound_name)
    
    url = paste0("https://massbank.eu/MassBank/Result.jsp?compound=",compound_name,
                 "&formula=",formula,
                 "&mz=",Exact_Mass,
                 "&tol=",mz_error,
                 "&op1=and&op2=and&type=quick&searchType=keyword&sortKey=not&sortAction=1&pageNo=1&exec=",
                 "&inst_grp=ESI&inst=CE-ESI-TOF&inst=ESI-ITFT&inst=ESI-ITTOF&inst=ESI-QIT&inst=ESI-QQ&inst=ESI-QTOF&inst=ESI-TOF&inst=LC-ESI-FT&inst=LC-ESI-IT&inst=LC-ESI-ITFT&inst=LC-ESI-ITTOF&inst=LC-ESI-Q&inst=LC-ESI-QFT&inst=LC-ESI-QIT&inst=LC-ESI-QQ&inst=LC-ESI-QQQ&inst=LC-ESI-QTOF&inst=LC-ESI-TOF",
                 "&ms=",MS_Type,
                 "&ion=",ion)
  }else{
    # search by InChIKey
    url = paste0("https://massbank.eu/MassBank/Result.jsp?inchikey=",InChIKey,
                 "&type=inchikey&searchType=inchikey",
                 "&op1=and&op2=and&sortKey=not&sortAction=1&pageNo=1&exec=&inst_grp=ESI&inst=CE-ESI-TOF&inst=ESI-ITFT&inst=ESI-ITTOF&inst=ESI-QIT&inst=ESI-QQ&inst=ESI-QTOF&inst=ESI-TOF&inst=LC-ESI-FT&inst=LC-ESI-IT&inst=LC-ESI-ITFT&inst=LC-ESI-ITTOF&inst=LC-ESI-Q&inst=LC-ESI-QFT&inst=LC-ESI-QIT&inst=LC-ESI-QQ&inst=LC-ESI-QQQ&inst=LC-ESI-QTOF&inst=LC-ESI-TOF",
                 "&ms=",MS_Type,
                 "&ion=",ion)
  }
  web <- try(rvest::read_html(url),silent = TRUE)
  if(all(class(web) %in% "try-error")){
    message("Get url maybe go wrong or check your input!")
    return(NA)
  }
  ms2_link <- rvest::html_nodes(web, "a") %>% rvest::html_attr("href")
  ms2_url <- grep("RecordDisplay",ms2_link,value = TRUE)
  if(length(ms2_url)==0){
    message("This compound do not have MS/MS!")
    return(NA)
  }
  ms2_urls = paste0("https://massbank.eu/MassBank/",ms2_url)
  message(paste0("The number of ms/ms:",length(ms2_urls),"\nGetting MS/MS..."))
  ms2_spectra <- lapply(ms2_urls, function(x){

    html_document <- readLines(paste0(x, ".xml"), warn = FALSE)
    # get ms1 infomation
    CH = grep("CH$",html_document,value = TRUE,fixed = TRUE)
    AC = grep("AC$",html_document,value = TRUE,fixed = TRUE)
    MS = grep("MS$",html_document,value = TRUE,fixed = TRUE)
    ms1_xml <- c(CH,AC,MS)
    entries <-  lapply(ms1_xml, function(x){
      
      str1 <- strsplit(x,"$",fixed = TRUE)[[1]][2]
      key <- strsplit(str1,":</b> ",fixed = TRUE)[[1]][1]
      value <- strsplit(str1,":</b> ",fixed = TRUE)[[1]][2]
      value <- strsplit(value,"<br>",fixed = TRUE)[[1]][1]
      data.frame(key,value)
      
    }) %>% do.call(rbind,.)
    #-----------------------#
    Name = entries %>% filter(key == "NAME")
    Name = paste0(Name[,2],collapse = ";")
    Compound_class = entries %>% filter(key == "COMPOUND_CLASS")
    ExactMass = entries %>% filter(key == "EXACT_MASS")
    SMILES = entries %>% filter(key == "SMILES")
    Instrument_type = entries %>% filter(key == "INSTRUMENT_TYPE")
    Instrument = entries %>% filter(key == "INSTRUMENT")
    chromatography = entries %>% filter(key == "CHROMATOGRAPHY")
    chromatography = paste0(chromatography[,2],collapse = ";")
    
    mass_spectrometry = entries %>% filter(key == "MASS_SPECTROMETRY")
    mass_spectrometry <- lapply(mass_spectrometry[,2], function(x){
      # 使用sub函数替换第一个空白为分号（;），然后进行分割  
      str <- strsplit(sub("^([^ ]*) (.*)$", "\\1;\\2", x), ";")[[1]] 
      # str = strsplit(x," ")[[1]]
      value = data.frame(str[2])
      colnames(value) = str[1]
      value 
    }) %>% do.call(cbind,.)
    
    FOCUSED_ION = entries %>% filter(key == "FOCUSED_ION")
    FOCUSED_ION <- lapply(FOCUSED_ION[,2], function(x){
      str = strsplit(x," ")[[1]]
      value = data.frame(str[2])
      colnames(value) = str[1]
      value 
    }) %>% do.call(cbind,.)
    #-----------------------#
    ms1_info <- data.table::data.table(Name = Name,
                                       Compound_class = Compound_class[,2],
                                       ExactMass = as.numeric(ExactMass[,2]),
                                       Instrument_type = Instrument_type[,2],
                                       Instrument = Instrument[,2],
                                       chromatography = chromatography,
                                       mass_spectrometry,
                                       FOCUSED_ION)
    
    ms1_info <- data.table::data.table(Name = ms1_info$Name,
                                       Compound_class = ifelse(is.null(ms1_info$Compound_class),
                                                               NA,ms1_info$Compound_class),
                                       ExactMass = ifelse(is.null(ms1_info$ExactMass),
                                                          NA,ms1_info$ExactMass),
                                       PrecursorMZ = ifelse(is.null(ms1_info$`PRECURSOR_M/Z`),
                                                            NA,ms1_info$`PRECURSOR_M/Z`),
                                       Precursor_type = ifelse(is.null(ms1_info$PRECURSOR_TYPE),
                                                               NA,ms1_info$PRECURSOR_TYPE),
                                       Ion_mode = ifelse(is.null(ms1_info$ION_MODE),
                                                         NA,ms1_info$ION_MODE),
                                       Ionization = ifelse(is.null(ms1_info$IONIZATION),
                                                           NA,ms1_info$IONIZATION),
                                       Fragmentation_mode = ifelse(is.null(ms1_info$FRAGMENTATION_MODE),
                                                                   NA,ms1_info$FRAGMENTATION_MODE),
                                       COLLISION_ENERGY = ifelse(is.null(ms1_info$COLLISION_ENERGY),
                                                                 NA,ms1_info$COLLISION_ENERGY),
                                       Instrument_type = ifelse(is.null(ms1_info$Instrument_type),
                                                                NA,ms1_info$Instrument_type),
                                       Instrument = ifelse(is.null(ms1_info$Instrument),
                                                           NA,ms1_info$Instrument),
                                       chromatography = ifelse(is.null(ms1_info$chromatography),
                                                               NA,ms1_info$chromatography),
                                       MS_type = ifelse(is.null(ms1_info$MS_TYPE),
                                                        NA,ms1_info$MS_TYPE),                                       
                                       RESOLUTION = ifelse(is.null(ms1_info$RESOLUTION),
                                                           NA,ms1_info$RESOLUTION),
                                       BASE_PEAK = ifelse(is.null(ms1_info$BASE_PEAK),
                                                          NA,ms1_info$BASE_PEAK))
    #get ms2
    end_idx = which(html_document == "//")-1
    strat_idx = which(html_document == "<b>PK$PEAK:</b> m/z int. rel.int.<br>")+1
    ms2_xml = html_document[strat_idx:end_idx]
    ms2 <- lapply(ms2_xml, function(str){
      # 提取字符串中的数字
      extracted <- str_extract_all(str, "\\d+(\\.\\d+)?(\\s+\\d+)?")[[1]] 
      ms2 = data.frame(mz = as.numeric(extracted[1]),
                       intensity = as.numeric(extracted[2]))
    }) %>% do.call(rbind,.)
    compound = ms1_info %>% mutate(ms2 = list(ms2))
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
