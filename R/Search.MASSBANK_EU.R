#'@title Search massbank_eu Database
#'
#'@author JiaMin ZHU
#'@description
#'
#'@param compound_name description
#'@param formula
#'@param InChIKey
#'@param Exact_Mass
#'@param mz_error
#'@param Instrument_Type
#'@param plot_ms2_match
#'@param ms2_matrix
#'
#'@example 
#' use name seacrh : search.massbank_eu(compound_name = "Mabuterol",
#'                                      formula = NULL,
#'                                      Instrument_Type = "LC-ESI-QTOF")

search.massbank_eu <- function(compound_name = "Dihydrotestosterone",
                               formula = NULL,
                               InChIKey = NULL,
                               Exact_Mass = NULL,
                               mz_error = 0.3,
                               MS_Type = "MS2",
                               Instrument_Type = "ESI",
                               ms2.mz.tol = 0.02,
                               mz.tol.type.ms2 = "Da",
                               min.frag.int = 5,
                               plot_ms2_match = FALSE,
                               ms2_matrix = NULL,
                               author = "ZJM"){
  
  if(is.null(InChIKey)){
    url = paste0("https://massbank.eu/MassBank/Result.jsp?compound=",compound_name,
                 "&formula=",formula,
                 "&mz=",Exact_Mass,
                 "&tol=",mz_error,
                 "&op1=and&op2=and&type=quick&searchType=keyword&sortKey=not&sortAction=1&pageNo=1&exec=",
                 "&inst_grp=ESI&inst=CE-ESI-TOF&inst=ESI-ITFT&inst=ESI-ITTOF&inst=ESI-QIT&inst=ESI-QQ&inst=ESI-QTOF&inst=ESI-TOF&inst=LC-ESI-FT&inst=LC-ESI-IT&inst=LC-ESI-ITFT&inst=LC-ESI-ITTOF&inst=LC-ESI-Q&inst=LC-ESI-QFT&inst=LC-ESI-QIT&inst=LC-ESI-QQ&inst=LC-ESI-QQQ&inst=LC-ESI-QTOF&inst=LC-ESI-TOF",
                 "&ms=",MS_Type,
                 "&ion=0")
  }else{
    # search by InChIKey
    url = paste0("https://massbank.eu/MassBank/Result.jsp?inchikey=",InChIKey,
                 "&type=inchikey&searchType=inchikey",
                 "&op1=and&op2=and&sortKey=not&sortAction=1&pageNo=1&exec=&inst_grp=ESI&inst=CE-ESI-TOF&inst=ESI-ITFT&inst=ESI-ITTOF&inst=ESI-QIT&inst=ESI-QQ&inst=ESI-QTOF&inst=ESI-TOF&inst=LC-ESI-FT&inst=LC-ESI-IT&inst=LC-ESI-ITFT&inst=LC-ESI-ITTOF&inst=LC-ESI-Q&inst=LC-ESI-QFT&inst=LC-ESI-QIT&inst=LC-ESI-QQ&inst=LC-ESI-QQQ&inst=LC-ESI-QTOF&inst=LC-ESI-TOF&ms=MS2&ion=0")
  }
  web <- try(rvest::read_html(url))
  ms2_link <- rvest::html_nodes(web, "a") %>% rvest::html_attr("href")
  ms2_url <- grep("RecordDisplay",ms2_link,value = TRUE)
  if(length(ms2_url)==0){
    message("This compound do not have MS/MS!")
    return(NA)
  }
  ms2_urls = paste0("https://massbank.eu/MassBank/",ms2_url)
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
      str = strsplit(x," ")[[1]]
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