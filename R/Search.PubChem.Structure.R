#'@title Search.PubChem.Structure
#'
#'@author JiaMin ZHU
#'@description 2024/11/19
#'

Search.PubChem.Structure <- function(cids = 2244,
                                     similarity = c("2d","3d"),
                                     Threshold = 99){
  
  similarity <- paste0(similarity,collapse = "&")
  
  if(all(similarity == "2d")){
    message("Only searching 2d...")
    
    url0 <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/cid/"
    
    url <- paste0(url0,
                  cids,
                  "/property/MolecularWeight,MolecularFormula,RotatableBondCount/JSON?",
                  "Threshold=",Threshold)
    
  }else if(all(similarity == "3d")){
    message("Only searching 3d...")
    
    url0 <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_3d/cid/"
    
    url <- paste0(url0,
                  cids,
                  "/property/MolecularWeight,MolecularFormula,RotatableBondCount/JSON?",
                  "Threshold=",Threshold)
    
  }else if(similarity == "2d&3d"){
    # 2d
    message("Step1 :searching 2d...")
    
    url0 <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/cid/"
    
    url <- paste0(url0,
                  cids,
                  "/property/MolecularWeight,MolecularFormula,RotatableBondCount/JSON?",
                  "Threshold=",Threshold)
    
    response <- try(httr::GET(url))
    
    if (response$status_code == 200) {
      
      json_data <- httr::content(response, as = "text",encoding = "UTF-8")  
      parsed_data <- jsonlite::fromJSON(json_data)
      similarity_compounds_2d <- parsed_data$PropertyTable$Properties %>% as.data.frame()
      
      similarity_compounds_2d$label <- "2d"
        
      message(paste0("Number of compouds:",nrow(similarity_compounds_2d)))
      
      
    }else{
      similarity_compounds_2d <- NULL
      message("Number of compouds: 0")
    }
    
    # 3d
    message("Step2 :searching 3d...")
    
    url0 <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_3d/cid/"
    
    url <- paste0(url0,
                  cids,
                  "/property/MolecularWeight,MolecularFormula,RotatableBondCount/JSON?",
                  "Threshold=",Threshold)
    
    response <- try(httr::GET(url))
    
    if (response$status_code == 200) {
      
      json_data <- httr::content(response, as = "text",encoding = "UTF-8")  
      parsed_data <- jsonlite::fromJSON(json_data)
      similarity_compounds_3d <- parsed_data$PropertyTable$Properties %>% as.data.frame()

      similarity_compounds_3d$label <- "3d"
      message(paste0("Number of compouds:",nrow(similarity_compounds_3d)))
      
    }else{
      
      similarity_compounds_3d <- NULL
      message("Number of compouds: 0")
      
    }
    t <- rbind(similarity_compounds_2d,similarity_compounds_3d)
    
    if(nrow(t)>0){
  
      t$query.id = cids
      
    }

    return(t)
    
  }
  
  response <- try(httr::GET(url))
  
  if (response$status_code == 200) {
    
    json_data <- httr::content(response, as = "text",encoding = "UTF-8")  
    parsed_data <- jsonlite::fromJSON(json_data)
    similarity_compounds <- parsed_data$PropertyTable$Properties %>% as.data.frame()
    
    similarity_compounds$query.id = cids
    message(paste0("Number of compouds:",nrow(similarity_compounds)))
    
  }else{
    
    similarity_compounds <- NULL
    message("Number of compouds: 0")
    
  }
  
  return(similarity_compounds)
}
