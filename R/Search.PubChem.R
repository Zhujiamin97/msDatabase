#'@title Search PubChem Database
#'
#'@author JiaMin ZHU
#'@description
#' PUG REST Tutorial (Documentation):
#' "https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest-tutorial#section=Input-Design-of-the-URL"
#'
#'@param compound_name 
#'@param cids IDs may also be specified in a comma-separated list : cids = c(100,101)
#'@param formula
#'@param smiles
#'@param InChIKey
#'@param Exact_Mass
#'@param mz_ppm
#'@param properties Name of the information to be scoped
#'           properties <- c("MolecularFormula", "MolecularWeight", 
#'                     "CanonicalSMILES", "IsomericSMILES", "InChI", "InChIKey", 
#'                     "IUPACName", "XLogP", "ExactMass", "MonoisotopicMass", 
#'                     "TPSA", "Complexity", "Charge", "HBondDonorCount", 
#'                     "HBondAcceptorCount", "RotatableBondCount", "HeavyAtomCount", 
#'                     "IsotopeAtomCount", "AtomStereoCount", "DefinedAtomStereoCount", 
#'                     "UndefinedAtomStereoCount", "BondStereoCount", "DefinedBondStereoCount", 
#'                     "UndefinedBondStereoCount", "CovalentUnitCount", 
#'                     "Volume3D", "XStericQuadrupole3D", "YStericQuadrupole3D", 
#'                     "ZStericQuadrupole3D", "FeatureCount3D", "FeatureAcceptorCount3D", 
#'                     "FeatureDonorCount3D", "FeatureAnionCount3D", "FeatureCationCount3D", 
#'                     "FeatureRingCount3D", "FeatureHydrophobeCount3D", 
#'                     "ConformerModelRMSD3D", "EffectiveRotorCount3D", 
#'                     "ConformerCount3D", "Fingerprint2D")
#'
#'@example 
#' use name seacrh : search.pubchem(cids = c(100,101)，
#'                                  properties = c("MolecularFormula", "MolecularWeight"))

search.pubchem <- function(compound_name = NULL,
                           CAS = "102-65-8",
                           cids = c(100,101),
                           formula = NULL,
                           smiles = NULL,
                           InChIKey = NULL,
                           Exact_Mass = NULL,
                           mz_ppm = 5,
                           properties = NULL,
                           author = "ZJM"){
  
  ## any type :fist step is get cids
  
  # search by Exact_Mass
  if(is.numeric(Exact_Mass)&is.numeric(mz_ppm)){
    
    mz_diff <- mz_ppm*10^(-6)*Exact_Mass
    mass_range = paste0(Exact_Mass-mz_diff,":",Exact_Mass+mz_diff)
    url = paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pccompound",
                 "&retmax=10000&term=",
                 mass_range,"%5Bexactmass%5D")
    web <- try(rvest::read_html(url))
    if(any(class(web) %in% "try-error")){
      message("please check your network!!")
      return(compoundlists = NA)
    }
    # If the connection time is exceeded ,return NA
    news <- web %>% rvest::html_nodes("body esearchresult idlist id")
    if(length(news)>0){
      cids = lapply(1:length(news),function(x){
        as.numeric(news[x] %>% rvest::html_text())
      }) %>% unlist()
      cids = c(paste0(cids, collapse = ","))
    }else{
      return(NA)
    }
  }else if(!is.null(InChIKey)){
    # search by InChIKey
    url_inchikey <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/",
                           InChIKey,"/cids/JSON")
    response <- try(httr::GET(url_inchikey))
    
    if (response$status_code == 200) {
      
      json_data <- content(response, as = "text",encoding = "UTF-8")  
      parsed_data <- jsonlite::fromJSON(json_data)
      cids = parsed_data$IdentifierList$CID
      cids = c(paste0(cids, collapse = ","))
      
    }else{  
      return(NA)
    }
  }else if(!is.null(smiles)){
    # search by InChIKey
    url_smiles <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/",
                           smiles,"/cids/JSON")
    response <- try(httr::GET(url_smiles))
    
    if (response$status_code == 200) {
      
      json_data <- content(response, as = "text",encoding = "UTF-8")  
      parsed_data <- jsonlite::fromJSON(json_data)
      cids = parsed_data$IdentifierList$CID
      cids = c(paste0(cids, collapse = ","))
      
    }else{  
      return(NA)
    }
  }else if(!is.null(compound_name)){
    # search by name
    url_name <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
                       compound_name,"/cids/JSON")
    response <- try(httr::GET(url_name))
    
    if (response$status_code == 200) {
      
      json_data <- content(response, as = "text",encoding = "UTF-8")  
      parsed_data <- jsonlite::fromJSON(json_data)
      cids = parsed_data$IdentifierList$CID
      cids = c(paste0(cids, collapse = ","))
      
    }else{  
      return(NA)
    }
  }else if(!is.null(cids)){
    # search by cids
    cids = c(paste0(cids, collapse = ","))
  }
  
  #--------------After getting the cids ,then start getting the compound information------------------------##
  
  prolog <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
  
  input <- "/compound/cid/"
  
  output <- "/JSON"
  
  if(is.null(properties)){
    properties <- c("MolecularFormula", "MolecularWeight", 
                    "CanonicalSMILES", "IsomericSMILES", "InChI", "InChIKey", 
                    "IUPACName", "XLogP", "ExactMass", "MonoisotopicMass", 
                    "TPSA", "Complexity", "Charge", "HBondDonorCount", 
                    "HBondAcceptorCount", "RotatableBondCount", "HeavyAtomCount", 
                    "IsotopeAtomCount", "AtomStereoCount", "DefinedAtomStereoCount", 
                    "UndefinedAtomStereoCount", "BondStereoCount", "DefinedBondStereoCount", 
                    "UndefinedBondStereoCount", "CovalentUnitCount", 
                    "Volume3D", "XStericQuadrupole3D", "YStericQuadrupole3D", 
                    "ZStericQuadrupole3D", "FeatureCount3D", "FeatureAcceptorCount3D", 
                    "FeatureDonorCount3D", "FeatureAnionCount3D", "FeatureCationCount3D", 
                    "FeatureRingCount3D", "FeatureHydrophobeCount3D", 
                    "ConformerModelRMSD3D", "EffectiveRotorCount3D", 
                    "ConformerCount3D", "Fingerprint2D")
  }
  
  properties <- paste0(properties, collapse = ",")
  
  url <- paste0(prolog,input,cids,"/property/",properties,output)
  
  response <- try(httr::GET(url))
  
  if (response$status_code == 200) {
    
    json_data <- content(response, as = "text",encoding = "UTF-8")  
    parsed_data <- jsonlite::fromJSON(json_data)
    compounds_data <- parsed_data$PropertyTable$Properties
    
  }else{  
      return(NA)
  }
  return(compounds_data)
}
