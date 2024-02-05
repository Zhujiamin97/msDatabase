#'@title Search Norman Database
#'
#'@author JiaMin ZHU
#'@description
#'
#'@param compound_name 
#'@param formula
#'@param InChIKey
#'@param Exact_Mass
#'@param mz_error
#'@param NormanID
#'@param CAS
#'@param ms2_matrix for example use: ms2_matrix_demo()
#'
#'@example 
#' use name seacrh : search.norman(compound_name = "Sulfaclozine",
#'                                 formula = NULL)

search.norman <- function(compound_name = "Sulfaclozine",
                          NormanID = "NS00000001",
                          CAS = "102-65-8",
                          CID = "66890",
                          formula = NULL,
                          InChIKey = NULL,
                          Exact_Mass = NULL,
                          mz_error = 0.3,
                          author = "ZJM"){
  # load norman csv
  try(data("Norman_SusDat"))
  if(nrow(Norman_SusDat)==0 |is.null(Norman_SusDat)){
    stop("No Norman database file in your dir!")
  }
  
  # first
  if(!is.null(CAS)){
    dat = Norman_SusDat %>% filter(CAS_RN_Dashboard == CAS)
    return(as.list(dat))
  }else if(!is.null(NormanID)){
    dat = Norman_SusDat %>% filter(Norman_SusDat_ID == NormanID)
    return(as.list(dat))
  }else if(!is.null(CID)){
    dat = Norman_SusDat %>% filter(PubChem_CID == CID)
    return(as.list(dat))
  }else if(!is.null(InChIKey)){
    dat = Norman_SusDat %>% filter(StdInChIKey == InChIKey)
    return(as.list(dat))
  }else if(!is.null(formula)){
    dat = Norman_SusDat %>% filter(Molecular_Formula == formula)
    return(as.list(dat))
  }else if(!is.null(Exact_Mass)){
    if(is.numeric(Exact_Mass)&is.numeric(mz_error)){
      dat = Norman_SusDat[which(as.numeric(Norman_SusDat$Monoiso_Mass) 
                                %between% c(Exact_Mass-mz_error,Exact_Mass+mz_error)),]
      return(as.list(dat))
    }else{
      stop("Exact_Mass & mz_error must be number")
    }
  }
  
}


