#' Run the P-model
#'
#' Runs the P-model and loads output in once.
#'
#' @param drivers A nested data frame with one row for each site and columns
#'  named according to the arguments of function `runread_pmodel_f_bysite`
#' @param par A named list of model parameters
#' @param makecheck A logical specifying whether checks are performed to verify
#'  forcings.
#' @param parallel A logical specifying whether simulations are to be
#'  parallelised (sending data from a certain number of sites to each core).
#'  Defaults to \code{FALSE}.
#' @param ncores An integer specifying the number of cores used for parallel
#'  computing (default = 2).
#'
#' @return A data frame (tibble) with one row for each site and outputs stored
#'  in the nested column \code{data}
#' @export
#'
#' @examples
#' \dontrun{
#'   mod <- runread_pmodel_f( 
#'    drivers,
#'    par,
#'    makecheck = TRUE,
#'    parallel = FALSE,
#'    ncores = 2 )
#' }

runread_pmodel_f <- function(
  drivers,
  par,
  makecheck = TRUE,
  parallel = FALSE,
  ncores = 1){
  
  if (parallel){

    cl <- multidplyr::new_cluster(n = ncores) %>%
      multidplyr::cluster_assign(par = par) %>%
      multidplyr::cluster_assign(makecheck = FALSE) %>%
      multidplyr::cluster_library(c("dplyr", "purrr", "rlang", "rsofun"))
    
    # distribute to to cores, making sure all data from
    # a specific site is sent to the same core
    df_out <- drivers %>%
      dplyr::group_by(id = row_number()) %>%
      tidyr::nest(
        input = c(
          sitename,
          params_siml,
          site_info,
          forcing,
          params_soil)
        ) %>%
      multidplyr::partition(cl) %>% 
      dplyr::mutate(data = purrr::map(input, 
       ~run_pmodel_f_bysite(
         sitename       = .x$sitename[[1]], 
         params_siml    = .x$params_siml[[1]], 
         site_info       = .x$site_info[[1]], 
         forcing        = .x$forcing[[1]], 
         params_soil = .x$params_soil[[1]], 
         par    = par, 
         makecheck      = makecheck )
      ))
    
     # collect the cluster data
     data <- df_out %>%
      dplyr::collect() %>%
      dplyr::ungroup() %>%
      dplyr::select( data )
     
     # meta-data
     meta_data <- df_out %>%
       dplyr::collect() %>%
       dplyr::ungroup() %>%
       dplyr::select( input ) %>%
       tidyr::unnest( cols = c( input )) %>%
       dplyr::select(sitename, site_info)
     
     # combine both data and meta-data
     # this implicitly assumes that the order
     # between the two functions above does
     # not alter! There is no way of checking
     # in the current setup
     df_out <- bind_cols(meta_data, data)
    
  } else {
    
    df_out <- drivers %>% 
      dplyr::mutate(
        data = purrr::pmap(.,
                           run_pmodel_f_bysite,
                           par = par,
                           makecheck = makecheck
        )) %>% 
      dplyr::select(sitename, site_info, data)
  }
  
  return(df_out)
}

#' Run LM3-PPA
#'
#' Runs the LM3-PPA model and loads output in once.
#'
#' @param drivers A nested data frame with one row for each site and columns
#'  named according to the  arguments of function `runread_pmodel_f_bysite()`
#' @param makecheck A logical specifying whether checks are performed to verify
#'  forcings.
#' @param parallel A logical specifying whether simulations are to be 
#'  parallelised (sending data from a certain number of sites to each core). 
#'  Defaults to \code{FALSE}.
#' @param ncores An integer specifying the number of cores used for parallel 
#' computing. Defaults to 2.
#'
#' @return A tibble with one row for each site and outputs stored 
#' in the nested column \code{data}
#' @export
#'
#' @examples 
#' \dontrun{
#'  mod <- runread_pmodel_f( df_drivers,
#'   params_modl, makecheck = TRUE, parallel = FALSE, ncores = 2 )
#' }

runread_lm3ppa_f <- function(
  drivers,
  makecheck = TRUE,
  parallel = FALSE,
  ncores = 2
  ){
  
  if (parallel){
    
    cl <- multidplyr::new_cluster(ncores) %>% 
      multidplyr::cluster_assign(makecheck = FALSE) %>% 
      multidplyr::cluster_library(c("dplyr", "purrr", "rlang", "rsofun"))
    
    ## distribute to to cores, making sure all data from a specific site is sent to the same core
    df_out <- drivers %>%
      dplyr::group_by( id = row_number() ) %>%
      tidyr::nest(input = c(sitename,
                            params_siml,
                            site_info,
                            forcing,
                            params_tile,
                            params_species,
                            params_soil,
                            init_cohort,
                            init_soil)) %>%
      multidplyr::partition(cl) %>% 
      dplyr::mutate(data = purrr::map( input, 
         ~run_lm3ppa_f_bysite(
           sitename       = .x$sitename[[1]], 
           params_siml    = .x$params_siml[[1]], 
           site_info       = .x$site_info[[1]], 
           forcing        = .x$forcing[[1]], 
           params_tile    = .x$params_tile[[1]], 
           params_species = .x$params_species[[1]], 
           params_soil    = .x$params_soil[[1]], 
           init_cohort    = .x$init_cohort[[1]], 
           init_soil      = .x$init_soil[[1]], 
           makecheck      = makecheck )
         
      )) %>% 
      dplyr::collect() %>%
      dplyr::ungroup() %>%
      dplyr::select( data )  %>% 
      tidyr::unnest( cols = c( data ))
    
  } else {
    
    df_out <- drivers %>% 
      dplyr::mutate(data = purrr::pmap(
        .,
        run_lm3ppa_f_bysite,
        makecheck = makecheck
      )) %>% 
      dplyr::select(sitename, data) 
    
  }
  
  return(df_out)
}