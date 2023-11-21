#' Running cross-sectional analysis for all timepoints
#'
#' @param data The dataframe with all covariates, targets and elements included in formula after having transformed the data
#' @param targets The dependant variable (phenotype) that is to be associated with
#' @param timepoints The timepoints which are to be analyzed and can be personlized to the timepoint variable
#' @param formula formula for providing the association without the dependant variable. Look at example for how to define the formula
#' @param fdr_correction Default to False but if FDR correction is required, the p-values are adjusted for the number of features
#'
#' @return Associations with the p-values, coefficients and variance explained for each timepoint and for all the features for each phenotype or target
#' @export
#'
#' @examples
run_cross_sectional_timepoints = function(data, targets,timepoints, formula, fdr_correction = FALSE){
  start_assoc = T
  for(timepoint in timepoints){

    temp_ = run_cross_sectional_association(data, targets,timepoint, formula, fdr_correction = FALSE)
    temp_$timepoint = timepoint

    if(start_assoc){

      assign('Association_Results', temp_)
      start_assoc = F

    }else{

      Association_Results = rbind(Association_Results,temp_)

    }
  }

  return(Association_Results)

}
