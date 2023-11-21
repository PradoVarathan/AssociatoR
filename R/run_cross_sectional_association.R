#' Cross Sectional Association
#'
#' @param data The dataframe with all covariates, targets and elements included in formula after having transformed the data
#' @param targets The dependant variable (phenotype) that is to be associated with
#' @param timepoint The timepoint which is to be analyzed and can be personlized to the timepoint variable
#' @param formula formula for providing the association without the dependant variable. Look at example for how to define the formula
#' @param fdr_correction Default to False but if FDR correction is required, the p-values are adjusted for the number of features
#'
#' @return Associations with the p-values, coefficients and variance explained for all the features for each phenotype or target
#' @export
#'
#' @examples
run_cross_sectional_association = function(data, targets,timepoint, formula, fdr_correction = FALSE){

  #Filtering data
  data = data %>% filter(timepoint == timepoint)
  columns_needed = strsplit(formula,' ')[[1]][sapply(strsplit(formula,' ')[[1]],function(str){!grepl("[^A-Za-z0-9_ ]", str)})]
  data = data[,c(columns_needed,targets)]

  ## Fitting line
  start_prep = T
  for(target in targets){

    target_formula = paste0(target,' ~ ',formula)
    lin_fit = lm(as.formula(target_formula),
                 data = data)


    temp_p_vals = summary(lin_fit)$coefficients[,4]
    temp_coef_vals = summary(lin_fit)$coefficients[,1]
    temp_main_adj_r2 = summary(lin_fit)$adj.r.squared
    factors = columns_needed[-1]
    temp_df = data.frame(matrix(ncol = length(factors), nrow = 0))
    colnames(temp_df) = sapply(factors, function(x){return(paste0(x,'_Var_Exp'))})

    # Obtaining variance explained
    for(i in 1:length(colnames(temp_df))){
      factors_temp = factors[-i]
      target_line_temp = paste0(target,' ~ ',paste(factors_temp,collapse = ' + '))
      lin_fit_temp = lm(as.formula(target_line_temp),
                        data = data)
      temp_part_adj_r2 = summary(lin_fit_temp)$adj.r.squared
      temp_df[1,i] = temp_main_adj_r2 - temp_part_adj_r2
    }

    if(start_prep){

      Association_dataframe = data.frame(matrix(ncol = (length(factors)*3)+1, nrow = 0))
      Association_colnames = c(sapply(factors, function(x){return(paste0(x,'_p_val'))}),
                               sapply(factors, function(x){return(paste0(x,'_Coef'))}),
                               sapply(factors, function(x){return(paste0(x,'_Var_Exp'))}))
      colnames(Association_dataframe) = Association_colnames

      for(f in factors){

        Association_dataframe[1,paste0(f,'_p_val')] = as.numeric(min(temp_p_vals[grep(f,names(temp_p_vals))]))
        Association_dataframe[1,paste0(f,'_Coef')] = as.numeric(min(temp_coef_vals[grep(f,names(temp_p_vals))]))
        Association_dataframe[1,paste0(f,'_Var_Exp')] = as.numeric(temp_df[paste0(f,'_Var_Exp')])

      }

      rownames(Association_dataframe) = c(target)
      start_prep = F

    }else{

      temp_Association_dataframe = data.frame(matrix(ncol = (length(factors)*3)+1, nrow = 0))
      colnames(temp_Association_dataframe) = Association_colnames

      for(f in factors){

        temp_Association_dataframe[1,paste0(f,'_p_val')] = as.numeric(min(temp_p_vals[grep(f,names(temp_p_vals))]))
        temp_Association_dataframe[1,paste0(f,'_Coef')] = as.numeric(min(temp_coef_vals[grep(f,names(temp_p_vals))]))
        temp_Association_dataframe[1,paste0(f,'_Var_Exp')] = as.numeric(temp_df[paste0(f,'_Var_Exp')])

      }
      rownames(temp_Association_dataframe) = c(target)
      Association_dataframe = rbind(Association_dataframe,temp_Association_dataframe)

    }

  }
  if(fdr_correction){

    for(f in factors){

      Association_dataframe[,paste0(f,'_p_val')] = p.adjust(Association_dataframe[,paste0(f,'_p_val')], method = 'fdr')

    }

  }

  return(Association_dataframe)
}
