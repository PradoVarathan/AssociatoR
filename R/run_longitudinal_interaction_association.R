#' Running linear mixed models for all timepoints
#'
#' @param data The dataframe with all covariates, targets and elements included in formula after having transformed the data
#' @param targets The dependant variable (phenotype) that is to be associated with
#' @param threshold_years Number of years as threshold, Has been default to 5
#' @param formula formula for providing the association without the dependant variable. Look at example for how to define the formula
#' @param plot_effect Default to True. Provides effect plots for the linear mixed models fit
#' @param plot_residuals Default to False Provides residual plots for the linear mixed models fit
#'
#' @return list of one dataframe that provides the p-value of the phenotype or target, and if plot_effect is TRUE, returns effect plot for each phenotype
#'     and if plot_residuals is TRUE, return residula plots for the fit
#' @export
#'
#' @examples
run_longitudinal_interaction_association = function(data, targets, predictor_main = "delta_age" ,type = "interaction",threshold_years = 5, formula, plot_effect = TRUE, plot_residuals = FALSE){

  #Filtering data
  data = data %>% filter(year <= threshold_years)
  columns_needed = unique(strsplit(formula,' ')[[1]][sapply(strsplit(formula,' ')[[1]],function(str){!grepl("[^A-Za-z0-9_ ]", str)})])
  columns_needed =  columns_needed[suppressWarnings(is.na(as.numeric(columns_needed)))]
  data = data[,c(columns_needed,targets)]
  target_factor = strsplit(formula,' ')[[1]][c(which(strsplit(formula,' ')[[1]] == '*')-1,which(strsplit(formula,' ')[[1]] == '*')+1)]
  data = na.omit(data)
  ## Fitting line
  start_prep = T
  effect_plots = list()
  residual_plots = list()
  out_data_frame  = list()
  beta_data_frame = list()
  t_data_frame = list()
  std_err_data_frame = list()
  for(target in targets){

    target_formula = paste(target,formula, sep = ' ~ ')

    og = lmer(as.formula(target_formula),data, REML = FALSE)
    if(type == "interaction"){
      og.n = lmer(as.formula(chartr("*","+",target_formula)),data, REML = FALSE)
      }
    else if(type == "none"){
      og.n = lmer(as.formula(gsub(paste0("[+] ",predictor_main," "),"",formula)),data, REML = FALSE)
      }
    temp_og = summary(og)
    k = anova(og,og.n)
    p_val = k$`Pr(>Chisq)`[2]
    out_data_frame[target] = p_val
    beta_data_frame[target] = temp_og$coefficients[nrow(temp_og$coefficients),1]
    t_data_frame[target] = temp_og$coefficients[nrow(temp_og$coefficients),3]
    std_err_data_frame[target] =  temp_og$coefficients[nrow(temp_og$coefficients),2]
    if(plot_effect){

      Quartiles <- quantile(data[,colnames(data) == predictor_main], probs = c(0, 0.25, 0.5, 0.75, 1))
      data$Quartiles <- cut(data[,colnames(data) == predictor_main], breaks = Quartiles, labels = c("Q1", "Q2", "Q3", "Q4"), include.lowest = TRUE)

      pred1=predict(og,re.form=NA,se.fit=TRUE,nsim=500)
      # save prediction results
      temp_model_data = data
      temp_model_data$boot_fit = pred1$fit

      # Only plot Q1 and Q4
      plot_Q1Q4 = temp_model_data %>% filter(Quartiles=="Q1"|Quartiles=="Q4")

      # Can change colors, linewidth, x labels, y labels or add title as needed
      p = ggplot(plot_Q1Q4,aes(x=year,y=boot_fit,col=Quartiles))+geom_smooth(aes(year,boot_fit),linewidth=1.5,method=lm,se=TRUE)+
        ylab(paste0(target,"_pred"))+xlab("Year")+ggtitle(paste0(" p-val: ",if(p_val > 0.00005){as.character(round(p_val,5))}else{p_val}))+
                                                            scale_color_manual(values=c("blue","red"))

      effect_plots[[target]] = p
    }
    if(plot_residuals){

      df = data.frame("residulas" = resid(og))
      p <- ggplot(df, aes(x=residulas)) +
        geom_histogram() + xlab('Residual Value') + ggtitle(target)
      residual_plots[[target]] = p

    }



  }

  out_list = list('p_values' = out_data_frame)
  if(plot_effect){
      if(plot_residuals){
    out_list = list('p_values' = out_data_frame,'effect_plots' =  effect_plots,'residual_plots' = residual_plots,'beta' = beta_data_frame,'t_values' = t_data_frame,'std_err_values' = std_err_data_frame)
  }else{
    out_list = list('p_values' = out_data_frame,'effect_plots' =  effect_plots,'beta' = beta_data_frame,'t_values' = t_data_frame,'std_err_values' = std_err_data_frame)
    }
  }else if(plot_residuals){
    out_list = list('p_values' = out_data_frame,'residual_plots' =  residual_plots,'beta' = beta_data_frame,'t_values' = t_data_frame,'std_err_values' = std_err_data_frame)
    }
  return(out_list)

}
