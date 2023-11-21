#' Making Data Counts
#'
#' @param data Dataset as a dataframe with the timepoint vairable as year column
#' @param target_cols The phenotype or column names for making the count
#' @param thresh Number of years as threshold, Has been default to 5
#' @param title The title of the phenotype for the graph
#'
#' @return Returns a barplot with counts for non-NA values inside target columns
#' @export
#'
#' @examples
#' making_data_counts(data, target_cols, 'Phenotype 1')
making_data_counts = function(data,target_cols, thresh = 5,title){

  start = TRUE
  for(item in target_cols){

    data_t = na.omit(data[,c(item,'year')])
    data_t = data_t %>% filter(year <= thresh)
    data_t = data_t[which(data_t$year%%1 == 0),]
    if(start){
      df = as.data.frame(table(data_t$year))
      df = df[which(df$Freq != 0),]
      df$Phenotype = item
      start = F
    }else{
      t = as.data.frame(table(data_t$year))
      t = t[which(t$Freq != 0),]
      t$Phenotype = item
      df = rbind(df,t)
    }

  }

  my_bar = ggplot(df,aes(x = Var1, y =Freq, fill = Phenotype)) +
    geom_bar(stat = "identity", position = "dodge") +ylab('Subjects') + ggtitle(title) +
    xlab(paste0('Time points upto ',as.character(thresh),' years'))

  return(my_bar)
}
