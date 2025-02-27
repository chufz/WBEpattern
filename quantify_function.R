################ S T E P  3 - QUANTIFICATION ###################################

## R E M A R K S ## 

# this script needs to be open in order to run the "quantify_" scripts 


quantify <- function(sample_frame, IS_sample_frame, cal_frame, levels, IS_cal_frame, IS_assignment, enr_factor=1, compounds="", output_pdf="./calib_curves.pdf"){
    if(ncol(sample_frame)!=ncol(cal_frame))message("cal_frame and sample_frame should have the same number of variables (compounds in columns)")
    if(nrow(cal_frame)!=nrow(IS_cal_frame))message("cal_frame and IS_cal_frame should have the same number of observables (samples in rows)")
    if(nrow(sample_frame)!=nrow(IS_sample_frame))message("sample_frame and IS_sample_frame should have the same number of observables (samples in rows)")
    if(ncol(IS_sample_frame)!=ncol(IS_cal_frame))message("IS_sample_frame and IS_cal_frame should have the same number of variables (compounds in columns")
    if(length(IS_assignment)!=ncol(sample_frame))message("IS_assignment should have the same number of compounds as in the data frames")
    if(length(levels)!=nrow(cal_frame))message("levels must be in the same number than the observables in cal_frame")
    #calculate response factor through models and plot cal curves
    pdf(file=output_pdf)
    rf <- list()
    for(i in 1:ncol(cal_frame)){
        # divide cals by IS 
        frame <- data.frame("y"=levels,"x"=cal_frame[,i] / IS_cal_frame[,IS_assignment[i]])
        frame <- frame[!is.na(frame$x),]
        if(length(frame$y[!is.na(frame$y)]) < 2){
            rf[[i]] = NA
            next
        }
        rf[[i]] <- lm(formula= y~0+x,weights = 1/x,data=frame)
        if(length(compounds) < 2){
            plot(frame$y~frame$x,xlab="Response Ratio",ylab="Conc. [ng/mL]",main=paste0("R2 = ", format(summary(rf[[i]])$adj.r.squared, digits = 2)))
        }else{
            plot(frame$y~frame$x,xlab="Response Ratio",ylab="Conc. [ng/mL]",main=paste0(compounds[i], " R2 = ", format(summary(rf[[i]])$adj.r.squared, digits = 2)))
        }
        
        abline(rf[[i]])
    }
    
    
    dev.off()
    
       # divide sample by IS
    r_sample_frame <- sample_frame
    for(i in 1:ncol(sample_frame)){
        r_sample_frame[,i] <- sample_frame[,i]/IS_sample_frame[,IS_assignment[i]]
    }
    
    # calculate concentrations
    conc_table <- r_sample_frame
    for(i in 1:ncol(sample_frame)){
        conc_table[,i] <- predict.lm(object=rf[[i]], newdata = data.frame(x=r_sample_frame[,i]), type="response")
    }
     
    # use enrichment factor
    conc_table <- conc_table * enr_factor
        
    return(conc_table)
}



linear <- function(sample_frame, IS_sample_frame, cal_frame, levels, IS_cal_frame, IS_assignment, enr_factor = 1, compounds = "") {
  # Checking dimensions
  if(ncol(sample_frame) != ncol(cal_frame)) message("cal_frame and sample_frame should have the same number of variables (compounds in columns)")
  if(nrow(cal_frame) != nrow(IS_cal_frame)) message("cal_frame and IS_cal_frame should have the same number of observables (samples in rows)")
  if(nrow(sample_frame) != nrow(IS_sample_frame)) message("sample_frame and IS_sample_frame should have the same number of observables (samples in rows)")
  if(ncol(IS_sample_frame) != ncol(IS_cal_frame)) message("IS_sample_frame and IS_cal_frame should have the same number of variables (compounds in columns)")
  if(length(IS_assignment) != ncol(sample_frame)) message("IS_assignment should have the same number of compounds as in the data frames")
  if(length(levels) != nrow(cal_frame)) message("levels must be in the same number than the observables in cal_frame")
  
  # Creating an empty list to store models and summaries
  rf2 <- list()
  adj_rsq <- numeric(ncol(cal_frame))
  
  for(i in 1:ncol(cal_frame)) {
    frame <- data.frame("y" = levels, "x" = cal_frame[, i] / IS_cal_frame[, IS_assignment[i]])
    frame <- frame[!is.na(frame$x),]
    
    if(length(frame$y[!is.na(frame$y)]) < 2) {
      # If there are insufficient data points, store NA
      summary_list[[i]] <- list(model = NA, summary = NA)
      next
    }
    
    # Fitting linear model
    rf2[[i]] <- lm(formula = y ~ 0 + x, weights = 1/x, data = frame)
    
    # Storing summary in summary_list
    adj_rsq[i] <- summary(rf2[[i]])$adj.r.squared
  }
  
  # Returning the list of models and summaries
  return(adj_rsq)
}



findBestIS <- function(sample_frame, IS_sample_frame, RT_comp, RT_IS){
    IS_assignment <- vector(length=ncol(sample_frame))
    for(i in 1:length(IS_assignment)){
        if(nrow(sample_frame)!=nrow(IS_sample_frame))message("Sample_frame and IS_sample_frame need to have the same number of observables (samples in rows)")
        if(ncol(sample_frame)!=length(RT_comp))message("Vector of RT of Compounds need to be the same length as compounds (columns) in sample_frame")
        if(ncol(IS_sample_frame)!=length(RT_IS))message("Vector of RT of IS need to be the same length as IS (columns) in IS_sample_frame")
        # chose IS that are found in all samples with a detection
        compound_detects <- IS_sample_frame[which(!is.na(sample_frame[,i])),]
        x <- apply(compound_detects, 2, function(x)sum(is.na(x)))
        if(sum(x==0)<1){
            message("For compound ", colnames(sample_frame)[i], " no IS showed sufficient detection rate, NA will be introduced.")
            # chose one with highest DR
            IS_assignment[i] <- which.min(x)
        }
        if(sum(x==0)==1){
            message("For compound ", colnames(sample_frame)[i], " only one IS showed sufficient detection rate.")
            IS_assignment[i] <- which(x==0)
        }
        if(sum(x==0)>1){
            message("For compound ", colnames(sample_frame)[i], " several IS possible for assignment, chosing closest RT.")
            # get rt deviations
            rt_dev <- abs(RT_IS - RT_comp[i])
            IS_assignment[i] <- which(x==0)[which.min(rt_dev[which(x==0)])]
        }
    }
    return(IS_assignment)
}

blankCorrection <- function(sample_frame, blank_frame, factor=5){
    if(ncol(sample_frame)!=ncol(blank_frame))message("Sample frame and blank frame need to have same number of compounds (columns)")
    mean_blanks <- apply(blank_frame,2,function(x)mean(x,na.rm=TRUE) )
    # remove Signal intensities below MDL
    for(i in 1:ncol(sample_frame)){
        if(!is.na(mean_blanks[i])){
            x <- which(sample_frame[,i] < factor*mean_blanks[i])
            sample_frame[x,i] <- NA
        }
    }
    return(sample_frame)
}

RSD <- function(x){(sd(x, na.rm=T)/mean(x, na.rm=T))*100}
DR <- function(x){(ncol(x) - apply(x, 1, function(x)sum(is.na(x))))/ ncol(x) *100}

