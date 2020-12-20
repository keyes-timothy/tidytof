# Write a modified FCS File (Code by Dmitry Tebaykin)
write_modified_FCS <- 
  function(
    tof_tibble, 
    fcs_name, 
    channel_descriptions = NULL, 
    reference_description = NULL, 
    old_description = NULL
  ) {
    if (is.matrix(fcs_exprs)) {
      # Don't write 'Leading_Push' to FCS files    
      fcs_exprs <- fcs_exprs[, which(colnames(fcs_exprs) != 'Leading_Push')]
      if(!is.matrix(fcs_exprs)){ 
        fcs_exprs <- t(as.matrix(fcs_exprs))
      }
      # Build metadata for FCS file
      pd <- c()  # 'params' phenoData
      dl <- list()  # 'description' list
      
      dl[['$DATATYPE']] <- 'F'
      
      if (!is.null(reference_description)) {
        if (!is.null(reference_description[['$DATE']])){
          dl[['$DATE']] <- reference_description[['$DATE']]
        }
        if (!is.null(reference_description[['$BTIM']])){
          dl[['$BTIM']] <- reference_description[['$BTIM']]
        }
        if (!is.null(reference_description[['$ETIM']])){
          dl[['$ETIM']] <- reference_description[['$ETIM']]
        }
        if (!is.null(reference_description[['$CYT']])){
          dl[['$CYT']] <- reference_description[['$CYT']]
        }
        if (!is.null(reference_description[['$CYTSN']])){
          dl[['$CYTSN']] <- reference_description[['$CYTSN']]
        }      
      }
      
      for (c in 1:ncol(fcs_exprs)) {
        c_name <- colnames(fcs_exprs)[c]
        c_desc <- colnames(fcs_exprs)[c]
        if (!is.null(channel_descriptions)){
          if (!is.na(channel_descriptions[c])){
            c_desc <- channel_descriptions[c]  
          }
        }
        
        c_min <- floor(min(0, min(fcs_exprs[, c])))  # Prevents flowCore from shifting range
        c_max <- ceiling(max(fcs_exprs[, c]))
        c_rng <- c_max - c_min + 1
        
        pl <- matrix(c(c_name, c_desc, c_rng, c_min, c_max), nrow = 1)
        colnames(pl) <- c('name', 'desc', 'range', 'minRange', 'maxRange')
        rownames(pl) <- paste('$P', c, sep = '') 
        pd <- rbind(pd, pl)
        
        if (!is.null(reference_description[[paste0('P', c, 'DISPLAY')]])){
          dl[[paste0('P', c, 'DISPLAY')]] <- reference_description[[paste0('P', c, 'DISPLAY')]]
        } 
        
        if (!is.null(reference_description[[paste0('$P', c, 'G')]])){
          dl[[paste0('$P', c, 'G')]] <- reference_description[[paste0('$P', c, 'G')]]
        } 
        
        if (!is.null(reference_description[[paste0('$P', c, 'R')]])){
          dl[[paste0('$P', c, 'R')]] <- reference_description[[paste0('$P', c, 'R')]]
        } else {
          dl[[paste('$P', c, 'R',sep = '')]] <- toString(c_rng); # Range
        }
        
        if (!is.null(reference_description[[paste0('$P', c, 'B')]])){
          dl[[paste0('$P', c, 'B')]] <- reference_description[[paste0('$P', c, 'B')]]
        } else {
          dl[[paste('$P', c, 'B', sep = '')]] <- '32';      # Number of bits
        }
        
        if (!is.null(reference_description[[paste0('$P', c, 'E')]])){
          dl[[paste0('$P', c, 'E')]] <- reference_description[[paste0('$P', c, 'E')]]
        } else { 
          dl[[paste('$P', c, 'E', sep = '')]] <- '0,0';      # Exponent
        }
        
        dl[[paste('$P', c, 'N', sep = '')]] <- c_name;	    # Name
        dl[[paste('$P', c, 'S', sep = '')]] <- c_desc;	    # Desc	
      }	
      
      if (!is.null(old_description)) {
        if (!is.null(old_description[['$CYT']])){
          dl[['$CYT']] <- old_description[['$CYT']]
        }		  
        if (!is.null(old_description[['$DATE']])){
          dl[['$DATE']] <- old_description[['$DATE']]
        }
        if (!is.null(old_description[['$BTIM']])){
          dl[['$BTIM']] <- old_description[['$BTIM']]
        }
        if (!is.null(old_description[['$ETIM']])){
          dl[['$ETIM']] <- old_description[['$ETIM']]
        }
      }
      
      fcs_exprs <- flowCore::flowFrame(fcs_exprs, 
                                       as(data.frame(pd), 'AnnotatedDataFrame'), 
                                       description = dl) 
    
    suppressWarnings(flowCore::write.FCS(fcs_exprs, fcs_name, what = 'numeric'))
    }
  }
