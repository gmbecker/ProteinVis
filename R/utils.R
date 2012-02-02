combineTMDomains = function(tmdf, poscols = c("start", "end"))
  {
    tmmat = matrix(NA, ncol = 2, nrow = nrow(tmdf))
    curnumrow = 1    
    for (i in 1:nrow(tmdf))
      {
        inserted = 0
        for(j in 1:curnumrow)
          {
            if(!is.na(tmmat[j, 1]))
              {
                if (any(tmdf[i ,poscols] >= tmmat[j, 1] & tmdf[i, poscols ] <= tmmat[j, 2]) )
                  {
                    tmmat[j , 1 ] = min( tmdf[ i , poscols[ 1 ] ] , tmmat[ j , 1 ] )
                    tmmat[j , 2 ] = max( tmdf[ i , poscols[ 2 ] ] , tmmat[ j , 2 ] )
                    inserted = 1
                    break()
                  }
              }
          }
        if (!inserted)
          {
            
            tmmat[curnumrow,] = as.numeric(tmdf[i,poscols])
            curnumrow = curnumrow + 1
          }
      }
    tmmat = tmmat[!is.na(tmmat[,1]), , drop=FALSE]
    tmmat
  }

calcPlotHeight = function(baseheight, type = "metaCount", pfam, categories)
  {
    if(!is.null(pfam) & nrow(pfam))
      {
    myPFIRanges = IRanges(start = pfam$start, end = pfam$end)
    bins = disjointBins(myPFIRanges)
  } else {
    bins = 1
  }
    denom = .20 + .25 + 1 #assume 1 pfam row and 5 categories
    if (type == "struct") 
      {
                                        #.20 + .25 + .25 #the .25 is for the hydro, which isn't in the above plot.
        categories = NULL
      }

    

    baseheight * (1 + .2 /denom * ( max(bins) - 1) + .2 / denom* (length(unique(categories)) - 5 + 2)) #2 fake categories

  }

spoofLevelsInDF = function(df, colname, newlevs, before = TRUE, force.factor = TRUE)
  {
    oldnrow = nrow(df)
    if(is.factor(df[[colname]]))
       oldlevs = levels(df[[colname]])
    else
      oldlevs = unique(df[[colname]]) #XXX this doesn't seem to give us the right order!!!
    df[[colname]] = as.character(df[[colname]])
    
     if(before)
      alllevs = unique(c(newlevs, oldlevs))
    else
      alllevs = unique(c(oldlevs, newlevs))
    
    for(i  in seq(along = newlevs))
      df[oldnrow + i, colname] = newlevs[i]

    if(force.factor)
      df[[colname]] = factor(df[[colname]], levels = alllevs)
    return(df)
  }
