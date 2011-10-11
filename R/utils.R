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

calcPlotHeight = function(baseheight, type = "metaCount", pfam, events, catname)
  {
    myPFIRanges = IRanges(start = pfam$start, end = pfam$end)
    bins = disjointBins(myPFIRanges)
    denom = if(type == "metaCount")
      .20 + .3 + .25 + 1 #assume 1 pfam row and 5 categories
    else
      .20 + .3 + .25

    
    baseheight * (1 + .2 /denom * ( max(bins) - 1) + .2 / denom* (length(unique(events[[catname]])) - 5))
  }
