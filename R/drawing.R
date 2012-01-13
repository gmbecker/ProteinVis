

axis.protstruct <-
    function(side, ...) 
{
    ylim <- current.panel.limits()$ylim
    xlim <- current.panel.limits()$xlim
    switch(side,
           left = {
               prettyY <- pretty(ylim)
               panel.axis(side = side, outside = TRUE,
                          at = prettyY, labels = prettyY)
           },
           bottom = {
             prettyX <- pretty(xlim)
             panel.axis(side = side, outside = TRUE,
                        at = prettyX, labels = prettyX)
           }
           
           )
  }

panel.protstruct <- function(x,y, subscripts, tmposcols, tm, sig, vertGuides, ... )
  {

    ylim = current.panel.limits()$ylim
    xlim = current.panel.limits()$xlim

    drawTM(tm, ylim)
    drawSigP(sig, ylim)
    #comes after tm/sigp stuff to deal with alpha issues
    grid.text("Hydrophobicity", unit(2, "mm"), unit(1, "npc") - unit(2, "mm"), just = c("left", "top"), gp = gpar(fontface="bold", cex=.9 ))

    pos = unit(2, "mm") + unit(1, "strwidth", "Hydrophobicity") + unit(2, "char")
    
    grid.text("Signal Peptide Prediction", pos, unit(1, "npc") - unit(2, "mm"), just = c("left", "top"), gp = gpar( cex=.9 ))
    pos = pos + unit(.8, "strwidth", "Signal Peptide Prediction") + unit(1, "char")
    grid.rect(pos, unit(1, "npc") - unit(2, "mm"), unit(1, "char"), unit(.9, "char"),gp = gpar(fill = rgb(180, 255, 180, max=255), alpha = .5, col = "white"), just = c("left", "top"))
    pos = pos + unit(1, "char") + unit(4, "char")
    
    grid.text("Transmembrane Domain Prediction", pos, unit(1, "npc") - unit(2, "mm"), just = c("left", "top"), gp = gpar( cex=.9 ))
    pos = pos + unit(.8, "strwidth", "Transmembrane Domain Prediction") 
    grid.rect(pos, unit(1, "npc") - unit(2, "mm"), unit(1, "char"), unit(.9, "char"),gp = gpar(fill = "grey75", alpha = .5, col = "white"), just = c("left", "top"))
    
              
    #This may get rid of the problem with hydrophobicity not plotting correctly. Not sure, I need to get the data from cory.
    if(length(x))
      panel.xyplot(x,y,...)
    TRUE
  }

drawTM = function(tm, ylim)
  {

    if(!is.null(tm) && nrow(tm) > 0 )
      {
        tmmat = combineTMDomains(tm, poscols = tmposcols)
        apply(tmmat, 1, function(x, ylim)
              {
                grid.polygon( c( x , rev( x ) ), rep( ylim , times = c(2, 2) ) , gp = gpar(fill = "grey75", alpha = .5, col = "white"), default.units = "native", draw = TRUE)
              }, ylim = ylim)
      }
    TRUE
  }

drawSigP = function(sig, ylim)
  {
    if(!is.null(sig) && nrow(sig) > 0)
      {
        apply(sig, 1, function(s, ylim)
              {
                sigst = as.numeric(s["start"])
                sigend = as.numeric(s["end"])
                sigx = c(sigst, sigend)
                grid.polygon(c(sigx, rev(sigx)), rep(ylim, times = c(2,2 )), gp = gpar(fill = rgb(180, 255, 180, max = 255), alpha = .5, col = "white"), default.units = "native", draw = TRUE)
                 }, ylim = ylim)
      }
   TRUE
  }

panel.psipred = function(x, y, subscripts, cutoff, strand, vertGuides, ...)
  {
    drawVertGuides(vertGuides, "grey80")
      
    drawPsipred(data.frame(start = x,helix = y, strand = strand), cutoff, current.panel.limits()$xlim)
  }

drawPsipred = function(dat, cutoff, xlim)
  {
    grid.text("Secondary Structure", unit(2, "mm"),unit(1, "npc") -  unit(2, "mm"), just = c("left", "top"), gp = gpar(fontface="bold", cex=.8 ))

    grid.text("strand", unit(2, "mm") + unit(1, "strwidth", data="Secondary Structure") + unit(6, "mm"),unit(1, "npc") -  unit(2, "mm"), just = c("left", "top"), gp = gpar(cex = .9))

    
    arrowpos = convertX( unit(2, "mm") + unit(1, "strwidth", data="Secondary Structure") + unit(6, "mm")+ unit(1, "strwidth", data="strand  "), "native", TRUE)
    ceny = convertY(unit(1, "npc") - unit(2, "mm") - unit(.5, "strheight", data="S"), "npc", TRUE)

    drawArrow( arrowpos, arrowpos + 15, height=.6*convertY(unit(1, "strheight", data="S"), "npc", TRUE), center.y = ceny, fill = "blue")

    grid.text("helix", unit(arrowpos + 15, "native") + unit(4, "mm"), unit(1, "npc") -  unit(2, "mm"), just = c("left", "top"), gp = gpar(cex = .9))

    coilpos = convertX( unit(arrowpos + 15, "native") + unit(4, "mm") + unit(1, "strwidth", data="helix  "), "native", TRUE)
    drawCoil(coilpos, coilpos + 15, height=.6*convertY(unit(1, "strheight", data="S"), "npc", TRUE), center.y = ceny, gp = gpar(fill = "red", alpha = .5), scorecol = "black")
    yline = .3
    height = .3
    grid.lines(unit(xlim, "native"), unit(yline, "npc"), gp = gpar(col = "grey50", lex = 1.5))
    datlag = dat[-1,]
    
    change.helix = which((datlag$helix >= cutoff) != (dat$helix[-nrow(dat)] >= cutoff))
    change.strand = which((datlag$strand >= cutoff) != (dat$strand[-nrow(dat)] >= cutoff))

    if(dat$helix[1] >= cutoff)
      change.helix = c(1, change.helix)

    if(dat$strand[1] >= cutoff)
      change.strand = c(1, change.strand)
    
    if(dat$helix[nrow(dat)] >= cutoff & !(nrow(dat) %in% change.helix))
      change.helix = c(change.helix, nrow(dat))
    
    if(dat$strand[nrow(dat)] >= cutoff & !(nrow(dat) %in% change.helix))
      change.strand = c( change.strand , nrow(dat))
    
    if(length(change.helix) >= 2)
      {
        for(i in seq(1, length(change.helix) - 1, by = 2))
          {
            if (change.helix[i+1] - change.helix[i] >= 10)
              drawCoil(change.helix[i], change.helix[i+1], height = height, center.y = yline, gp = gpar(fill = "red", alpha = .5), scorecol = "black")
          }
      }
    if(length(change.strand) >= 2)
      {
        
        for(j in seq(1, length(change.strand) - 1, by = 2))
          {
            if(change.strand[j+1] - change.strand[j] >= 10)
              
              drawArrow(change.strand[j], change.strand[j + 1], height = height, center.y = yline, fill = "blue")
          }
      }
    TRUE
  }


drawArrow = function(start, end, height, head.height = height*1.6, center.y , nativelims = current.panel.limits(), head.length = NULL, inlay.x, inlay.y, inlay.col, fill = "white", draw = TRUE, gp = gpar(fill = fill, alpha = .5))
  {
    xlim = nativelims$xlim
    ylim = c(0, 1) #nativelims$ylim
    rangex = xlim[2] - xlim[1]
    rangey = ylim[2] - ylim[1]

    tstart = (start - xlim[1]) / rangex
    tend = (end - xlim[1]) / rangex
    theight = height / rangey
    thead.height = head.height / rangey
    tcenter.y = (center.y - ylim[1]) / rangey

    if (is.null(head.length))
      thead.length = min(.02, .75*(tend - tstart))
    else
      thead.length = head.length / rangex

    xs = c(tstart, tend - thead.length, tend - thead.length, tend, tend - thead.length, tend - thead.length, tstart)
    ys = c(tcenter.y + theight/2, tcenter.y + theight/2, tcenter.y + thead.height/2, tcenter.y, tcenter.y - thead.height/2, tcenter.y - theight/2, tcenter.y - theight/2)
    #draw arrow
    pic = grid.polygon(xs, unit(ys, "npc"),default.unit = "npc",  gp = gp)
 
    pic
  }

drawCoil = function(start, end, height, center.y, nativelims = current.panel.limits(), inlay.x, inlay.y, inlay.lims, inlay.col, gp = gpar(fill = "white"), scorecol = "white")
  {
    
    xlim = nativelims$xlim
    ylim = c(0, 1) #nativelims$ylim
    rangex = xlim[2] - xlim[1]
    rangey = ylim[2] - ylim[1]

    tstart = (start - xlim[1]) / rangex
    tend = (end - xlim[1]) / rangex
    theight = height / rangey
    tcenter.y = (center.y - ylim[1]) / rangey
    #draw rect
    grid.rect(x = tstart, y = unit(tcenter.y - theight/2, "npc"), just = c("left", "bottom"), height = theight, width = tend - tstart, gp = gp, default.unit = "npc")
    #draw scoring
    #scoring is not perfect but will do for now
    if(tstart + .005 < tend)
      {
        myseq = seq(tstart + .005, tend, by = .005)
        grid.segments(myseq, tcenter.y + theight/2, myseq - .005 , tcenter.y - theight/2, gp = gpar(col = scorecol))
      }
    TRUE
  }

#x is the position on the protein, y is the pfam bin (pfam$bin)
panel.PFAM = function(x, y, end, subscripts, labs, vertGuides, ...)
  {
    drawVertGuides(vertGuides, "grey80")
    grid.text("PFAM Domains", unit(2, "mm"), unit(1, "npc") - unit(2, "mm"), just = c("left", "top"), gp = gpar(fontface="bold", cex=.8 ))
    #if(!is.null(pfamBins))
    if(any(!is.na(x)))
      drawPFAM(data.frame(start = x, end = end, labels = labs, bin = y), labcol = "labels", xlim = current.panel.limits()$xlim)
    else
      {
        grid.text("No PFAM data available")
      }
  }

drawPFAM = function(dat, poscolumns = c("start", "end"), labcol = "featureName", xlim)
  {
    nrow = length(levels(dat$bin))
    hseq = unique(dat$bin)#seq(min(dat$bin, na.rm=TRUE), max(dat$bins, na.rm=TRUE))
    hseq = hseq[!is.na(hseq)]
    grid.segments(y0 = unit(hseq, "native"), y1= unit(hseq, "native"), gp = gpar(col = "grey80"))
    
    dat = dat[apply(dat, 1, function(x) all(!is.na(x))),]
    
    #dat$bin = bins
    protlen = xlim[2] - xlim[1]
    apply(dat, 1, function(x,  poscolumns, nrow, protlen)
          {
            
            st = as.numeric(x[poscolumns[1]])
            end = as.numeric(x[poscolumns[2]])
     
            bin = as.numeric(x["bin"])
            #determine label and color
            fulln = x[labcol]

            fullnlen = convertWidth(unit(nchar(fulln), "char"), "native", valueOnly= TRUE)
            barlen = end - st
  
            lab  = fulln
            #grab the color from the color mapping dataframe
            col = pfColMap[as.character(pfColMap[,1]) == lab, 2]
            #check if the name is unrecognised
            if (!length(col))
              col = "#111111" #grey
            
            rot = 0
            #check if domain is very short eg NOTCH. If so, rotate the name
            if( (barlen) / protlen <= .04)
              {
                rot = 90
                lab = substr(lab, 1, 4)
                if(substr(lab, 4, 4) == "_")
                  lab = substr(lab, 1, 3)
                textypos = unit(nrow - bin + 1 , "native")
              } else {
                textypos =  unit(nrow - bin + 1 , "native") - unit(1.2, "char")
              }
                
            
            
            step = .6/nrow
            ypos = 1 - (.15 / nrow ) - step*bin
            #grid.lines(unit(c(st, end), "native"), unit(nrow - bin + 1, "native"), gp = gpar(col = col, lex = 20, alpha=.5))
            grid.rect(unit(st, "native"),
                      unit(nrow - bin + 1, "native"),
                      width = unit(end - st, "native"),
                      height = unit(.45, "native"),
                      just = c("left", "center"), gp = gpar(fill = col, alpha = .5)) 
            grid.text(lab, unit((st + end) /2, "native"), y =textypos, rot = rot, gp = gpar(cex = .67))

          },  poscolumns = poscolumns, nrow = nrow, protlen = protlen)

          TRUE
  }

if(FALSE)
  {
panel.metaCountSimple = function(x, y, subscripts, patientid, colpalette, ...)
  {
    if(missing(colpalette))
      cols =  rev(heat.colors(11) )
    else
      cols = colpalette

    patientid = patientid[subscripts]
    grid.rect(.5, .5, 1, 1, gp = gpar(fill = "grey95"))
    
    ylevs = unique(y)
    grid.segments(unit(0, "npc"), ylevs, unit(1, "npc"), ylevs, default.units="native", gp = gpar(col = "grey60"))
        thing = lapply(ylevs, function(ylev, x, y, patid)
      {
        inds = which(y == ylev)
        
        x = x[inds]
        patid = patid[inds]
        sampsize = length(unique(patid))
        counts = sort(table(x))
        non-zeros = which(sampsize != 0)
        props = counts / sampsize
        
        list(samplesize = sampsize, proportions = props, x = x, y  = ylev, counts = counts)
      }, x = x, y = y, patid = patientid)

    lapply(thing, function(myl)
           {
             counts = myl$counts
             p = myl$proportions
             y = myl$y
             n = myl$samplesize
             
             #browser()
             colpos = sapply(p, function(x) min(ceiling(x/.005), 11))
             grid.circle(as.numeric(names(counts)), y, default.units = "native", 5, gp = gpar(fill = cols[colpos]))
           })
    
  }

}

panel.metaCount = function(x, y, end, subscripts, patientid, scale.factor = 8, logscale = FALSE, logbase = exp(1), at.baseline = FALSE,colpalette = rev(brewer.pal(11, "RdBu")), legend.step = .005, levels = levels(y), vertGuides, lose1 = FALSE, sequence.counts = NULL, title, ...)
  {
    if(sum(!is.na(x)) == 0)
      return(TRUE)

    if(length(subscripts) != length(x))
      stop("Multiple panels are not currently supported. Please contact the maintainer if you need this functionality")

    #lose1 is true if we added a fake x value of 1 to the data.frame to get lattice to get the right plotting limits. It is not real data so we need to remove it.
    if(lose1)
      {
        len = length(x)
        x = x[ -len ]
        y = y[ -len ]
        end = end[ -len ]
        patientid = patientid[ -len ]
      }

    #if non-indels are indicated by end == NA, then change it so that start == end
    nas = which( is.na( end ) )
    end[ nas ] = x[ nas ]
    
    #find indels
    indels = which( x != end )
    if(length(indels))
      {
        indeldat = data.frame( start = x[indels], end = end[indels], category = y[indels])
                                          #now remove indels and plot snps as before
        x=x[-indels]
        y = y[-indels]
      }
    ylevs = unique(y)
    #get rid of the fake categories we added to make room for the legend they will not have ( as they did not have counts added to them.
    ylevs = ylevs[grepl("\\(", as.character(ylevs))]
    
    yseq = seq(min(as.integer(ylevs), na.rm=TRUE)  , max(as.integer(ylevs), na.rm=TRUE) )
   
    grid.rect(unit(.5, "npc"), unit(yseq, "native"), unit(1, "npc"), unit(1, "native"), gp = gpar(fill = c(rgb(217, 224, 235, max= 255), rgb(205, 205, 205, max = 255)), col = "grey95"))
 
        scaleseq = seq(min(as.integer(ylevs), na.rm=TRUE) - .5 , max(as.integer(ylevs), na.rm=TRUE) +.5)
 
    ygridseq = seq(min(as.integer(y), na.rm=TRUE) - .5, max(as.integer(y), na.rm=TRUE) + .5, by = 1)
    
    grid.segments(unit(0, "npc"), ygridseq  , unit(1, "npc"), ygridseq, gp = gpar(col = "white"), default.units = "native")
        
    ylim = current.panel.limits()$ylim
    gray.yrange = .5 - ylim[1]
    grid.rect(unit(.5, "npc"), unit(.5 - gray.yrange/2, "native"), unit(1, "npc"), unit(gray.yrange, "native"), gp = gpar(fill = "grey35"))
    drawVertGuides(vertGuides, "white")

    haveSeqCounts = !is.null(sequence.counts)
    thing = lapply(ylevs, function(ylev, x, y, patid, haveSC)
      {
        
        inds = which(y == ylev)
        x = x[inds]
        patid = patid[inds]
        if(!haveSC)
          sampsize = length(unique(patid))
        else
          {
            cat = gsub( " \\(.*", "", as.character(ylev))
            sampsize = sequence.counts$count[sequence.counts$category == cat]
          }

        counts = sort(table(x), decreasing = TRUE)
        if(sampsize == 0)
          props = rep(0, times = length(counts)) #XXX should this be NA, Inf, ...?
        else
          props = counts / sampsize
        list(samplesize = sampsize, proportions = props, x = x, y  = ylev, counts = counts)
      }, x = x, y = y, patid = patientid, haveSC = haveSeqCounts)

    lapply(thing, function(myl)
           {
             if(sum( !is.na(myl$x) ) > 0)
               {
                 counts = myl$counts
                 p = myl$proportions
                 y = myl$y
                 n = myl$samplesize
  

                 heights = calculateBarHeights(counts, logbase, scale.factor)

                 ypos = as.integer(y) - .5
                 
                 vjust = 0
                 colinds = sapply(p, function(x) min(ceiling(x / legend.step), 11))
                 
                 xpos = as.numeric(names(counts))
                 #grid.rect(xpos, ypos, vjust = vjust, default.units = "native", width = unit(.5, "mm"), heights, gp = gpar(fill = colpalette[colinds], col = colpalette[colinds]), )
                 grid.rect(xpos, ypos, just=c("center", "bottom"), default.units = "native", width = unit(.5, "mm"), heights, gp = gpar(fill = colpalette[colinds], col = colpalette[colinds]), )
                 #draw cross bar at the top of bars that hit the cap
                 
                 topbars = which(counts >= logbase ^ scale.factor)
                 if(length(topbars))
                   {
                     print("topbars required")
                     grid.rect(xpos[topbars], ypos + heights[topbars], width = unit(2, "mm"), height = unit(.25, "mm"), default.units = "native", gp = gpar(fill = colpalette[colinds[topbars]], col = colpalette[colinds[topbars]]))
                   }
               } #else
             #Logic for empty rows here
           })

        #now add indels

    if(length(indels))
      {
        by(indeldat, indeldat$category, function(x)
           {
             
             myIRange = IRanges(start = x$start, end = x$end)
             cov = as.numeric(coverage(myIRange)) #coverage returns coverage starting at 0, not at the beginning of our earliest observation!
             
             x.s = seq(min(x$start), max(x$end))
             cov = cov[x.s]
             covlag = c(0, cov[-length(cov)])
             covcat =cumsum( (cov > 0) & (covlag == 0) )  
             y = as.integer(x$category[1])
             nonzero = which(cov != 0)
             heights = rep(0, times=length(cov))
             heights[nonzero] = calculateBarHeights(cov[nonzero], logbase, scale.factor)
             if(FALSE)
               {
             heights = sapply(cov, function(x)
               {
                 if( x == 0)
                   0
                 else
                   {
                     calculateBarHeights(x, logbase, scale.factor)
                   }
               })
           }
             useseq = which(heights > 0)
  
             by(data.frame(x = x.s[useseq], ht = heights[useseq]), covcat[useseq], function(dat)
                {
                  drawOneIndel(dat$x, dat$ht, y = y)
                })
           } )
           
      }
    #add the legend!!!
    #avoid the lazy evaluation of the damned
    catht = convertY(unit(1, "native"), "mm", valueOnly = TRUE)
    makeFullLegend(c(5, 15, 45), colpalette, legend.step, scale.factor, logbase, title, catht, ylim = c(floor(ylim[2])  - 2 + .5, ylim[2]), yunit= "native")
    #makeFullLegend = function(values , colpalette, legend.step, scale.factor, log.base, title, one.cat.height,  xlim = c(0, 1), xunit = "npc", ylim = c(0, 1), yunit = "npc", draw.border = TRUE, indelcol ="#00AA00")
    
    TRUE
  }

drawOneIndel = function(xs, heights, y)
  {
      indelseq = c(y - .5 + heights, rep(y - .5, times = length(xs)))
    
    grid.polygon(x = unit( c(xs, rev(xs) ), "native"), y = unit(  indelseq, "native"), gp = gpar(stroke=NULL, fill="#00AA00", alpha=.5) )
    
    TRUE
  }
 
proteinStructPlot = function(pfam, structPred, hydro, transMem, sigP, xlim, tmposcol = c("start", "end"), main = NULL, pfamLabels = "featureName", draw = FALSE)
  {

    pfam = fixPFAM( pfam , pfamLabels)
       
    plots = makeStructPlots(pfam, structPred, hydro, transMem, sigP, xlim, tmposcol, pfamLabels = pfamLabels)
    hydroPlot = plots$hydro

    structPredPlot = plots$structPred
    pfamPlot = plots$pfam
    
    cplot = update(
      c(pfamPlot,
        structPredPlot,
        hydroPlot,
        x.same = TRUE),
      layout = c(1, 3),
      xlim = xlim,
      main = main,
      par.settings = list(layout.heights = list(panel = c(.55, .55, 1.5))) #only one pfam row now!
      )
    if(draw)
      print(cplot)
    else
      cplot

  }

makeStructPlots = function(pfam, structPred, hydro, transMem, sigP, xlim, tmposcol = c("start", "end"), main = NULL, pfamLabels = "featureName", cutoff = if(max(structPred$helix > 1)) 6 else .6, vertGuides)
  {

    if(dim(hydro)[2] == 0 | is.null(hydro))
      hydro = data.frame(start = numeric(), featureValue = numeric())

    hydroPlot = xyplot(featureValue ~ start, type = "l", col = "black", data = hydro, tm = transMem, sig = sigP, xlim = xlim, panel = panel.protstruct, axis = axis.combined, tmposcol = tmposcol, ylab= NULL, vertGuides = vertGuides) 

    if(dim(structPred)[2] == 0 | is.null(structPred))
      structPred = data.frame(start = numeric(), helix = numeric(), strand = numeric() )
    
    structPredPlot = xyplot(helix ~ start, data = structPred, strand = structPred$strand, cutoff= cutoff, panel= panel.psipred, ylim = c(0, 1), ylab = NULL, xlab = NULL, vertGuides = vertGuides)

    #if(dim(pfam)[2] == 0 | is.null(pfam))
    if(is.null(pfam))
      {
        pfam = data.frame(start = numeric(), end = numeric(), pfamLabels = character(), bin = numeric())
        labs = character()
        ylim = 1:2
        xlim = 1:10
      }
    else
      {
        #We pad the pfam domains with an extra bin for spacing
        newlev = max(as.integer( pfam$bin ), na.rm = TRUE ) + 1
        levels( pfam$bin ) = c( levels( pfam$bin ) , newlev )
        pfam[ nrow( pfam ) + 1, "bin" ] = newlev
        labs = pfam[ , pfamLabels ]
        ylim = range( as.integer( pfam$bin ) )
        xlim = range( pfam$start )
      }
    pfamPlot = xyplot(bin~start, end = pfam$end,  data= pfam, labs = labs, panel = panel.PFAM, ylab = NULL, xlab = "Amino Acid Position", vertGuides = vertGuides, ylim = ylim, xlim = xlim)
    
    list(hydro = hydroPlot, structPred = structPredPlot, pfam = pfamPlot)
  }

metaCountStructPlot = function(events,catname = "PRIMARY_TISSUE", requiredCats = NULL, position = c("protpos", "protposend"),  pfam, pfamLabels = "featureName",structPred, hydro, transMem, sigP, xlim, tmposcol = c("start", "end"), main = NULL, simple = FALSE, at.baseline = TRUE, logscale = TRUE, logbase = 1.506, scale.factor = 10, colpalette = rev(brewer.pal(11, "RdYlBu")), legend.step = .01, sampleID, key , subtitle = "Amino Acid Position", draw = FALSE, vertGuides = 10, sequence.counts = NULL)
  {

    pfam = fixPFAM(pfam, pfamLabels)
    
    plots = makeStructPlots(pfam, structPred, hydro, transMem, sigP, xlim, tmposcol, pfamLabels = pfamLabels, vertGuides = vertGuides)
    
    hydroPlot = plots$hydro

    structPredPlot = plots$structPred
    pfamPlot = plots$pfam
 
    if(!is.null(requiredCats))
      {
#        obscats = unique(as.character(events[[catname]]))
#        additcats = obscats[which( !( obscats %in% requiredCats ) )]
        #XXX add fake categories last!!!!
 #       cats = unique(c(requiredCats, additcats, "fake1XXX", "fake2XXX", "fake3XXX"))
  #      charvals = as.character(events[[catname]])
        
   #     events[[catname]] = factor(charvals, levels = cats)

   #     for(reqcat in c( requiredCats, "fake1XXX", "fake2XXX", "fake3XXX") )
    #      events[nrow(events) + 1 , catname] = reqcat
        events = spoofLevelsInDF(events, catname, requiredCats)
      }


                                        #deal with missing categories (even though there really shouldn't be any!!!!)
    #this will interfere with matching the passed in sequence counts. We weren't using it anyway, so it's disabled for now
    if(FALSE)
      {
        y = events[[catname]]
        levs = levels(y)
        missingCat = which(is.na(y))
        tmpcat = as.character(y)
        tmpcat[missingCat] = "UnCategorized"
        events[[catname]] = factor(tmpcat, levels = c("UnCategorized", levels(events[[catname]])))
      }
    #if there are no non-na x values lattice refuses to draw properly
    if( all( is.na( events[[ position[ 1 ] ]] ) ) )
      {
        lose1 = TRUE
        events[nrow(events) + 1, position[1] ] = 1
      } else {
        lose1 = FALSE
      }

            #sampleID may be passed in as the name of a column, but we need to have the actual data.
    if(length(sampleID) != 1 & length(sampleID) != nrow(events) )
      stop("sampleID must be either the name of a column in the events data.frame or a vector with the same number of elements that events has rows.")
    if(is(sampleID, "character") & length(sampleID) == 1)
      sampleID = events[,sampleID]
# counts = tapply(c(rep(NA, times=length(requiredCats)), sampleID),


    numcats = length(requiredCats)
        #add counts to category names
    if (!is.factor(events[[catname]]))
      events[[catname]] = factor(events[[catname]])
    if(is.null(sequence.counts))
      scounts = rep(NA, times = numcats)
    else
      {
        #XXX If we are demanding this information be the same, why do we need it in 2 places??
        if(length(setdiff(as.character(sequence.counts$category), requiredCats)))
          stop("Categories in sequence.counts and requiredCats do not match")
        ordinds = unlist(sapply(as.character(sequence.counts$category), function(x, cats)
          {
            which(x == cats)
          }, cats = levels(events[[catname]])))
        scounts = sequence.counts$count[ordinds]
      }
   
    mutcounts = tapply( sampleID,
      events[[catname]] ,
      function(x) sum(!is.na(x)))
    
    levels(events[[catname]]) = paste(levels(events[[catname]]), " (", mutcounts, " , ", scounts, ")", sep = "")

    events = spoofLevelsInDF(events, catname, c("fake1", "fake2"), before = FALSE)

    countPlot = xyplot(as.formula(paste(catname, "~", position[1])), end = events$end, data = events, panel = panel.metaCount,  patientid = sampleID, at.baseline = at.baseline, logscale = logscale, scale.factor = scale.factor, logbase = logbase, colpalette = colpalette, legend.step = legend.step, vertGuides = vertGuides, lose1 = lose1, sequence.counts = sequence.counts, title=main)

    #leg = makeColorLegend(colpalette, scale.factor, legend.step)
    cat.names = levels(events[[catname]])
    combPlot = combinePlots(countPlot, pfamPlot, structPredPlot,  cat.names, main, subtitle, xlim , colpalette)
      
    if(draw)
      print(combPlot)
    else
      combPlot
  }

combinePlots = function(countPlot, pfamPlot, structPlot,  cat.names, main, subtitle, xlim, col.palette)
  {
  
    combPlot = c(  #hydroPlot,
      structPlot,
      pfamPlot,
      countPlot,
      x.same=TRUE, y.same=NA)#, merge.legends=TRUE)

    leftpad = max(nchar(cat.names))*.67
  
    panelLayout = c(  #.3,
      .25,
      .20*2, #*max(pfamBins, 2),
        .20*length(cat.names))
    combPlot = update(combPlot,
                                        #layout = c(1, 4),
      layout = c(1, 3),
      xlim = xlim,
      axis = axis.combined,
      par.settings = list(
        layout.heights = list(
          panel = panelLayout,
          key.top = .20,
          xlab.top = 0,
          axis.top = 0),
        layout.widths = list(
          right.padding = 5,
          left.padding = leftpad)),
      #legend = list(top = list( fun = makeColorLegend, args = list( colpalette = col.palette))),    
      #main = main,
      xlab = subtitle,
      ylab.right = list(label = "variant counts", vjust = 0, rot = -90,
        y =  1 - panelLayout[3] / ( 2 * sum(panelLayout) ))
      )
    
    return(combPlot)
  }

fixPFAM = function(pfam,  labels)
  {

    if(!is.null(pfam) & sum(dim(pfam)[1]))
      {
        pfamIRange = IRanges(start = pfam$start, end = pfam$end, names = pfam[ , labels ])
        pfamBins = disjointBins(pfamIRange)
        pfam$bin = factor(pfamBins, levels = 0:(max(pfamBins + 1)))
        
      } else {
        pfamBins = NULL
        pfam = NULL
      }
    pfam
  }


makeFullLegend = function(values , colpalette, legend.step, scale.factor, log.base, title, one.cat.height,  xlim = c(0, 1), xunit = "npc", ylim = c(0, 1), yunit = "npc", draw.border = TRUE, indelcol ="#00AA00")
  {
    totxwid = diff(xlim)
    totywid = diff(ylim)
    legvp = viewport(unit(xlim[1] + totxwid/2, xunit), unit(ylim[1] + totywid / 2, yunit), width = unit(totxwid, xunit), height = unit(totywid, yunit), xscale = xlim, yscale = ylim)
    pushViewport(legvp)
    allgrobs = list()
    if(draw.border)
      allgrobs[[length(allgrobs) + 1]] = grid.rect( draw=TRUE, gp = gpar(fill= "#FFFFFF"))  #box the entire area

    #figure out where the two boxes go.
    
    
    xwid2 = .190
    xwid1 =  .15
    xpos2 = convertX( unit(1, "npc") - unit(2, "mm") -unit(xwid2/2, "npc"), "npc", valueOnly = TRUE)
    xpos1 = xpos2 - xwid2 / 2 - xwid1 / 2

    height = unit(1, "npc") - unit(4, "mm")

    allgrobs[[length(allgrobs) + 1]] = grid.rect(c(xpos1, xpos2), y = unit(.5, "npc"), width = c(xwid1, xwid2), height = height, draw = TRUE)
     
    yscale = c(convertY(unit(2, "mm"), "native", valueOnly = TRUE), convertY(unit(1, "npc") - unit(2, "mm"), "native", valueOnly = TRUE))
    vp1 = viewport(x = xpos1, y = unit(.5, "npc"), width = xwid1, height = height, yscale = yscale)
    vp2 = viewport(x = xpos2, y = unit(.5, "npc"), width = xwid2, height = height, yscale = yscale)

    #draw colorscale and text
    pushViewport(vp2)
    allgrobs = colorScaleBox(colpalette, legend.step, allgrobs, indelcol)
    popViewport(1) #vp2

    #draw height scale and text
    pushViewport(vp1)
    allgrobs = heightScaleBox(values, scale.factor, one.cat.height, log.base, allgrobs)
    popViewport(1) #vp1
    
    allgrobs[[length(allgrobs) + 1]] = grid.text(title, unit(2, "mm"), gp = gpar(fontsize = 24, fontfamily = "Helvetica"), #, fontface = "bold"),
              just = "left")
    popViewport(1) #legvp
    do.call(gList, allgrobs)
  }

heightScaleBox = function(values, scale.factor, one.cat.height, log.base, groblist = list())
  {
    #put lines at positions 1/5, 2/5, 3/5, 4/5 npc
    positions = (1:4) / 5
    
    if(length(values) ==3) #didn't specify max twice
      values[4] = floor(log.base ^ scale.factor)

    heights = calculateBarHeights(values, log.base, scale.factor, one.cat.height)

    groblist[[length(groblist) + 1]] = grid.segments(unit(positions, "npc"), .25, unit(positions, "npc"),unit(.25, "npc") +  unit(heights, "mm"))
    yTpos = unit(.25, "npc") + unit(max(heights), "mm")
    Twidth = convertX(unit(2, "mm"), "npc", valueOnly= TRUE)
    groblist[[length(groblist) + 1]] = grid.lines((4/5  + c(-1, 1)*Twidth), yTpos)
    htlabs = as.character(values)
    htlabs[4] = paste(htlabs[4], "+", sep="")
    groblist[[length(groblist) + 1]] = grid.text(htlabs, positions, .125, gp = gpar(cex=.7))

    groblist[[length(groblist) + 1]] = grid.text("variant counts\n(capped log scale)", x = unit(1, "mm"), y = unit(1, "npc") - unit(1, "mm"), gp = gpar(cex=.7), just = c("left", "top"))
    groblist
  }

calculateBarHeights = function(counts, logBase, denom, totHeight= 1)
  {
    if(!is.null(logBase)) #we are in log scale
      heights = sapply(counts, function(x) .05 + .875 * min( log( x, base = logBase), denom ) / denom)
    else
      heights = sapply(counts, function(x) .9 * min( x, denom ) / denom)

    heights * totHeight
  }
    
colorScaleBox = function(colpalette, legend.step, groblist = list(), indelcol)
  {
                                        #top text
    groblist[[length(groblist) + 1]]  = grid.text("variant rel. frequency \n(by position and category)", x = .05, y = .95, just = c("left", "top"), gp = gpar(cex = .7))
    
                                        #color scale and labels
    maxpct = ( length(colpalette) -1 ) * legend.step
    if (maxpct < 1)
      maxpct = 100 * maxpct
    maxpct = floor(maxpct) #get an integer
    maxpctlab = paste(">", maxpct, "%", sep="")
    
    pos = unit(2, "mm")
    totalColSpace = convertX(unit(1, "npc") - unit(5, "mm") - unit(.7, "strwidth", "0%") - unit(.7, "strwidth", maxpctlab), "npc", valueOnly = TRUE)
    groblist[[length(groblist) + 1]] = grid.text("0%", x = pos, y = 3/6, just = "left", gp = gpar(cex=.7))
    pos = pos + unit(.7, "strwidth", "0%") + unit(1, "mm")
    
    xpos = convertX(pos, "npc", valueOnly = TRUE)
    allxpos = xpos + seq(0, by = totalColSpace /length(colpalette), length.out = length(colpalette))
    groblist[[length(groblist) + 1]] =
      grid.rect(x = unit(allxpos, "npc"), y = 3/6, hjust = 0, width = totalColSpace / length(colpalette), height = 2/9, gp = gpar(fill = colpalette, col=colpalette))
    pos = pos + unit(totalColSpace, "npc") +  unit(1, "mm")

    groblist[[length(groblist) + 1]] = grid.text(maxpctlab, x = pos, y = 3/6, just = "left", gp = gpar(cex=.7))
                                        #( 0 : length( colpalette ) ) * totalColSpace / length(colpalette), "npc") 
    
    #indel legend
    groblist[[length(groblist) + 1]] = grid.text("indel:", x = unit(2, "mm"), y = 1/6, just = "left", gp=gpar(cex=.7))

    indelxpos = unit(2, "mm") + unit(.7, "strwidth", "indel:") + unit(.5, "mm")
    groblist[[length(groblist) + 1]] = grid.rect(x = indelxpos, y = 1/6,  width = .2, height = 2/9, gp = gpar(col = indelcol, fill = indelcol, alpha = .5), just = "left")

    groblist
  }


makeColorLegend = function(colpalette)
  {

    #lyt = trellis.currentLayout()
    
    ncols = length(colpalette) + 3
    boxwidth = 1 / (2 * (ncols + 3) )
    lab = grid.text("Location Percent Mutated", y = .85, x = boxwidth, gp = gpar(cex = .9), draw = FALSE, just = "left")

    colboxes = lapply(seq(along = colpalette), function(pos, cols, width)
      {
        grid.rect(x = (pos+1)*width , y = .30, width = width, height = .40, gp = gpar(col = cols[pos], fill = cols[pos]), draw = FALSE)
      },  cols = colpalette, width = boxwidth )


    length(colboxes) = ncols
    n = length(colpalette)
    colboxes[[n+1]]= grid.text("0%", x = unit(1*boxwidth, "npc"), y = .3, draw = FALSE, just = c("right", "center"), gp = gpar(cex = .8))
    colboxes[[n+2]]= grid.text(">10%", x = (n + 2)*boxwidth, y = .3, just = c("left", "center"), draw = FALSE, gp = gpar(cex = .8))
    colboxes[[n+3]] = grid.rect(x = (n + 7) * boxwidth, y = .3, width = boxwidth, height = .4, gp = gpar(stroke=NULL, fill="#00AA00", alpha=.5), draw = FALSE )
    colboxes[[n+4]] = grid.text("indel", x = (n+8) * boxwidth, y = .3, just = c("left", "center"), draw = FALSE, gp = gpar(cex = .8))
    colboxes[[n+5]] = lab
    
    mychildren = do.call("gList", c(colboxes))# , txtlabs ))
    gTree(height = unit(1, "inches"), children = mychildren)
  }

axis.combined = function(side, ...)
  {
    #packnum = get("packet.number", sys.frame(3))
    packnum = packet.number()
    switch(side,
           left = {
              #This is  a super-hack but it's the only way I have figured out to do it.
             #It is likely that this will fail if we ever try to do a conditional form of these plots
             
             if(!is.null(packnum))
               {
                     if(packnum == 3)
                       {
                         
                         args = list(...)
                         #get rid of fake labels!
                         labs = args$components$left$labels$labels
                         labs[!grepl("(", labs, fixed=TRUE)] = "."
                         args$components$left$labels$labels = labs
                         #get rid of fake ticks!
                         length(args$components$left$ticks$at) = length(args$components$left$ticks$at) - 3
                         args$rot = 0
                         args$components$left$tck = c(0, 0)
                         do.call("axis.default", c(side, args))
                       }                 
               }
           },
            bottom = {
              axis.default(side , ...)
            }
         
            )
   }

drawVertGuides = function(num, col = "black")
  {

    if(num)
      {
        myseq = seq(0, 1, length.out = num)
        
        grid.segments(x0 = unit(myseq, "npc"),
                      x1 = unit(myseq, "npc"), gp = gpar(col = col))
      }
    TRUE
  }
