

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

    if(!is.null(tm) && nrow(tm) > 0 )
      {
        tmmat = combineTMDomains(tm, poscols = tmposcols)
        apply(tmmat, 1, function(x, ylim)
              {
                grid.polygon( c( x , rev( x ) ), rep( ylim , times = c(2, 2) ) , gp = gpar(fill = "grey75", alpha = .5, col = "white"), default.units = "native", draw = TRUE)
                 }, ylim = ylim)
      }

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
    #comes after tm/sigp stuff to deal with alpha issues
    grid.text("Hydrophobicity", unit(2, "mm"), unit(1, "npc") - unit(2, "mm"), just = c("left", "top"), gp = gpar(fontface="bold", cex=.9 ))

   #browser()
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


panel.psipred = function(x, y, subscripts, cutoff, strand, vertGuides, ...)
  {
    drawVertGuides(vertGuides, "grey80")
      
    drawPsipred(data.frame(start = x,helix = y, strand = strand), cutoff, current.panel.limits()$xlim)
  }

drawPsipred = function(dat, cutoff, xlim)
  {
    #grid.text("Secondary Structure", unit(1, "char"), .75, just = "left")
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


panel.PFAM = function(x, y, end, subscripts, labs,  pfamBins, vertGuides, ...)
  {
    drawVertGuides(vertGuides, "grey80")
    grid.text("PFAM Domains", unit(2, "mm"), unit(1, "npc") - unit(2, "mm"), just = c("left", "top"), gp = gpar(fontface="bold", cex=.8 ))
    if(!is.null(pfamBins))
      drawPFAM(data.frame(start = x, end = end, labels = labs, bin = y), labcol = "labels", bins = pfamBins, xlim = current.panel.limits()$xlim)
    else
      {
        grid.text("No PFAM data available")

      }
  }

drawPFAM = function(dat, poscolumns = c("start", "end"), labcol = "featureName", bins, xlim)
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
            #check if domain is very short eg NOTCH
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
                
            
            #.45 for one row ...
            step = .6/nrow
            ypos = 1 - (.15 / nrow ) - step*bin
            #grid.lines(unit(c(st, end), "native"), unit(nrow - bin + 1, "native"), gp = gpar(col = col, lex = 20, alpha=.5))
            grid.rect(unit(st, "native"),
                      unit(nrow - bin + 1, "native"),
                      width = unit(end - st, "native"),
                      #height = unit(10, "points"),
                      height = unit(.45, "native"),
                      just = c("left", "center"), gp = gpar(fill = col, alpha = .5)) 
            grid.text(lab, unit((st + end) /2, "native"), y =textypos, rot = rot, gp = gpar(cex = .67))

          },  poscolumns = poscolumns, nrow = nrow, protlen = protlen)

          TRUE
  }

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
        props = counts/sampsize
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

panel.metaCount = function(x, y, end, subscripts, patientid, scale.factor = 8, logscale = FALSE, logbase = exp(1), at.baseline = FALSE,colpalette = rev(brewer.pal(11, "RdBu")), legend.step = .005, levels = levels(y), indel.overlay = FALSE, vertGuides, lose1 = FALSE, sequence.counts = NULL, ...)
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
    
    #patientid = patientid[subscripts]
    #grid.rect(.5, .5, 1, 1, gp = gpar(fill = "grey95"))
    
    #end = end[subscripts]
    
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
    yseq = seq(min(as.integer(ylevs))  , max(as.integer(ylevs)) )
   
    grid.rect(unit(.5, "npc"), unit(yseq, "native"), unit(1, "npc"), unit(1, "native"), gp = gpar(fill = c(rgb(217, 224, 235, max= 255), rgb(205, 205, 205, max = 255)), col = "grey95"))
    if(indel.overlay)
      {
        scaleseq = seq(min(as.integer(ylevs)) - .5 , max(as.integer(ylevs)) +.5)
        hdenom = scale.factor
        #grid.segments(unit(0, "npc"), scaleseq, unit(1, "npc"), scaleseq, default.units = "native", gp = gpar(col = c("black"), lex = 1.5))
    }
    else
      {
        hdenom = 2 * scale.factor
        grid.segments(unit(0, "npc"), yseq, unit(1, "npc"), yseq, default.units = "native", gp = gpar(col = c("black"), lex = 1.5))
    }

    ygridseq = seq(min(as.integer(y), na.rm=TRUE) - .5, max(as.integer(y), na.rm=TRUE) + .5, by = 1)
    
    grid.segments(unit(0, "npc"), ygridseq  , unit(1, "npc"), ygridseq, gp = gpar(col = "white"), default.units = "native")
    #grid.segments(unit(seq(1, max(x, na.rm=TRUE), length.out = 10), "native"), unit(0, "npc"), unit(seq(1, max(x, na.rm=TRUE), length.out = 10 ), "native"), unit(1, "npc"), gp = gpar(col = "white"))
    
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
          sampsize = sequence.counts$count[as.character(sequence.counts$category) == as.character(ylev)]
        counts = sort(table(x), decreasing = TRUE)
        props = counts/sampsize
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
                                        #browser()
    
                 if(logscale)
                   heights = sapply(counts, function(x) .05 + .9 * min( log( x, base = logbase), scale.factor ) / hdenom)
                 else
                   heights = sapply(counts, function(x) .9 * min( x, scale.factor ) / hdenom)

                 if(indel.overlay)
                   ypos = as.integer(y) - .5
                 else
                   ypos = as.integer(y)

                 vjust = 0
                 colinds = sapply(p, function(x) min(ceiling(x / legend.step), 11))

             
                 grid.rect(as.numeric(names(counts)), ypos, vjust = vjust, default.units = "native", width = unit(.5, "mm"), heights, gp = gpar(fill = colpalette[colinds], col = colpalette[colinds]), )
               } #else
             #Logic for empty rows here
           })

        #now add indels

    if(length(indels))
      {
        indelIRange = IRanges(start = indeldat$start, end = indeldat$end)
        by(indeldat, indeldat$category, function(x)
           {
             
             myIRange = IRanges(start = x$start, end = x$end)
             cov = as.numeric(coverage(myIRange)) #coverage returns coverage starting at 0, not at the beginning of our earliest observation!
             x.s = seq(min(x$start), max(x$end))
             cov = cov[x.s]
             y = as.integer(x$category[1])
             heights = sapply(cov, function(x)
               {
                 if( x == 0)
                   0
                 else
                   {
                     if(logscale)
                       .05 + .9 * min( log( x, base = logbase), scale.factor ) / hdenom
                 else
                   .9 * min( x, scale.factor ) / hdenom
                   }
               })

             if(indel.overlay)
               indelseq = c(y - .5 + heights, rep(y - .5, times = length(x.s)))
             else
               indelseq = c(y - heights, rep(y, times = length(x.s)))
             
              grid.polygon(x = unit( c(x.s, rev(x.s) ), "native"), y = unit(  indelseq, "native"), gp = gpar(stroke=NULL, fill="#00AA00", alpha=.5) )
           } )
      }
        
    TRUE
  }
 
proteinStructPlot = function(pfam, structPred, hydro, transMem, sigP, xlim, tmposcol = c("start", "end"), main = NULL, pfamLabels = "featureName", draw = FALSE)
  {
    
    pfamIRange = IRanges(start = pfam$start, end = pfam$end, names = pfam[,pfamLabels])
    pfamBins = disjointBins(pfamIRange)
    pfam$bin = factor(pfamBins, levels = 0:(max(pfamBins) + 1))
      
    plots = makeStructPlots(pfam, structPred, hydro, transMem, sigP, xlim, tmposcol, pfamLabels = pfamLabels, pfamBins = pfamBins)
    hydroPlot = plots$hydro

    structPredPlot = plots$structPred
    pfamPlot = plots$pfam

    
    cplot = update(c(pfamPlot,
             structPredPlot,
             hydroPlot,
             x.same = TRUE),
           layout = c(1, 3),
           xlim = xlim,
           main = main,
           par.settings = list(layout.heights = list(panel = c(.3 + .25*max(pfamBins), .55, 1.5)))
           )
    if(draw)
      print(cplot)
    else
      cplot

  }

makeStructPlots = function(pfam, structPred, hydro, transMem, sigP, xlim, tmposcol = c("start", "end"), main = NULL, pfamLabels = "featureName", cutoff = if(max(structPred$helix > 1)) 6 else .6, pfamBins, vertGuides)
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
#        pfam = data.frame(start = c(NA, NA, 1, 2), end = rep(NA, 4), pfamLabels = rep(NA, 4), bin = c(1, 2, NA,NA))
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
    pfamPlot = xyplot(bin~start, end = pfam$end,  data= pfam, labs = labs, panel = panel.PFAM, ylab = NULL, xlab = "Amino Acid Position", pfamBins = pfamBins, vertGuides = vertGuides, ylim = ylim, xlim = xlim)
    
    list(hydro = hydroPlot, structPred = structPredPlot, pfam = pfamPlot)
  }

metaCountStructPlot = function(events,catname = "PRIMARY_TISSUE", requiredCats = NULL, position = c("protpos", "protposend"),  pfam, pfamLabels = "featureName",structPred, hydro, transMem, sigP, xlim, tmposcol = c("start", "end"), main = NULL, simple = FALSE, at.baseline = TRUE, logscale = TRUE, logbase = 1.5, scale.factor = 10, colpalette = rev(brewer.pal(11, "RdYlBu")), legend.step = .01, sampleID, key , subtitle = "Amino Acid Position", draw = FALSE, indel.overlay = TRUE, vertGuides = 10, sequence.counts = NULL)
  {
    
    if(!is.null(pfam) & sum(dim(pfam)[1]))
      {
        pfamIRange = IRanges(start = pfam$start, end = pfam$end, names = pfam[,pfamLabels])
        pfamBins = disjointBins(pfamIRange)
        pfam$bin = factor(pfamBins, levels = 0:(max(pfamBins + 1)))
        
      } else {
        pfamBins = NULL
        pfam = NULL
      }
    plots = makeStructPlots(pfam, structPred, hydro, transMem, sigP, xlim, tmposcol, pfamLabels = pfamLabels, pfamBins = pfamBins, vertGuides = vertGuides)
    
    hydroPlot = plots$hydro

    structPredPlot = plots$structPred
    pfamPlot = plots$pfam


 
    if(!is.null(requiredCats))
      {
        obscats = unique(as.character(events[[catname]]))
        additcats = obscats[which( !( obscats %in% requiredCats ) )]
        cats = unique(c(requiredCats, additcats))
        charvals = as.character(events[[catname]])
        
        events[[catname]] = factor(charvals, levels = cats)

        for(reqcat in requiredCats)
          events[nrow(events) + 1 , catname] = reqcat
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


    numcats = length(unique(events[[catname]]))
        #add counts to category names
    if (!is.factor(events[[catname]]))
      events[[catname]] = factor(events[[catname]])
    if(is.null(sequence.counts))
      scounts = rep(NA, times = numcats)
    else
      {
        if(nrow(sequence.counts) != numcats )
          stop("Number of categories in sequence.counts does not match number of categories in data")
        ordinds = sapply(as.character(sequence.counts$category), function(x, cats)
          {
            which(x == cats)
          }, cats = levels(events[[catname]]))
        scounts = sequence.counts$count[ordinds]
      }
   
          mutcounts = tapply( sampleID,
      events[[catname]] ,
      function(x) sum(!is.na(x)))
    

    
    levels(events[[catname]]) = paste(levels(events[[catname]]), " (", mutcounts, " / ", scounts, ")", sep = "")
    

    countPlot = xyplot(as.formula(paste(catname, "~", position[1])), end = events$end, data = events, panel = panel.metaCount,  patientid = sampleID, at.baseline = at.baseline, logscale = logscale, scale.factor = scale.factor, logbase = logbase, colpalette = colpalette, legend.step = legend.step,  indel.overlay = indel.overlay, vertGuides = vertGuides, lose1 = lose1, sequence.counts = sequence.counts)


    leg = makeColorLegend(colpalette, scale.factor, legend.step)
    
    combPlot = c(  #hydroPlot,
      structPredPlot,
      pfamPlot,
      countPlot,
      x.same=TRUE, y.same=NA)#, merge.legends=TRUE)

    leftpad = max(nchar(as.character(unique(events[[catname]] ))))*.67
    panelLayout = c(  #.3,
      .25,
      .20*max(pfamBins, 2),
      .20*length(unique(events[[catname]])))
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
     legend = list(top = list(fun=leg)),
     
       
     main = main, xlab = subtitle,
     ylab.right = list(label = "mutation counts", vjust = -1, rot = -90,
     y =  1 - panelLayout[3] / ( 2 * sum(panelLayout) ))
   )  
   # c(combPlot, draw.key(key), layout = c(2, 1))
    
    if(draw)
      print(combPlot)
    else
      combPlot
  }


makeColorLegend = function(colpalette, scalefactor, step)
  {
    ncols = length(colpalette) + 3
    boxwidth = 1 / (ncols + 3)
    lab = grid.text("Location Percent Mutated", y = .85, x = boxwidth, gp = gpar(cex = .9), draw = FALSE, just = "left")

    colboxes = lapply(seq(along = colpalette), function(pos, cols, width)
      {
        grid.rect(x = (pos+1)*width , y = .30, width = width, height = .40, gp = gpar(col = cols[pos], fill = cols[pos]), draw = FALSE)
      },  cols = colpalette, width = boxwidth )

    if(FALSE)
      {
    
        txtlabs  = lapply(seq(along = colpalette), function(pos, cols, width)
          {
            grid.text(paste(pos, "%", sep = ""), x = (pos + 1)*width, y = .3, draw = FALSE)
          }, cols = colpalette, width = boxwidth)
      }
    length(colboxes) = ncols
    n = length(colpalette)
    colboxes[[n+1]]= grid.text("0%", x = unit(1*boxwidth, "npc"), y = .3, draw = FALSE, just = c("right", "center"), gp = gpar(cex = .8))
    colboxes[[n+2]]= grid.text(">10%", x = (n + 2)*boxwidth, y = .3, just = c("left", "center"), draw = FALSE, gp = gpar(cex = .8))
    colboxes[[n+3]] = grid.rect(x = (n + 3.5) * boxwidth, y = .3, width = boxwidth, height = .4, gp = gpar(stroke=NULL, fill="#00AA00", alpha=.5), draw = FALSE )
    colboxes[[n+4]] = grid.text("indel", x = (n+4.5) * boxwidth, y = .3, just = c("left", "center"), draw = FALSE, gp = gpar(cex = .8))
    colboxes[[n+5]] = lab
    
    mychildren = do.call("gList", c(colboxes))# , txtlabs ))
    gTree(height = unit(1, "inches"), children = mychildren)
  }



axis.combined = function(side, ...)
  {
    packnum = get("packet.number", sys.frame(3))
    switch(side,
           left = {
              #This is  a super-hack but it's the only way I have figured out to do it.
             #It is likely that this will fail if we ever try to do a conditional form of these plots
             
             if(!is.null(packnum))
               {
                     if(packnum == 3)
                       {
                         args = list(...)
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
