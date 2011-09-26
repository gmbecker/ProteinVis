

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

panel.protstruct <- function(x,y, subscripts, tmposcols, tm, sig,... )
  {

    ylim = current.panel.limits()$ylim
    xlim = current.panel.limits()$xlim

    if(!is.null(tm) && nrow(tm) > 0 )
      {
        tmmat = combineTMDomains(tm, poscols = tmposcols)
        apply(tmmat, 1, function(x, ylim)
              {
                grid.polygon( c( x , rev( x ) ), rep( ylim , times = c(2, 2) ) , gp = gpar(fill = "grey75", alpha = .5, col = "white"), default.units = "native", draw = TRUE)
                grid.text("TM",unit(sum(x)/2, "native"), unit(.05, "npc"), gp = gpar(cex = .75))
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
                grid.text("SigP",unit(sum(sigx)/2, "native"), unit(.05, "npc"), gp = gpar(cex = .75))
              }, ylim = ylim)
      }
    #comes after tm/sigp stuff to deal with alpha issues
    grid.text("Hydrophobicity", unit(2, "mm"), unit(1, "npc") - unit(2, "mm"), just = c("left", "top"), gp = gpar(fontface="bold", cex=.9 ))
    #This may get rid of the problem with hydrophobicity not plotting correctly. Not sure, I need to get the data from cory.
    if(length(x))
      panel.xyplot(x,y,...)
    TRUE
  }


panel.psipred = function(x, y, subscripts, cutoff, strand, ...)
  {
    drawPsipred(data.frame(start = x,helix = y, strand = strand), cutoff, current.panel.limits()$xlim)
  }

drawPsipred = function(dat, cutoff, xlim)
  {
    #grid.text("Secondary Structure", unit(1, "char"), .75, just = "left")
    grid.text("Secondary Structure", unit(2, "mm"),unit(1, "npc") -  unit(2, "mm"), just = c("left", "top"), gp = gpar(fontface="bold", cex=.9 ))

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


panel.PFAM = function(x, y, end, subscripts, labs,  pfamBins, ...)
  {
    drawPFAM(data.frame(start = x, end = end, labels = labs, bin = y), labcol = "labels", bins = pfamBins, xlim = current.panel.limits()$xlim)
  }

drawPFAM = function(dat, poscolumns = c("start", "end"), labcol = "featureName", bins, xlim)
  {
    nrow = max(bins)
    
    grid.text("PFAM Domains", unit(2, "mm"), unit(1, "npc") - unit(2, "mm"), just = c("left", "top"), gp = gpar(fontface="bold", cex=.9 ))
    dat$bin = bins
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
              } 
            
            #.45 for one row ...
            step = .6/nrow
            ypos = 1 - (.15 / nrow ) - step*bin
            #grid.lines(unit(c(st, end), "native"), unit(nrow - bin + 1, "native"), gp = gpar(col = col, lex = 20, alpha=.5))
            grid.rect(unit(st, "native"), unit(nrow - bin + 1, "native"), width = unit(end - st, "native"), height = unit(15, "points"), just = c("left", "center"), gp = gpar(fill = col, alpha = .5)) 
            grid.text(lab, unit((st + end) /2, "native"), y = unit(nrow - bin + 1 , "native") - unit(1.2, "char") , rot = rot)

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

panel.metaCount = function(x, y, end, subscripts, patientid, scale.factor = 8, logscale = FALSE, logbase = exp(1), at.baseline = FALSE,colpalette = rev(brewer.pal(11, "RdBu")), legend.step = .005, levels = levels(y), ...)
  {
    if(sum(!is.na(x)) == 0)
      return(TRUE)

    patientid = patientid[subscripts]
    grid.rect(.5, .5, 1, 1, gp = gpar(fill = "grey95"))
 
    end = end[subscripts]
    
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
    
    scaleseq = seq(min(as.integer(ylevs)) - .5, max(as.integer(ylevs)) + .5, by = 1/2)
    grid.segments(unit(0, "npc"), scaleseq, unit(1, "npc"), scaleseq, default.units = "native", gp = gpar(col = c("black", "grey50"), lex = c(1.5, 1)))
    ylim = current.panel.limits()$ylim
    gray.yrange = .5 - ylim[1]
    grid.rect(unit(.5, "npc"), unit(.5 - gray.yrange/2, "native"), unit(1, "npc"), unit(gray.yrange, "native"), gp = gpar(fill = "grey35"))
    thing = lapply(ylevs, function(ylev, x, y, patid)
      {
        inds = which(y == ylev)
        
        x = x[inds]
        patid = patid[inds]
        sampsize = length(unique(patid))
        counts = sort(table(x), decreasing = TRUE)
        props = counts/sampsize
        list(samplesize = sampsize, proportions = props, x = x, y  = ylev, counts = counts)
      }, x = x, y = y, patid = patientid)

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
                   heights = sapply(counts, function(x) .05 + .9 * min( log( x, base = logbase), scale.factor ) / ( 2 * scale.factor))
                 else
                   heights = sapply(counts, function(x) .9 * min( x, scale.factor ) / (2 * scale.factor ))
                 if (at.baseline)
                   {
                     ypos = as.integer(y) - .5
                     vjust = 0
                   }
                 else
                   {
                     ypos = y
                     vjust = NULL
                   }
                 colinds = sapply(p, function(x) min(ceiling(x / legend.step), 11))

             
                 grid.rect(as.numeric(names(counts)), y, vjust = vjust, default.units = "native", 5, heights, gp = gpar(fill = colpalette[colinds]), )
               } else
             grid.text("NO DATA", unit(.5, "npc"), unit(as.integer(myl$y) , "native"))
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
                       .05 + .9 * min( log( x, base = logbase), scale.factor ) / ( 2 * scale.factor)
                 else
                   .9 * min( x, scale.factor ) / ( 2 * scale.factor)
                   }
               })
              grid.polygon(x = unit( c(x.s, rev(x.s) ), "native"), y = unit(c( y -  heights, rep( y, times = length(x.s) ) ), "native"), gp = gpar(stroke=NULL, fill="#00AA00", alpha=.5) )
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

makeStructPlots = function(pfam, structPred, hydro, transMem, sigP, xlim, tmposcol = c("start", "end"), main = NULL, pfamLabels = "featureName", cutoff = if(max(structPred$helix > 1)) 6 else .6, pfamBins)
  {

    hydroPlot = xyplot(featureValue ~ start, type = "l", col = "black", data = hydro, tm = transMem, sig = sigP, xlim = xlim, panel = panel.protstruct, axis = axis.protstruct, tmposcol = tmposcol)

    structPredPlot = xyplot(helix ~ start, data = structPred, strand = structPred$strand, cutoff= cutoff, panel= panel.psipred, ylim = c(0, 1), ylab = NULL, xlab = NULL)

    pfamPlot = xyplot(bin~start, end = pfam$end,  data= pfam, labs = pfam[,pfamLabels], panel = panel.PFAM, ylab = NULL, xlab = "Amino Acid Position", axis = axis.combined, pfamBins = pfamBins)
    list(hydro = hydroPlot, structPred = structPredPlot, pfam = pfamPlot)
  }

metaCountStructPlot = function(events,catname = "PRIMARY_TISSUE", requiredCats = NULL, position = c("protpos", "protposend"),  pfam, pfamLabels = "featureName",structPred, hydro, transMem, sigP, xlim, tmposcol = c("start", "end"), main = NULL, simple = FALSE, at.baseline = TRUE, logscale = TRUE, logbase = 1.5, scale.factor = 10, colpalette = rev(brewer.pal(11, "RdYlBu")), legend.step = .01, sampleID, key , subtitle, draw = FALSE)
  {
    pfamIRange = IRanges(start = pfam$start, end = pfam$end, names = pfam[,pfamLabels])
    pfamBins = disjointBins(pfamIRange)
    pfam$bin = factor(pfamBins, levels = 0:(max(pfamBins + 1)))
    plots = makeStructPlots(pfam, structPred, hydro, transMem, sigP, xlim, tmposcol, pfamLabels = pfamLabels, pfamBins = pfamBins)
    
    hydroPlot = plots$hydro

    structPredPlot = plots$structPred
    pfamPlot = plots$pfam

    if(!is.null(requiredCats))
      {
        #add fake empty observations for each required category
        colnm = names(events)
        colcl = sapply(events, mode)
        emptyrow = rep(NA, times = length(colnm))
        names(emptyrow) = colnm
        for(i in seq(along = requiredCats))
          {
            emptyrow[catname] = requiredCats[i]
            events = rbind(emptyrow, events)
          }
        #check if this is all of them
        obscats = unique(as.character(events[[catname]]))
        additcats = obscats[which( !( obscats %in% requiredCats ) )]
        cats = c(requiredCats, additcats)

        #reorder factors class comes out as c("ordered", "factor")
        events[[catname]] = factor(as.character(events[[catname]]), levels = cats)
        for(i in seq(along = colcl))
          {
            if(!any( class( events[[ colnm[ i ] ]] ) == "factor" ) )
              mode(events[[colnm[i]]] ) = colcl[i]
          }
      }

        #deal with missing categories (even though there really shouldn't be any!!!!)
    y = events[[catname]]
    levs = levels(y)
    missingCat = which(is.na(y))
    tmpcat = as.character(y)
    tmpcat[missingCat] = "UnCategorized"
    events[[catname]] = factor(tmpcat, levels = c("UnCategorized", levels(events[[catname]])))

    #sampleID is passed in as the name of a column, but we need to have the actual data.
    if(is(sampleID, "character") & length(sampleID) == 1)
      sampleID = events[,sampleID]
    
    countPlot = xyplot(as.formula(paste(catname, "~", position)), end = events$end, data = events, panel = panel.metaCount,  patientid = sampleID, at.baseline = at.baseline, logscale = logscale, scale.factor = scale.factor, logbase = logbase, colpalette = colpalette, legend.step = legend.step)

    combPlot = c(pfamPlot, structPredPlot, hydroPlot, countPlot, x.same=TRUE, y.same=NA, merge.legends=TRUE)
    
   combPlot = update(combPlot,
       layout = c(1, 4),
       ylab = list(c( "Hydrophobicity", "", ""),
       y = c( .35, NA, NA)),
       xlim = xlim,
       par.settings = list(layout.heights = list(panel = c(.3 + .25*max(pfamBins), .55, 1.5, 2)), layout.widths = list(right.padding = 10)),
     
      #    par.settings = list(layout.heights = list(panel = c(unit(3, "cm"), unit(2, "cm"), unit(8, "cm") , unit(10, "cm")))),
       main = main, sub = subtitle )
   # c(combPlot, draw.key(key), layout = c(2, 1))

    if(!is.null(key))
      {
  
        
      }
        
    if(draw)
      print(combPlot)
    else
      combPlot

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
                 if(packnum > 2)
                   {
                     if(packnum == 4)
                       {
                         args = list(...)
                         args$rot=0;
                         do.call("axis.default", c(side, args))
                       }
                     else
                       axis.default(side, ...)
                   }
               }
           },
            bottom = {
              axis.default(side , ...)
            },
            right = {
              if(packnum == 4)
                {
                  #even more superhacky!
                  trell = get("x", sys.frame(2))
                  myargs = trell$panel.args[[4]] # 4 is for panel.metacount
                  ys = myargs$y
                  xs = myargs$x
                  yinds = as.integer(ys[ seq( which( !is.na( xs ) )[1] , length(ys) ) ])
                  yvals = levels(ys)[yinds]
                  patid = myargs$patientid
                  mutcounts = by( yinds , yinds , length)
                  patcounts = as.numeric(by( patid , yinds , function(x) length(unique(x))))
                  panel.axis(side = side, outside = TRUE,
                             at = as.numeric(names(mutcounts)), labels = paste(as.numeric(mutcounts) , patcounts, sep = " / ")) 
                }
            }
            )
   }


