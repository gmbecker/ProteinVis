

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

    if(!is.null(tm))
      {
        tmmat = combineTMDomains(tm, poscols = tmposcols)
        apply(tmmat, 1, function(x, ylim)
              {
                grid.polygon( c( x , rev( x ) ), rep( ylim , times = c(2, 2) ) , gp = gpar(fill = "grey75", alpha = .5, col = "white"), default.units = "native", draw = TRUE)
                grid.text("TM",unit(sum(x)/2, "native"), unit(.05, "npc"), gp = gpar(cex = .75))
              }, ylim = ylim)
      }

    if(!is.null(sig))
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
    
    panel.xyplot(x,y,...)
    TRUE
  }


panel.psipred = function(x, y, subscripts, cutoff, strand, ...)
  {
    drawPsipred(data.frame(start = x,helix = y, strand = strand), cutoff, current.panel.limits()$xlim)
  }

drawPsipred = function(dat, cutoff, xlim)
  {
    grid.text("Secondary Structure", unit(1, "char"), .75, just = "left")
    grid.lines(unit(xlim, "native"), unit(.25, "npc"), gp = gpar(col = "grey50", lex = 1.5))
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
          drawCoil(change.helix[i], change.helix[i+1], height = .35, center.y = .25, gp = gpar(fill = "red", alpha = .5), scorecol = "black")
      }
  }
    if(length(change.strand) >= 2)
      {
        
    for(j in seq(1, length(change.strand) - 1, by = 2))
      {
        if(change.strand[j+1] - change.strand[j] >= 10)
          
          drawArrow(change.strand[j], change.strand[j + 1], height = .35, center.y = .25, fill = "blue")
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

#this is a bit hacky, it expects to be passed the starting and ending of the ranges as x and y, respectively. Would it be better to pass in the end as an arg and have y be the labels?

panel.PFAM = function(x, y, end, subscripts, labs, colors, pfamBins, ...)
  {
    drawPFAM(data.frame(start = x, end = end, labels = labs, bin = y), labcol = "labels", colors = colors, bins = pfamBins)
  }

drawPFAM = function(dat, colors, poscolumns = c("start", "end"), labcol = "featureName", bins)
  {
    nrow = max(bins)
    
    grid.text("PFAM Domains", unit(1, "char"), 1-.2/nrow, just = "left")
    dat$bin = bins
    apply(dat, 1, function(x, colors, poscolumns, nrow)
          {
            
            st = as.numeric(x[poscolumns[1]])
            end = as.numeric(x[poscolumns[2]])
     
            bin = as.numeric(x["bin"])
            #determine label and color
            fulln = x[labcol]

            fullnlen = convertWidth(unit(nchar(fulln), "char"), "native", valueOnly= TRUE)
            barlen = end - st
  
            lab  = fulln
            rot = 0
#            col = x["col"]
            col = pfColMap[as.character(pfColMap[,1]) == lab, 2]
            if (!length(col))
              col = "#111111" #grey
            #.45 for one row ...
            step = .6/nrow
            ypos = 1 - (.15 / nrow ) - step*bin
            grid.lines(unit(c(st, end), "native"), unit(nrow - bin + 1, "native"), gp = gpar(col = col, lex = 10))
            grid.text(lab, unit((st + end) /2, "native"), y = unit(nrow - bin + 1 , "native") - unit(1, "char") , rot = rot)

          }, colors = colors, poscolumns = poscolumns, nrow = nrow)

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

panel.metaCount = function(x, y, subscripts, patientid, scale.factor = 8, logscale = FALSE, logbase = exp(1), at.baseline = FALSE,colpalette = rev(brewer.pal(11, "RdBu")), legend.step = .005, ...)
  {
    if(sum(!is.na(x)) == 0)
      return(TRUE)
    library(gridSVG)
    patientid = patientid[subscripts]
    grid.rect(.5, .5, 1, 1, gp = gpar(fill = "grey95"))
    
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
             counts = myl$counts
             p = myl$proportions
             y = myl$y
             n = myl$samplesize
                                        #browser()
             if(logscale)
               heights = sapply(counts, function(x) .05 + .9 * min( log( x, base = logbase), scale.factor ) / scale.factor)
             else
               heights = sapply(counts, function(x) .9 * min( x, scale.factor ) / scale.factor )
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
             grid.rect(as.numeric(names(counts)), ypos, vjust = vjust, default.units = "native", 5, heights, gp = gpar(fill = colpalette[colinds]), )
           })
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
           par.settings = list(layout.heights = list(panel = c(.3 + .25*max(pfamBins), .4, 1.5, 2)))
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

metaCountStructPlot = function(events,catname = "PRIMARY_TISSUE", position = c("protpos", "protposend"),  pfam, pfamLabels = "featureName",structPred, hydro, transMem, sigP, xlim, tmposcol = c("start", "end"), main = NULL, simple = FALSE, at.baseline = TRUE, logscale = TRUE, logbase = 1.5, scale.factor = 10, colpalette = rev(brewer.pal(11, "RdYlBu")), legend.step = .01, sampleID, key , subtitle, draw = FALSE)
  {
    pfamIRange = IRanges(start = pfam$start, end = pfam$end, names = pfam[,pfamLabels])
    pfamBins = disjointBins(pfamIRange)
    pfam$bin = factor(pfamBins, levels = 0:(max(pfamBins + 1)))
    plots = makeStructPlots(pfam, structPred, hydro, transMem, sigP, xlim, tmposcol, pfamLabels = pfamLabels, pfamBins = pfamBins)
    
    hydroPlot = plots$hydro

    structPredPlot = plots$structPred
    pfamPlot = plots$pfam

    countPlot = xyplot(as.formula(paste(catname, "~", position)), data = events, panel = panel.metaCount,  patientid = sampleID, at.baseline = at.baseline, logscale = logscale, scale.factor = scale.factor, logbase = logbase, colpalette = colpalette, legend.step = legend.step)

    combPlot = c(pfamPlot, structPredPlot, hydroPlot, countPlot, x.same=TRUE, y.same=NA, merge.legends=TRUE)
    
   combPlot = update(combPlot,
       layout = c(1, 4),
       ylab = list(c( "Hydrophobicity", "", ""),
       y = c( .35, NA, NA)),
       xlim = xlim,
       par.settings = list(layout.heights = list(panel = c(.3 + .25*max(pfamBins), .4, 1.5, 2))),
      #    par.settings = list(layout.heights = list(panel = c(unit(3, "cm"), unit(2, "cm"), unit(8, "cm") , unit(10, "cm")))),
       main = main, sub = subtitle )
   # c(combPlot, draw.key(key), layout = c(2, 1))
    if(draw)
      print(combPlot)
    else
      combPlot

  }

axis.combined = function(side, ...)
  {
     switch(side,
           left = {
              #This is  a super-hack but it's the only way I have figured out to do it.
             #It is likely that this will fail if we ever try to do a conditional form of these plots
             packnum = get("packet.number", sys.frame(3))
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
            }
            )
   }



