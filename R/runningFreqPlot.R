running.freq = function(x, tiles = 4e5, info){

  listgr = vector("list", info$part_nr)
  names(listgr) = names(x[[1]])
  for(chr in names(listgr)){
    myt = tileGenome( tilewidth = tiles,
                      cut.last.tile.in.chrom=TRUE, seqlengths = seqlengths(x[[1]][[chr]])[chr])
    for(samp in names(x)){
      hits = findOverlaps(x[[samp]][[chr]], myt)
      tp = table(hits@to)
      mcols(myt)[,samp ] = 0
      mcols(myt)[[samp]][as.numeric(names(tp))] = tp
    }
    myt$co.freq = apply(mcols(myt), 1, function(x) sum(x)/(tiles/1e3))
    listgr[[chr]] = myt
  }
  return(listgr)
}

# tiles = 4e5
plotFreqgen = function(myx, tiles, file, groups = NULL, info, verbose = TRUE){
  if(sum(c("Gviz") %in% rownames(installed.packages())) != 2) stop("To generate this plot you need to have installed Gviz Bioconductor package.\n
                                                                   https://bioconductor.org/packages/release/bioc/html/Gviz.html")

  if(verbose) cat("Plotting Gen frequencies.\n")
  total.running = running.freq(myx, tiles = tiles, info = info)
  if(is.null(groups)){
    pdf(file)
    for(chr in info$part_names){
      # cat("Inside chromosome ", chr, "\n")
      FinalRes_chr = granges(total.running[[chr]])
      FinalRes_chr$co.freq = total.running[[chr]]$co.freq
      # FinalRes_chr$wt = wt.running[[chr]]$co.freq

      dat = FinalRes_chr[,c("co.freq")]

      datgviz = Gviz::DataTrack(dat, type = "l", name = "COs per MB", col = c("#fc8d62"))
      axistrack = Gviz::GenomeAxisTrack()

      mySize = c(.3,1)
      # png(paste("~/Documents/PhDthesis/Chapter4/Figures/runningCo/MutvsWt-", chr, ".png", sep=""))
      Gviz::plotTracks(c(axistrack,datgviz),  groups = colnames(mcols(dat)),
                 cex.feature=0.7, background.title="darkgrey", lwd=2,
                 from=1,
                 to= seqlengths(dat)[chr] ,
                 # to = 1733000,
                 main =chr,
                 sizes=mySize,
                 showFeatureId=FALSE,
                 fontcolor.feature="black", cex.feature=0.7, background.title="darkgrey",
                 showId=TRUE)
    }
    # cat("Finishing loop. \n")
    dev.off()
    # cat("Closing the file.\n")
  } else{
    cat("I enter else\n")
    if(!is.factor(groups)) groups = as.factor(groups)
    dat = lapply(levels(groups), function(group){
      myg = total.running[groups == group]
    })
    for(chr in info$part_names){
      FinalRes_chr = granges(total.running[[chr]])
      for(group in levels(groups)){
        FinalRes_chr[[group]] = dat[[group]][[chr]]$co.freq
      }

      # FinalRes_chr$wt = wt.running[[chr]]$co.freq

      # dat = FinalRes_chr[,c("co.freq")]
      colors = c("#fc8d62", "#8da0cb","#b35806","#e08214","#fdb863","#fee0b6","#f7f7f7","#d8daeb","#b2abd2","#8073ac","#542788")

      datgviz = Gviz::DataTrack(dat, type = "l", name = "CO freq", col = colors[1:length(levels(groups))])
      axistrack = Gviz::GenomeAxisTrack()
      pdf(file)

      mySize = c(.3,1)
      # png(paste("~/Documents/PhDthesis/Chapter4/Figures/runningCo/MutvsWt-", chr, ".png", sep=""))
      Gviz::plotTracks(c(axistrack,datgviz),  groups = colnames(mcols(dat)),
                 cex.feature=0.7, background.title="darkgrey", lwd=2,
                 from=1,
                 to= seqlengths(dat)[chr] ,
                 # to = 1733000,
                 main =chr,
                 sizes=mySize,
                 showFeatureId=FALSE,
                 fontcolor.feature="black", cex.feature=0.7, background.title="darkgrey",
                 showId=TRUE)
    }
    dev.off()
  }
  # cat("IS her the proble?")
  #
  # cat("Yes it fucking is!")
}

