# Create expDesign data frame ---------------------------------------------

path = system.file("extdata",  package = "Rpackages")
files = list.files(path, full.names = T)

expDesign = data.frame(files = files, name = as.character(sapply(files, function(samp)  unlist(strsplit(unlist(strsplit(samp, split = ".", fixed = T))[2], split = "/"))[4])))

data("ATseqlengths")


myDat = DataSetImportFromtxt(experimentDesign = expDesign,
                             bin.length = 200,
                             min.samples = 2,
                             quant = .2,
                             min.counts = 5,
                             seqlengths = ATseqlengths)


# Check samples quality ---------------------------------------------------

myAnalysis = analysis.samples(myDat)

myAnalysis$lib.size

barplot(myAnalysis$lib.size, col = (!myAnalysis$passThresh) +1, las = 2)

