

#############################################################
#####         low level utility routines       ##############
#############################################################
  
namevalue <- function(nam, val) {
  paste(nam, "=\"", val, "\"", sep="") 
  }

numbered.id <- function(before, g, after) {
  paste("id=\"", before, g, after, "\"", sep="") 
  }
  
id <- function(id) { 
  paste("id=\"", id, "\"", sep="")
  }
  
  
idref <- function(element, id) {
  paste("<", element, " idref=\"", id, "\" />", sep="") 
  }  

numbered.idref <- function(element, before, g, after) {
  paste("<", element, " idref=\"", before, g, after, "\" />", sep="") 
  }  
  

catln <- function(indent, ...) {
  if (indent == 0) {
    cat(..., "\n", sep="")
    } else {
    tabs <- paste(rep("\t", indent), sep="", collapse="")
    cat(tabs, ..., "\n") 
    }
  }
 
commentln <- function(txt) {
cat("<!-- ", txt, " -->\n"); 
}

  
#############################################################
##########   functions for real data           ##############
#############################################################
 

noftetraploids.fromTable <- function(taxatable) {
noftetraploids <- 0
spp <- unique(taxatable[,"species"])
for (sp in 1:length(spp)) {
  wsp <- which(taxatable[,"species"] == spp[sp])
  ABs <- taxatable[wsp,"genome"]
  if (length(grep("B", ABs)) > 0) {
    noftetraploids <- noftetraploids + 1
    }
  }
noftetraploids
}


nofdiploids.fromTable <- function(taxatable) {
nofdiploids <- 0
spp <- unique(taxatable[,"species"])
for (sp in 1:length(spp)) {
  wsp <- which(taxatable[,"species"] == spp[sp])
  ABs <- taxatable[wsp,"genome"]
  if (length(grep("B", ABs)) == 0) {
    nofdiploids <- nofdiploids + 1
    }
  }
nofdiploids
}


make.taxa.and.alignments.fromRealData <- function(taxatable, alignments, alignmentnames) {
commentln("Taxa. Each taxon is a sequence (A or B) from an individual")
commentln("in a species.")
  
ntaxaIds <- nrow(taxatable)
catln(1, "<taxa id=\"taxa\">")
for (sq in 1:ntaxaIds) {
  catln(2, "<taxon",  id(taxatable[sq,"ID"]), "/>")
  }
catln(1, "</taxa>")

nofalignments <- length(alignments)
for (g in 1:nofalignments) {
  commentln(paste("Alignment for gene", g, "from nex file", alignmentnames[g]))
  catln(1, "<alignment", id(paste("alignment", g, sep="")), "dataType=\"nucleotide\">")
  algn <- alignments[[g]]
  seqlen <- length(algn[[1]])
  anames <-  names(algn)
  for (sq in 1:ntaxaIds) {
    catln(2, "<sequence>")
    ID <- taxatable[sq,"ID"]
    catln(3, idref("taxon", ID))
    gotseq <- FALSE
    for (n in 1:length(anames)) {
      if (anames[n] == ID) {
        cat("\t\t\t", toupper(algn[[n]]), "\n", sep="")
         gotseq <- TRUE
        }
      }
    if (!gotseq) {
      cat("\t\t\t", rep("?", seqlen), "\n", sep="")
      }
    catln(2, "</sequence>")
    }
  catln(1, "</alignment>")
  }
}


 

make.xml.realdata <- function(taxatable, alignments, beastAlloppDTxmlinfo) {
           
dpath <- beastAlloppDTxmlinfo$data.dpath           
nofalignments <- length(beastAlloppDTxmlinfo$alignmentnames)
sink(paste(dpath, beastAlloppDTxmlinfo$beastXMLfilename, sep=""))
make.beastxml.start("Made by R code AlloppDT, for Allopolyploid model for diploids and tetraploids", 
                     paste("Input data from", dpath))
make.taxa.and.alignments.fromRealData(taxatable, alignments, beastAlloppDTxmlinfo$alignmentnames)
commentln("Patterns")
for (g in 1:nofalignments) {
  make.patterns(g)
  }

nofdiploids <- nofdiploids.fromTable(taxatable)
noftetraploids <- noftetraploids.fromTable(taxatable) 
totalnseqs <- nofalignments * (nofdiploids+noftetraploids)
# this totalnseqs is a estimate. It is used only for operator weights  
make.beastxml.after.patterns(nofdiploids, noftetraploids, nofalignments, totalnseqs, 
                  taxastruct=NULL, taxatable,
                  popmean.eventrate.relativerates.prior=c(-9,1, 4.6,2, 0,1),
                  beastAlloppDTxmlinfo$sampledgtrees.fpathbase,
                  beastAlloppDTxmlinfo$sampledmultrees.fpath,
                  beastAlloppDTxmlinfo$sampledparams.fpath,
                  beastAlloppDTxmlinfo$DBUGTUNE.fpath)
sink(NULL)
}




make.beastxml.AlloppDT.forRealData <- function(beastAlloppDTxmlinfo) {
dpath <- beastAlloppDTxmlinfo$data.dpath
nofalignments <- length(beastAlloppDTxmlinfo$alignmentnames)
alignments <- NULL
for (a in 1:nofalignments) {
  nexfpath <- paste(dpath, alignmentnames[a], ".nex", sep="")
  aln <- read.nexus.data(nexfpath)
  alignments <- c(alignments, list(aln))
  }  
taxatable <- read.table(paste(dpath, beastAlloppDTxmlinfo$fpath.taxatable, sep=""), header=TRUE, colClasses="character")
stopifnot(length(unique(taxatable[,"ID"])) == length(taxatable[,"ID"]))
ABs <- sort(unique(taxatable[,"genome"]))
stopifnot(length(ABs) == 2)
stopifnot(ABs[1] == "A"  &&  ABs[2] == "B")
make.xml.realdata(taxatable, alignments, beastAlloppDTxmlinfo)
}




 
  
  
#############################################################
##########           models and priors         ##############
#############################################################
  

make.patterns <- function(g) {
  catln(1, "<patterns", numbered.id("", g, ".patterns"), namevalue("from", 1), ">")
  catln(2, numbered.idref("alignment", "alignment", g, ""))
  catln(1, "</patterns>")
  }
  
make.constantSize <- function() { 
 	catln(1, "<constantSize",  id("constant"), namevalue("units", "substitutions"), ">")
  catln(2, "<populationSize>")
  catln(3, "<parameter",  id("constant.popSize"),  namevalue("value", 0.016),
                 namevalue("lower", 0.0), namevalue("upper", "Infinity"), "/>")
  catln(2, "</populationSize>")
  catln(1, "</constantSize>")
 }
  
make.coalescentTree <- function(g, roothgt) {
  catln(1, "<coalescentTree ", numbered.id("", g, ".startingTree"),
          namevalue("rootHeight", roothgt), ">")
  catln(2, idref("taxa", "taxa"))
  catln(2, idref("constantSize", "constant"))
  catln(1, "</coalescentTree>")
  }


make.treeModel <- function(g) {
  catln(1, "<treeModel", numbered.id("", g, ".treeModel"), ">")
  catln(2, numbered.idref("coalescentTree", "", g, ".startingTree"))
  
  catln(2, "<rootHeight>")
  catln(3, "<parameter", numbered.id("", g, ".treeModel.rootHeight"), "/>")
  catln(2, "</rootHeight>")
  
  catln(2, "<nodeHeights internalNodes=\"true\" >")
  catln(3, "<parameter", numbered.id("", g, ".treeModel.internalNodeHeights"), "/>")
  catln(2, "</nodeHeights>")

  catln(2, "<nodeHeights internalNodes=\"true\" rootNode=\"true\" >")
  catln(3, "<parameter", numbered.id("", g, ".treeModel.allInternalNodeHeights"), "/>")
  catln(2, "</nodeHeights>")
  
  catln(1, "</treeModel>")
  }



make.strictClockBranchRate <- function(g) {
  catln(1, "<strictClockBranchRates", numbered.id("", g, ".branchRates"), ">")
  catln(2, "<rate>")
  catln(3, "<parameter", numbered.id("", g, ".clock.rate"), 
           namevalue("value", 1.0), namevalue("lower", 0.0),
           namevalue("upper", "Infinity"), "/>") 
  catln(2, "</rate>")
  catln(1, "</strictClockBranchRates>")
  }


  
make.HKYModel <- function(g) {
  catln(1, "<HKYModel", numbered.id("", g, ".hky"), ">")
  catln(2, "<frequencies>")
  catln(3, "<frequencyModel", namevalue("dataType", "nucleotide"), ">")
  catln(4, "<frequencies>")
  catln(5, "<parameter", numbered.id("", g, ".frequencies"), 
                 namevalue("value", "0.25 0.25 0.25 0.25"), "/>")
  catln(4, "</frequencies>")
  catln(3, "</frequencyModel>")
  catln(2, "</frequencies>")
  catln(2, "<kappa>")
  catln(3, "<parameter", numbered.id("", g, ".kappa"), namevalue("value", 2.0),
              namevalue("lower", 0), namevalue("upper", "Infinity"), "/>")
  catln(2, "</kappa>")
  catln(1, "</HKYModel>")
  }


make.siteModel <- function(g) {
  catln(1, "<siteModel", numbered.id("", g, ".siteModel"), ">")
  catln(2, "<substitutionModel>")
  catln(3, numbered.idref("HKYModel", "", g, ".hky"))
  catln(2, "</substitutionModel>")
  catln(1, "</siteModel>")
  }


make.treeLikelihood <- function(g) {
  catln(1, "<treeLikelihood", numbered.id("", g, ".treeLikelihood"), 
                namevalue("useAmbiguities", "false"), ">")
  catln(2, numbered.idref("patterns", "", g, ".patterns"))
  catln(2, numbered.idref("treeModel", "", g, ".treeModel"))
  catln(2, numbered.idref("siteModel", "", g, ".siteModel"))
  catln(2, numbered.idref("strictClockBranchRates", "", g, ".branchRates"))
  catln(1, "</treeLikelihood>")
  }




make.species.list.fromStruct <- function(taxastruct) {
for (sp in 1:length(taxastruct)) {
  spname <- taxastruct[[sp]]$sp
  ploidylevel <- 2 * length(taxastruct[[sp]]$sqs)
  catln(2, "<apsp", id(spname), namevalue("ploidylevel", ploidylevel), ">")
  for (iv in 1:length(taxastruct[[sp]]$indivs)) {
    ivname <- paste(taxastruct[[sp]]$indivs[iv], spname, sep="")
    catln(3, "<individual", id(ivname), ">")
    for (sq in 1:length(taxastruct[[sp]]$sqs)) {
      txname <- paste(ivname, taxastruct[[sp]]$sqs[sq], sep="")
      catln(4, idref("taxon", txname)) 
      }
    catln(3, "</individual>")
    }
  catln(2, "</apsp>")
  }
}
  



make.species.list.fromTable <- function(taxatable) {
spp <- unique(taxatable[,"species"])
for (sp in 1:length(spp)) {
  wsp <- which(taxatable[,"species"] == spp[sp])
  ploidylevel <- 2 * length(unique(taxatable[wsp,"genome"]))
  stopifnot(ploidylevel==2  ||  ploidylevel==4)
  catln(2, "<apsp", id(spp[sp]), namevalue("ploidylevel", ploidylevel), ">")
  indivs <- unique(taxatable[wsp,"individual"])
  for (iv in 1:length(indivs)) {
    wiv <- which((taxatable[,"species"] == spp[sp])  &  (taxatable[,"individual"] == indivs[iv]))
    stopifnot(2 * length(wiv) == ploidylevel)
    catln(3, "<individual", id(paste(spp[sp], indivs[iv], sep="")), ">")
    for (sq in 1:length(wiv)) {
      catln(4, idref("taxon", taxatable[wiv[sq],"ID"])) 
      }
    catln(3, "</individual>")
    }
  catln(2, "</apsp>")
  }
}

  
make.species.list <- function(taxastruct, taxatable) {  
if (is.null(taxastruct)) {
  make.species.list.fromTable(taxatable)        
  } else {
  make.species.list.fromStruct(taxastruct)        
  }
}

  
make.genetree.list <- function(nofgenes) {
catln(2, "<geneTrees", id("geneTrees"), ">")
for (g in 1:nofgenes) {
  catln(3, numbered.idref("treeModel", "", g, ".treeModel"))
  }
catln(2, "</geneTrees>")
}


make.alloppSpecies <- function(nofgenes, taxastruct, taxatable, mingenenodeheight) {
  catln(1, "<alloppSpecies", id("alloppSpecies"), 
               namevalue("minGeneNodeHeight", mingenenodeheight), ">")
  make.species.list(taxastruct, taxatable)
  make.genetree.list(nofgenes)
  catln(1, "</alloppSpecies>")
  } 
  
  

make.alloppSpeciesNetwork <- function(initpopval) { 
  catln(1, "<alloppSpeciesNetwork", id("apspnetwork"), "oneHybridization=\"false\" diploidRootIsRoot=\"true\" >")
  catln(2, idref("alloppSpecies", "alloppSpecies")) 
  catln(2, "<tipPopulations", namevalue("value", initpopval), ">")
  catln(3, "<parameter", id("apspnetwork.tipPopSizes"), "/>")
  catln(2, "</tipPopulations>") 
  catln(2, "<rootPopulations", namevalue("value", initpopval), ">")
  catln(3, "<parameter", id("apspnetwork.rootPopSizes"), "/>")
  catln(2, "</rootPopulations>")
  commentln("hybridPopulations is not a normal parameter, since it changes dimension.")
  commentln("It can be logged (unused values are set to 0). It needs special operators.")
  catln(2, "<hybridPopulations", namevalue("value", initpopval), ">")
  catln(3, "<parameter", id("apspnetwork.hybridPopSizes"), "/>")
  catln(2, "</hybridPopulations>")
  catln(1, "</alloppSpeciesNetwork>")
  }



    
make.alloppNumHybsStatistic <- function() {
  catln(1, "<alloppNumHybsStatistic", id("allopp.numHybs"), ">") 
  catln(2, "<apspNetwork>")
  catln(3, idref("alloppSpeciesNetwork", "apspnetwork"))
  catln(2, "</apspNetwork>")
  catln(1, "</alloppNumHybsStatistic>")

}


make.alloppNetworkPriorModel <- function() { 
  catln(1, "<alloppNetworkPriorModel", id("network.prior.model"), 
                   namevalue("units", "substitutions"), ">")
  catln(2, "<eventRate>")
  catln(3, "<parameter", id("apspnetwork.prior.eventRate"), 
             namevalue("value", 100), namevalue("lower", 0.0), namevalue("upper", "Infinity"),"/>")
  catln(2, "</eventRate>")

  catln(2, "<populationScalingFactor>")
  catln(3, "<parameter", id("population.scaling.factor"), 
             namevalue("value", 0.2), namevalue("lower", 0.0), namevalue("upper", "Infinity"),"/>")
  catln(2, "</populationScalingFactor>")
  
  catln(2, "<tipPopulationDistribution>")
  catln(3, "<gammaDistributionModel>")
  catln(4, "<shape>")
  catln(5, "4")
  catln(4, "</shape>")
  catln(4, "<scale>")
  catln(5, idref("parameter", "population.scaling.factor"))
  catln(4, "</scale>")
  catln(3, "</gammaDistributionModel>")
  catln(2, "</tipPopulationDistribution>")
  
  catln(2, "<rootPopulationDistribution>")
  catln(3, "<gammaDistributionModel>")
  catln(4, "<shape>")
  catln(5, "2")
  catln(4, "</shape>")
  catln(4, "<scale>")
  catln(5, idref("parameter", "population.scaling.factor"))
  catln(4, "</scale>")
  catln(3, "</gammaDistributionModel>")
  catln(2, "</rootPopulationDistribution>")

  catln(2, "<hybridPopulationDistribution>")
  catln(3, "<gammaDistributionModel>")
  catln(4, "<shape>")
  catln(5, "1")
  catln(4, "</shape>")
  catln(4, "<scale>")
  catln(5, idref("parameter", "population.scaling.factor"))
  catln(4, "</scale>")
  catln(3, "</gammaDistributionModel>")
  catln(2, "</hybridPopulationDistribution>")
  
  catln(1, "</alloppNetworkPriorModel>")
  }


make.apspNetworkPrior <- function() { 
  catln(1, "<apspNetworkPrior", id("apspnetwork.prior"), ">")
  catln(2, "<model>")
  catln(3, idref("alloppNetworkPriorModel", "network.prior.model"))
  catln(2, "</model>")
  catln(2, "<apspNetwork>")
  catln(3, idref("alloppSpeciesNetwork", "apspnetwork"))
  catln(2, "</apspNetwork>")
  catln(1, "</apspNetworkPrior>")
  }  
  
  
  
make.apspCoalescent <- function() { 
  catln(1, "<apspCoalescent", id("apsp.coalescent"), ">")
  catln(2, idref("alloppSpecies", "alloppSpecies"))
  catln(2, idref("alloppSpeciesNetwork", "apspnetwork"))
  catln(1, "</apspCoalescent>")    
  }




make.gammaDistributionModel <- function(shape, initpopval, first) {
  catln(3, "<gammaDistributionModel>")
  catln(4, "<shape>")
  catln(5, shape)
  catln(4, "</shape>")
  catln(4, "<scale>")
  if (first) {
    catln(5, "<parameter", id("species.popMean"), namevalue("value", initpopval), 
               namevalue("lower", 0.0), namevalue("upper", "Infinity"),"/>")
    } else {
    catln(5, idref("parameter", "species.popMean"))
    }
  catln(4, "</scale>")
  catln(3, "</gammaDistributionModel>")
  }
 
                                         

#############################################################
##########            MCMC operators           ##############
#############################################################
 


make.general.numbered.scaleOperator <- function(g, sfactor, weight, before, after) {
  catln(2, "<scaleOperator", namevalue("scaleFactor", sfactor),
         namevalue("weight", weight), ">")
  catln(3, numbered.idref("parameter", before, g, after))
  catln(2, "</scaleOperator>")
  }		
		
make.general.scaleOperator <- function(sfactor, weight, param) {
  catln(2, "<scaleOperator", namevalue("scaleFactor", sfactor),
         namevalue("weight", weight), ">")
  catln(3, idref("parameter", param))
  catln(2, "</scaleOperator>")
  }			
		
		
make.scaleOperator.kappa <- function(g, weight) {
  make.general.numbered.scaleOperator(g, 0.75, weight, "", ".kappa")
  }
  

make.scaleOperator.clockrate <- function(g, weight) {
  make.general.numbered.scaleOperator(g, 0.75, weight, "", ".clock.rate")
  }  
	
make.deltaExchange.frequencies <- function(g, weight) {
  catln(2, "<deltaExchange", namevalue("delta", 0.01),
         namevalue("weight", weight), ">")
  catln(3, numbered.idref("parameter", "", g, ".frequencies"))
  catln(2, "</deltaExchange>")  
  }
  

make.subtreeSlide.treeModel <- function(g, weight) {
  catln(2, "<subtreeSlide", namevalue("size", .002), 
       namevalue("gaussian", "true"), namevalue("weight", weight), ">")
  catln(3, numbered.idref("treeModel", "", g, ".treeModel"))
  catln(2, "</subtreeSlide>") 
  }
  

  
make.narrowExchange.treeModel <- function(g, weight) {
  catln(2, "<narrowExchange", namevalue("weight", weight), ">")
  catln(3, numbered.idref("treeModel", "", g, ".treeModel"))
  catln(2, "</narrowExchange>") 
  }
  
   
make.wideExchange.treeModel <- function(g, weight) {
  catln(2, "<wideExchange", namevalue("weight", weight), ">")
  catln(3, numbered.idref("treeModel", "", g, ".treeModel"))
  catln(2, "</wideExchange>") 
  }
  
   
make.wilsonBalding.treeModel <- function(g, weight) {
  catln(2, "<wilsonBalding", namevalue("weight", weight), ">")
  catln(3, numbered.idref("treeModel", "", g, ".treeModel"))
  catln(2, "</wilsonBalding>") 
  
  }
  
  
make.scaleOperator.treeModel.rootHeight <- function(g, weight) {
  make.general.numbered.scaleOperator(g, 0.75, weight, "", ".treeModel.rootHeight")
  }
  

 
make.uniformOperator.treeModel.internalNodeHeights <- function(g, weight) {
  catln(2, "<uniformOperator", namevalue("weight", weight), ">")
  catln(3, numbered.idref("parameter", "", g, ".treeModel.internalNodeHeights"))
  catln(2, "</uniformOperator>")
  }
  
   
   
make.upDownOperator.clockrate.allInternalNodeHeights <- function(g, weight) {
  catln(2, "<upDownOperator", namevalue("scaleFactor", 0.75), 
           namevalue("weight", weight), ">")
  catln(3, "<up>")
  if (g>1) {
    catln(4, numbered.idref("parameter", "", g, ".clock.rate")) 
    }  
  catln(3, "</up>")
  catln(3, "<down>")
  catln(4, numbered.idref("parameter", "", g, ".treeModel.allInternalNodeHeights"))
  catln(3, "</down>")
  catln(2, "</upDownOperator>")
  } 
  
  
      
make.network.upDownOperator.clocks.pops.heights <- function(nofgenes, weight) {
  catln(2, "<upDownOperator", namevalue("scaleFactor", 0.75), 
           namevalue("weight", weight), ">")
           
  catln(3, "<up>")
  for (g in 1:nofgenes) {
    if (g>1) {
      catln(4, numbered.idref("parameter", "", g, ".clock.rate")) 
    }
  }
  catln(4, idref("parameter", "apspnetwork.prior.eventRate")) 
  catln(3, "</up>")
  
  catln(3, "<down>")
  catln(4, idref("alloppSpeciesNetwork", "apspnetwork"))
  catln(4, idref("parameter" , "population.scaling.factor"))

  for (g in 1:nofgenes) {
    catln(4, numbered.idref("parameter", "", g, ".treeModel.allInternalNodeHeights"))
  }
  catln(3, "</down>")
  catln(2, "</upDownOperator>")
  }
  



make.scaleOperator.species.popSF <- function(weight) {
  make.general.scaleOperator(0.9, weight, "population.scaling.factor")
  }
  
  
make.scaleOperator.apspnetwork.tipPopSizes <- function(weight) {
  make.general.scaleOperator(0.5, weight, "apspnetwork.tipPopSizes")
  }
  
make.scaleOperator.apspnetwork.rootPopSizes <- function(weight) {
  make.general.scaleOperator(0.5, weight, "apspnetwork.rootPopSizes")
  }

make.scaleOperator.apspnetwork.prior.eventRate <- function(weight) {
  make.general.scaleOperator(0.5, weight, "apspnetwork.prior.eventRate")
  }
  
make.hybPopSizesScaleOperator <- function(weight) {
  catln(2, "<hybPopSizesScaleOperator",  namevalue("scaleFactor", "0.5"), namevalue("weight", weight), ">")
  catln(3, idref("alloppSpecies", "alloppSpecies"))
  catln(3, idref("alloppSpeciesNetwork", "apspnetwork"))
  catln(2, "</hybPopSizesScaleOperator>")  
  }
    
make.networkNodeReHeightOperator <- function(weight) {
  catln(2, "<networkNodeReHeight",  namevalue("weight", weight), ">")
  catln(3, idref("alloppSpecies", "alloppSpecies"))
  catln(3, idref("alloppSpeciesNetwork", "apspnetwork"))
  catln(2, "</networkNodeReHeight>")   
  }

  
  
make.sequenceReassignmentOperator <- function(weight) {
  catln(2, "<sequenceReassignment",  namevalue("weight", weight), ">")
  catln(3, idref("alloppSpecies", "alloppSpecies"))
  catln(3, idref("alloppSpeciesNetwork", "apspnetwork"))
  catln(2, "</sequenceReassignment>")   
  }


make.changeNumHybridizationsOperator <- function(weight) {
  catln(2, "<changeNumHybridizations",  namevalue("weight", weight), ">")
  catln(3, idref("alloppSpecies", "alloppSpecies"))
  catln(3, idref("alloppSpeciesNetwork", "apspnetwork"))
  catln(2, "</changeNumHybridizations>")   
  }     

  
      

