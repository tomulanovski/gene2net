
#################### This is 'exported' ########################################

make.beastxml.from.sim.network <- function(network, nofgenes, taxa, popmean.eventrate.relativerates.prior, 
                   seqs.fpathbase, sampledgtrees.fpathbase,
                   sampledmultrees.fpath, sampledparams.fpath, DBUGTUNE.fpath) {
  nofdiploids <- network$nofdips
  noftetraploids <- network$noftets
  taxastruct <- taxa.array.to.structure(taxa)

  make.beastxml.start("Simulated data for allopolyploids made by seqgen and R code.",
                      "Model for arbitrary numbers of diploids, tetraploids, hybridizations.")
  
  totalnofseqs <- make.taxa.and.alignments(nofgenes, seqs.fpathbase)
  commentln("Patterns")
  for (g in 1:nofgenes) {
    make.patterns(g)
    }

  make.beastxml.after.patterns(nofdiploids, noftetraploids, nofgenes, totalnofseqs,
                   taxastruct, taxatable=NULL, popmean.eventrate.relativerates.prior,
                   sampledgtrees.fpathbase, sampledmultrees.fpath, sampledparams.fpath,
                   DBUGTUNE.fpath)
}
  
  
  
################################################################################
####################### rest is 'private' ###################################### 



make.beastxml.start <- function(comment1, comment2) {
  catln(0, "<?xml version=\"1.0\" standalone=\"yes\"?>")
  commentln(comment1)
  commentln(comment2)
  catln(0, "<beast>")
}




make.beastxml.after.patterns <- function(nofdiploids, noftetraploids, nofgenes, totalnofseqs,
                   taxastruct, taxatable, popmean.eventrate.relativerates.prior,
                   sampledgtrees.fpathbase, sampledmultrees.fpath, sampledparams.fpath,
                   DBUGTUNE.fpath) {
                   
  stopifnot( !(is.null(taxastruct)  &&  is.null(taxatable)) )
  stopifnot( !(!is.null(taxastruct)  &&  !is.null(taxatable)) )   
  
  commentln("Constant size population parameter used to generate starting trees.")
  make.constantSize()
  commentln("Randomly generated starting trees for genes.")
  for (g in 1:nofgenes) {
    make.coalescentTree(g, 0.2)
    }
  commentln("Tree models for gene trees.")
  for (g in 1:nofgenes) {
    make.treeModel(g)
    }
  commentln("Branch models for gene trees.")
  for (g in 1:nofgenes) {
    make.strictClockBranchRate(g)
    }
  commentln("HKY substitution models for gene trees.")
  for (g in 1:nofgenes) {
    make.HKYModel(g)
    }  
  commentln("Site models for gene trees.")
  for (g in 1:nofgenes) {
    make.siteModel(g)
    }   
  commentln("Likelihoods for gene trees.")
  for (g in 1:nofgenes) {
    make.treeLikelihood(g)
    } 
    
  commentln("Species network model.")    
  initpopval <- 0.2
  commentln("Assignments of sequences to individuals and individuals to species.") 
  
  make.alloppSpecies(nofgenes, taxastruct, taxatable, 0.2)
  commentln("Species network, like species tree in *BEAST")    
  make.alloppSpeciesNetwork(initpopval)
  commentln("Statistic for the number of hybridizations")    
  make.alloppNumHybsStatistic()
  commentln("Prior model for species network (Yule-like model with one parameter, population prior)")
  make.alloppNetworkPriorModel()
  commentln("Prior probability for species network")
  make.apspNetworkPrior()
  commentln("Coalescent for species network: probability that gene trees fit")
  commentln("into the network.")
  make.apspCoalescent()

  make.operators(nofgenes, totalnofseqs)
  
  make.mcmc(nofgenes, popmean.eventrate.relativerates.prior,
          sampledgtrees.fpathbase, sampledmultrees.fpath, 
          sampledparams.fpath, DBUGTUNE.fpath)

  catln(0, "</beast>")
}

                   


 
  
  
taxa.array.to.structure <- function(taxa) {
  struct <- NULL
  spp <- collect.unique.string.parts(taxa, 3, 3)
  for (i in 1:length(spp)) {
    taxa.of.spp <- NULL
    for (j in 1:length(taxa)) {
      if (substr(taxa[j], 3, 3) == spp[i]) {
        taxa.of.spp <- c(taxa.of.spp, taxa[j])
      }
    }
    struct <- c(struct, list(taxa.of.spp.to.structure(spp[i], taxa.of.spp))) 
    }
  struct
  }


taxa.of.spp.to.structure <- function(sp, tofs) {
  indivs <- collect.unique.string.parts(tofs, 1, 2)
  sqs <- collect.unique.string.parts(tofs, 4, 4)
  list(sp=sp, indivs=indivs, sqs=sqs)
  } 


collect.unique.string.parts <- function(strings, beg, last) {
  parts <- substr(strings[1], beg, last)
  if (length(strings) > 1) {
    for (i in 2:length(strings)) {
      part <- substr(strings[i], beg, last)
      got = FALSE
      for (j in 1:length(parts)) {
        if (part == parts[j]) { got = TRUE }
        }
      if (!got) {
        parts <- c(parts, part)
        }
      }
    }
  parts
  } 





  
# Read in the sequences files and join to make some XML
make.taxa.and.alignments <- function(nofgenes, seqs.fpathbase) {
commentln("Taxa. Each taxon is a sequence (labelled A,B,...) from")
commentln("an individual (01,02,03,...) in a species (a,b,c...).")
totalnofseqs <- 0
  for (g in 1:nofgenes) {
    numastext <- as.character(g)
    textlines <- readLines(paste(seqs.fpathbase, numastext, ".txt", sep=""))
    n <- length(textlines)
    totalnofseqs <- totalnofseqs + n
    tokens <- strsplit(textlines, split=" +")
  
    if (g == 1) {
      taxa <- rep("", n-1)
      for (i in 2:n) {
        taxa[i-1] <- tokens[[i]][1]
        }
      taxa <- sort(taxa)
      
      cat("\t<taxa id=\"taxa\">\n")
      for (i in 2:n) {
        cat("\t\t<taxon id=\"", taxa[i-1], "\"/>\n", sep="")
        }
      cat("\t</taxa>\n")
      } else {
      # error check would be good
      }
    commentln(paste("Simulated alignment for gene", g, "of", nofgenes, "genes."))
    cat("\t<alignment id=\"alignment", numastext, "\" dataType=\"nucleotide\">\n", sep="")
    for (i in 2:n) {
      taxon <- tokens[[i]][1]
      seqnce <- tokens[[i]][2]
      cat("\t\t<sequence>\n", sep="")
      cat("\t\t\t<taxon idref=\"", taxon, "\" />\n", sep="")
      cat(seqnce, "\n", sep="")
      cat("\t\t</sequence>\n", sep="")
      }
    cat("\t</alignment>\n")
    }
  totalnofseqs
  }





gene.tree.weight <- function(wt, nofgenes, totalnofseqs) {
relativeweight <- totalnofseqs / (100*nofgenes)
as.integer(round(wt+wt*relativeweight))
}


make.operators <- function(nofgenes, totalnofseqs) { 
  commentln("Operators for gene trees.")
  catln(1, "<operators", id("operators"), ">")
  for (g in 1:nofgenes) {
    make.scaleOperator.kappa(g, weight=gene.tree.weight(1, nofgenes, totalnofseqs))
    }    
  for (g in 1:nofgenes) {
    make.deltaExchange.frequencies(g, weight=gene.tree.weight(1, nofgenes, totalnofseqs))
    }
  for (g in 1:nofgenes) {
    if (g > 1) {
      make.scaleOperator.clockrate(g, weight=gene.tree.weight(1, nofgenes, totalnofseqs))
      }
    }
  for (g in 1:nofgenes) {
    make.subtreeSlide.treeModel(g, weight=gene.tree.weight(15, nofgenes, totalnofseqs)) 
    }
  for (g in 1:nofgenes) {
    make.narrowExchange.treeModel(g, weight=gene.tree.weight(15, nofgenes, totalnofseqs)) 
    }    
  for (g in 1:nofgenes) {
    make.wideExchange.treeModel(g, weight=gene.tree.weight(3, nofgenes, totalnofseqs)) 
    }    
  for (g in 1:nofgenes) {
    make.wilsonBalding.treeModel(g, weight=gene.tree.weight(3, nofgenes, totalnofseqs)) 
    }    
  for (g in 1:nofgenes) {
    make.scaleOperator.treeModel.rootHeight(g, weight=gene.tree.weight(3, nofgenes, totalnofseqs)) 
    }    
  for (g in 1:nofgenes) {
    make.uniformOperator.treeModel.internalNodeHeights(g, weight=gene.tree.weight(30, nofgenes, totalnofseqs)) 
    }    
  # this moves the inode heights plus the clock rates relative 
  # to the first gene - note g==1 case done differently
  for (g in 1:nofgenes) {
    make.upDownOperator.clockrate.allInternalNodeHeights(g, weight=gene.tree.weight(30, nofgenes, totalnofseqs)) 
    }    

  commentln("Operator for both network and gene trees,stretches/squeezes everything.")
  make.network.upDownOperator.clocks.pops.heights(nofgenes, weight=30)
  commentln("Operators for network.")
  
  make.scaleOperator.species.popSF(weight=5+nofgenes)
  make.scaleOperator.apspnetwork.tipPopSizes(weight=10+2*nofgenes) 
  make.scaleOperator.apspnetwork.rootPopSizes(weight=10+2*nofgenes) 
  make.scaleOperator.apspnetwork.prior.eventRate(weight=30)
  make.hybPopSizesScaleOperator(weight=10+2*nofgenes)
  make.networkNodeReHeightOperator(weight=30+15*nofgenes)
  make.sequenceReassignmentOperator(weight=30+30*nofgenes)   
  make.changeNumHybridizationsOperator(weight=30+15*nofgenes)
    
  catln(1, "</operators>")
  }
  

 
  
make.logs <- function(nofgenes,
                      sampledgtrees.fpathbase, sampledmultrees.fpath,
                      sampledparams.fpath, DBUGTUNE.fpath) {
  commentln("Log to screen.")
  catln(2, "<log", id("screenLog"), namevalue("logEvery", beast.screen.logevery), ">") 
  catln(3, "<column", namevalue("label", "Posterior"), namevalue("dp", 4),
                       namevalue("width", 12), ">")
  catln(3, idref("posterior", "posterior"))
  catln(3, "</column>")
  catln(3, "<column", namevalue("label", "Prior"), namevalue("dp", 4),
                       namevalue("width", 12), ">")
  catln(3, idref("prior", "prior"))
  catln(3, "</column>")
  catln(3, "<column", namevalue("label", "Likelihood"), namevalue("dp", 4),
                       namevalue("width", 12), ">")
  catln(3, idref("likelihood", "likelihood"))
  catln(3, "</column>")
  catln(3, "<column", namevalue("label", "PopSF"), namevalue("sf", 6),
                       namevalue("width", 12), ">")
  catln(3, idref("parameter", "population.scaling.factor"))
  catln(3, "</column>")
  catln(2, "</log>")
  
  commentln("Main log to file.")
  catln(2, "<log", id("fileLog"), namevalue("logEvery", beast.params.logevery), 
            namevalue("fileName", sampledparams.fpath), ">") 
  catln(3, idref("posterior", "posterior"))
  catln(3, idref("prior", "prior"))
  catln(3, idref("likelihood", "likelihood"))

  catln(3, idref("alloppSpecies", "alloppSpecies"))
  catln(3, idref("apspCoalescent", "apsp.coalescent"))
  catln(3, idref("apspNetworkPrior", "apspnetwork.prior"))
  catln(3, idref("alloppNumHybsStatistic", "allopp.numHybs"))
  catln(3, idref("parameter", "population.scaling.factor"))
  catln(3, idref("parameter", "apspnetwork.tipPopSizes"))
  catln(3, idref("parameter", "apspnetwork.rootPopSizes"))
  catln(3, idref("parameter", "apspnetwork.hybridPopSizes"))
  catln(3, idref("parameter", "apspnetwork.prior.eventRate"))

  for (g in 1:nofgenes) {
    catln(3, numbered.idref("parameter", "", g, ".treeModel.rootHeight"))
    }
  for (g in 1:nofgenes) {
    catln(3, numbered.idref("parameter", "", g, ".kappa"))
    }
  for (g in 1:nofgenes) {
    catln(3, numbered.idref("parameter", "", g, ".frequencies"))
    }
  for (g in 1:nofgenes) {
    catln(3, numbered.idref("parameter", "", g, ".clock.rate"))
    }
  for (g in 1:nofgenes) {
    catln(3, numbered.idref("treeLikelihood", "", g, ".treeLikelihood"))
    }
  catln(2, "</log>")
  
  commentln("Gene tree log files.")
  for (g in 1:nofgenes) {
    fname <- paste(sampledgtrees.fpathbase, g, ".txt", sep="")
    catln(2, "<logTree", numbered.id("", g, ".treeFileLog"), namevalue("logEvery", beast.gtrees.logevery), 
              namevalue("fileName", fname), 
              namevalue("sortTranslationTable", "true"),
              namevalue("nexusFormat", "true"), ">") 
    catln(3, numbered.idref("treeModel", "", g, ".treeModel"))
    catln(2, "</logTree>")
    }
  
  commentln("MUL-tree log file.")
  catln(2, "<logTree", id("multreeFileLog"), namevalue("logEvery", beast.multree.logevery), 
         namevalue("fileName", sampledmultrees.fpath), 
         namevalue("sortTranslationTable", "true"),
         namevalue("nexusFormat", "true"), ">") 
  catln(3, idref("alloppSpeciesNetwork", "apspnetwork"))  
  catln(2, "</logTree>")    

  if (!is.null(DBUGTUNE.fpath)) {
    commentln("Debugging and tuning log.")
    catln(2, "<log", id("alloppDBUGTUNE"), namevalue("logEvery", beast.dbugtune.logevery), 
                      namevalue("fileName", DBUGTUNE.fpath), ">")
    catln(3, idref("alloppSpeciesNetwork", "apspnetwork"))  
    catln(2, "</log>")
    }
  } 


  
 
make.mcmc <- function(nofgenes, popmean.eventrate.relativerates.prior, sampledgtrees.fpathbase,
                   sampledmultrees.fpath, sampledparams.fpath, DBUGTUNE.fpath) { 
  catln(1, "<mcmc", id("mcmc"), namevalue("chainLength", beast.chain.length),
                       namevalue("autoOptimize", "true") , ">")
  catln(2, "<posterior", id("posterior"), namevalue("threads", 1), ">")
  catln(3, "<prior", id("prior"), ">") 
  catln(4, idref("apspCoalescent", "apsp.coalescent"))
  catln(4, idref("apspNetworkPrior", "apspnetwork.prior"))
  for (g in 1:nofgenes) {
    catln(4, "<logNormalPrior", namevalue("mean", 1.0), namevalue("stdev", 1.25),
           namevalue("offset", 0),  namevalue("meanInRealSpace", "false"), ">")
    catln(5, numbered.idref("parameter", "", g, ".kappa"))
    catln(4, "</logNormalPrior>")
    }
  # popmean hyperprior
  popmean.meanlog <- popmean.eventrate.relativerates.prior[1]      
  popmean.sdlog <- popmean.eventrate.relativerates.prior[2]      
  catln(4, "<logNormalPrior", namevalue("mean", popmean.meanlog), namevalue("stdev", popmean.sdlog),
         namevalue("offset", 0),  namevalue("meanInRealSpace", "false"), ">")
  catln(5, idref("parameter", "population.scaling.factor"))
  catln(4, "</logNormalPrior>")
  
  # speciation/hybrate prior
  eventrate.meanlog <- popmean.eventrate.relativerates.prior[3]      
  eventrate.sdlog <- popmean.eventrate.relativerates.prior[4]      
  catln(4, "<logNormalPrior", namevalue("mean", eventrate.meanlog), namevalue("stdev", eventrate.sdlog),
         namevalue("offset", 0),  namevalue("meanInRealSpace", "false"), ">")
  catln(5, idref("parameter", "apspnetwork.prior.eventRate"))
  catln(4, "</logNormalPrior>")
  
  # relative clock rates prior for genes 2,3,..
  relativerates.meanlog <- popmean.eventrate.relativerates.prior[5]      
  relativerates.sdlog <- popmean.eventrate.relativerates.prior[6]      
  for (g in 1:nofgenes) {
    if (g > 1) {
      catln(4, "<logNormalPrior", namevalue("mean", relativerates.meanlog), namevalue("stdev", relativerates.sdlog),
         namevalue("offset", 0),  namevalue("meanInRealSpace", "false"), ">")
      catln(5, numbered.idref("parameter", "", g, ".clock.rate"))
      catln(4, "</logNormalPrior>")
      }
    }
  catln(3, "</prior>")
  
  catln(3, "<likelihood", id("likelihood"), ">")
  for (g in 1:nofgenes) {
    catln(4, numbered.idref("treeLikelihood", "", g, ".treeLikelihood"))
    }
  catln(3, "</likelihood>")
  catln(2, "</posterior>")
  
  catln(2, idref("operators", "operators"))
  
  make.logs(nofgenes, sampledgtrees.fpathbase, sampledmultrees.fpath, 
                      sampledparams.fpath, DBUGTUNE.fpath)

  catln(1, "</mcmc>")
  }
  
  





