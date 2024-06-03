library(hisse)
library(diversitree)
library(tictoc)

#set working directory and options
	setwd("C:/a/geohisse")
	filename <- "Micrathena_GeoHiSSE"
	pdf_w <- 10
	pdf_h <- 20
	state_colors <- c("black", "white", "grey") #endemic1, endemic2 and widespread
	proportion_sampled_species <- c(1,1,1)
	number_models <- 10

#reads trees
	single_tree <- read.nexus("1_log12_MCCT_medianheights.nex") #file with the MCCT tree
	multiple_trees <- read.nexus("10_pruned_trees.nex") #file with multiple trees
	analysis_multiple_trees <- TRUE

#reads data
	data <- read.csv("1_all_species_elevation.txt", header = T) #0 is widespread area, 1 is highland (area 0 sensu GeoHiSSE), 2 is lowland (area 1 sensu GeoHiSSE)
	species <- data$Species

sink(paste(filename, "_logfile.txt", sep=""), append=FALSE, split=TRUE)
class(multiple_trees) <- "multiPhylo"

if(analysis_multiple_trees) {replication_number <- length(multiple_trees)} else {replication_number <- 1}

#prepares output tables
all_aic_weights <-  as.data.frame(matrix(NA, ncol = number_models+1, nrow = 0))
names(all_aic_weights) <- c("tree", paste("model_", 1:number_models, sep=""))
all_tau_estimates <-  as.data.frame(matrix(NA, ncol = 9, nrow = 0))
names_summary_table <- c("tree", "AIC","AICw","tau00A","tau11A", "tau01A","tau00B","tau11B", "tau01B")
names(all_tau_estimates) <- names_summary_table

# runs models for each tree
for(x in 1:replication_number){
start_time <- Sys.time()

#reads tree
	if(analysis_multiple_trees==FALSE) {phy <- keep.tip(single_tree, species)} else {phy <- keep.tip(multiple_trees[[x]], species)} #prunes tree to match data  
	if(analysis_multiple_trees) {output_filename <- paste("tree",x,"/",filename,sep="")} else {output_filename <- paste("MCCtree/",filename,sep="")}
	
#creates output directory
if(analysis_multiple_trees==FALSE) {
		if(file.exists("MCCtree")==FALSE) dir.create("MCCtree")
	}
if(file.exists(paste("tree",x,sep=""))==FALSE) dir.create(paste("tree",x,sep=""))

tic()

## Model 1 - No range-dependent diversification
	turnover <- c(1,1,0) #diversification rates are the same in the single-area ranges, and 0 in the widespread range
	eps <- c(1,1) #extinction rates are the same in the single-area ranges
	trans.rate <- TransMatMakerGeoHiSSE(
						hidden.traits=0, #no hidden area
						include.jumps=FALSE,
						separate.extirpation=FALSE) 
	mod1 <- GeoHiSSE(phy = phy, data = data, f=proportion_sampled_species,
                 turnover=turnover, eps=eps,
                 hidden.states=FALSE, trans.rate=trans.rate,
                 turnover.upper=100, trans.upper=10)
toc()
tic()

## Model 2. GeoHiSSE model with 1 hidden trait, no range-dependent diversification
	# Note below how parameters vary among hidden classes but are the same within each hidden class
	turnover <- c(1,1,0,2,2,0) #diversification rates are the same within non-hidden and hidden categories, but differ between hidden and non-hidden; diversification in widespread ranges removed from model
	eps <- c(1,1,2,2) #extinction rates are the same in the single-area ranges
	trans.rate <- TransMatMakerGeoHiSSE(
						hidden.traits=1,
						include.jumps=FALSE,
						separate.extirpation=FALSE) #one hidden area
	mod2 <- GeoHiSSE(phy = phy, data = data, f=proportion_sampled_species,
                   turnover=turnover, eps=eps,
                   hidden.states=TRUE, trans.rate=trans.rate,
                   turnover.upper=100, trans.upper=10)
toc()
tic()
 
## Model 3. GeoHiSSE model with 4 hidden traits, no range-dependent diversification.
	turnover <- c(1,1,0,2,2,0,3,3,0,4,4,0,5,5,0) #diversification rates are the same within non-hidden and hidden categories, but differ between hidden and non-hidden; diversification in widespread ranges removed from model
	eps <- c(1,1,2,2,3,3,4,4,5,5) #extinction rates are the same in the single-area ranges
	trans.rate <- TransMatMakerGeoHiSSE(
						hidden.traits=4, #four hidden areas
						include.jumps=FALSE,
						separate.extirpation=FALSE) 
	mod3 <- GeoHiSSE(phy = phy, data = data, f=proportion_sampled_species,
                   turnover=turnover, eps=eps,
                   hidden.states=TRUE, trans.rate=trans.rate,
                   turnover.upper=100, trans.upper=10)
toc()
tic()

## Model 4. Canonical GeoSSE model, range effect on diversification
	turnover <- c(1,2,3) #each single-range area and the widespread area have different diversification rates
	eps <- c(1,2) #extirpation rates differ among areas
	trans.rate <- TransMatMakerGeoHiSSE(
						hidden.traits=0, #no hidden area
						include.jumps=FALSE,
						separate.extirpation=FALSE) 
	mod4 <- GeoHiSSE(phy = phy, data = data, f=proportion_sampled_species,
                   turnover=turnover, eps=eps,
                   hidden.states=FALSE, trans.rate=trans.rate,
                   turnover.upper=100, trans.upper=10)
toc()
tic()

## Model 5. GeoHiSSE model with 1 hidden trait, range-dependent diversification.
	turnover <- c(1,2,3,4,5,6) #diversification rates are different for each of the single-area and widespread ranges, and vary between non-hidden and hidden areas
	eps <- c(1,2,3,4) #extirpation rates differ among areas
	trans.rate <- TransMatMakerGeoHiSSE(
						hidden.traits=1, #one hidden area
						include.jumps=FALSE,
						separate.extirpation=FALSE) 
	mod5 <- GeoHiSSE(phy = phy, data = data, f=proportion_sampled_species,
                   turnover=turnover, eps=eps,
                   hidden.states=TRUE, trans.rate=trans.rate,
                   turnover.upper=100, trans.upper=10)
toc()
tic()

## Model 6. GeoHiSSE model with 4 hidden traits, range-dependent diversification.
	turnover <- c(1:15) #diversification rates are different for each of the single-area and widespread ranges, and vary between non-hidden and hidden areas
	eps <- c(1:10) #extinction rates are the same in the single-area ranges
	trans.rate <- TransMatMakerGeoHiSSE(
						hidden.traits=4, #four hidden areas
						include.jumps=FALSE,
						separate.extirpation=FALSE) 
	mod6 <- GeoHiSSE(phy = phy, data = data, f=proportion_sampled_species,
                   turnover=turnover, eps=eps,
                   hidden.states=TRUE, trans.rate=trans.rate,
                   turnover.upper=100, trans.upper=10)				   
toc()
tic()


## Model 7. MuSSE-like model with no hidden trait, no cladogenetic effects, no range-dependent diversification.
	turnover <- c(1,1,0)
	eps <- c(1,1)
	trans.rate <- TransMatMakerGeoHiSSE(
						hidden.traits=0,
						include.jumps=FALSE,
                        separate.extirpation = FALSE)
	mod7 <- GeoHiSSE(phy = phy, data = data, f=proportion_sampled_species,
                   turnover=turnover, eps=eps,
                   hidden.states=FALSE, trans.rate=trans.rate,
                   turnover.upper=100, trans.upper=10, sann=FALSE,
                   assume.cladogenetic = FALSE)
toc()
tic()

## Model 8. MuSSE-like model with 1 hidden trait, no cladogenetic effects, no range-dependent diversification.
	turnover <- c(1,1,0,2,2,0)
	eps <- c(1,1,2,2)
	trans.rate <- TransMatMakerGeoHiSSE(
							hidden.traits=1,
							include.jumps=FALSE,
                            separate.extirpation = FALSE)
	mod8 <- GeoHiSSE(phy = phy, data = data, f=proportion_sampled_species,
                   turnover=turnover, eps=eps,
                   hidden.states=TRUE, trans.rate=trans.rate,
                   turnover.upper=100, trans.upper=10, sann=FALSE,
                   assume.cladogenetic = FALSE)
toc()
tic()

## Model 9. MuSSE-like model with no hidden trait, no cladogenetic effects, range effect on diversification
	turnover <- c(1,2,0)
	eps <- c(1,2)
	trans.rate <- TransMatMakerGeoHiSSE(
						hidden.traits=0,
						include.jumps=FALSE,
                        separate.extirpation = FALSE)
	mod9 <- GeoHiSSE(phy = phy, data = data, f=proportion_sampled_species,
                   turnover=turnover, eps=eps,
                   hidden.states=FALSE, trans.rate=trans.rate,
                   turnover.upper=100, trans.upper=10, sann=FALSE,
                   assume.cladogenetic = FALSE)
toc()
tic()

## Model 10. MuSSE-like model with 1 hidden trait, no cladogenetic effects, range effect on diversification
	turnover <- c(1,2,0,3,4,0)
	eps <- c(1,2,3,4)
	trans.rate <- TransMatMakerGeoHiSSE(
						hidden.traits=1,
						include.jumps=FALSE,
                        separate.extirpation = FALSE)
	mod10 <- GeoHiSSE(phy = phy, data = data, f=proportion_sampled_species,
                   turnover=turnover, eps=eps,
                   hidden.states=TRUE, trans.rate=trans.rate,
                   turnover.upper=100, trans.upper=10, sann=FALSE,
                   assume.cladogenetic = FALSE)
toc()
tic()

#saves
	results_GeoHiSSE <- list(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10)
	saveRDS(results_GeoHiSSE, file=paste(output_filename,"_model_results.Rdata", sep =""))

#gets AIC weights of the models
  aicweights <- GetAICWeights(list(model1 = mod1, model2 = mod2, model3 = mod3, model4 = mod4, model5 =mod5, model6=mod6, model7 =mod7, model8 =mod8, model9=mod9, model10=mod10), criterion="AIC")
  write.csv(aicweights, paste(output_filename,"_aicweights.csv", sep="_"))
  all_aic_weights <- rbind(all_aic_weights, c(x, aicweights))

#marginal reconstruction for each of the models in the set
  recon.mod1 <- MarginReconGeoSSE(phy = mod1$phy, data = mod1$data, f = mod1$f,
                                  pars = mod1$solution, hidden.states = 1,
                                  root.type = mod1$root.type, root.p = mod1$root.p,
                                  AIC = mod1$AIC, n.cores = 1)
  recon.mod2 <- MarginReconGeoSSE(phy = mod2$phy, data = mod2$data, f = mod2$f,
                                  pars = mod2$solution, hidden.states = 2,
                                  root.type = mod2$root.type, root.p = mod2$root.p,
                                  AIC = mod2$AIC, n.cores = 1)
  recon.mod3 <- MarginReconGeoSSE(phy = mod3$phy, data = mod3$data, f = mod3$f,
                                  pars = mod3$solution, hidden.states = 5,
                                  root.type = mod3$root.type, root.p = mod3$root.p,
                                  AIC = mod3$AIC, n.cores = 1)
  recon.mod4 <- MarginReconGeoSSE(phy = mod4$phy, data = mod4$data, f = mod4$f,
                                  pars = mod4$solution, hidden.states = 1,
                                  root.type = mod4$root.type, root.p = mod4$root.p,
                                  AIC = mod4$AIC, n.cores = 1)
  recon.mod5 <- MarginReconGeoSSE(phy = mod5$phy, data = mod5$data, f = mod5$f,
                                  pars = mod5$solution, hidden.states = 2,
                                  root.type = mod5$root.type, root.p = mod5$root.p,
                                  AIC = mod5$AIC, n.cores = 1)
  recon.mod6 <- MarginReconGeoSSE(phy = mod6$phy, data = mod6$data, f = mod6$f,
                                  pars = mod6$solution, hidden.states = 5,
                                  root.type = mod6$root.type, root.p = mod6$root.p,
                                  AIC = mod6$AIC, n.cores = 1)
  recon.mod7 <- MarginReconGeoSSE(phy = mod7$phy, data = mod7$data, f = mod7$f,
                                  pars = mod7$solution, hidden.states = 1,
                                  root.type = mod7$root.type, root.p = mod7$root.p,
                                  AIC = mod7$AIC, n.cores = 1)
  recon.mod8 <- MarginReconGeoSSE(phy = mod8$phy, data = mod8$data, f = mod8$f,
                                  pars = mod8$solution, hidden.states = 2,
                                  root.type = mod8$root.type, root.p = mod8$root.p,
                                  AIC = mod8$AIC, n.cores = 1)									  
  recon.mod9 <- MarginReconGeoSSE(phy = mod9$phy, data = mod9$data, f = mod9$f,
                                  pars = mod9$solution, hidden.states = 1,
                                  root.type = mod9$root.type, root.p = mod9$root.p,
                                  AIC = mod9$AIC, n.cores = 1)	
  recon.mod10 <- MarginReconGeoSSE(phy = mod10$phy, data = mod10$data, f = mod10$f,
                                  pars = mod10$solution, hidden.states = 2,
                                  root.type = mod10$root.type, root.p = mod10$root.p,
                                  AIC = mod10$AIC, n.cores = 1)	
toc()
tic()
								  
#saves
	recon.models <- list(recon.mod1, recon.mod2, recon.mod3, recon.mod4, recon.mod5, recon.mod6, recon.mod7, recon.mod8, recon.mod9, recon.mod10)
	saveRDS(recon.models, file=paste(output_filename,"_marginal_reconstructions.Rdata", sep =""))

#model average
  model.ave.rates <- GetModelAveRates(x = recon.models, type = "tips")	

#plot results
	#each model
		for(i in 1:length(recon.models)){
			pdf(file=paste(output_filename,"_mod",i,"_netdiv.pdf", sep =""), width=pdf_w, height=pdf_h)
			plot.geohisse.states(x = recon.models[[i]], rate.param = "net.div", type= "phylogram", show.tip.label = T, legend = T, fsize= 0.5, state.colors = state_colors)
			while (!is.null(dev.list())) {dev.off()}
		}
	
	#model average
		pdf(file=paste(output_filename,"_model_average_netdiv.pdf", sep =""), width=pdf_w, height=pdf_h)
		plot.geohisse.states(x = recon.models, rate.param = "net.div", type= "phylogram", show.tip.label = T, legend = T, fsize= 0.5, state.colors = state_colors)
		while (!is.null(dev.list())) {dev.off()}

#get AIC and number of free parameters
	aic_models <- rep("-", 10)
	for(a in 1:length(results_GeoHiSSE)){
		aic_models[a] <- results_GeoHiSSE[[a]]$AIC
	}

#get net diversification for each area/model
	tau_mod1 <- mod1$solution[c("tau00A","tau11A", "tau01A","tau00B","tau11B", "tau01B")]
	tau_mod2 <- mod2$solution[c("tau00A","tau11A", "tau01A","tau00B","tau11B", "tau01B")]
	tau_mod3 <- mod3$solution[c("tau00A","tau11A", "tau01A","tau00B","tau11B", "tau01B")]
	tau_mod4 <- mod4$solution[c("tau00A","tau11A", "tau01A","tau00B","tau11B", "tau01B")]
	tau_mod5 <- mod5$solution[c("tau00A","tau11A", "tau01A","tau00B","tau11B", "tau01B")]
	tau_mod6 <- mod6$solution[c("tau00A","tau11A", "tau01A","tau00B","tau11B", "tau01B")]
	tau_mod7 <- mod7$solution[c("tau00A","tau11A", "tau01A","tau00B","tau11B", "tau01B")]
	tau_mod8 <- mod8$solution[c("tau00A","tau11A", "tau01A","tau00B","tau11B", "tau01B")]
	tau_mod9 <- mod9$solution[c("tau00A","tau11A", "tau01A","tau00B","tau11B", "tau01B")]
	tau_mod10 <- mod10$solution[c("tau00A","tau11A", "tau01A","tau00B","tau11B", "tau01B")]
	
	tau_all_models <- as.data.frame(rbind(tau_mod1, tau_mod2, tau_mod3, tau_mod4, tau_mod5, tau_mod6, tau_mod7, tau_mod8, tau_mod9, tau_mod10))
	if(analysis_multiple_trees) {tree_number <- rep(x, number_models)} else {tree_number <- rep("MCCT", length(number_models))}
	tau_all_models <- cbind(tree_number, aic_models, aicweights, tau_all_models)
	names(tau_all_models) <- names_summary_table
	write.csv(tau_all_models, paste(output_filename,"_tau_all_models.csv",sep=""))
	all_tau_estimates <- rbind(all_tau_estimates,tau_all_models)

print(paste("Completed tree", x))
print("Total time elapsed:")
Sys.time() - start_time 
toc()

}

write.csv(all_tau_estimates, paste(filename,"_all_tau_estimates.csv",sep=""))
write.csv(all_aic_weights, paste(filename,"_all_aic_weights.csv",sep=""))
write.csv(model.ave.rates, paste(filename,"_model_ave_rates.csv",sep=""))

sink(file = NULL)