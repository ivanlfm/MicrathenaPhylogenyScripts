#################################################################################
### BioGeoBEARS                                                               ###
### Estimates ancestral ranges under the DEC, DIVALIKE and BAYAREALIKE models ###
### using different time stratifications and dispersal matrices               ###
### and compares fit of different models using AICc.                          ###
### Modified by Ivan. L. F. Magalhaes										  ###
### from original example script by Nicholas J. Matzke:              		  ###
### http://phylo.wikidot.com/biogeobears#script                               ###
#################################################################################

# loads packages
	library(rexpokit)
	library(cladoRcpp)
	library(BioGeoBEARS)
		#if you don't have BioGeoBEARS, install it from github (you need Rtools and devtools installed):
		#install_github(repo="nmatzke/BioGeoBEARS")
	library(GenSA)
	library(FD)
	library(snow)
	library(parallel)
	library(phytools)
	library(tictoc)
	extdata_dir = np(system.file("extdata", package="BioGeoBEARS")) # BioGeoBEARS needs this
	scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

###############################
### INPUT FILES & FILENAMES ###
###############################

	# sets working directory, logfile, number of cores to use
		setwd("d:/ILFM/all_species_6_areas")
		analysis_name <- "Micrathena_all_6areas"
		output_dir <- "results_model_test"
		out_name <- paste(output_dir, "/", analysis_name, sep="")
		dir.create(output_dir)
		sink(paste(out_name, "_logfile.txt", sep=""), append=FALSE, split=TRUE) 
		n_cores <- 8
		pdf_w <- 10 #dimensions for pdf output
		pdf_h <- 18 
	
	#set this to FALSE to load results from previous runs, TRUE to run the analyses
		run_analysis = TRUE
			# use start_from_model to add new models (e.g. add a 4th model to a script that had only 3 models), run only new models and re-run model comparison
			# run_analysis must be set to FALSE
			# warning: this will set run_analysis to TRUE after the (start_from_model)th model, so
			# if start_model is set to e.g. 1, it WILL re-run all models starting from model 1
			# so if you only want to re-load previous models, better comment this line...
				start_from_model = 1
	
	# define the names of each of the models you intend to run. You have to set them in the script below	
	# these will be the names of output files (plots, result objects, etc)
		filename_vector <- c(
			paste(out_name, "_01unconstrained_DEC", sep=""), #model 1
			paste(out_name, "_02unconstrained_DVL", sep=""), #model 2
			paste(out_name, "_03unconstrained_BAL", sep=""), #model 3... etc
			paste(out_name, "_04unconstrained_DECJ", sep=""),
			paste(out_name, "_05unconstrained_DVLJ", sep=""),
			paste(out_name, "_06unconstrained_BALJ", sep=""),
			paste(out_name, "_07dispmat_DEC", sep=""),
			paste(out_name, "_08dispmat_DVL", sep=""),
			paste(out_name, "_09dispmat_BAL", sep=""),
			paste(out_name, "_10dispmat_DECJ", sep=""),
			paste(out_name, "_11dispmat_DVLJ", sep=""),
			paste(out_name, "_12dispmat_BALJ", sep="")	
		)
		
	#creates a vector to store the inputs for the analysis and their names
		number_analyses <- length(filename_vector) # the number of models you intend to set
		BGB_inputs <- vector(mode = "list", length = number_analyses)
		names(BGB_inputs) <- filename_vector
	
	# reads distributions, give in Phylip format [first row: Ntaxa Nareas (A B C...)]
	# other rows:
	# SpeciesA 100
	# SpeciesB 010
	# SpeciesC 011
		geogfn = "1_all_species_6.txt"
		tipranges = getranges_from_LagrangePHYLIP(geogfn)
		# gets the species names from the distribution file
		# species without range data will be dropped from the input tree later, so check carefully
			species <- rownames(tipranges@df)
		
	# reads tree - nexus format, needs to be fully resolved, have no very short branches
	# ages should be in millions of years, values between 0 and 1000
	# preferrably a 'target' tree -- consensus or maximum clade credibility
	# trfn = tree file name. BioGeoBEARS needs the complete path
		#reads the tree, prunes terminals to keep only species with geographical data, saves a copy
			unpruned_tree <- read.nexus("1_log12_MCCT_medianheights.nex")
			pruned_tree <- ladderize(keep.tip(unpruned_tree, tip = species),right=TRUE)
			trfn <-  paste(out_name,"pruned_input_tree.nwk", sep= "")
			write.tree(pruned_tree, file = trfn)
			tr <- read.tree(trfn)
			plot(tr) # let's plot it just because we love cladograms
			axisPhylo() # especially if they are dated!
		
	# number of areas in your dataset
		numareas <- length(tipranges@df[1,])
	# maximum range size observed in the tip ranges
		max_observed_range_size <- max(rowSums(dfnums_to_numeric(tipranges@df)))
	# defines max range size
		# by default, the number of areas; this can slow down the analysis, especially if numareas is greater than 9
		# you can set max_range_size manually to any value between max_observed_range_size and numareas
			max_range_size <- numareas
		
	# the total number of possible ranges depends both on the number of areas and max number of areas a species can occupy
	# should be less than 1000 to run analysis in under a day, less than 1500 to run in under a week, less than 2500 to run at all
		numstates_from_numareas(numareas=numareas, maxareas=max_range_size, include_null_range=TRUE)
		
	#this will make model testing easier later by adjusting the size of the table and collecting the names and results
		analysis_count <- 1
		model_count <- 1
		model_names <- vector(mode = "list", length = 0)
		all_analyses <- vector(mode = "list", length = 0)
	
######################
######################
### SET UP MODELS  ###
######################
######################

	# below, define the parameters and inputs for each of the models you want to test
	# modify the code below for your own purposes	
	
##########################
### BioGeoBEARS object ###
##########################

	# defines a BioGeoBEARS object, inputs trees and ranges, sets various parameters for the program
	# these parameters will be used for downstream models
	# other models will be constructed as modifications of this object
		BioGeoBEARS_run_object_DEC = define_BioGeoBEARS_run() 		# Initialize a default model
		BioGeoBEARS_run_object_DEC$trfn = trfn 						# Give BioGeoBEARS the location of the phylogeny Newick file
		BioGeoBEARS_run_object_DEC$geogfn = geogfn 					# Give BioGeoBEARS the location of the geography text file
		BioGeoBEARS_run_object_DEC$max_range_size = max_range_size 	# Input the maximum range size
		BioGeoBEARS_run_object_DEC$min_branchlength = 0.000001 		# Min to treat tip as a direct ancestor (no speciation event)
		BioGeoBEARS_run_object_DEC$include_null_range = TRUE 		# set to FALSE for e.g. DEC* model, DEC*+J, etc.
		BioGeoBEARS_run_object_DEC$on_NaN_error = -1e50 			# returns very low lnL if parameters produce NaN error (underflow check)
		BioGeoBEARS_run_object_DEC$speedup = TRUE          			# shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
		BioGeoBEARS_run_object_DEC$use_optimx = "GenSA"   			# if FALSE, use optim() instead of optimx(); you can also use "GenSA" that is slower but searchers the parameter space more thoroughly
		BioGeoBEARS_run_object_DEC$num_cores_to_use <- n_cores 		# use more cores to speed it up; this requires library(parallel) and/or library(snow)
		BioGeoBEARS_run_object_DEC$force_sparse = FALSE				# force_sparse=TRUE causes pathology & isn't much faster at this scale
		BioGeoBEARS_run_object_DEC$return_condlikes_table = TRUE 	# Good default settings to get ancestral states
		BioGeoBEARS_run_object_DEC$calc_TTL_loglike_from_condlikes_table = TRUE # Good default settings to get ancestral states
		BioGeoBEARS_run_object_DEC$calc_ancprobs = TRUE    			# get ancestral states from optim run

###############################
### UNCONSTRAINED ANALYSES  ###
###############################
	
	#################
	### DEC model ###
	#################
		# this is identical to Lagrange's original Dispersal-Extinction-Cladogenesis model
		
		# the object set above runs the DEC model by default, so just store the run object in the BGB_inputs vector
			BGB_inputs[[model_count]] <- BioGeoBEARS_run_object_DEC
			filename_vector[[model_count]]
			check_BioGeoBEARS_run(BGB_inputs[[model_count]]) # checks the object and parameters
			model_count <- model_count+1
			
	######################
	### DIVALIKE model ###
	######################
		# this is a modified version of Dispersal-Extinction-Cladogenesis model, but with DIVA's inheritance mode at cladogenetic events
		
		#copies the run object of DEC and only modifies relevant parameters
			BioGeoBEARS_run_object_DVL = BioGeoBEARS_run_object_DEC
		
		# Set up DIVALIKE model
		# Remove subset-sympatry
			BioGeoBEARS_run_object_DVL$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
			BioGeoBEARS_run_object_DVL$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
			BioGeoBEARS_run_object_DVL$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
			BioGeoBEARS_run_object_DVL$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
			BioGeoBEARS_run_object_DVL$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
			BioGeoBEARS_run_object_DVL$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
			BioGeoBEARS_run_object_DVL$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
		# Allow classic, widespread vicariance; all events equiprobable
			BioGeoBEARS_run_object_DVL$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
			BioGeoBEARS_run_object_DVL$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
			BioGeoBEARS_run_object_DVL$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
		
		# store the run object in the BGB_inputs vector
			BGB_inputs[[model_count]] <- BioGeoBEARS_run_object_DVL
			filename_vector[[model_count]]
			check_BioGeoBEARS_run(BGB_inputs[[model_count]]) # checks the object and parameters
			model_count <- model_count+1
				
	#########################
	### BAYAREALIKE model ###
	#########################
		# this is a modified version of Dispersal-Extinction-Cladogenesis model, but with descendants inheriting exactly the same range as the ancestors
		
		# copies the run object of DEC and only modifies relevant parameters
			BioGeoBEARS_run_object_BAL = BioGeoBEARS_run_object_DEC
			
		# Set up BAYAREALIKE model
		# No subset sympatry
			BioGeoBEARS_run_object_BAL$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
			BioGeoBEARS_run_object_BAL$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
			BioGeoBEARS_run_object_BAL$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
		# No vicariance
			BioGeoBEARS_run_object_BAL$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
			BioGeoBEARS_run_object_BAL$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
			BioGeoBEARS_run_object_BAL$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
		# Adjust linkage between parameters
			BioGeoBEARS_run_object_BAL$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
			BioGeoBEARS_run_object_BAL$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
			BioGeoBEARS_run_object_BAL$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
		# Only sympatric/range-copying (y) events allowed, and with exact copying (both descendants always the same size as the ancestor)
			BioGeoBEARS_run_object_BAL$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
			BioGeoBEARS_run_object_BAL$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
			BioGeoBEARS_run_object_BAL$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
		
		# store the run object in the BGB_inputs vector
			BGB_inputs[[model_count]] <- BioGeoBEARS_run_object_BAL
			filename_vector[[model_count]]
			check_BioGeoBEARS_run(BGB_inputs[[model_count]]) # checks the object and parameters
			model_count <- model_count+1
			
	#############################################
	### MODELS WITH FOUNDER-EVENT SPECIATION  ###
	#############################################
		# beware of overfitting!! see Ree & SanMartín 2018: https://doi.org/10.1111/jbi.13173
		# copies the run objects above and adds jump dispersal (founder-event speciation)
		
		### DEC+J
			BioGeoBEARS_run_object_DECJ <- BioGeoBEARS_run_object_DEC
			BioGeoBEARS_run_object_DECJ$BioGeoBEARS_model_object@params_table["j","type"] = "free"
			BioGeoBEARS_run_object_DECJ$BioGeoBEARS_model_object@params_table["j","init"] = 0.0001
			BioGeoBEARS_run_object_DECJ$BioGeoBEARS_model_object@params_table["j","est"] = 0.0001
			
			# store the run object in the BGB_inputs vector
				BGB_inputs[[model_count]] <- BioGeoBEARS_run_object_DECJ
				filename_vector[[model_count]]
				check_BioGeoBEARS_run(BGB_inputs[[model_count]]) # checks the object and parameters
				model_count <- model_count+1
		
		### DIVALIKE+J
			BioGeoBEARS_run_object_DVLJ <- BioGeoBEARS_run_object_DVL
			BioGeoBEARS_run_object_DVLJ$BioGeoBEARS_model_object@params_table["j","type"] = "free"
			BioGeoBEARS_run_object_DVLJ$BioGeoBEARS_model_object@params_table["j","init"] = 0.0001
			BioGeoBEARS_run_object_DVLJ$BioGeoBEARS_model_object@params_table["j","est"] = 0.0001
			# Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
			BioGeoBEARS_run_object_DVLJ$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
			BioGeoBEARS_run_object_DVLJ$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999
			
			# store the run object in the BGB_inputs vector
				BGB_inputs[[model_count]] <- BioGeoBEARS_run_object_DVLJ
				filename_vector[[model_count]]
				check_BioGeoBEARS_run(BGB_inputs[[model_count]]) # checks the object and parameters
				model_count <- model_count+1	
				
		### BAYAREALIKE+J
			BioGeoBEARS_run_object_BALJ <- BioGeoBEARS_run_object_BAL
			BioGeoBEARS_run_object_BALJ$BioGeoBEARS_model_object@params_table["j","type"] = "free"
			BioGeoBEARS_run_object_BALJ$BioGeoBEARS_model_object@params_table["j","init"] = 0.0001
			BioGeoBEARS_run_object_BALJ$BioGeoBEARS_model_object@params_table["j","est"] = 0.0001
			# Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
			BioGeoBEARS_run_object_BALJ$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
			BioGeoBEARS_run_object_BALJ$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999
			
			# store the run object in the BGB_inputs vector
				BGB_inputs[[model_count]] <- BioGeoBEARS_run_object_BALJ
				filename_vector[[model_count]]
				check_BioGeoBEARS_run(BGB_inputs[[model_count]]) # checks the object and parameters
				model_count <- model_count+1
				
	###################################################
	### ANALYSES WITH MODIFIED DISPERSAL MATRIX + w	###
	###################################################
			
		###### defines a dispersal matrix ######
				dispersal_matrix = "2_dispersalmatrix_adjacency.txt"
				
			# adds dispersal matrix to each object (DEC, DIVALIKE and BAYAREALIKE and + J variants)
				#DEC + matrix
					BioGeoBEARS_run_object_DEC_matrix <- BioGeoBEARS_run_object_DEC
					BioGeoBEARS_run_object_DEC_matrix$BioGeoBEARS_model_object@params_table["w","type"] <- "free"
					BioGeoBEARS_run_object_DEC_matrix$dispersal_multipliers_fn = dispersal_matrix
					BioGeoBEARS_run_object_DEC_matrix = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object_DEC_matrix)
				
					# store the run object in the BGB_inputs vector
						BGB_inputs[[model_count]] <- BioGeoBEARS_run_object_DEC_matrix
						filename_vector[[model_count]]
						check_BioGeoBEARS_run(BGB_inputs[[model_count]]) # checks the object and parameters
						model_count <- model_count+1

				#DVL + matrix
					BioGeoBEARS_run_object_DVL_matrix <- BioGeoBEARS_run_object_DVL
					BioGeoBEARS_run_object_DVL_matrix$BioGeoBEARS_model_object@params_table["w","type"] <- "free"
					BioGeoBEARS_run_object_DVL_matrix$dispersal_multipliers_fn = dispersal_matrix
					BioGeoBEARS_run_object_DVL_matrix = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object_DVL_matrix)
				
					# store the run object in the BGB_inputs vector
						BGB_inputs[[model_count]] <- BioGeoBEARS_run_object_DVL_matrix
						filename_vector[[model_count]]
						check_BioGeoBEARS_run(BGB_inputs[[model_count]]) # checks the object and parameters
						model_count <- model_count+1						
				
				#BAL + matrix
					BioGeoBEARS_run_object_BAL_matrix <- BioGeoBEARS_run_object_BAL
					BioGeoBEARS_run_object_BAL_matrix$BioGeoBEARS_model_object@params_table["w","type"] <- "free"
					BioGeoBEARS_run_object_BAL_matrix$dispersal_multipliers_fn = dispersal_matrix
					BioGeoBEARS_run_object_BAL_matrix = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object_BAL_matrix)
				
					# store the run object in the BGB_inputs vector
						BGB_inputs[[model_count]] <- BioGeoBEARS_run_object_BAL_matrix
						filename_vector[[model_count]]
						check_BioGeoBEARS_run(BGB_inputs[[model_count]]) # checks the object and parameters
						model_count <- model_count+1
										
				#DECJ + matrix
					BioGeoBEARS_run_object_DECJ_matrix <- BioGeoBEARS_run_object_DECJ
					BioGeoBEARS_run_object_DECJ_matrix$BioGeoBEARS_model_object@params_table["w","type"] <- "free"
					BioGeoBEARS_run_object_DECJ_matrix$dispersal_multipliers_fn = dispersal_matrix
					BioGeoBEARS_run_object_DECJ_matrix = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object_DECJ_matrix)
				
					# store the run object in the BGB_inputs vector
						BGB_inputs[[model_count]] <- BioGeoBEARS_run_object_DECJ_matrix
						filename_vector[[model_count]]
						check_BioGeoBEARS_run(BGB_inputs[[model_count]]) # checks the object and parameters
						model_count <- model_count+1	
										
				#DVLJ + matrix
					BioGeoBEARS_run_object_DVLJ_matrix <- BioGeoBEARS_run_object_DVLJ
					BioGeoBEARS_run_object_DVLJ_matrix$BioGeoBEARS_model_object@params_table["w","type"] <- "free"
					BioGeoBEARS_run_object_DVLJ_matrix$dispersal_multipliers_fn = dispersal_matrix
					BioGeoBEARS_run_object_DVLJ_matrix = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object_DVLJ_matrix)
				
					# store the run object in the BGB_inputs vector
						BGB_inputs[[model_count]] <- BioGeoBEARS_run_object_DVLJ_matrix
						filename_vector[[model_count]]
						check_BioGeoBEARS_run(BGB_inputs[[model_count]]) # checks the object and parameters
						model_count <- model_count+1			
						
				#BALJ + matrix
					BioGeoBEARS_run_object_BALJ_matrix <- BioGeoBEARS_run_object_BALJ
					BioGeoBEARS_run_object_BALJ_matrix$BioGeoBEARS_model_object@params_table["w","type"] <- "free"
					BioGeoBEARS_run_object_BALJ_matrix$dispersal_multipliers_fn = dispersal_matrix
					BioGeoBEARS_run_object_BALJ_matrix = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object_BALJ_matrix)
				
					# store the run object in the BGB_inputs vector
						BGB_inputs[[model_count]] <- BioGeoBEARS_run_object_BALJ_matrix
						filename_vector[[model_count]]
						check_BioGeoBEARS_run(BGB_inputs[[model_count]]) # checks the object and parameters
						model_count <- model_count+1	
					
########################
#### RUNS ALL MODELS ###
########################

	# just in case, let's re-check that the BGB objects are in good shape:
		lapply(BGB_inputs, check_BioGeoBEARS_run)
		
	# for each of the inputs stored in BGB_inputs,
	# runs analyses or loads previous results
	
	# by default this runs all analysis stored in BGB_inputs
	# change the value of run_how_many_analyses if you want to run just the first e.g. 3 models
		run_how_many_analyses = length(BGB_inputs)
	
	while (!is.null(dev.list())) {dev.off()} #clear the plotting space
	for (i in 1:run_how_many_analyses){
		if(check_BioGeoBEARS_run(BGB_inputs[[i]])){ #skips the analysis if the BGB object is not good to go
			filename = filename_vector[[i]] #gets the name for the output files
			if (i==start_from_model) {run_analysis=TRUE} # by default it starts from model 1 (e.g. runs all the models). Change start_from_model if you have added more models later and want to run model comparison between them and the ones you had already run
			if (run_analysis) { #runs analysis...
				tic()
				print(paste("Model ", i, ": ", filename_vector[[i]], sep=""))
				BGB_res = bears_optim_run(BGB_inputs[[i]])
				save(BGB_res, file=paste(filename,".Rdata", sep =""))
				pdf(paste(filename,"_pie.pdf", sep =""), width=pdf_w, height=pdf_h)
					plot_BioGeoBEARS_results(BGB_res, filename, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
				while (!is.null(dev.list())) {dev.off()}
				pdf(paste(filename,"_text.pdf", sep =""), width=pdf_w, height=pdf_h)
					plot_BioGeoBEARS_results(BGB_res, filename, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
				while (!is.null(dev.list())) {dev.off()}
				jpeg(paste(filename,"_pie.jpg", sep =""), w= pdf_w*100, h= pdf_h*100, pointsize = 18, quality = 90)
						plot_BioGeoBEARS_results(BGB_res, filename, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
				while (!is.null(dev.list())) {dev.off()}
				toc()
			} else { #or loads results.
				load(paste(filename,".Rdata", sep =""))
			}
		
			model_names[[analysis_count]] <- filename
			all_analyses[[analysis_count]] <- BGB_res
			analysis_count <- analysis_count+1
		}
	}

#################################
#### TESTS MODEL FIT WITH AIC ###
#################################

	# set up empty tables to store the statistical results
		n_models <- length(all_analyses)
		n_taxa = length(tr$tip.label)
		model_comparison <- matrix(NA, ncol = 8, nrow = n_models)
		rownames(model_comparison) <-  model_names
		colnames(model_comparison) <-  c("LnL", "k", "AICc", "AICw", "d", "e", "j", "w")
	
	# defines a function to collect relevant results from the output of BioGeoBEARS	
		fetch_BGB_results <- function(BioGeoBEARS_result_object, n_taxa) {
			# gets log-likelihood of the data given this model
				LnL <- get_LnL_from_BioGeoBEARS_results_object(BioGeoBEARS_result_object)
			#gets number of free parameters the model
				n_parameters <- sum(BioGeoBEARS_result_object$inputs$BioGeoBEARS_model_object@params_table[,1]=="free")
			# calculates AICc
				aicc <- (-2*LnL)+((2*n_parameters)+(((2*n_parameters)*(n_parameters+1))/(n_taxa-n_parameters-1)))
			# gets estimated rates of dispersal and extinction
				dispersal <- BioGeoBEARS_result_object$outputs@params_table["d","est"]
				extinction <- BioGeoBEARS_result_object$outputs@params_table["e","est"]
				j <- BioGeoBEARS_result_object$outputs@params_table["j","est"]
				w <- BioGeoBEARS_result_object$outputs@params_table["w","est"]
			return(c(LnL, n_parameters, aicc, dispersal, extinction, j, w))
		}
	
	# apply the above function to extract results from all models, and organize results in a table
		aic_models <- sapply(all_analyses,fetch_BGB_results, n_taxa = n_taxa)
		model_comparison[,1] <- aic_models[1,]
		model_comparison[,2] <- aic_models[2,]
		model_comparison[,3] <- aic_models[3,]
		model_comparison[,4] <- aic.w(aic_models[3,])
		model_comparison[,5] <- aic_models[4,]
		model_comparison[,6] <- aic_models[5,]
		model_comparison[,7] <- aic_models[6,]
		model_comparison[,8] <- aic_models[7,]
	
	
	# outputs table with model comparison
		model_comparison
		write.csv(model_comparison, paste(out_name, "_model_comparison.csv", sep=""))
	
	# collects models without a +J free parameter
		non_J_models <- vector(length=0)
		names_non_J_models <- vector(length=0)
		for (i in 1:length(all_analyses)){
			if(all_analyses[[i]]$inputs$BioGeoBEARS_model_object@params_table["j","type"]=="fixed"){
				non_J_models <- c(non_J_models, i)
				names_non_J_models <- c(names_non_J_models, model_names[[i]])
			}
		}
	
	#recalculates AICw without +J models
		model_comparison_noJ <- matrix(NA, ncol = 8, nrow = length(non_J_models))
		rownames(model_comparison_noJ) <- names_non_J_models
		colnames(model_comparison_noJ) <-  c("LnL", "k", "AICc", "AICw", "d", "e", "j", "w")
		for (i in 1:length(non_J_models)){
			model_comparison_noJ[i,] <- model_comparison[non_J_models[[i]],]
		}
		model_comparison_noJ[,4] <- aic.w(model_comparison_noJ[,3])
		
	# outputs table with model comparison
		model_comparison_noJ
		write.csv(model_comparison_noJ, paste(out_name, "_model_comparison_noJ.csv", sep=""))
	
	# saves all inputs and all outputs
		save(BGB_inputs, file=paste(out_name, "_all_inputs.Rdata", sep =""))
		save(all_analyses, file=paste(out_name, "_all_results.Rdata", sep =""))

sink(file = NULL)
#finished