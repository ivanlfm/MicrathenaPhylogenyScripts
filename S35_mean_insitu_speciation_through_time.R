library(castor)
library(ape)
library(plyr)
library(dplyr)
library(ggplot2)

setwd("C:/a/BGB/BGB_multiple_trees")

#defines areas of interest -- we will obtain estimates of in-situ speciation from them
#also define size of time bins
	all_areas <- c("E", "M", "I", "L", "D", "C") #list all areas in your analysis
	plot_transitions <- c("E", "M") #list only the areas for which you want to plot speciation through time
	size_time_bin <- 1
	only_in_situ_speciation <- TRUE 
	
#reads the same trees used in the analysis
	trees <- read.nexus("trees/100_pruned_trees.nex")
	BSM_transitions <- read.csv("All_transitions_by_node_BSM.csv")

###########################
#gets number of lineages of each tree in each time bin
###########################

	matrix_all_trees <- matrix(NA, ncol = 3, nrow = 0)
	for(i in 1:length(trees)){
		#defines maximum age of tree to collect LTTs
			maximum_age <- as.integer(max(node.depth.edgelength(trees[[i]])))
			times <- seq(from=0, to=maximum_age, by=size_time_bin)
			results <- count_lineages_through_time(trees[[i]], times = times)
		
			Tree = rep(results$Ntimes+1, x=i)
			EndTime = c(rev(results$times)[1]+1, rev(results$times))
			N_lineages = c(1,results$lineages)
		
			results_this_tree <- as.data.frame(cbind(Tree, EndTime, N_lineages))
		
			matrix_all_trees <- rbind(matrix_all_trees,results_this_tree)
		}

	write.csv(matrix_all_trees, file = "LTT_all_trees.csv")

###########################
#post-process the BSM table
###########################

	#adds time bins
	clado_time_bin <- round_any(BSM_transitions$ageStart, size_time_bin, f=floor) #I use floor to round up to smaller integer - so the time bin matchs the LTT matrix
	# if we want anagenetic transitions, we need to get transition times from BSM_transitions$transitionTime, but it is tricky (numbers are strings, and mixed with "-"...)
	# perhaps to ad-hoc subtable only with those anagenetic transitions and coerce them into numbers
	BSM_transitions <- cbind(BSM_transitions, clado_time_bin)

	write.csv(BSM_transitions, file = paste("All_transitions_by_node_BSM_modified.csv",sep=""), row.names = FALSE, fileEncoding="UTF-8")

##########################################################
#creates table with number of events per tree+BSM+age bin
##########################################################

	number_events_per_bin_per_replicate <- matrix(NA, ncol = 5, nrow = 0)
	colnames(number_events_per_bin_per_replicate) <- c("Tree", "EndTime", "N_lineages", "Tree-BSM-Age", "BSM")

	number_BSM <- max(BSM_transitions$BSM)

	#repeats (by the number of BSM replicates) the matrix with number of lineages in each time bin
	for(i in 1:number_BSM){
		matrix_this_rep <- matrix_all_trees
		BSM <- i
		matrix_this_rep <- cbind(matrix_this_rep, BSM)
		number_events_per_bin_per_replicate <- rbind(number_events_per_bin_per_replicate, matrix_this_rep)
	}

	#rearranges by tree and BSM
	number_events_per_bin_per_replicate <- arrange(number_events_per_bin_per_replicate, number_events_per_bin_per_replicate$Tree, number_events_per_bin_per_replicate$BSM)
	number_rows <- length(number_events_per_bin_per_replicate[,1])

	#creates columns for the transitions we are interested in
	#relative columns are for calculating the ratio between the number of events and number of lineages
	events_interest <- matrix(0, ncol = (2*length(all_areas)+2), nrow = number_rows)
	relative_names <- gsub(" ", "", paste(all_areas, "_relative"))
	colnames(events_interest) <- c(all_areas, "Total", relative_names, "Total_relative")
	number_events_per_bin_per_replicate <- cbind(number_events_per_bin_per_replicate, events_interest)
	Tree_previous=0
	
	#populates the table with the corresponding events
	print("This might take a while...:")
	if(only_in_situ_speciation){
				for(b in 1:number_rows){
					Tree <- number_events_per_bin_per_replicate$Tree[b]
					if (Tree != Tree_previous) print(paste((Tree_previous/max(number_events_per_bin_per_replicate$Tree))*100,"%..."))
					Tree_previous=Tree
					BSM <- number_events_per_bin_per_replicate$BSM[b]
					timebin <- number_events_per_bin_per_replicate$EndTime[b]
					number_events_total <- 0
					for(a in 1:length(all_areas)){
						area <- all_areas[a]
						number_events <- length(which(BSM_transitions$replicate==Tree & BSM_transitions$BSMreplicate==BSM & BSM_transitions$clado_time_bin==timebin & BSM_transitions$transitionCladog == paste(area,"->",area)))
						number_events_per_bin_per_replicate[b,(which(colnames(number_events_per_bin_per_replicate)==all_areas[a]))] <- number_events
						number_events_per_bin_per_replicate[b,(which(colnames(number_events_per_bin_per_replicate)==relative_names[a]))] <- number_events / number_events_per_bin_per_replicate$N_lineages[b]
						number_events_total <- number_events_total + number_events
					}
					number_events_per_bin_per_replicate$Total[b] <- number_events_total
					number_events_per_bin_per_replicate$Total_relative[b] <- number_events_total / number_events_per_bin_per_replicate$N_lineages[b]
				}	
			} else {
				for(b in 1:number_rows){
					Tree <- number_events_per_bin_per_replicate$Tree[b]
					if (Tree != Tree_previous) print(paste((Tree_previous/max(number_events_per_bin_per_replicate$Tree))*100,"%..."))
					Tree_previous=Tree
					BSM <- number_events_per_bin_per_replicate$BSM[b]
					timebin <- number_events_per_bin_per_replicate$EndTime[b]
					number_events_total <- 0
					for(a in 1:length(all_areas)){
						area <- all_areas[a]
						number_events <- length(which(BSM_transitions$replicate==Tree & BSM_transitions$BSMreplicate==BSM & BSM_transitions$clado_time_bin==timebin & BSM_transitions$rangeStart == area))
						number_events_per_bin_per_replicate[b,(which(colnames(number_events_per_bin_per_replicate)==all_areas[a]))] <- number_events
						number_events_per_bin_per_replicate[b,(which(colnames(number_events_per_bin_per_replicate)==relative_names[a]))] <- number_events / number_events_per_bin_per_replicate$N_lineages[b]
						number_events_total <- number_events_total + number_events
					}
					number_events_per_bin_per_replicate$Total[b] <- number_events_total
					number_events_per_bin_per_replicate$Total_relative[b] <- number_events_total / number_events_per_bin_per_replicate$N_lineages[b]
				}
			}
					
				write.csv(number_events_per_bin_per_replicate, "number_events_per_bin_per_replicate.csv")
	
####################################################
#summarizes and plots number of events through time
####################################################
	
	oldest_age <- max(number_events_per_bin_per_replicate$EndTime)
	summaries_events <- matrix(NA, ncol = (2*length(all_areas))+4, nrow = oldest_age)
	colnames(summaries_events) <- c("EndTimeBin",gsub(" ", "", paste(all_areas,"_mean")), gsub(" ", "",paste(all_areas,"_SD")), "Total_mean", "Total_SD", "Average_mean")
	summaries_events <- as.data.frame(summaries_events)
	summaries_events$EndTimeBin <- rep(0:(-1*(oldest_age-1)))
	
	for(x in 0:(oldest_age)){
		for(a in 1:length(all_areas)){
			relevant_events <- number_events_per_bin_per_replicate[
									which(number_events_per_bin_per_replicate$EndTime==x),
									which(colnames(number_events_per_bin_per_replicate)==relative_names[a])
								]
			mean <- mean(relevant_events)
			sd <- sd(relevant_events)
			summaries_events[x,which(colnames(summaries_events)==gsub(" ", "", paste(all_areas[a],"_mean")))] <- mean
			summaries_events[x,which(colnames(summaries_events)==gsub(" ", "", paste(all_areas[a],"_SD")))] <- sd
		}
		total <- number_events_per_bin_per_replicate[
									which(number_events_per_bin_per_replicate$EndTime==x),
									which(colnames(number_events_per_bin_per_replicate)=="Total_relative")
								]
								summaries_events[x,which(colnames(summaries_events)==gsub(" ", "", paste(all_areas[a],"_mean")))] <- mean
		summaries_events[x,which(colnames(summaries_events)=="Total_mean")] <- mean(total)
		summaries_events[x,which(colnames(summaries_events)=="Total_SD")] <- sd(total)			
		summaries_events[x,which(colnames(summaries_events)=="Average_mean")] <- mean(total)/length(all_areas)
	}
	
	write.csv(summaries_events, "summaries_events.csv")

###################################################
#plot plot plot plot plot plot plot plot plot plot
###################################################

	ylim = max(summaries_events)
	if (ylim>1) ylim =1

		var1 = plot_transitions[1]
		var2 = plot_transitions[2]
		var3 = "Average of all areas"
		
		plot<-
			ggplot(summaries_events, aes(x = EndTimeBin)) +
			xlab("Age (Myr)") + 
			ylab("Relative in-situ speciation") +
			scale_x_continuous(limits = c(round(-1*oldest_age),0))+
			scale_y_continuous(limits = c(0, ylim))+
			geom_errorbar(aes (ymin = pmax(0, E_mean - E_SD),  ymax = pmin(ylim, E_mean + E_SD), color = "Variable 1"), width = 0, size=5, alpha= 0.2) +
			geom_errorbar(aes (ymin = pmax(0, M_mean - M_SD), ymax = pmin(ylim, M_mean + M_SD), color = "Variable 2"), width = 0, size=5, alpha= 0.2) +
			geom_smooth(se = F, aes(y = E_mean,color = "Variable 1")) +
			geom_smooth(se = F, aes(y = M_mean,color = "Variable 2")) +
			geom_point(aes(y = E_mean, fill = "Variable 1"), size = 4, shape = 24, alpha=0.9) +
			geom_point(aes(y = M_mean, fill = "Variable 2"), size = 4, shape = 21, alpha=0.9) +
			geom_point(aes(y = Average_mean, fill = "Variable 3"), size = 2, shape = 21, alpha=0.8) +
			scale_fill_manual(values = c("Variable 1"="gold", "Variable 2" = "blue", "Variable 3" = "black")) +
			scale_color_manual(values = c("Variable 1"="gold", "Variable 2" = "blue", "Variable 3" = "black")) +
			theme_bw()+
			guides(color = FALSE)+
			guides(fill = guide_legend(name="Area", override.aes = list(labels=c("H","L","Average"),shape = c(24, 21, 21), size = c(4, 4, 2), colors=c("gold", "blue", "black"))))

		plot
		
		ggsave("speciation_through_time.pdf", plot, width = 10, height = 7, units = "in", device = "pdf")