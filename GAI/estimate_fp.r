est_fp_spline = function(yr_data, min_visit = NULL, min_pos = NULL, add_anchor = FALSE, tp_col = "WEEKNO", season_lim = c(1,26), verbose = FALSE, maxit = NULL, site_effect = FALSE, gam_func = "gam"){
	# Data checks and argument checks
	# Check that data contains the expect columns
	if(any(!c("SITE","COUNT",tp_col) %in% names(yr_data))){
		stop("yr_data must contain columns titled SITE,COUNT,and ",tp_col)
	}
	# Check that yr_data only contains data for 1 year and 1 species
	if(length(unique(yr_data$YEAR)) > 1 ){
		stop("yr_data contains data for more than 1 year")
	}
	if(length(unique(yr_data$SPECIES)) > 1){
		stop("yr_data contains data for more than 1 species")
	}
	# Check that there is data
	if(nrow(yr_data) == 0){
		warning("Data frame passed to flight period GAM contains no data (0 rows) returning NULL")
		return(NULL)
	}
	# If min_pos and min_visit suuplied then check that min_pos >= min_visit
	if(!is.null(min_visit) & !is.null(min_pos)){
		if(min_pos > min_visit){
			stop("Minimum number of positive counts cannot be less than the minimum number of visits")
		}
	}
	
	# Clean the yr_data
		yr_data = clean_data(yr_data, tp_col = tp_col, season_lim = season_lim)
	# Double check that there is still data if not the set nm_data to NULL and return
		if(nrow(yr_data) == 0){
			nm_data = NULL
			# Output progress statement
			if(verbose)
				cat("WARNING: No positve counts for this species/year cobmination (exiting and returning NULL)\n")
			return(nm_data)
		}
	# Attempt to load mgcv package and return error if not found
	temp = require(mgcv)
	if(!temp)
		stop("Required package 'mgcv' not found")
	
	# Determine the number of visits and positive counts for each site
		visit_inds = which(yr_data$COUNT >= 0) # A visit should always have a 0 or a positive count
		n_visit = sapply(tapply(yr_data[visit_inds,tp_col],yr_data$SITE[visit_inds], unique),length)
		pos_inds = which(yr_data$COUNT > 0)
		n_pos = sapply(tapply(yr_data[pos_inds,tp_col],yr_data$SITE[pos_inds], unique),length)
		# Cleanup
		rm(visit_inds, pos_inds)
	
	# If min_visit and/or min_pos not null then find sites passing critera
	if(!is.null(min_pos) & !is.null(min_visit)){
		pass_sites = sort(intersect(names(n_pos)[n_pos >= min_pos], names(n_visit)[n_visit >= min_visit]))
		# Output progress statement
		if(verbose)
			cat("Excluding sites that have less than ",min_visit," visits or ",min_pos," positive counts\n", sep="")
	} else if(!is.null(min_pos)){
		pass_sites = names(n_pos)[n_pos >= min_pos]
		# Output progress statement
		if(verbose)
			cat("Excluding sites that have less ",min_pos," positive counts\n", sep="")
	} else if(!is.null(min_visit)){
		pass_sites = names(n_visit)[n_visit >= min_visit]
		# Output progress statement
		if(verbose)
			cat("Excluding sites that have less than ",min_visit," visits\n", sep="")
	}
	
	# Determine which sites are to be excluded based on min_visit and min_pos criteria and then exclude them
	if(exists("pass_sites")){
		# Initialise excluded site vector as all site codes
			exc_sites = unique(yr_data$SITE)
		# Remove any sites codes that are in pass_sites
			exc_sites = sort(exc_sites[which(!exc_sites %in% pass_sites)])
		# Output progress statement (if verbose = TRUE)
		if(verbose){
			if(length(exc_sites) > 0){
				cat("  Excluding ", length(exc_sites), " of ",length(exc_sites) + length(pass_sites)," sites:\n    ", paste(exc_sites,collapse=", "),"\n",sep="")
			} else {
				cat("  No sites excluded as none found that fall below specified thresholds\n", sep="")
			}
		}
		# Now drop these sites from yr_data
			gam_inds = which(yr_data$SITE %in% pass_sites)
			yr_data = yr_data[gam_inds,]
			rm(pass_sites)
	}
	
	# Check that if there is still data
	if(nrow(yr_data) >= 0){
		# Are anchor points to be added (need to be added for figuring out which indices to use for analysis)
			if(add_anchor){
				# Determine anchors days/weeks 
					# If tp is dayno then add anchors of 5 consecutive days on each side of the range where the first day of the anchor points is +/- 10 days from range limits
					if(toupper(tp_col) == "DAYNO"){
						anc_tp = c(season_lim[1]-14:10,season_lim[2]+10:14)
					}
					# If tp is weekno then add anchors to the 5 consecutive weeks on each side of the range
					if(toupper(tp_col) == "WEEKNO"){
						anc_tp = c(season_lim[1]-1:5,season_lim[2]+1:5)
					}
				# Create an anchor dataset (for only pass sites) - only needs columns that are used in fitting GAM
					site_codes = unique(yr_data$SITE)
					anc_data = data.frame(SITE = rep(as.numeric(site_codes), each = 10), ANC_TP = rep(anc_tp, length(site_codes)), COUNT = 0)
					# Rename the ANC_TP columns
						names(anc_data)[2] = tp_col
					rm(site_codes, anc_tp)
			}
		# Fit GAM
			# Output progress statement
			if(verbose)
				cat("Fitting GAM\n",sep="")
			# Determine no of sites in GAM dataset
				n_sites = length(unique(yr_data$SITE))
			# if maxit (maximum number of iterations for GAM convergence) not supplied then use default
			if(is.null(maxit)){
				maxit = gam.control()$maxit
			}
			# Setup model formula using the correct tp column and including a site term if there is more than 1 site
				mod_form = as.formula(paste0("COUNT ~ s(",tp_col,",bs =\"cr\")",ifelse(n_sites > 1 & site_effect, "+factor(SITE)","")))
				# If there is only 1 site and verbose is true then print statement to say that site term has been removed
				if(verbose & n_sites == 1)
					cat("  Only one site so removing site term from GAM model\n", sep="")
			# Now try fitting model using the bam function
			if(add_anchor){
				gam_obj = try(do.call(gam_func,list(formula=mod_form, data = rbind(yr_data[,c("SITE",tp_col,"COUNT")], anc_data), family = poisson(link="log"), control = list(maxit = maxit))), silent = TRUE)
			} else {
				gam_obj = try(do.call(gam_func,list(formula = mod_form, data = yr_data[,c("SITE",tp_col,"COUNT")], family = poisson(link="log"), control = list(maxit = maxit))), silent = TRUE)
			}
			
			# If it hasn't failed check that it has converged again
			if(!"try-error" %in% class(gam_obj)){
				if(gam_obj$iter == maxit | !gam_obj$converged){
					class(gam_obj) = "try-error"
					# Output progress statement
					if(verbose)
						cat("  GAM failed to converge\n",sep="")
				} else {
					# Output progress statement
					if(verbose)
						cat("  GAM fitting successful (",gam_obj$iter," / ",maxit," iterations)\n", sep="")
				}
			}
			
			if(!"try-error" %in% class(gam_obj)){
				# Output progress statement
				if(verbose)
					cat("Estimate NM values\n")
				# If model includes site term the should predict values for each site, normalise individually for each sites after which the values for each site should be the same.
				# However to save time/resources can instead predict for 1st site and normalise that site to give same values but without predicting rest of sites
				# The above also means that the same method can be used regardless of whether the model includes a site term as there is only a single set of WEEKNO/DAYNO values and the SITE column with the 1st site code will only be used if the model includes a site term
					# Setup data to be used to predict
						pred_data = data.frame(SITE = yr_data$SITE[1], tp_col = season_lim[1]:season_lim[2])
						# Rename tp_col colummn
						names(pred_data)[grepl("tp_col",names(pred_data))] = tp_col
						# Output progress statement
						if(verbose)
							cat("  Predicting values from fitted model\n")
				# Now predicted values from model
					pred_data$FITTED = predict(gam_obj,newdata = pred_data, type="response")
				# Output progress statement
				if(verbose)
					cat("  Standardising values to produce NM values\n")
				# Normalise data at each site so that each site the NM values follow the fitted but sum to 1
					pred_data$NM = pred_data$FITTED/sum(pred_data$FITTED)
				# Save to nm_data (NEED TO ADD SPECIES AND YEAR?) and remove pred_data
					nm_data = data.frame(SPECIES = yr_data$SPECIES[1],YEAR = yr_data$YEAR[1],pred_data[,c(tp_col,"NM")])
				# Output progress statement
				if(verbose)
					cat("NM values calculated\n")
					
				
			} else {
				nm_data = NULL
				# Output progress
				if(verbose)
					cat("ERROR: GAM fitting failed (exiting and returning NULL)\n")
			}
		
	} else {
		nm_data = NULL
		# Output progress statement
		if(verbose)
			cat("WARNING: No data for this species/year cobmination (exiting and returning NULL)\n")
	}
	
	return(nm_data)
	
}

est_fp_spline_wrpr = function(yr, sp_data, other_args){
	# Find indices that match current year
		yr_data = sp_data[which(sp_data$YEAR == yr),]
	# Add yr_data to arugments passed using other_args list
		fun_args = c(list(yr_data = yr_data), other_args)
	# Now estimate the flight period for this year using est_fp_spline funciton
		cur_fp = f_do_call("est_fp_spline",fun_args)
}

# A small wrapper function to estimate all yearly flight periods for a species (only useful if manually running some code)
est_fp_spline_all = function(sp_data, ...){
	temp = by(sp_data,sp_data$YEAR,FUN= est_fp_spline, ...)
	nm_data = do.call("rbind",temp)
	row.names(nm_data) = NULL
	rm(temp)
	return(nm_data)
}

est_fp_spline_clust = function(i){
	# Store start_time
		start_time = Sys.time()
	# Deterime whether clust_params object already exists in workspace
	if(!exists("fp_est_clust_params")){
		# Look in current working directory so need to ensure working directory is set prior to this or clust_params exported to each node
		# Load all clust_params
			clust_params = readRDS("fp_est_clust_params.rds")
		# Extract elements for ith run from object and assign as objects in current environment
			list2env(clust_params[i,], environment())
		# remove clust params
			rm(clust_params)
	}
	if(!exists("fp_est_fixed_params")){
		# Look in current working directory so need to ensure working directory is set prior to this or fixed_params exported to each node
			fixed_params = readRDS("fp_est_fixed_params.rds")
		# Assign all the objects in the fixed_params file as objects in the current environment
			list2env(fixed_params, environment())
	}
	# Determine if data_filter was defined in fixed_params.rds if not create and set to null
	if(!exists("data_filter")){
			data_filter = NULL
	}
	# Get the count data
		sp_data = read_data_files(sp = species, yrs = NULL, data_type="count", dpath = paths["data"], data_filter = data_filter)
	# If yrs is not null then restrict data to only the years given in years
		if(!is.null(yrs)){
			sp_data = sp_data[which(sp_data$YEAR %in% yrs),]
		}
	# Determine what years have a positve count (no point attempting to estimate flight period for years with no positive count
		pos_yrs = sort(unique(sp_data$YEAR[which(sp_data$COUNT > 0)]))
		# if no positve years found then exit and return NULL
		if(length(pos_yrs) == 0){
			return(NULL)
		}
		
	# Build a data.frame containing the years/bootstrap combinations that this node will run based on the current batch size specified in batch_nboot (should be defined and passed to the function via fixed_params.rds)
		# If batch_nboot not found then create it and set it to NULL (in which case no bootstraps will be performed and only the flight periods for the non-bootstrapped data will be analysed
		if(!exists("batch_nboot")){
				batch_nboot = NULL
		} else {
			# If batch_nboot has been specified as 0 then set it to NULL (these are equivalent and both only analyses the real non-bootstrapped data)
			if(batch_nboot == 0){
				batch_nboot = NULL
			}
		}
		# Determine bootstrap sets that will be run on this node and build data.frame to hold year/bootstrap combinations
		if(is.null(batch_nboot)){
			# If batch_nboot is null then only need to run real data (iboot = 0) for each year
				node_boots = 0
			# Estimate flight period using est_fp_spline_wrpr wrapper function which call est_fp_spline with correct arguments
				sp_nm = do.call("rbind",lapply(pos_yrs,FUN = est_fp_spline_wrpr, sp_data = sp_data, other_args = other_args))
		} else {
			# When iboot_start == 0 then need to run 1 set on real analysis and then batch_nboot number of bootstraps (real results will be saved seperately so it makes more sense if each boostrap output file contains the same number of results sets rather than having the first contain 1 less
			if(iboot_start == 0){
				node_boots = seq(from = iboot_start, by = 1, length.out = batch_nboot+1)
			} else {
				node_boots = seq(from = iboot_start, by = 1, length.out = batch_nboot)
			}
			# Check that the bootstraps required from the current node doesn't exceed the max number of bootstraps requested if so truncate
			if(is.null(nboot)){
				nboot = 0
			}
			if(max(node_boots) > nboot){
				node_boots = node_boots[which(node_boots <= nboot)]
			}
			# Read in boot_sample file containing all site codes for each boot_sample (checking that the file exists if not set boot_samp to NULL)
			if(file.exists(file.path(paths["boot_samp"],paste0("Sp ",species,"_bootsamps.rds")))){
				boot_samp = readRDS(file.path(paths["boot_samp"],paste0("Sp ",species,"_bootsamps.rds")))
				
				# Only keep those for current bootstraps
					boot_samp = boot_samp[,c("site",paste("iboot",node_boots[which(node_boots > 0)],sep="_"))]
			} else {
				boot_samp = NULL
			}
				
			# Look in current workspace and build list of arguments that match arguments in est_fp_spline function
				other_args = mget(ls()[which(ls() %in% names(formals(est_fp_spline)))])
				
			# Setup object to hold the output data
				boot_nm = vector("list",length(node_boots[which(node_boots > 0)]))
			# Loop through bootstraps each time resampling the data and then estimating the flight period
			for(cur_boot in node_boots){
				if(cur_boot == 0){
					# Estimate flight period using est_fp_spline_wrpr wrapper function which call est_fp_spline with correct arguments
						sp_nm = do.call("rbind",lapply(pos_yrs,FUN = est_fp_spline_wrpr, sp_data = sp_data, other_args = other_args))
					# 
				} else {
					if(!is.null(boot_samp)){
						boot_inds = which(node_boots[which(node_boots > 0)] == cur_boot)
						# Produce bootstrap dataset using species codes in boot_samp for current bootstrap, renaming columns in boot_samp so that site becomes SITE (new site codes, integers 1:nsites) and iboot_x becomes ORG_SITE 
							boot_data = boot_sample(sp_data,boot_samp = setNames(boot_samp[,c("site",paste("iboot",cur_boot,sep="_"))],c("SITE","ORG_SITE")))
						# Determine which years are positive for the bootstrap data sample
							boot_pos_yrs = sort(unique(boot_data$YEAR[which(boot_data$COUNT > 0)]))
						if(length(boot_pos_yrs) > 0){
							# Estimate flight period using est_fp_spline_wrpr wrapper function which call est_fp_spline with correct arguments
								boot_nm[[boot_inds]] = do.call("rbind",lapply(boot_pos_yrs,FUN = est_fp_spline_wrpr, sp_data = boot_data, other_args = other_args))
							# Add column with current bootstrap number to output but only if the returned object isn't null
							if(!is.null(boot_nm[[boot_inds]])){
								boot_nm[[boot_inds]][,"iboot"] = cur_boot
							}
						}
					}	
				}
			}
			
			# Collapse boot_nm list to form a single dataset
				boot_nm = do.call("rbind",boot_nm)
		}
		
		# Save outputs to file
			# Flight periods from real data
			if(exists("sp_nm")){
				saveRDS(sp_nm, file=file.path(paths["nm"],paste0("Sp ",species, "_nm.rds")))
				# if sp_nm exists then also return it (done for now just to make 
			}
			
			# Flight periods from 
			if(exists("boot_nm")){
				if(ifelse(is.null(boot_nm),0,nrow(boot_nm)) > 0){
					f_name = paste0("Sp ",species,"_iboot ", paste(range(node_boots[which(node_boots > 0)]), collapse=" to "),"_nm.rds")
					saveRDS(boot_nm, file=file.path(paths["boot_nm"],f_name))
				}
			}
}

est_fp_spline_clust_old = function(i){
	# Store start_time
		start_time = Sys.time()
	# Deterime whether clust_params object already exists in workspace
	if(!exists("fp_est_clust_params")){
		# Look in current working directory so need to ensure working directory is set prior to this or clust_params exported to each node
		# Load all clust_params
			clust_params = readRDS("fp_est_clust_params.rds")
		# Extract elements for ith run from object and assign as objects in current environment
			list2env(clust_params[i,], environment())
		# remove clust params
			rm(clust_params)
	}
	if(!exists("fp_est_fixed_params")){
		# Look in current working directory so need to ensure working directory is set prior to this or fixed_params exported to each node
			fixed_params = readRDS("fp_est_fixed_params.rds")
		# Assign all the objects in the fixed_params file as objects in the current environment
			list2env(fixed_params, environment())
	}
		
	# Setup log file and redirect normal output (i.e. cat/print) and also other for warnings/errors messages) for the species using sink to direct the output there
		# Block sending output to log file if the function is flagged for debugging
		if(!(isdebugged(est_fp_spline) | isdebugged(est_fp_spline_clust))){
			# Set options so warnings appear at time they are generated
				options(warn = 1)
			# Setup a file connection in write mode
				fcon_log = file(file.path(paths["log"],paste0("Sp ",species,"_Yr ", year, ifelse(exists("iboot"), paste0("_iboot ",iboot),""),"_fp_est_log.txt")), open = "w")
			# Use sink to redirect standard output and error messages to this file
				sink(file = fcon_log, type = "output", split = TRUE)
				sink(file = fcon_log, type = "message", append = TRUE)
			# On exist close the redirects via sink and then close the file connection
			on.exit({
				sink(NULL, type="output")
				sink(NULL, type="message")
				log_path = summary(fcon_log)$description
				close(fcon_log)
				if(!is.null(sp_nm)){
					file.remove(log_path)
				} else {
					file.rename(log_path, gsub("[.]txt$","_complete.txt",log_path))
				}
			})
		}
		
	# Print header
		cat(rep("-",80),"\nSTARTING FIRST STAGE OF GAI ANALYSIS\nFLIGHT PERIOD ESTIMATION FOR SPECIES = ",species, " YEAR = ", year, ifelse(exists("iboot"), paste0(" BOOTSTRAP = ",iboot),""), "\nANALYSIS STARTED: ",format(start_time,"%d/%m/%Y %T"),"\n",rep("-",80),"\n\n",sep="")
	
	# Get the data
		# Check if there is a data filter (if not add object and set to null		
		if(!exists("data_filter")){
			data_filter = NULL
		}
		yr_data = read_data_files(sp = species, data_type="count", yrs = year, dpath = paths["data"], data_filter = data_filter)
		if(is.null(yr_data)){
			cat("  WARNING: no count data found for given species/year combination")
		}
	# Deterine if year data has any positive counts
		yr_tot = sum(yr_data$COUNT)
		if(yr_tot <= 0){
			cat("  WARNING: no postive counts found for this year\n")
			yr_data = NULL
		}
	# Setup variable to hold path under which results will be saved
		out_path = file.path(paths["nm"],paste0("Sp ",species, "_Yr ",year,"_nm.rds"))
	# Determine if this is to be a bootstrap run or not
	if(exists("iboot")){
		if(iboot != 0){
			# Read in boot_sample file containing all site codes for each boot_sample
				boot_samp = readRDS(file.path(paths["boot_samp"],paste0("Sp ",species,"_bootsamps.rds")))
				# Only keep those for current bootstrap
					boot_samp = boot_samp[,c("SITE",paste("iboot",iboot,sep="_"))]
				# Rename iboot site codes to ORG_SITE
					names(boot_samp)[2] = "ORG_SITE"
			# Apply the bootstrap sample codes to the yr data to create a bootstrapped version of the year data
				yr_data = boot_sample(yr_data, boot_samp)
			# Modify the output file path so bootstrap results kept seperate
				out_path = file.path(paths["boot_nm"],paste0("Sp ",species, "_Yr ",year,"_iboot ",iboot,"_nm.rds"))
		}
	}
	# If data found then continue
	if(!is.null(yr_data)){
		if(nrow(yr_data) > 0){
			# Get list of formals for fp fitting function
				fun_args = names(formals("est_fp_spline"))
				# Reduce to only those that exist in workspace
					fun_args = fun_args[fun_args %in% ls()]
			# Attempt to fit the spline using mget to create a list of the arguments that are in thw workspace
				sp_nm = f_do_call("est_fp_spline", mget(fun_args))
			# If sp_nm fitting return values and iboot exists then add iboot value to sp_nm data.frame
			if(!is.null(sp_nm)){
				if(exists("iboot")){
					sp_nm[,"iboot"] = iboot
				}
			}
			# Save the data for the species
				saveRDS(sp_nm, file = out_path)
		} else {
			sp_nm = NULL
			cat("  WARNING: no yr_data remained after cleaning\n")
		}
	} else {
		sp_nm = NULL
	}
	
	# Print summary information and footer to log file
		# Store end time
			end_time = Sys.time()	
		# Calculate time taken
			run_time = end_time - start_time
		# Determine memory usage (use gc on unix but memory.size on windows)
			if(.Platform$OS == "unix"){
				mem_raw = gc()
				mem_use = apply(mem_raw[,c(2,4)],2,sum)
			} else {
				mem_use = c(memory.size(), memory.size(max = TRUE))
			}
		# Print footer to log file and/or screen
			cat("\n",rep("-",80),"\nANALYSIS COMPLETED: ",format(end_time,"%d/%m/%Y %T"),"\nTIME TAKEN = ",format(run_time),"\nMEMORY LIMITS = ", mem_use[1],"Mb (",mem_use[2],"Mb Max)\n",rep("-",80),"\n",sep="")
			
	return(sp_nm)
}

fit_fp_splines = function(dpath_data, dpath_out = NULL, spp = NULL, yrs = NULL, parallel = TRUE, parallel_pkg = "parallel", n_cpus = NULL, nboot = NULL, batch_nboot = 10, iboot_initial = 0, ...){
	# Check that dpath_data exits
	if(!file.exists(dpath_data)){
		stop("Supplied directory for data files could not be found")
	} else {	
		# Normalise dpath_data to ensure file path is absolute not relative
			dpath_data = normalizePath(dpath_data, winslash = "/")
	}
	# Check the path that contains the species data to ensure that it contains count data files
	f_list = list.files(dpath_data, pattern = "count[.]rds$")
	if(length(f_list) == 0){
		stop("No data files found in supplied directory")
	}
	# Determine whether the species data uses seperate files for each year or if all the data for a species is in a single file
	if(all(!grepl("_Yr [[:digit:]]{4}",f_list))){
		yr_files = FALSE
	} else {
		yr_files = TRUE
	}
		
	# if dpath_out not supplied then setup directory to hold data (will be created in next step if it doesn't already exist
	if(is.null(dpath_out)){
		dpath_out = file.path(getwd(),paste0("Results_spline_NM_",format(Sys.time(),"%Y_%m_%d_%H%M%S")))
	} else {
		# Check whether output folder already exists (if not create it)
		if(!file.exists(dpath_out)){
			temp = dir.create(dpath_out)
			if(!temp){
				if(!file.exists(dirname(dpath_out))){
					stop("ERROR: could not create output directory. Issue with path prior to output directory, check full path to ensure all directories (with exception of output directory) exist and are properly named")
				} else {
					stop("ERROR: could not create output directory. Issue creating output directory, check name supplied in path for output directory is valid")
				}
			}
		}
		# Normalise input path to ensure paths are absolute and not relative
		dpath_out = normalizePath(dpath_out, winslash = "/")
	}
	# Combine dpaths into a paths vector (and add paths for NM results and Logs
		paths = c(data = dpath_data, results = dpath_out, nm = file.path(dpath_out,"NM"), log = file.path(dpath_out,"Logs"))
		# Add a path for storing bootstrap files
		if(!is.null(nboot)){
			paths = c(paths, boot = file.path(dpath_out,"Boot"), boot_samp = file.path(dpath_out,"Boot","samp"), boot_nm = file.path(dpath_out,"Boot","NM"))
		}
		# Check to see whether these directories already exist
		d_test = file.exists(paths)
		if(any(!d_test)){
			for(i in which(!d_test)){
				stopifnot(dir.create(paths[i]))
			}
		}
		# Cleanup old pathways
		rm(dpath_data,dpath_out,d_test)
	# Change working directory to results folder
		# Determine original working directory (will change back to it on exist)
			org_wd = getwd()
			setwd(paths["results"])
			on.exit(setwd(org_wd))
	# Save call
		temp = match.call()
		capture.output(print(temp),file="fp_est_call.txt")
		rm(temp)
	# Determine species from files
		# If data files are for each species/year then use names otherwise will have to look for all_sp_yr_combs.rds or read files and build it 
			spp_codes = sort(as.numeric(unique(gsub("^Sp ([[:digit:]]{1,})(_Yr ([[:digit:]]{4}))?_.*$","\\1",f_list))))
		# If spp or yrs is not null then filter combs to only keep combs with supplied spp and/or yrs
		if( !is.null(spp) ){
			# print warning of any species codes given in spp which did not have data file(s)
			miss_inds = which(!spp %in% spp_codes)
			if(length(miss_inds) > 0){
				cat("WARNING: the following species codes supplied to spp argument did not have data files\n\t",paste(spp[miss_inds],collapse=","))
			}
			spp_codes = spp_codes[which(spp_codes %in% spp)]
		}
		# If bootstrapping required determine bootstrap numbers at which each node will start and then construct data.frame of all species/bootstrap start combinations
		if(!is.null(nboot)){
			# Determine bootstraps numbers at which each node will start (starting at iboot_initial and progressing in batch_nboot jumps until nboot)
				boot_seq = seq(from = 1, to = nboot, by = batch_nboot)
				# Exception when iboot_initial == 0 (the default) where first run will actually have batch_nboot+1 runs (1 for the real data + then batch_nboot bootstrap runs)
				if(iboot_initial == 0){
					boot_seq[1] = 0
				}
				
			# Construct data.frame of combinations of species code and bootstrap starts
				clust_combs = expand.grid(species = spp_codes, iboot_start = boot_seq, KEEP.OUT.ATTRS = FALSE)
			
		} else {
			# If bootstrapping not required then simply build data.frame of species codes setting iboot_start = 0 for them all
				clust_combs = data.frame(species = spp_codes, iboot_start = 0)
		}
		# TODO - ADD OPTION OF CALCULATING FP SEPERATELY FOR DIFFERENT SUBSETS OF DATA VIA A GROUPING COLUMM
		# Save clust_combs as RDS file
			saveRDS(clust_combs, file=file.path(paths["results"],"fp_est_clust_params.rds"))
	# Setup object to store the variables that will be fixed across all runs (e.g. file paths, etc)
		# Determine what arguments were passed to function via dots
			dot_args = list(...)
			dot_argnames = names(dot_args)
		# Determine formals for function (use base function as clust wrapper has only 1 arugment i
			form_args = c(names(formals("est_fp_spline")),"data_filter")
			# Exlcude yr_data which is loaded by est_fp_spline_clust using species, year arguments already in clust_params
			form_args = form_args[-grepl("yr_data",form_args)]
		# Create list of fun_args from dot_args that match formals plus any other args that need adding (i.e. paths)
			fun_args = c(dot_args[which(dot_argnames %in% form_args)],list(paths = paths, yrs = yrs, nboot = nboot, batch_nboot = batch_nboot))
		# Save the other fixed variables to an RDS file
			saveRDS(fun_args, file=file.path(paths["results"], "fp_est_fixed_params.rds"))
	
	# Now run the est_fp_spline_clust function for all sets in the list
	if(parallel){
		# Setup a vector containing all objects functions that will need exported to the cluster nodes
		exp_objs = c("read_data_files","filter_data","clean_data","est_fp_spline","est_fp_spline_clust","est_fp_spline_wrpr","f_do_call","create_boot_samples","boot_sample")
		if(parallel_pkg == "parallel"){
			# Load parallel library
			require(parallel)
			if(is.null(n_cpus)){
				# Determine number of cpus to use if not specified
					n_l_cores = detectCores()
					n_p_cores = detectCores(logical = FALSE)
					n_cpus = n_l_cores - (n_l_cores/n_p_cores)
			}
			# Initialise cluster (and setup on.exit to close the cluster on exit)
				cl = makeCluster(n_cpus)
				on.exit(stopCluster(cl), add = TRUE)
			# Load required packages onto nodes
				clusterEvalQ(cl,library(mgcv))
			# Export required objects/functions to clusters
				clusterExport(cl, varlist = exp_objs)
			# Use cluster version of apply to produce bootstrap sample codes for each species
			if(!is.null(nboot)){
				parLapplyLB(cl,spp_codes,create_boot_samples, yrs = yrs, dpath_data = paths["data"], dpath_out = paths["boot_samp"], nboot = nboot)
			}
			# Use cluster version of apply to run fitting for all combinations
				parLapplyLB(cl,1:nrow(clust_combs),est_fp_spline_clust)
		}
		if(parallel_pkg == "snowfall"){
			# Load snowfall library
			require(snowfall)
			# Initialise cluster
			if(.Platform$OS == "unix"){
				hosts = as.character(read.table(Sys.getenv('PBS_NODEFILE'),header=FALSE)[,1]) # read the nodes to use
				sfSetMaxCPUs(length(hosts)) # ensure that snowfall can cope with this many hosts
				sfInit(parallel=TRUE,type="MPI",cpus=length(hosts),useRscript=TRUE) # initialise the connection
			} else {
				if(is.null(n_cpus)){
					# Determine number of cpus to use if not specified
						n_l_cores = parallel::detectCores()
						n_p_cores = parallel::detectCores(logical = FALSE)
						n_cpus = n_l_cores - (n_l_cores/n_p_cores)
				}
				sfInit(parallel = TRUE, type="SOCK", cpus = n_cpus, useRscript=TRUE)
			}
			# Ensure that cluster is stopped when this function exits
				on.exit(sfStop(), add = TRUE)
			# Load required libraries on nodes
				sfLibrary(mgcv)
			# Export required variables and functions to nodes
				sfExport(list = exp_objs)
			# Create bootstrap samples of all sites codes for each species
			if(!is.null(nboot)){
					sfClusterApplyLB(spp_codes,create_boot_samples, yrs = yrs, dpath_data = paths["data"], dpath_out = paths["boot_samp"], nboot = nboot)
			}
			# Use load balancing cluster version of apply to run fitting for all combinations
				sfClusterApplyLB(1:nrow(clust_combs), est_fp_spline_clust)
		}
	} else {
		# Create bootstrap samples of all sites codes for each species
		if(!is.null(nboot)){
				lapply(spp_codes,create_boot_samples, yrs = yrs, dpath_data = paths["data"], dpath_out = paths["boot_samp"], nboot = nboot)
		}
		# If not running in parallel use standar lapply to run analysis for all combinations
			lapply(1:nrow(clust_combs),est_fp_spline_clust)
		cat("Flight Periods Estimation Complete\n\n")
	}
}

UKBMS_fit_fp_splines = function(dpath_data = "~/UKBMS data", dpath_out = NULL, spp = NULL, yrs = NULL, extract_data = FALSE, parallel = TRUE, parallel_pkg = "snowfall", n_cpus = NULL, nboot = NULL, boot_only = FALSE, channel = NULL, i_start = 1, ... ){
	# If extract_data is TRUE then the function is to extract the data from the UKBMS and save it to the default location
	if(extract_data){
		# First remove all species data from that folder
			# Get list of all files that start Sp and end .rds (to ensure no other files can be deleted by accident)
				f_list = list.files(dpath_data, "^Sp .*[.]rds$", full.names = TRUE)
			# Now remove all these files
			if(length(f_list) > 0){
				file.remove(f_list)
			}
		# Now Extract data from UKBMS database and create rds files
		temp = f_do_call("create_UKBMS_files",c(list(spp_list = spp, years = yrs, dpath = dpath_data, channel = channel), list(...)))
		#if("verbose" %in% names(list(...))){
		#	temp = create_UKBMS_files(spp_list = spp, years = yrs, dpath = dpath_data, channel = channel, verbose = list(...)[["verbose"]])
		#} else {
		#	temp = create_UKBMS_files(spp_list = spp, years = yrs, dpath = dpath_data, channel = channel)
		#}
	}
	
	# Fit splines
	nm_vals = fit_fp_splines(dpath_data = dpath_data, dpath_out = dpath_out, spp = spp, yrs = yrs, parallel = parallel, parallel_pkg = parallel_pkg, n_cpus = n_cpus, nboot = nboot, boot_only = boot_only, i_start = i_start, ...)
	
	return(nm_vals)
}