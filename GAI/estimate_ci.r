# Function to find the nearest value to a value(s) from a vector of possible values
find_nearest = function(x, vec){
	# Check x and vec are numeric
	# Convert x to numeric if factor or character (assuming these can be safely converted to numbers)
	if(class(x) == "factor"){
		x = as.numeric(as.character(x))
	} else if(class(x) == "character"){
		x = as.numeric(x)
	}
	# Convert vec to numeric if factor or character (assuming these can be safely converted to numbers)
	if(class(vec) == "factor"){
		vec = as.numeric(as.character(vec))
	} else if(class(x) == "character"){
		vec = as.numeric(vec)
	}
	
	ret_obj = rep(NA, length(x))
	for(i in 1:length(x)){
		ret_obj[i] = vec[which.min(abs(vec-x[i]))]
	}
	return(ret_obj)
}

miss_nm_yrs = function(sp_data, nm_data){
	data_yrs = sort(unique(sp_data$YEAR))
	nm_yrs = sort(unique(nm_data$YEAR))
	miss_yrs = data_yrs[which(!data_yrs %in% nm_yrs)]
	if(length(miss_yrs) > 0){
		near_yrs = find_nearest(miss_yrs,nm_yrs)
		out_obj = data.frame(MISSING = miss_yrs, USE = near_yrs)
	} else {
		out_obj = NULL
	}
		return(out_obj)
}

nm_fill_miss = function(sp_data, si_data, nm_data){
	# Determine if the NM data is missing any years for which there is data (excluding site index only data) and if so also return which is the closest year with data
		miss_yrs = miss_nm_yrs(sp_data, nm_data)
	# If missing years found then create NM values for these years by using the nearest year if not then do nothing and return nm_data unchanged
	if(!is.null(miss_yrs)){
		# Print warning to console
		cat("NM values missing for years present in actual data\n",paste(sprintf("  No NM values for %s so using values for %s",miss_yrs[,1], miss_yrs[,2]),collspe="\n"),sep="")
		# Get indices of values in NM which are to be used to replace missing years
			add_inds = lapply(miss_yrs$USE,function(x,y){ which(y  ==  x) },nm_data$YEAR)
		# Extract the NM data for years that are to replace missing years
			add_nm = nm_data[unlist(add_inds),]
		# Alter the year in the replacement to be the year(s) that were missing 
			add_nm$YEAR = rep(miss_yrs$MISSING,sapply(add_inds,length))
		# Add the replacement data into the nm_data
			nm_data = rbind(nm_data, add_nm)
	}
	return(nm_data)
}

gai_si = function(sp_data, nm_data, site_yr_sum, tp_col = "WEEKNO"){
	# Determine number of tp vals
		all_yrs = sort(unique(site_yr_sum$YEAR))
		all_tps = sort(unique(nm_data[,tp_col]))
	# Setup data to hold sindex values
		si_list = vector("list",length(all_yrs))
	# Loop through by year
	for(i_yr in 1:length(all_yrs)){
		# Get subset of data from site_yr_sum for current year
			temp = site_yr_sum[which(site_yr_sum$YEAR == all_yrs[i_yr]),]
		# Get nm_data for current year
			yr_nm = nm_data[which(nm_data$YEAR == all_yrs[i_yr]),]
			# Check that nm is order by weekno
			if(any(yr_nm[,tp_col] != all_tps)){
				yr_nm = yr_nm[order(yr_nm[,tp_col]),]
			}
		# Create an expanded version of data to include a value for each dayno/weekno for each site/year entry (adding extra columns that will be needed later)
			exp_yr_data = data.frame(temp[rep(1:nrow(temp),each = length(all_tps)),c("SITE","N_EST")], TP_COL = rep(all_tps,nrow(temp)), COUNT = NA, NM = rep(yr_nm$NM,nrow(temp)), IMP = NA)
			# Rename tp col
				names(exp_yr_data)[grepl("TP_COL",names(exp_yr_data))] = tp_col
			# Remove temp
				rm(temp)
			# Set row names
				row.names(exp_yr_data) = apply(exp_yr_data[,c("SITE",tp_col)],1,paste, collapse="_")
		# Add count data
			yr_inds = which(sp_data$YEAR == all_yrs[i_yr])
			exp_yr_data[apply(sp_data[yr_inds,c("SITE",tp_col)],1,paste,collapse="_"),"COUNT"] = sp_data$COUNT[yr_inds]
		# Calculate imputted (count where exists or if NA then N_EST * nm)
			exp_yr_data$IMP = ifelse(!is.na(exp_yr_data$COUNT),exp_yr_data$COUNT, exp_yr_data$N_EST * exp_yr_data$NM)
		# Calculate the site index values for the current year
			if("ORG_SITE" %in% names(si_data)){
				si_list[[i_yr]] = data.frame(YEAR = all_yrs[i_yr], trap_index(exp_yr_data,data_col="IMP",time_col="WEEKNO",by_cols = c("SITE","ORG_SITE")))
			} else {
				si_list[[i_yr]] = data.frame(YEAR = all_yrs[i_yr], trap_index(exp_yr_data,data_col="IMP",time_col="WEEKNO",by_cols = "SITE"))
			}
		# Cleanup
			rm(exp_yr_data, yr_inds, yr_nm)
	}
	# Combine to data.frame
		si_out = do.call("rbind",si_list)
	# Return si values
		return(si_out)
}

pred_SI_glm = function(glm_obj, site_yr_sum){
	# Get vectors of unique sites and years in site_yr_sum data.frame
		sites = unique(site_yr_sum$SITE)
		years = unique(site_yr_sum$YEAR)
	# Build an expanded version of site_yr_sum table that includes all combinations of site and year
		full_site_yr_sum = expand.grid(SITE = sites,YEAR = years)
	# If original site_yr_sum included ORG_SITE column then it is a bootstrap and the original sites codes need to also be retained
	if("ORG_SITE" %in% names(site_yr_sum)){
		combs = unique(site_yr_sum[,c("SITE","ORG_SITE")])
		full_site_yr_sum = merge(full_site_yr_sum,combs, by = "SITE", sort = FALSE, all.x = TRUE)
	}
	# Predict the N_EST values for each site/year combination back from the fitted model
		full_site_yr_sum$PRED_N_EST = predict(glm_obj, newdata = full_site_yr_sum, type="response")
	# Create a PRED_SINDEX column by rounding predicted N_EST values
		full_site_yr_sum$PRED_SINDEX = round(full_site_yr_sum$PRED_N_EST)
	# Merge the full expanded version of site_yr_sum with the original version to have a copy that has real N_EST/SINDEX values as well as predicted ones
		out_site_yr_sum = merge(site_yr_sum, full_site_yr_sum, all = TRUE, sort = TRUE)
	# Return new site_yr_sum
	return(out_site_yr_sum)
}


gai_coll_ind = function(sp_data, nm_data, si_data = NULL, tp_col = "WEEKNO", methodCI = "GLM_N", sindex_fname = NULL, fill_miss_nm = TRUE, min_sites = NULL, min_sites_pos = NULL, zero_fp_thres = NULL, season_lim = c(1,26), iboot = NULL){
	# Check data
		# Only data for 1 species
		if(length(unique(sp_data$SPECIES)) > 1){
			stop("sp_data contains data for more than 1 species")
		}
		if(length(unique(nm_data$SPECIES)) > 1){
			stop("nm_data contains data for more than 1 species")
		}
		# Data has required/execpted columns
		req_cols = c("SPECIES","YEAR","SITE","COUNT",tp_col)
		if(any(!req_cols %in% names(sp_data))){
			stop("sp_data must contain the following column: ",paste(req_cols, collapse=", "))
		}
		req_cols = c("SPECIES","YEAR","NM",tp_col)
		if(any(!req_cols %in% names(nm_data))){
			stop("nm_data must contain the following column: ",paste(req_cols, collapse=", "))
		}
		req_cols = c("GLM_N","GLM_T","SIMPLE_N")
		if(!methodCI %in% req_cols ){
			stop("unrecognised value supplied for methodCI, please use one of the following: ",paste(req_cols, collapse=", "))
		}
		# Cleanup
		rm(req_cols)
		
	# Define year summary function that will be used later
	calc_yr_summ = function(site_yr_sum){
		# Calculate summary for each year
		summ_out = do.call("rbind", by(site_yr_sum, site_yr_sum$YEAR, 
			function(x){
				ret_obj = data.frame(YEAR = unique(x$YEAR), TOTAL = sum(x$TOTAL), NSITES = length(unique(x$SITE)), NSITES_POS = length(unique(x$SITE[which(x$TOTAL > 0)])) )
				return(ret_obj)
			})
		)
		row.names(summ_out) = NULL
		return(summ_out)
	}
	
	# Clean species data
		if("ORG_SITE" %in% names(sp_data)){
			sp_data = clean_data(sp_data, tp_col = tp_col, season_lim = season_lim, add_cols = "ORG_SITE")
		} else {
			sp_data = clean_data(sp_data, tp_col = tp_col, season_lim = season_lim)
		}
		
	# If fill_miss_nm is TRUE then look for missing years and if so use the nm values for closest year
	if(fill_miss_nm){
		nm_data = nm_fill_miss(sp_data = sp_data, si_data = si_data, nm_data = nm_data)
	}
		
	# Look for any years with no positive count
		# if(!is.null(si_data)){
			# yr_sum = tapply(c(sp_data$COUNT,si_data$SINDEX), c(sp_data$YEAR,si_data$YEAR), sum, na.rm = TRUE)
		# } else {
			# yr_sum = tapply(sp_data$COUNT, sp_data$YEAR, sum, na.rm = TRUE)
		# }
		# if(any(yr_sum == 0)){
			# bad_yrs = names(yr_sum)[which(yr_sum == 0)]
			# cat("Warning: data includes the following years where there is no postive count (years will be excluded):\n")
			# cat("  ",paste(bad_yrs,collapse=", "),"\n")
			# # Exclude from sp_data
			# rm_inds = which(sp_data$YEAR %in% bad_yrs)
			# if(length(rm_inds) > 0){
				# sp_data = sp_data[-rm_inds,]
			# }
			# # Exclude from nm_data
			# rm_inds = which(nm_data$YEAR %in% bad_yrs)
			# if(length(rm_inds) > 0){
				# nm_data = nm_data[-rm_inds,]
			# }
			# # Exclude from si_data (if present)
			# if(!is.null(si_data)){
				# rm_inds = which(si_data$YEAR %in% bad_yrs)
				# if(length(rm_inds) > 0){
					# si_data = si_data[-rm_inds,]
				# }
			# }
		# }
		
	# Estimate N_i based on concentrated likelihood approach (common across all methods)
		# Add row.names to nm_data so that the row names of each is the cocatination of the year and tp_col seperated by an underscore
			row.names(nm_data) = apply(nm_data[,c("YEAR",tp_col)],1,paste, collapse="_")
		# Setup data.frame to hold site/year total counts and sum of FP corresponding to visits at same time populate with total counts (note names of columns in list elements in aggregate will be row.names in produced data.frame)
			site_yr_sum = data.frame(aggregate(list(TOTAL = sp_data$COUNT), sp_data[,c("SITE","YEAR")],sum, na.rm = TRUE),FP_SUM = NA, N_EST = NA, N_POS = NA, CNT_TRAP = NA, INDEX_ONLY = NA)
			
			# If sp_data contains org_site column then run is a bootstrap and the org_site info should be added to site_yr_sum (but keeping same order, hence not using merge!)
			if("ORG_SITE" %in% names(sp_data)){
				combs = unique(sp_data[,c("SITE","ORG_SITE")])
				row.names(combs) = combs$SITE
				site_yr_sum[,"ORG_SITE"] = combs[as.character(site_yr_sum$SITE),"ORG_SITE"] 
			}
			
		# Get sum of fp corresponding to weeks with counts
			temp_rn = apply(sp_data[,c("YEAR",tp_col)],1, paste,collapse="_")
			site_yr_sum[,"FP_SUM"] = aggregate(list(FP_SUM = nm_data[temp_rn,"NM"]), list(SITE = sp_data$SITE, YEAR = sp_data$YEAR),sum, na.rm = TRUE)[,"FP_SUM"]
		# Determine number of visits with positive count for each species year combination
			site_yr_sum[,"N_POS"] = aggregate(list(N_POS = sp_data$COUNT), list(SITE = sp_data$SITE, YEAR = sp_data$YEAR),function(x){ x = x[which(x > 0)]; length(x)})[,"N_POS"]
		# Determine the trapazoidal area of the raw count data
			site_yr_sum[,"CNT_TRAP"] = trap_index(sp_data,data_col="COUNT",time_col=tp_col,by_cols = c("SITE","YEAR"))[,"SINDEX"]
			
			# Cleanup
				rm(temp_rn)
				
		# Calculate N_i
		# If zero_fp_thres specified then find all all cases where the are no positive counts and the estimated flight period coverage is less than the threshold then exclude from N_i calculation
		if(!is.null(zero_fp_thres)){
			thres_fail = which(site_yr_sum$TOTAL == 0 & site_yr_sum$FP_SUM < zero_fp_thres)
			# Calculate N_i
			if(length(thres_fail) > 0){
				site_yr_sum$N_EST[-thres_fail] = site_yr_sum$TOTAL[-thres_fail] / site_yr_sum$FP_SUM[-thres_fail]
			} else {
				site_yr_sum$N_EST = site_yr_sum$TOTAL / site_yr_sum$FP_SUM
			}
		} else {
			# Calculate N_i
				site_yr_sum$N_EST = site_yr_sum$TOTAL / site_yr_sum$FP_SUM
		}
		# Show progress
			cat("N_EST values calculated\n")
	
	
	if(methodCI == "GLM_N" | methodCI == "GLM_T"){
		if(methodCI == "GLM_N"){
			# GLM_N method
			#-------------
			# Show progress
				cat("Fitting GAI - GLM_N method\n")
			# If there is site index only data then add it on to site_yr_sum data
			if(ifelse(!is.null(si_data),nrow(si_data),0) > 0){
				if("ORG_SITE" %in% names(si_data)){
					site_yr_sum = rbind(site_yr_sum, data.frame(SITE = si_data$SITE, YEAR = si_data$YEAR, ORG_SITE = si_data$ORG_SITE, TOTAL = si_data$SINDEX, FP_SUM = 1, N_EST = si_data$SINDEX, N_POS = NA, CNT_TRAP = NA, INDEX_ONLY = 1))
				} else {
					site_yr_sum = rbind(site_yr_sum, data.frame(SITE = si_data$SITE, YEAR = si_data$YEAR, TOTAL = si_data$SINDEX, FP_SUM = 1, N_EST = si_data$SINDEX, N_POS = NA, CNT_TRAP = NA, INDEX_ONLY = 1))
				}
				# Show progress
				cat("Added Site Index only data to N_EST values estimated from counts\n")
			}
			# Calculate year summary values
				yr_summary = calc_yr_summ(site_yr_sum)
			# Determine if min_sites and/or min_sites_pos specified and if so then limit data to only years that pass thresholds
			if(!is.null(min_sites) | !is.null(min_sites_pos)){
				if(!is.null(min_sites) & !is.null(min_sites_pos)){
					inc_years = yr_summary$YEAR[which(yr_summary$NSITES >= min_sites & yr_summary$NSITES_POS >= min_sites_pos)]
				} else if(!is.null(min_sites) & is.null(min_sites_pos)){
					inc_years = yr_summary$YEAR[which(yr_summary$NSITES >= min_sites)]
				} else {
					inc_years = yr_summary$YEAR[which(yr_summary$NSITES_POS >= min_sites_pos)]
				}
				# Print statement showing that years have been removed
					# Detemrine which years will be excluded
						exc_years = yr_summary$YEAR[which(!yr_summary$YEAR %in% inc_years)]
					# Print output
					if(length(exc_years) > 0){
						cat("Warning: the following years have been excluded due to being below the minimum site or minimum positive site thresholds:\n  ",paste(exc_years, collapse=", "), "\n", sep="")
					}
			} else {
				inc_years = yr_summary$YEAR
			}
			
			# Determine which rows need to be excluded
				exc_inds = which(!site_yr_sum$YEAR %in% inc_years | is.na(site_yr_sum$N_EST))
			# Fit Poisson GLM model to N_EST values in site_yr_sum
				if(length(exc_inds) > 0){
					glm_obj = try(speedglm(round(N_EST) ~ factor(YEAR) + factor(SITE) - 1, data = site_yr_sum[-exc_inds,], family = poisson(), weights = site_yr_sum$FP_SUM[-exc_inds]), silent = TRUE)
				} else {
					glm_obj = try(speedglm(round(N_EST) ~ factor(YEAR) + factor(SITE) - 1, data = site_yr_sum, family = poisson(), weights = site_yr_sum$FP_SUM), silent = TRUE)
				}
				if("try-error" %in% class(glm_obj)){
					# Show progress
						cat("WARNING: fitting GAI to N_EST using speedglm failed trying standard glm(see error below)\n",geterrmessage(),"\n")
					# Try using standard glm
					if(length(exc_inds) > 0){
						glm_obj = try(glm(round(N_EST) ~ factor(YEAR) + factor(SITE) - 1, data = site_yr_sum[-exc_inds,], family = poisson(), weights = site_yr_sum$FP_SUM[-exc_inds]), silent = TRUE)
					} else {
						glm_obj = try(glm(round(N_EST) ~ factor(YEAR) + factor(SITE) - 1, data = site_yr_sum, family = poisson(), weights = site_yr_sum$FP_SUM), silent = TRUE)
					}
				}
		} else {
			# GLM_T method
			#-------------
			# Show progress
					cat("Fitting GAI - GLM_T method\n")
			# Build expanded dataset (so every site/year comb has every tp value), imput missing counts and then calculate site indices
				sindex = gai_si(sp_data, nm_data, site_yr_sum, tp_col = tp_col)
			# Show progress
				cat("Estimate site index values from count data (via imputting)\n")
			# Add in site index only data
			if(ifelse(!is.null(si_data),nrow(si_data),0) > 0){
				if("ORG_SITE" %in% names(si_data)){
					sindex = rbind(sindex,si_data[,c("YEAR","SITE","ORG_SITE","SINDEX")])
					site_yr_sum = rbind(site_yr_sum, data.frame(SITE = si_data$SITE, YEAR = si_data$YEAR, ORG_SITE = si_data$ORG_SITE, TOTAL = si_data$SINDEX, FP_SUM = 1, N_EST = si_data$SINDEX, N_POS = NA, CNT_TRAP = NA, INDEX_ONLY = 1))
				} else {
					sindex = rbind(sindex,si_data[,c("YEAR","SITE","SINDEX")])
					site_yr_sum = rbind(site_yr_sum, data.frame(SITE = si_data$SITE, YEAR = si_data$YEAR, TOTAL = si_data$SINDEX, FP_SUM = 1, N_EST = si_data$SINDEX, N_POS = NA, CNT_TRAP = NA, INDEX_ONLY = 1))
				}
				# Show progress
				cat("Added Site Index only data to sindex values estimated from counts\n")
			}
			# Merge sindex values from glm into site_yr_sum
				site_yr_sum = merge(site_yr_sum, sindex, by=c("SITE","YEAR"), all = TRUE)
			# Calculate year summary values
				yr_summary = calc_yr_summ(site_yr_sum)
			# Determine if min_sites and/or min_sites_pos specified and if so then limit data to only years that pass thresholds
			if(!is.null(min_sites) | !is.null(min_sites_pos)){
				if(!is.null(min_sites) & !is.null(min_sites_pos)){
					inc_years = yr_summary$YEAR[which(yr_summary$NSITES >= min_sites & yr_summary$NSITES_POS >= min_sites_pos)]
				} else if(!is.null(min_sites) & is.null(min_sites_pos)){
					inc_years = yr_summary$YEAR[which(yr_summary$NSITES >= min_sites)]
				} else {
					inc_years = yr_summary$YEAR[which(yr_summary$NSITES_POS >= min_sites_pos)]
				}
				# Print statement showing that years have been removed
					# Detemrine which years will be excluded
						exc_years = yr_summary$YEAR[which(!yr_summary$YEAR %in% inc_years)]
					# Print output
					if(length(exc_years) > 0){
						cat("Warning: the following years have been excluded due to being below the minimum site or minimum positive site thresholds:\n  ",paste(exc_years, collapse=", "), "\n", sep="")
					}
			} else {
				inc_years = yr_summary$YEAR
			}
			# Fit glm with sindex as response
			# Determine which rows need to be excluded (if any)
				exc_inds = which(!site_yr_sum$YEAR %in% inc_years | is.na(site_yr_sum$N_EST))
			if(length(exc_inds) > 0){
				glm_obj = try(speedglm(round(SINDEX) ~ factor(YEAR) + factor(SITE) - 1, data = site_yr_sum[-exc_inds,], weights = site_yr_sum$FP_SUM[-exc_inds], family = poisson()), silent = TRUE)
			} else {
				glm_obj = try(speedglm(round(SINDEX) ~ factor(YEAR) + factor(SITE) - 1, data = site_yr_sum,  weights = site_yr_sum$FP_SUM, family = poisson()), silent = TRUE)
			}
			if("try-error" %in% class(glm_obj)){
				# Show progress
						cat("WARNING: fitting GAI to SINDEX using speedglm failed trying standard glm (see error below)\n",geterrmessage(),"\n")
					# Try using standard glm
				if(length(exc_inds) > 0){
					glm_obj = try(glm(round(SINDEX) ~ factor(YEAR) + factor(SITE) - 1, data = site_yr_sum[-exc_inds,], weights = site_yr_sum$FP_SUM, family = poisson()), silent = TRUE)
				} else {
					glm_obj = try(glm(round(SINDEX) ~ factor(YEAR) + factor(SITE) - 1, data = site_yr_sum, weights = site_yr_sum$FP_SUM, family = poisson()), silent = TRUE)
				}
			}
			
			# Cleanup
			rm(sindex)
		}
			# For both GLM_N and GLM_T methods once glm fitted then process of getting collated index values is the same
			if(!"try-error" %in% class(glm_obj)){
				# Predict site indices back from model (allowing indices to be given for years/site combinations which did not have data originally) - this is only done if the sindex values are being saved (e.g. a filename/path provides) otherwise there is no point!
					# Add SINDEX column to site_yr_sum data
						site_yr_sum[,"SINDEX"] = round(site_yr_sum$N_EST)
					# Predict site indices
						site_yr_sum = pred_SI_glm(glm_obj, site_yr_sum)
					# Calculate yearly totals from predicted site indices
						pred_yr_tots = tapply(site_yr_sum$PRED_N_EST, site_yr_sum$YEAR, sum, na.rm = TRUE)
				if(!is.null(sindex_fname)){
						# Progress statement
						cat("N_EST values predicted from fitted GLM\n")
					# Save to file
					if(grepl("[.]csv[.]gz$",sindex_fname)){
						# When passed a file name that ends .csv.gz then output to be saved to gz csv file to reduce file size but still allow the file to be appended to
						# First need to determine whether data is being saved to a new file or appended to existing file
						if(!file.exists(sindex_fname)){
							f_append = FALSE
						} else {
							f_append = TRUE
						}
						# Now need to build temporary dataset to hold data to be saved
						if(!"ORG_SITE" %in% names(site_yr_sum)){
							temp_sindex = data.frame(SPECIES = sp_data$SPECIES[1], site_yr_sum[,c("SITE","YEAR")], ORG_SITE = NA, site_yr_sum[,which(!names(site_yr_sum) %in% c("SITE","YEAR"))], iboot = ifelse(is.null(iboot),NA, iboot))
						} else {
							temp_sindex = data.frame(SPECIES = sp_data$SPECIES[1], site_yr_sum, iboot = ifelse(is.null(iboot),NA, iboot))
						}
						# Open connection to gz file in write or append mode (depending on value of f_append
							f_con = gzfile(sindex_fname, open = ifelse(TRUE,"a","w"))
						# Write data to connection (turning on/off column names based on !f_append)
							write.table(temp_sindex, file = f_con, col.names = !f_append, row.names = FALSE, quote = TRUE, sep=",", na = "")
						# Close connection
							close(f_con)	
						# Remove temp file
							rm(temp_sindex)
					} else {
						if(is.null(iboot)){
							saveRDS(data.frame(SPECIES = sp_data$SPECIES[1], site_yr_sum), file=sindex_fname)
						} else {
							saveRDS(data.frame(SPECIES = sp_data$SPECIES[1], site_yr_sum, iboot = iboot), file=sindex_fname)
						}
					}
					
					# Show progress
						# first need to determine base directory of site index file path is given
						base_dir = dirname(sindex_fname)
						cat("Site index data saved to file '",ifelse(base_dir == ".",sindex_fname, basename(sindex_fname)),"' in directory:\n   ", ifelse(base_dir == ".",getwd(), base_dir),"\n", sep="")
				}
				# Extract year coefs
					# Get all coefficients and standard erroors
						all_coefs = summary(glm_obj)$coef
					# If all_coefs is a matrix rather than data.frame (when using glm rather than speedglm) then convert to data.frame
					if(class(all_coefs) == "matrix"){
						all_coefs = as.data.frame(all_coefs)
					}
					# Extract year coefs from all coefs
						yr_coefs = all_coefs[grepl("YEAR",row.names(all_coefs)),c("Estimate","Std. Error")]
						# Get that estimate is numeric and not factor
						if(class(yr_coefs$Estimate) == "factor"){
							yr_coefs$Estimate = as.numeric(as.character(yr_coefs$Estimate))
							# if yr_coef$Estimate a factor then so is all_coefs and also need to covert this so that they mean site coef can be calculated
							all_coefs$Estimate = as.numeric(as.character(all_coefs$Estimate))
						}
						# Set row names to be the year extracted from the coef name
						row.names(yr_coefs) = gsub("^(factor\\())?YEAR(\\))?([[:digit:]]{,4})$","\\3",row.names(yr_coefs))
					# Calculate a version of the coefs that is on the orginal count scale
						# Determine mean of  all site coefs
							if(class(all_coefs$Estimate) == "factor"){
								mean_site = mean(c(0,as.numeric(as.character(all_coefs$Estimate[grepl("SITE",row.names(all_coefs))]))), na.rm = TRUE)
							} else {
								mean_site = mean(c(0,all_coefs$Estimate[grepl("SITE",row.names(all_coefs))]), na.rm = TRUE)
							}
							n_sites = length(c(0,all_coefs$Estimate[grepl("SITE",row.names(all_coefs))]))
						# Add the mean site coefficient to each year coef (and convert from log scale to normal scale)
							adj_yr_coefs = yr_coefs
							adj_yr_coefs$Estimate = adj_yr_coefs$Estimate + mean_site
							adj_yr_coefs = exp(adj_yr_coefs)
							adj_yr_coefs[,"total_cnt"] = adj_yr_coefs$Estimate * n_sites
					# rm all _coefs
						rm(all_coefs)
				# Convert year coefs and std errors from being on natural log scale to log10 scale by dividing by log(10)
					yr_coefs = yr_coefs/log(10)
				# Show progress
					cat("Year coefficients extracted from fitted GLM (and converted to log10)\n")
				# Create data.frame to hold results
					all_yrs = inc_years
					ci_data = data.frame(SPECIES = sp_data$SPECIES[1], yr_summary[which(yr_summary$YEAR %in% inc_years),], TR0OBS = NA, TRMOBS = NA, TR0SE = NA, YR_COEF_CNT = NA, TRCOUNT = NA)
					# Clean up
						rm(all_yrs)
					# Set row.names to year
						row.names(ci_data) = ci_data$YEAR
				# Save original coefs
					ci_data[row.names(yr_coefs),"YR_COEFS"] = yr_coefs$Estimate
				# Rescale so that first year = 2 
					ci_data[row.names(yr_coefs),"TR0OBS"] = 2 + (yr_coefs$Estimate-yr_coefs$Estimate[1])
				# Convert index from natural log to log10 (by divinding by log(10)) and then rescale so that mean of series = 2
					ci_data[row.names(yr_coefs),"TRMOBS"] = 2 + (yr_coefs$Estimate-mean(yr_coefs$Estimate))
				# Add std error to table
					ci_data[row.names(yr_coefs),"TR0SE"] = yr_coefs[,"Std. Error"]
				# Add adjusted version of year coefs on original count scale for a typical site
					ci_data[row.names(adj_yr_coefs),c("YR_COEF_CNT","TRCOUNT")] = adj_yr_coefs[,c("Estimate","total_cnt")]
				# Add yearly total predicted N_EST values to output object
					ci_data[names(pred_yr_tots),"TR_PRED_N_EST"] = pred_yr_tots
				# Show progress
					cat("Year coefficients rescaled to provided collated indices\n")
			} else {
				ci_data = NULL
				# Show progress and warning
				cat("WARNING: GLM fitting using glm function failed returning NULL (see error below)\n",geterrmessage(),"\n")
			}
	} else if(methodCI == "SIMPLE_N"){
		# SIMPLE_N method
		#-------------
		# Show progress
			cat("Fitting GAI - SIMPLE_N method\n")
		# Add in site index only data (if supplied)
		if(ifelse(!is.null(si_data),nrow(si_data),0) > 0){
			if("ORG_SITE" %in% names(si_data)){
				site_yr_sum = rbind(site_yr_sum, data.frame(SITE = si_data$SITE, YEAR = si_data$YEAR, ORG_SITE = si_data$ORG_SITE, TOTAL = si_data$SINDEX, FP_SUM = 1, N_EST = si_data$SINDEX, N_POS = NA, CNT_TRAP = NA, INDEX_ONLY = 1))
			} else {
				site_yr_sum = rbind(site_yr_sum, data.frame(SITE = si_data$SITE, YEAR = si_data$YEAR, TOTAL = si_data$SINDEX, FP_SUM = 1, N_EST = si_data$SINDEX, N_POS = NA, CNT_TRAP = NA, INDEX_ONLY = 1))
			}
			# Show progress
				cat("Added Site Index only data to N_EST values estimated from counts\n")
		}
		# Save sindex N_EST & site index only data as sindex values (if filename/path supplied as argument)
		if(!is.null(sindex_fname)){
			saveRDS(data.frame(SPECIES = sp_data$SPECIES[1], site_yr_sum, SINDEX = round(site_yr_sum$N_EST), INDEX_TYPE = "GAI - SIMPLE N"), file=sindex_fname)
			# Show progress
				cat("N_EST values saved to file as sindex values\n")
		}
		# Determine if min_sites and/or min_sites_pos specified and if so then limit data to only years that pass thresholds
			if(!is.null(min_sites) | !is.null(min_sites_pos)){
				yr_summary = calc_yr_summ(site_yr_sum)
				if(!is.null(min_sites) & !is.null(min_sites_pos)){
					inc_years = yr_summary$YEAR[which(yr_summary$NSITES >= min_sites & yr_summary$NSITES_POS >= min_sites_pos)]
				} else if(!is.null(min_sites) & is.null(min_sites_pos)){
					inc_years = yr_summary$YEAR[which(yr_summary$NSITES >= min_sites)]
				} else {
					inc_years = yr_summary$YEAR[which(yr_summary$NSITES_POS >= min_sites_pos)]
				}
				# Print statement showing that years have been removed
					# Detemrine which years will be excluded
						exc_years = yr_summary$YEAR[which(!yr_summary$YEAR %in% inc_years)]
					# Print output
					cat("The following years have been excluded based on min_sites and/or min_sites_pos:\n  ",paste(exc_years, collapse=", "))
			} else {
				inc_years = yr_summary$YEAR
			}
		# Calculate yearly mean N_EST and its standard error along with 
			inc_inds = which(site_yr_sum$YEAR %in% inc_years)
			mn_nest = tapply(site_yr_sum$N_EST[inc_inds], site_yr_sum$YEAR[inc_inds], mean, na.rm = TRUE)
			n_sites = sapply(tapply(site_yr_sum$SITE, site_yr_sum$YEAR, unique),length)
			se_nest = tapply(site_yr_sum$N_EST[inc_inds], site_yr_sum$YEAR[inc_inds], sd, na.rm = TRUE) / sqrt(n_sites)
		# Show progress
			cat("Yearly mean N_EST calculated\n")
		# Setup data.frame to hold data and fill with scaled data
			if(!is.null(si_data)){
				all_yrs = sort(unique(c(sp_data$YEAR, si_data$YEAR)))
			} else {
				all_yrs = sort(unique(sp_data$YEAR))
			}
			ci_data = data.frame(SPECIES = sp_data$SPECIES[1], YEAR = all_yrs,NSITES = NA, NSITES_POS = NA, TR0OBS = NA, TRMOBS = NA, TR0SE = NA)
				# Clean up
					rm(all_yrs)
				# Set row.names to year
					row.names(ci_data) = ci_data$YEAR
			# Add n_sites
				ci_data[names(n_sites),"NSITES"] = n_sites
			# Add number of positive sites
				pos_inds = which(site_yr_sum$TOTAL > 0)
				n_sites = sapply(tapply(site_yr_sum$SITE[pos_inds], site_yr_sum$YEAR[pos_inds], unique),length)
				# Add to ci results
				ci_data[names(n_sites),"NSITES_POS"] = n_sites
				# cleanup 
				rm(n_sites, pos_inds)
			# Calculate index by scaling means so that they are centred on mean of series and scale by divinding by standard deviation
				temp = scale(mn_nest)
			# Now rescale to index where the first year = 2
				ci_data[row.names(temp),"TR0OBS"] = 2 + (temp-temp[1]) / log(10)
			# Now rescale to index where mean of series = 2 (temp should already be centered on zero so can just add 2)
				ci_data[row.names(temp),"TRMOBS"] = 2 + temp / log(10)
			# Add standard error centring so that the mean is zero and then rescale by the same scale as the mean was rescaled
				ci_data[names(se_nest),"TR0SE"] = scale(se_nest,center = FALSE, scale = attr(temp,"scaled:scale")) / log(10)
			# Show progress
				cat("Yearly mean N_EST VALUES rescaled to provided collated indices\n")
	}
	
	# Remove row.nams from ci_data
		row.names(ci_data) = NULL
	# Return ci_data
		return(ci_data)
}

gai_coll_ind_wrpr = function(iboot, sp_data, nm_data, other_args, si_data = NULL, boot_samp = NULL){
	# Create bootstrap version of data if required
	if(iboot != 0){
		# If boot_samp is null then can't create bootstrap sample of data so skip and return NULL for nm_data
		if(!is.null(boot_samp)){
			# Get current bootstrap site codes from boot_samp
				cur_boot = setNames(boot_samp[,c("site",paste("iboot",iboot,sep="_"))],c("SITE","ORG_SITE"))
			# bootstrap count data
				sp_data = boot_sample(sp_data, boot_samp = cur_boot)
			# 
			# bootstrap site index data if not null
			if(!is.null(si_data)){
				si_data = boot_sample(si_data, boot_samp = cur_boot)
			}
			# keep only nm_data for current bootstrap
				nm_data = nm_data[nm_data$iboot == iboot,]
		} else {
			nm_data = NULL
		}
	} else {
		# Check whether nm_data contains a column iboot
		if("iboot" %in% names(nm_data)){
			# If so then check it only contains data for iboot for current value of iboot
			nm_data = nm_data[nm_data$iboot == iboot,]
		}	
	}
	# Check that nm_data isn't a zero row data.frame if it is then set to NULL so in next block it won't attempt to estimate collated index for this bootstrap
	if(!is.null(nm_data)){
		if(nrow(nm_data) == 0){
			nm_data = NULL
		}
	}
	# If nm_data is not null then estimate collated index otherwise return NULL from function
	if(!is.null(nm_data)){
		# Estimate collated index
			cur_ci = f_do_call("gai_coll_ind",c(list(sp_data = sp_data, si_data = si_data, nm_data = nm_data, iboot = iboot), other_args))
		# add iboot to cur_ci so that when collated results from individual bootstraps are disguishable
		if(!is.null(cur_ci)){
			cur_ci[,"iboot"] = iboot
		}
	} else {
		cur_ci = NULL
	}
	# Now return results
		return(cur_ci)
}

read_boot_nm = function(sp, yrs, iboot_range, dpath){
	# Get list of files in dpath that are for current species
		f_list = list.files(dpath, pattern = paste0("Sp ",sp,"_.*[.]rds"))
		
	if(length(f_list) > 0){
		# Determine min/max bootstraps in each file
			# Extract bootstrap numbers from file name
				f_boots = gsub(paste0("Sp ",sp,"_iboot ([[:digit:]]{1,}) to ([[:digit:]]{1,})_nm[.]rds"),"\\1 \\2", f_list)
			# Seperate out min and max and then bind into a matrix
				f_boots = do.call("rbind",strsplit(f_boots,split = " "))
				# convert matrix from text to numeric using mode function
					mode(f_boots) = "numeric"
		# Determine which files are  need to cover the current bootstrap range specified by iboot_range
			inds = which( ( iboot_range[1] >= f_boots[,1] & iboot_range[1] <= f_boots[,2] ) | ( iboot_range[2] >= f_boots[,1] & iboot_range[2] <= f_boots[,2] ) )
		# Check that at least 1 files found that matched boot range specified
		if(length(inds) > 0){
			# Now load in files containing the nm values required for the current iboot_range
				boot_nm = do.call("rbind",lapply(file.path(dpath, f_list[inds]), readRDS))
			# Check which bootstrap were in the files loaded as they may contain unneccesary bootstraps if bootstrap batch sizes were/are different between 1st and 2nd stage
				boot_lims = range(boot_nm$iboot)
				# Only keep nm values for required iboot_range (if iboot_range is not same as boot_nm)
				if(any(boot_lims !=  iboot_range)){
					boot_nm = boot_nm[which(boot_nm$iboot >= iboot_range[1] & boot_nm$iboot <= iboot_range[2]),]
				}
		} else {
			boot_nm = NULL
		}
	} else {
		boot_nm = NULL
	}
	# Return the required bootstrap	flight periods
	return(boot_nm)
}

gai_coll_ind_clust = function(i){
	# Store start_time
		start_time = Sys.time()
	# Load in parameter values file and then assign values for current run as object in environment 
		# Cluster parameters/combinations
			# Look in current working directory so need to ensure working directory is set prior to this or clust_params exported to each node
			# Load all clust_params
				clust_params = readRDS("gai_clust_params.rds")
			# Extract elements for ith run from object and assign as objects in current environment
				list2env(clust_params[i,], environment())
			# remove clust params
				rm(clust_params)
		# Fixed/constant parameters
			# Look in current working directory so need to ensure working directory is set prior to this or fixed_params exported to each node
				fixed_params = readRDS("gai_fixed_params.rds")
			# Assign all the objects in the fixed_params file as objects in the current environment
				list2env(fixed_params, environment())
	# Determine if data_filter was defined in fixed_params.rds if not create and set to null
	if(!exists("data_filter")){
			data_filter = NULL
	}

	
	# Get the count data
		sp_data = read_data_files(sp = species, yrs = yrs, dpath = paths["data"], data_type="count", data_filter = data_filter)
		# Show progress
		#if(verbose)
			#cat("  Species count data loaded (",nrow(sp_data), " rows)\n",sep="")
	# Get the sindex only data
		si_data = read_data_files(sp = species, yrs = yrs, dpath = paths["data"], data_type="sindex", data_filter = data_filter)
		# Show progress
		if(ifelse(is.null(si_data),0, nrow(si_data)) > 0){
			#if(verbose)
			#cat("  Site Index only data loaded (",nrow(si_data), " rows)\n",sep="")
		}
		
	# Determine which bootstraps are going to be run in the current batch
		# If batch_nboot not found then create it and set it to NULL (in which case no bootstraps will be performed and only the flight periods for the non-bootstrapped data will be analysed
		if(!exists("batch_nboot")){
				batch_nboot = NULL
		} else {
			# If batch_nboot has been specified as 0 then set it to NULL (these are equivalent and both only analyses the real non-bootstrapped data)
			if(batch_nboot == 0){
				batch_nboot = NULL
			}
		}
		
	# Get the normalised flight period data (if bootstrap then get this from the bootstrap NM directory
	if(is.null(batch_nboot)){
			# If batch_nboot is null then only need to run real data (iboot = 0) for each year
				node_boots = 0
			# if no bootstrapping then set boot_samp to null
				boot_samp = NULL
		} else {
			# Determine which bootstraps are going to be run on the current node
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
			
		# Check whether species has a boot sample file (may not if species didn't have > 1 site in the data for the time period in question)
		if(file.exists(file.path(paths["boot_samp"],paste0("Sp ",species,"_bootsamps.rds")))){
			# Read in boot_sample file containing all site codes for each boot_sample
				boot_samp = readRDS(file.path(paths["boot_samp"],paste0("Sp ",species,"_bootsamps.rds")))
				# Only keep those for current bootstraps
					boot_samp = boot_samp[,c("site",paste("iboot",node_boots[which(node_boots > 0)],sep="_"))]
		} else {
			boot_samp = NULL
		}
	}
		
	# Determine if sindex values to be saved and if construct the file name
	if(ifelse(!exists("save_sindex"),FALSE, save_sindex)){
		if(is.null(batch_nboot) | nboot == 0){
			sindex_fname = file.path(paths["sindex"],paste0("Sp ",species,"_sindex.rds"))
		} else {
			sindex_fname = file.path(paths["boot_sindex"],paste0("Sp ",species,"_iboot ",paste(range(node_boots),collapse=" to "),"_sindex.csv.gz"))
		}
		# If file already exists then remove it
		if(file.exists(sindex_fname)){
			invisible(file.remove(sindex_fname))
		}
	}
	
	# Look in current workspace and build list of arguments that match arguments in est_fp_spline function
		# get list of gai_coll_ind arguments
			fun_args = names(formals(gai_coll_ind))
			# Exclude sp_data, si_data and nm_data from this list (as this will be directly added)
				fun_args = fun_args[!fun_args %in% c("sp_data","si_data","nm_data")]
		# Now look for these arguments in the workspace and collect into a list
			other_args = mget(ls()[which(ls() %in% fun_args)])
	

			
	# Read in flight periods for these bootsrap samples (from files created in stage 1)
	if(is.null(batch_nboot)){
		# If no bootrrapping then simply load real flight period estimates
		nm_data = read_data_files(sp = species, yrs = yrs, dpath = paths["nm"], data_type="nm")
	} else {
		if(iboot_start == 0){
			# If nodes is to run first set of bootsrapping then it also need the real flight periods for the first run
			# Read in real flight period estimates
				sp_nm = read_data_files(sp = species, yrs = yrs, dpath = paths["nm"], data_type="nm")
				# Add iboot column to allow real data to be stacked with bootsrap results
				if(!is.null(sp_nm)){
					sp_nm[,"iboot"] = 0
				}
			# Now add sp_nm flight periods to those for bootstraps
				nm_data = rbind(sp_nm, read_boot_nm(sp = species, yrs = yrs, iboot_range =  range(node_boots[which(node_boots > 0)]),dpath = paths["boot_nm"]))
			# now remove sp_nm
				rm(sp_nm)
		} else {
			# If node is running later bootstraps then it can ignore real flight periods and only load relevant ones for current bootstraps
				nm_data = read_boot_nm(sp = species, yrs = yrs, iboot_range =  range(node_boots),dpath = paths["boot_nm"])
		}
	}
	# Continue if any NM data returned from previous block, if not then stop and return NULL from function
	if(!is.null(nm_data)){
		# Check NM data (if nm_data found!)
			# Check if tp_col and/or season_lim exist in workspace otherwise use default from gai_coll_ind funciton
				temp_tp_col = if(exists("tp_col")){ tp_col } else { eval(formals(gai_coll_ind)$tp_col) }
				temp_season_lim = if(exists("season_lim")){ season_lim } else { eval(formals(gai_coll_ind)$season_lim) }
			# Now check nm and remove NM which fail criteria (currently set as defaults in check_nm function should be opened up later to be passable as arguments from higher level functions!!)
				nm_data = check_nm(sp_data = sp_data, sp_nm = nm_data, tp_col = temp_tp_col, season_lim = temp_season_lim )
			# If the nm check ends up removing all of the rows of data from NM data then exit the function without calculating the collated index
			if(nrow(nm_data) == 0){
				return(NULL)
			}
	
		# If save_sindex exists/was supplied then use this value otherwise set to FALSE
		if(!exists("save_sindex")){
			save_sindex = FALSE
		}
		# Estimate flight period using est_fp_spline_wrpr wrapper function which call est_fp_spline with correct arguments after subsetting required datasets and/or recreating bootstrap samples of data 
			sp_ci = do.call("rbind",lapply(node_boots,FUN = gai_coll_ind_wrpr, sp_data = sp_data, nm_data = nm_data, si_data = si_data, boot_samp = boot_samp, other_args = other_args))
		# Save results
		if(is.null(batch_nboot) | nboot == 0){
			# Save real results (exclude iboot column)
				saveRDS(sp_ci[,!grepl("iboot",names(sp_ci))],file=file.path(paths["coll_ind"],paste0("Sp ",species,"_coll_ind.rds")))
		} else {
			if(iboot_start == 0){
				# Save real results (exclude iboot column)
					real_inds = which(sp_ci$iboot == 0)
				if(length(real_inds) > 0){
					saveRDS(sp_ci[real_inds,!grepl("iboot",names(sp_ci))],file=file.path(paths["coll_ind"],paste0("Sp ",species,"_coll_ind.rds")))
				}
				# Construct file name for bootstrap file
					boot_fname = paste0("Sp ",species,"_iboot ", paste(range(node_boots[which(node_boots > 0)]), collapse=" to "),"_coll_ind.rds")
				# Save bootstrap results (excluding iboot == 0 run)
					saveRDS(sp_ci[which(sp_ci$iboot > 0),],file=file.path(paths["boot_coll_ind"],boot_fname))
			} else {
				# Construct file name for bootstrap file
					boot_fname = paste0("Sp ",species,"_iboot ", paste(range(node_boots), collapse=" to "),"_coll_ind.rds")
				# Save bootstrap results (save all as should be no iboot == 0)
					saveRDS(sp_ci,file=file.path(paths["boot_coll_ind"],boot_fname))
			}
		}
	} else {
		sp_ci = NULL
	}
	return(sp_ci)	
}

gai_coll_ind_clust_old = function(i){
	# Store start_time
		start_time = Sys.time()
	# Load in parameter values file and then assign values for current run as object in environment 
		# Deterime whether clust_params object already exists in workspace
		if(!exists("clust_params")){
			# Look in current working directory so need to ensure working directory is set prior to this or clust_params exported to each node
			# Load all clust_params
				clust_params = readRDS("gai_clust_params.rds")
			# Extract elements for ith run from object and assign as objects in current environment
				list2env(clust_params[i,], environment())
			# remove clust params
				rm(clust_params)
		}
		# Look for and load in the fixed/constant parameters used for runs
		if(!exists("fixed_params")){
			# Look in current working directory so need to ensure working directory is set prior to this or fixed_params exported to each node
				fixed_params = readRDS("gai_fixed_params.rds")
			# Assign all the objects in the fixed_params file as objects in the current environment
				list2env(fixed_params, environment())
		}
			
	# Setup log file and redirect normal output (i.e. cat/print) and also other for warnings/errors messages) for the species using sink to direct the output there
		# Block sending output to log file if the function is flagged for debugging
		if(!(isdebugged(gai_coll_ind) | isdebugged(gai_coll_ind_clust))){
			# Set options so warnings appear at time they are generated
				options(warn = 1)
			# Setup a file connection in write mode
				fcon_log = file(file.path(paths["log"],paste0("Sp ",species,ifelse(exists("iboot"), paste0("_iboot ",iboot),""),"_coll_ind_log.txt")), open = "w")
			# Use sink to redirect standard output and error messages to this file
				sink(file = fcon_log, type = "output", split = TRUE)
				sink(file = fcon_log, type = "message", append = TRUE)
			# On exist close the redirects via sink and then close the file connection
			on.exit({
				sink(NULL, type="output")
				sink(NULL, type="message")
				close(fcon_log)
			})
		}
		
	# Print header to log file and/or screen
		cat(rep("-",80),"\nSTARTING SECOND STAGE OF GAI ANALYSIS\nCOLLATED INDEX ESTIMATION FOR SPECIES = ",species,ifelse(exists("iboot"), paste0(" BOOTSTRAP = ",iboot),""), "\nANALYSIS STARTED: ",format(start_time,"%d/%m/%Y %T"),"\n",rep("-",80),"\n\n",sep="")

	# Show progress
		cat("Loading data\n")
	# Check if there is a data filter (if not add object and set to null		
		if(!exists("data_filter")){
			data_filter = NULL
		}
	# Get the count data
		sp_data = read_data_files(sp = species, yrs = yrs, dpath = paths["data"], data_type="count", data_filter = data_filter)
		# Show progress
		cat("  Species count data loaded (",nrow(sp_data), " rows)\n",sep="")
	# Get the sindex only data
		si_data = read_data_files(sp = species, yrs = yrs, dpath = paths["data"], data_type="sindex", data_filter = data_filter)
		# Show progress
		if(ifelse(is.null(si_data),0, nrow(si_data)) > 0){
			cat("  Site Index only data loaded (",nrow(si_data), " rows)\n",sep="")
		}
	# Get the normalised flight period data (if bootstrap then get this from the bootstrap NM directory
		if(ifelse(exists("iboot"),iboot,0) == 0){
			nm_data = read_data_files(sp = species, yrs = yrs, dpath = paths["nm"], data_type="nm")
		} else {
			nm_data = read_data_files(sp = species, yrs = yrs, dpath = paths["boot_nm"], iboot = iboot, data_type="nm")
		}
		# Show progress
		cat("  Flight Period data loaded (",nrow(nm_data), " rows)\n")
		
	# Setup variable to hold path under which results will be saved
		out_path = file.path(paths["coll_ind"],paste0("Sp ",species,"_coll_ind.rds"))

	# Determine if the site index values are to be stored and if so setup the filename/path of the file to which the data is to be stored
	if(ifelse(exists("save_sindex"),save_sindex,FALSE)){
		if(ifelse(exists("iboot"),iboot,0) == 0){
			# current run is not for a bootstrap sample
			sindex_fname = file.path(paths["sindex"],paste0("Sp ",species,"_sindex.rds"))
		} else {
			# Current run is for a bootstrap sample
			sindex_fname = file.path(paths["boot_sindex"],paste0("Sp ",species,"_iboot ",iboot,"_sindex.rds"))
		}
	}
	# Determine if this is to be a bootstrap run or not
	if(ifelse(exists("iboot"),iboot,0) != 0){
		# Show progress
			cat("  Bootstrap sample sites info loaded\n")
		# Read in boot_sample file containing all site codes for each boot_sample
			boot_samp = readRDS(file.path(paths["boot_samp"],paste0("Sp ",species,"_bootsamps.rds")))
		# Only keep those for current bootstrap
			boot_samp = boot_samp[,c("SITE",paste("iboot",iboot,sep="_"))]
		# Rename iboot site codes to ORG_SITE
			names(boot_samp)[2] = "ORG_SITE"
		# Apply the bootstrap sample codes to the sp data to create a bootstrapped version of the year data
			sp_data = boot_sample(sp_data, boot_samp)
		# Apply the bootstrap sample codes to the site index only data to create a bootstrapped version of that data
			si_data = boot_sample(si_data, boot_samp)
		# Modify the output file path so bootstrap results kept seperate
			out_path = file.path(paths["boot_coll_ind"],paste0("Sp ",species, "_iboot ",iboot,"_coll_ind.rds"))
		# Show progress
			cat("Bootstrapped versions of data created\n")
	}
	if(!is.null(sp_data) & !is.null(nm_data)){
		if(nrow(sp_data) > 0 & nrow(nm_data) > 0){
			# Get list of formals for fp fitting function
				fun_args = names(formals("gai_coll_ind"))
				# Reduce to only those that exist in workspace
					fun_args = fun_args[fun_args %in% ls()]
			# Attempt to fit the GAI using mget to create a list of the arguments that are in thw workspace
				coll_ind = f_do_call("gai_coll_ind", mget(fun_args))
			# If sp_nm fitting return values and iboot exists then add iboot value to sp_nm data.frame
			if(!is.null(coll_ind)){
				# Add iboot column to the data.frame is iboot exists
				if(exists("iboot")){
					coll_ind[,"iboot"] = iboot
				}
				# Save the data for the species
					saveRDS(coll_ind, file = out_path)
				# Show progress
					cat("Collated Index values saved to file\n")
			}
		} else {
			if(nrow(sp_data) == 0)
				cat("ERROR: Species data file contains no rows (aborting GAI fitting)\n")
			if(nrow(nm_data) == 0)
				cat("ERROR: Normalised flight period data file contains no rows (aborting GAI fitting)\n")
			coll_ind = NULL
		}
	} else {
		if(is.null(sp_data))
			cat("ERROR: No species data found (aborting GAI fitting)\n")
		if(is.null(nm_data))
			cat("ERROR: No normalised flight period data found (aborting GAI fitting)\n")
		coll_ind = NULL
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
		
	# Return the collated indices
		return(coll_ind)
}

fit_gai = function(dpath_data, dpath_fp_res = NULL, spp = NULL, yrs = NULL, parallel = TRUE, parallel_pkg = "parallel", n_cpus = NULL, nboot = NULL, batch_nboot = 10, iboot_initial = 0, ...){
	# Check that dpath_data exits
	if(!file.exists(dpath_data)){
		stop("Directory supplied to dpath_data does not exist")
	} else {
		# Normalise dpath_data to ensure file path is absolute not relative
			dpath_data = normalizePath(dpath_data, winslash = "/")
	}
	# Determine which species and years have data files in the dpath_data folder (assumes they follow the template
	if(length(list.files(dpath_data)) == 0)
		stop("No data files were found in directory supplied to dpath_data")
	# If dpath_fp_res is null then assume the working directory is where the NM values to be used are store and also the directory where the results from the GAI will be stored (in sub-directories)
	if(is.null(dpath_fp_res)){
		dpath_fp_res = getwd()
	} else {
		# If a path was supplied then check it exists and if so normalise the path to ensure that it is absolute and not relative
		if(!file.exists(dpath_fp_res)){
			stop("Directory supplied to dpath_fp_res does not exists")
		} else {
			dpath_fp_res = normalizePath(dpath_fp_res, winslash = "/")
		}
	}
	# Check that dpath_fp_res includes a folder title NM that contains the normalised flight period values
	if(!file.exists(file.path(dpath_fp_res,"NM")) | !length(list.files(file.path(dpath_fp_res,"NM"))) > 0)
		stop("Could not find normalised flight period results in the directory supplied to dpath_fp_res (should be path to directory containing output from flight period estimation including nm values in a sub-directory title NM)")
	# Setup a variable to hold all the different path information
		paths = c(data = dpath_data, results = dpath_fp_res, nm = file.path(dpath_fp_res,"NM"), log = file.path(dpath_fp_res,"Logs"), sindex = file.path(dpath_fp_res,"sindex"), coll_ind = file.path(dpath_fp_res,"coll_ind"))
		# Add a paths for bootstrap files
		if(!is.null(nboot)){
			paths = c(paths, boot = file.path(dpath_fp_res,"Boot"), boot_samp = file.path(dpath_fp_res,"Boot","samp"), boot_nm = file.path(dpath_fp_res,"Boot","NM"), boot_sindex = file.path(dpath_fp_res,"Boot","sindex"), boot_coll_ind = file.path(dpath_fp_res, "Boot","coll_ind"))
			# if nboot is not null then check that the 
		}
	# Check to see whether the directories needed to store data for this run exist (only a subset of these paths)
		req_paths = paths[grepl("^boot$|sindex|coll_ind",names(paths))]
		d_test = file.exists(req_paths)
		# If any of the directories for this stage of analysis don't exist then create them
		if(any(!d_test)){
			for(i in which(!d_test)){
				stopifnot(dir.create(req_paths[i]))
			}
		}
		# Cleanup old pathways
		rm(dpath_data,dpath_fp_res,d_test, req_paths)
		
	# Determine original working directory (will change back to it on exit)
		org_wd = getwd()
	# Change working directory to results folder
		setwd(paths["results"])
		on.exit(setwd(org_wd))
	# Save call
		temp = match.call()
		capture.output(print(temp),file="gai_call.txt")
		rm(temp)
	# Get list of species codes (if not supplied) or check species codes against data files
		spp_codes = sort(as.numeric(unique(gsub("^Sp ([[:digit:]]{1,})_.*$","\\1",list.files(paths["data"], pattern = "^Sp ")))))
		if( !is.null(spp) ){
			# print warning of any species codes given in spp which did not have data file(s)
			miss_inds = which(!spp %in% spp_codes)
			if(length(miss_inds) > 0){
				cat("WARNING: the following species codes supplied to spp argument did not have data files\n\t",paste(spp[miss_inds],collapse=","))
			}
			spp_codes = spp_codes[which(spp_codes %in% spp)]
		}
	# Create dataframe that will be used by nodes to determinewhich species bootstrap combinations they need to start on
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
	
	# Save run combinations to RDS file
		saveRDS(clust_combs, file=file.path(paths["results"],"gai_clust_params.rds"))
		# changes
			
	# Determine which other arguments need to be passed to the the cluster wrapper function (these arguments will be constant over all runs)
		# Determine what arguments were passed to function via dots
			dot_args = list(...)
			dot_argnames = names(dot_args)
		# Determine formals for function (use base function as clust wrapper has only 1 arugment i
			form_args = c(names(formals("gai_coll_ind")),"data_filter")
		# If save_sindex object passed to function ensure that it is passed to gai_coll_ind_clust function by adding it to form_args
			if("save_sindex" %in% dot_argnames){
				form_args = c(form_args,"save_sindex")
			}
		# Create list of fun_args from dot_args that match formals plus any other args that need adding (i.e. paths)
			fun_args = c(dot_args[which(dot_argnames %in% form_args)],list(paths = paths, yrs = yrs, nboot = nboot, batch_nboot = batch_nboot))
		# Save additional constant/fixed run parameters to RDS file
			saveRDS(fun_args, file=file.path(paths["results"], "gai_fixed_params.rds"))
		
	# Now run the est_fp_spline_clust function for all sets in the list
	if(parallel){
		# Setup a vector containing all objects functions that will need exported to the cluster nodes
		exp_objs = c("read_data_files","filter_data","clean_data","gai_coll_ind","gai_coll_ind_clust","gai_coll_ind_wrpr","read_boot_nm","f_do_call","boot_sample","find_nearest","miss_nm_yrs","nm_fill_miss","trap_index","trap_area","pred_SI_glm","check_nm")
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
				clusterEvalQ(cl,library(speedglm))
			# Export required objects/functions to clusters
				clusterExport(cl, varlist = exp_objs)
			# Use cluster version of apply to run fitting for all combinations
				ci_vals = do.call("rbind",parLapplyLB(cl,1:nrow(clust_combs),gai_coll_ind_clust))
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
				sfLibrary(speedglm)
			# Export required variables and functions to nodes
				sfExport(list = exp_objs)
			# Use load balancing cluster version of apply to run fitting for all combinations
				ci_vals = do.call("rbind",sfClusterApplyLB(1:nrow(clust_combs), gai_coll_ind_clust))
		}
	} else {
		require(speedglm)
		ci_vals = do.call("rbind",lapply(1:nrow(clust_combs),gai_coll_ind_clust))
	
	}
	
	# If bootstrapping return values with confidence intervals
	if(!is.null(nboot)){
		ci_vals = boot_conf(boot_data = ci_vals, stat_col = c("TR0OBS","TRMOBS","YR_COEF_CNT","TRCOUNT","TR_PRED_N_EST"), by_cols = c("SPECIES","YEAR"))
		for(sp in spp){
			ci_inds = which(ci_vals$SPECIES == sp)
			write.table(ci_vals[ci_inds,], file = file.path(paths["results"],paste0("Sp ",sp,"_coll_ind_conf.csv")),sep=",", row.names = FALSE, col.names = TRUE, na="", quote=TRUE)
		}
	}
	return(ci_vals)	
		
}

restart_gai = function(dpath_results, i_start = NULL, i_end = NULL, i_vec = NULL, parallel = TRUE, parallel_pkg = "snowfall",n_cpus = NULL){
	# Check dpath_results exists
	if(!file.exists(dpath_results))
		stop("Directory supplied to dpath_results does not exist")
	# Check that it contains a clust parameter file
	if(!file.exists(file.path(dpath_results,"gai_clust_params.rds")))
		stop("Supplied directory does not contain a gai_clust_params RDS file")
	# Determine original working directory (will change back to it on exit)
		org_wd = getwd()
	# Change working directory to results folder
		setwd(dpath_results)
		on.exit(setwd(org_wd))
	# Load the clust parameters file contain the variables/parameters for each run
		clust_params = readRDS("gai_clust_params.rds")
	# Determine the i values for which the gai is to be rerun/restarted
		if(is.null(i_start) & is.null(i_end) & is.null(i_vec)){
			run_i = 1:length(clust_params)
		} else if(!is.null(i_vec)){
			run_i = i_vec
		} else {
			run_i = ifelse(!is.null(i_start),i_start,1):ifelse(!is.null(i_end),i_end,length(clust_params))
		}
	# Setup cluster if parallel = TRUE
	if(parallel){
		# Setup a vector containing all objects functions that will need exported to the cluster nodes
		exp_objs = c("read_data_files","clean_data","gai_coll_ind","gai_coll_ind_clust","f_do_call","boot_sample","find_nearest","miss_nm_yrs","nm_fill_miss","trap_index","trap_area")
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
				clusterEvalQ(cl,library(speedglm))
			# Export required objects/functions to clusters
				clusterExport(cl, varlist = exp_objs)
			# Use cluster version of apply to run fitting for all combinations
				ci_vals = do.call("rbind",parLapplyLB(cl,run_i,gai_coll_ind_clust))
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
				sfLibrary(speedglm)
			# Export required variables and functions to nodes
				sfExport(list = exp_objs)
			# Use load balancing cluster version of apply to run fitting for all combinations
				ci_vals = do.call("rbind",sfClusterApplyLB(run_i, gai_coll_ind_clust))
		}
	} else {
		require(speedglm)
		ci_vals = do.call("rbind",lapply(run_i,gai_coll_ind_clust))
	}
	# Return results
	return(ci_vals)
}