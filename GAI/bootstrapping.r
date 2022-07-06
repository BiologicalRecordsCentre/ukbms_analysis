create_boot_samples = function(sp, yrs = NULL, dpath_data, dpath_out = NULL, nboot = 100){
	# Read in raw data
		# Read in count data for species
			sp_counts = read_data_files(sp, yrs = yrs, data_type="count",dpath = dpath_data)
		# Read in any sindex only data
			sp_sindex = read_data_files(sp, yrs = yrs, data_type="sindex",dpath = dpath_data)
	# Determine all unique site codes across both datasets
		site_codes = unique(c(sp_counts$SITE, sp_sindex$SITE))
	# Determine the number of sites
		n_sites = length(site_codes)
	if(n_sites > 1){
		# Create nboot number of bootstrap samples of the site codes
			temp_samps = matrix(sample(site_codes, n_sites*nboot, replace = TRUE), nrow = n_sites, ncol = nboot, dimnames = list(NULL,paste("iboot",1:nboot,sep="_")))
		# Add new site codes to start and store as data.frame
			all_samp = data.frame(site = 1:n_sites, temp_samps)
		# Save single bootstrap file for each species (rather than 1 per species/bootstrap as was previously done)
			saveRDS(all_samp, file = file.path(dpath_out,paste0("Sp ",sp,"_bootsamps.rds")))
	} else {
		all_samp = NULL
	}
	# Return invisbly
		invisible(all_samp)
}

boot_sample = function(sp_data, boot_samp){
	if(!is.null(sp_data)){
		if(nrow(sp_data) > 0){
			# Get orignal order of columns in sp_data and add ORG_SITE TO end
				col_names = c(names(sp_data),"ORG_SITE")
			# Rename the site column in the sp data to ORG_SITE
				names(sp_data)[grepl("SITE",names(sp_data))] = "ORG_SITE"
			# Merge the datasets TO ad the new_sites column to the 
				out_data = merge(sp_data, boot_samp, by = "ORG_SITE", sort = FALSE)
			# reorder columsn to match original order
				out_data = out_data[,col_names]			
		} else {
			out_data = NULL
		}
	} else {
		out_data = NULL
	}
	# Return boot samp version of sp_data
		return(out_data)
}

boot_conf = function(boot_data, stat_col, by_cols = c("SPECIES","YEAR"), probs = c(0.025, 0.975)){
	# Does the data contain the real values or only bootstrap values
	real_inds = which(boot_data$iboot == 0)
	if(length(real_inds) > 0){
		real_data = boot_data[real_inds,!grepl("iboot",names(boot_data))]
	} else {
		real_data = NULL
	}
	# Which indices correspond to bootstrap samples
		boot_data = boot_data[which(boot_data$iboot > 0),]
	# If real data is not null then realign bootstrap data so that the mean of each bootstrap CI series aligns with mean of real data series covering same time period (this is to control for bootsrap samples that cover only a subset of the years covered by the full data)
	if(!is.null(real_data)){
		boot_data = align_boot_ci(boot_data = boot_data, real_data = real_data)
	}
	# Add row.names to real_data where row.names = by_cols seperate by _
		row.names(real_data) = apply(real_data[,by_cols],1,paste, collapse="_")
		# Add new columns to hold output data for each stat col
			new_cols = paste(rep(c("LOW","UPP","NBOOT"),length(stat_col)),rep(stat_col,each = 3), sep="_")
			real_data[,new_cols] = NA
	# Loop through cur_stat and calculate conf ints	
	for(cur_stat in stat_col){
		# Determine the confidence intervals
			temp_ci = aggregate(boot_data[,cur_stat], boot_data[,by_cols], quantile, probs = probs, simplify = TRUE, type=8, na.rm = TRUE) 
			temp_n = aggregate(boot_data[,cur_stat], boot_data[,by_cols], function(x){x = na.omit(x); length(x)}, simplify = TRUE)
		# Seperate x column into two columns
			temp_ci[,c(paste("LOW",cur_stat,sep="_"),paste("UPP",cur_stat,sep="_"))] = temp_ci[,"x"]
		# Now remove x column
			temp_ci = temp_ci[,!grepl("^x$",names(temp_ci))]
		# Add in the n_boot
			temp_ci[,paste0("NBOOT_",cur_stat)] = temp_n$x
		# Add key statistics from temp_ci to real_data
			# get columns names of columns to be added (any but by_cols which are already in real_data)
				add_cols = names(temp_ci)[!names(temp_ci) %in% by_cols]
			# Now add these columns to the correct rows of real_data
				real_data[apply(temp_ci[,by_cols],1,paste, collapse="_") ,add_cols] = temp_ci[,add_cols]
	}
	# Return results
		return(real_data)
}


align_ci_base = function(bootsamp_res, real_res){
	# Determine year limits for real data
		real_range = range(real_res$YEAR)
	# Determine boot sample year limits
		boot_range = range(bootsamp_res$YEAR)
	# if either value from boot_range is not matched by the value in real_range then need to realign the data
	if(any(!boot_range == real_range)){
		# Calculate mean of real data TR0OBS corresponding to boot_range
			sub_means = colMeans(real_res[which(real_res$YEAR >= boot_range[1] & real_res$YEAR <= boot_range[2]),c("TR0OBS","TRMOBS")])
			boot_means = colMeans(bootsamp_res[,c("TR0OBS","TRMOBS")])
		# Rescale bootsamp CI values so that boot_means = sub_means
			bootsamp_res[,"TR0OBS"] = bootsamp_res[,"TR0OBS"] + (sub_means["TR0OBS"] - boot_means["TR0OBS"])
			bootsamp_res[,"TRMOBS"] = bootsamp_res[,"TRMOBS"] + (sub_means["TRMOBS"] - boot_means["TRMOBS"])
	}
	# Return the aligned bootsamp_res
		return(bootsamp_res)
}


align_boot_ci = function(boot_data, real_data = NULL, ci_col = c("TR0OBS","TRMOBS")){
	# If no real data supplied then needs to be included in boot_data with iboot value of 0
	if(is.null(real_data)){
		inds = which(boot_data$iboot == 0)
		if(length(inds) > 0){
			# Extract real data from boot_data
			real_data = boot_data[inds,]
			boot_data = boot_data[-inds,]
		} else {
			stop("Must supply the real results either directly through real_data argument or via boot_data where they can have an iboot value of 0")
		}
	}
	# Apply align_ci_base function to all bootstraps
	aligned_boot = do.call("rbind",by(boot_data, boot_data$iboot, align_ci_base, real_res = real_data))
	# Return aligned boot_data
	return(aligned_boot)
}