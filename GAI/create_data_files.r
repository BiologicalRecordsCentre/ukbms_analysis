create_data_files = function(sp_data, dpath = NULL, yr_files = TRUE, verbose = FALSE){
	# Get species code
		sp_code = unique(sp_data$SPECIES)
	# If more than 1 species then stop
		stopifnot(length(sp_code) == 1)
	# Determine data type (COUNT, SINDEX or NM)
		col_ind = which(names(sp_data) %in% c("COUNT","SINDEX","NM"))
	if(length(col_ind) == 1){
		data_type = tolower(names(sp_data)[col_ind])
	} else {
		stop("sp_data must included one of the following columns; COUNT, SINDEX, NM")
	}
	# Datermine which yeara are in the dataset
		all_yrs = unique(sp_data$YEAR)
	# If dpath NULL then use working directory
	if(is.null(dpath)){
		dpath = getwd()
	}
	# Show progress
	if(verbose)
		cat("Creating data file(s) for Species ",sp_code, "\n  File(s) in directory: ",dpath,"\n", sep="")
	# If seperate files for each year
	if(yr_files){
		# Loop through years
		for(cur_yr in all_yrs){
			# Show progress
			if(verbose)
				cat("  Creating file for ",cur_yr," data\n",sep="")
			# Build file name
			f_name = sprintf("Sp %s_Yr %s_%s.rds",as.character(sp_code),as.character(cur_yr), data_type)
			# Save data for current year to file
			saveRDS(sp_data[which(sp_data$YEAR == cur_yr),], file.path(dpath, f_name))
		}
	} else {
		# Build file name
		f_name = sprintf("Sp %s_%s.rds",as.character(sp_code),data_type)
		# Save data for species to file
		saveRDS(sp_data, file.path(dpath, f_name))
	}
	
	# Show progress
	if(verbose)
		cat("File(s) created\n\n",sep="")
}

file_sp_yr_combs = function(f_list){
	# Setup list to hold results from each file
		all_combs = vector("list",length(f_list))
	# Loop through files
	for(i_f in 1:length(f_list)){
		temp = readRDS(f_list[i_f])
		all_combs[[i_f]] = unique(temp[,c("SPECIES","YEAR")])
	}
	# collapse list to data.frame
		sp_yr_out = do.call("rbind",all_combs)
	# Order by species then year
		sp_yr_out = sp_yr_out[order(sp_yr_out$SPECIES,sp_yr_out$YEAR),]
	# Reset row.names
		row.names(sp_yr_out) = NULL
	# Save species year combs as a RDS file (so they don't have to be calculated again)
		saveRDS(sp_yr_out, file=file.path(dirname(f_list[1]), "all_sp_yr_combs.rds"))
	# Invisibly return sp_yr_out
		invisible(sp_yr_out)
}


filter_data = function(dataset, filter_arg){
	if(class(filter_arg) == "character"){
		# Stop if character contains a vector of character strings
		stopifnot(length(filter_arg) == 1)
		filt_code = filter_arg
	} else if(class(filter_arg) == "list"){
		# Check all names in data filter are columns in the dataset
			stopifnot(all(names(filter_arg) %in% names(dataset)))
		# Determine length of filter_arg list
			n_filt = length(filter_arg)
		# Setup object to store filter elements
			df_code = vector("character",n_filt)
		# Loop through length of data_filter and build code that will be used to subset			
		for(i_df in 1:n_filt){
			# Get length of current filter elements
				cur_len = length(filter_arg[[i_df]])
				stopifnot(cur_len >= 1)
			if(cur_len == 1){
				df_code[[i_df]] = paste(names(filter_arg)[i_df],shQuote(filter_arg[[i_df]]), sep=" == ")
			} else if(cur_len > 1){
				df_code[[i_df]] = paste0(names(filter_arg)[i_df]," %in% c(", paste(shQuote(filter_arg[[i_df]]), collapse=","),")")
			}
		}
		# Now build all components together to make full filter expression
			filt_code = paste(df_code, collapse=" & ")
	
	} else {
		stop("filter_arg needs to be a character string or a list")
	}
	# Use completed filter expression to subset data
		dataset = subset(dataset,eval(parse(text = filt_code)))
	# Return filtered data
		return(dataset)
}

read_data_files = function(sp, data_type, yrs = NULL, iboot = NULL, dpath = NULL, data_filter = NULL){
	# If dpath is null then set to working directory
	if(is.null(dpath)){
		dpath = getwd()
	}
	# Find all data type files for species (which may or may not have years in it
		f_pat = f_pat = paste0("Sp ",as.character(sp),"(_Yr [[:digit:]]{4})?_",ifelse(!is.null(iboot),paste0("iboot ",iboot,"_"),""),data_type,"[.]rds")
	# Get list of files in dpath that match pattern
		f_list = list.files(dpath, pattern = f_pat)
	# Determine if species data is in a single file or split across seperate files for each year
	if(all(grepl("_Yr ([[:digit:]]{4})",f_list))){
		# Species year files
			yr_files = TRUE
		# If more than 1 year supplied to yrs then reduce f_list to only keep those that are in years
		if(length(yrs) >= 1){
			# Extract years from file names
				f_yr = gsub(".*_Yr ([[:digit:]]{4}).*$","\\1",f_list)
			# Find which files are within the years specified
				inds = which(f_yr %in% yrs)
			# Keep only those which are within the years listed
				f_list = f_list[inds]
			# Cleanup
				rm(f_yr,inds)
		}
	} else {
		yr_files = FALSE
	}
	
	if(length(f_list) == 1){
		sp_data = readRDS(file.path(dpath, f_list[1]))
	} else if(length(f_list) == 0){
		sp_data = NULL
	} else {
		sp_data = do.call("rbind",lapply(file.path(dpath,f_list),readRDS))
	}
	
	# If yrs supplied and yr_files = FALSE then limit data to only the required years
	if(!is.null(yrs) & !yr_files){
		sp_data = sp_data[which(sp_data$YEAR %in% yrs),]
		if(nrow(sp_data) == 0){
			sp_data = NULL
		}
	}
	
	# If values supplied to data_filter then subset data
	if(!is.null(data_filter) & !is.null(sp_data)){
		sp_data = filter_data(sp_data, filter_arg = data_filter)
	}
	
	# Reset row names
		row.names(sp_data) = NULL
	# Return
	return(sp_data)
}

clean_data = function(sp_data, tp_col = "WEEKNO", add_cols = NULL, season_lim = c(1,26)){
	# Remove any negative counts (shouldn't be any but better safe!)
	rm_inds = which(sp_data$COUNT < 0)
	if(length(rm_inds) > 0){
		sp_data = sp_data[-rm_inds,]
	}
	# Create vector of data by which the sp_data is to be reduced
	by_cols = c("SPECIES","SITE","YEAR",tp_col, add_cols)
	# Now determine the max value for each unique combination of the by_cols
	temp = aggregate(list(COUNT = sp_data[,"COUNT"]), sp_data[,by_cols], FUN = max, na.rm = TRUE)
	temp = temp[order(temp$YEAR, temp$SITE, temp[,tp_col]),]
	row.names(temp) = NULL
	# Remove any counts outwith season
	if(!is.null(season_lim)){
		stopifnot(season_lim[1] < season_lim[2])
		rm_inds = which(temp[,tp_col] < season_lim[1] | temp[,tp_col] > season_lim[2])
		if(length(rm_inds) > 0){
			temp = temp[-rm_inds,]
		}
	}
	# Remove any sites which never have a positive count within season
	site_totals = tapply(temp$COUNT, temp$SITE, sum)
	rm_sites = names(which(site_totals == 0))
	if(length(rm_sites) > 0){
		rm_inds = which(temp$SITE %in% rm_sites)
		temp = temp[-rm_inds,]
	}
	return(temp)
}

create_UKBMS_files = function(spp_list = NULL, years = NULL, dpath = NULL, channel = NULL, dataset = "BOTH",  yr_files = TRUE, add_where = NULL, add_out_cols = NULL, verbose = FALSE){
	require(RODBC)
	if(verbose)
		cat("Extracting data from UKBMS database and creating data files:\n")
	if(is.null(dpath)){
		# If no dpath supplied then use a folder Visit Data within working directory (creating it where necessary)
		dpath = file.path(getwd(), "Visit_Data")
		if(!file.exists(dpath)){
			dir.create(dpath)
		}
	}
	if(is.null(channel)){
		# If no channel passed then open new connection and then close when funciton finishes
		channel = setupDBcon()
		on.exit(odbcClose(channel))
	}
	# If spp_list not supplied then search species table to get full list of species for which analysis is conducted
	if(is.null(spp_list)){
		spp_list = sqlQuery(channel, "select bmscode from ukbms_species where sp_for_anal = 1 order by bmscode")[,1]
	}
	# Loop through species and extract data from UKBMS database
	for(sp in spp_list){
		# Show progress
		if(verbose)
			cat("  Species ", sp,"\n")
		# Extract count data
			sp_data = extract_visit_totals(sp, years = years, channel = channel, dataset = dataset, add_where = add_where, add_out_cols = add_out_cols)
		# Create count data files
		if(nrow(sp_data) > 0){
			create_data_files(sp_data, dpath = dpath, yr_files = yr_files)
		}
		# Cleanup
			rm(sp_data)
		# Extract site index only data
			si_data = extract_sindex_only(sp, years = years, channel = channel, add_where = add_where, add_out_cols = add_out_cols)
		# Create sindex only data files
		if(nrow(si_data) > 0){
			create_data_files(si_data, dpath = dpath, yr_files = yr_files)
		}	
	}
	
	return(dpath)
} 
