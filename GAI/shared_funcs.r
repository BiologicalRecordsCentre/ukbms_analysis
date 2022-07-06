setupDBcon = function(dsn = "BMS", uid = "BMS", pwd = NULL){
	# If using RStudio on Windows then the Microsoft ODBC password dialog box that would normally appear is prevented, instead a normal text box will be used to collect password from user at run time rather than writing it into the script
	if(.Platform$GUI == "RStudio" & is.null(pwd)){
		bms_pass = winDialogString("Enter BMS password", "") # Force a windows dialog box to appear for user to type password
		ret_obj = odbcConnect(uid, pwd = bms_pass)
	} else {
		if(!is.null(pwd)){
			ret_obj = odbcConnect(dsn, uid, pwd = pwd)
		} else {
			ret_obj = odbcConnect(dsn, uid)
		}
	}
	return(ret_obj)
}

extract_visit_totals = function(spp = NULL, years = NULL, dataset = "BOTH", add_where = NULL, add_out_cols = NULL, channel = NULL){
	if(!is.null(add_where)){
		add_where = paste(" and",add_where)
	}
	# Construct where conditions
	if(!is.null(spp)){
		if(length(spp) == 1){
			add_where = paste(add_where, " and sv.species = ", spp, sep="")
		} else {
			add_where = paste(add_where, " and sv.species in (", paste(spp, collapse = ","),")", sep="")
		}
	} 
	if(!is.null(years)){
		if(length(years) == 1){
			add_where = paste(add_where, " and v.year = ", years, sep="")
		} else {
			# If all years are sequential with no gap/jump then use between otherwise use in statement
			if(all(diff(years)) == 1){
				add_where = paste(add_where, " and v.year between ",min(years)," and ", max(years), sep="")
			} else {
				add_where = paste(add_where, " and v.year in (", paste(years, collapse=","),")", sep="")
			}
		}
	}
	
	# Construct base query
		temp_sql = paste(
			"select 
				sv.species, 
				sv.siteno as site, 
				v.year, 
				v.month, 
				v.day, 
				v.weekno,
				date_to_weekno(v.day, v.month, v.year, i_output => 'days') as dayno,
				round(avg(nvl(vt.visit_total, 0))) as count", ifelse(!is.null(add_out_cols),paste0(", ",add_out_cols),""),
			" from ukbms_visit v
				join
					(
					select vt.species,vt.siteno,min(vt.year) as first_yr
					from ukbms_visit_total vt
						join ukbms_visit v on vt.siteno = v.transect_code and vt.year = v.year and vt.month = v.month and vt.day = v.day
					where v.non_crit is null group by vt.species, vt.siteno
					) sv on v.transect_code = sv.siteno and v.year >= sv.first_yr 
				left join ukbms_visit_total vt on sv.species = vt.species and sv.siteno = vt.siteno and v.year = vt.year and v.month = vt.month and v.day = vt.day
				left join ukbms_site st on sv.siteno = st.siteno
				left join ukbms_species sp on sv.species = sp.bmscode
			where V.NON_CRIT is null", add_where,  
			" group by sv.species, sv.siteno, v.year, v.month, v.day, v.weekno, date_to_weekno(v.day, v.month, v.year, i_output => 'days')", ifelse(!is.null(add_out_cols),paste0(", ",add_out_cols),""),
			" order by species, site, year, dayno", sep="")
	
	# Extract species data for spp
	if(dataset %in% c("UKBMS","BOTH")){
		spp_data = sqlQuery(channel, temp_sql, stringsAsFactors = FALSE)
	}
		
	# Modify query to extract WCS data if required
	if(dataset %in% c("WCS","BOTH")){
		# Modify sql
			# Table names
				temp_sql = gsub("ukbms_visit","ukbms_wcbs_visit", temp_sql)
			# Transect code
				temp_sql = gsub("transect_code","siteno", temp_sql)
		
		# Extract data
		if(dataset == "WCS"){
			spp_data = sqlQuery(channel, temp_sql, stringsAsFactors = FALSE)
		} else {
			
			spp_data = rbind(spp_data, sqlQuery(channel, temp_sql, stringsAsFactors = FALSE))
		}
		
	}
			  
	# Return data
	return(spp_data)
}

extract_sindex_only = function(spp = NULL, years = NULL, add_where = NULL, add_out_cols = NULL, channel = NULL){
	# If connection via channel has not been supplied then open a new connection for this function
	if(is.null(channel)){
		channel = setupDBcon()
		on.exit(odbcClose(channel)) # If the connection was opened by this function then close it when the function exits
	}
	# Construct where conditions
		if(!is.null(add_where)){
			add_where = paste(" and",add_where)
		}
		if(!is.null(spp)){
			if(length(spp) == 1){
				add_where = paste(add_where, " and si.species = ", spp, sep="")
			} else {
				add_where = paste(add_where, " and si.species in (", paste(spp, collapse = ","),")", sep="")
			}
		} 
		if(!is.null(years)){
			if(length(years) == 1){
				add_where = paste(add_where, " and si.year = ", years, sep="")
			} else {
				if(all(diff(years)) == 1){
					add_where = paste(add_where, " and si.year between ",min(years)," and ", max(years), sep="")
				} else {
					add_where = paste(add_where, " and si.year in (", paste(years, collapse=","),")", sep="")
				}
			}
		}
	# Build SQL query
		temp_sql = paste(
			"select si.species, si.year, si.site, round(avg(sindex)) as sindex",
			ifelse(!is.null(add_out_cols),paste0(", ",add_out_cols),""),
			" from ukbms_sindex_only si
				join ukbms_site st on si.site = st.siteno
				join ukbms_species sp on si.species = sp.bmscode
			where si.sindex is not null and si.sindex >= 0 and brood = 0", add_where,  
			" group by si.species, si.year, si.site", 
			ifelse(!is.null(add_out_cols),paste0(", ",add_out_cols),""),
			" order by si.species, si.year, site", sep="")
	# Extract data
		si_data = sqlQuery(channel, temp_sql, stringsAsFactors = FALSE)
	
	# Return data
	return(si_data)
}


# Function used to calculate the trapezoidal area of a set of x,y values
trap_area = function(x,y = NULL){
	# If y is null and x has multiple columns then set y to x[,2] and x to x[,1]
		if(is.null(y)){
			if(length(dim(x)) == 2){
				y = x[,2]
				x = x[,1]
			} else {
				stop("ERROR: need to either specify both x and y or supply a two column data.frame/matrix to x")
			}
		}
		
	# Check x and y are same length
		if(length(x) != length(y)){
			stop("ERROR: x and y need to be the same length")
		}
	
	# Need to exclude any pairs that are NA for either x or y
		rm_inds = which(is.na(x) | is.na(y))
		if(length(rm_inds) > 0){
			x = x[-rm_inds]
			y = y[-rm_inds]
		}
		
		
	# Determine values of trapezoids under curve
		# Get inds
			inds = 1:(length(x)-1)
		# Determine area using trapezoidal rule Area = ( (b1 + b2)/2 ) * h where b1 and b2 are lengths of bases (the parallel sides) and h is the height (the perpendicular distance between two bases)
			areas = ( ( y[inds] + y[inds+1] ) / 2 ) * diff(x)
	
	# total area is sum of all trapezoid areas
	tot_area = sum(areas)
	

	# Return total area
	return(tot_area)
}

# Modified version of trap_index function from 2stage_stage2.r (modified to work with a single by column and also to remove the sindex values being divided by 7
trap_index = function(sp_data, data_col = "IMP", time_col = "DAYNO", by_cols = c("SPECIES","SITE","YEAR")){
	
	# Build output data.frame (use drop = FALSE to stop 1 column data.frame being forced to a vector)
		out_obj = unique(sp_data[,by_cols, drop=FALSE])
		# Set row.names to be equal to collapsing of output rows (will be unique, you need them to make uploading values back to data.frame will be easier)
		if(length(by_cols) == 1){
			row.names(out_obj) = out_obj[,by_cols]
		} else {
			row.names(out_obj) = apply(out_obj, 1, paste, collapse = "_")
		}
			
	# Determine which rows are site index only data as these will need treated separately
		#si_inds = which(sp_data$SITE <= 0)
		si_inds = NULL
		# Create vector of row.names from sp_data
		if(length(by_cols) == 1){
			temp_rn = as.character(sp_data[,by_cols])
		} else {
			temp_rn = apply(sp_data[,by_cols], 1, paste, collapse="_")
		}
		
	# Using this row.names from out_obj above as index in by function to loop through values all unique combs of by_colss and fit trap_area to data
		if(length(si_inds) > 0){
			ind_dat = by(sp_data[-si_inds,c(time_col,data_col)], temp_rn[-si_inds], trap_area)
		} else {
			ind_dat = by(sp_data[,c(time_col,data_col)], temp_rn, trap_area)
		}
	
	# Add this data to output object
		out_obj[names(ind_dat), "SINDEX"] = round(ind_dat, 1)
	
	# Now copy site index only data across (if any)
		if(length(si_inds) > 0){
			out_obj[temp_rn[si_inds],"SINDEX"] = sp_data$COUNT[si_inds]
		}
	
	# Set row.names to defaults
		row.names(out_obj) = NULL
		
	# Return output object
	return(out_obj)
}

f_do_call = function(fun, fun_args){
	# Check fun is a function that exisits
		if(!exists(fun,mode = "function")){
			stop("name passed to fun must be a valid function name")
		}
		
	# Check arguments in dots argument and compare against the formula for the specified function
		# Investigate dot arguments
			fun_argnames = names(fun_args)
		# Determine which dot arguments apply to the function supplied to fun
			match_args = intersect( fun_argnames, names(formals(fun)))
			fun_args = modifyList(as.list(formals(fun)[match_args]), fun_args[match_args])
	
	# Now call function fun supplying arguments
		do.call(fun, fun_args)
}

check_nm = function(sp_data,sp_nm, max_nm = 0.7, min_overlap = 0.3, tp_col = "WEEKNO", season_lim = c(1,26)){
	# Create small wrapper function used later to check nm across all year/bootstrap combinations
	nm_overlap = function(cur_nm, cnt_nm){
		prop_over = (2 - sum(abs(cur_nm - cnt_nm)))/2
		max_nm = max(cur_nm)
		ret_obj = c(prop_over = prop_over, max_nm = max_nm)
		return(ret_obj)
	}
	# Get indices of sp_data that are withing flight period
		inds = which(sp_data$WEEKNO >= season_lim[1] & sp_data$WEEKNO <= season_lim[2])
	# Estimate a normalised flight period directly from counts by calculating weekly sums and then dividing by number of unique sites in data (not of that week, as this can cause issues with dodgy records outside of normal flight period for species e.g. if only 1 record ever comes from week 52 and it is a moderate count it can end up being the highest value if you take the mean for each week. Likewise the total sum can be biased by 1 single large count
		wk_sum = tapply(sp_data$COUNT[inds], sp_data[inds,tp_col], sum, na.rm = TRUE)
		tot_n_sites = length(unique(sp_data$SITE[inds]))
		wk_avg = wk_sum/tot_n_sites
	# Setup object to hold normalised count flight period that spans entire season (may be weeks/days without counts)
		cnt_nm = setNames(rep(0,length(season_lim[1]:season_lim[2])),season_lim[1]:season_lim[2])
	# Normalise weekly mean count
		cnt_nm[names(wk_avg)] = wk_avg / sum(wk_avg)
	# figure out which grouping columns to work over when applying nm checks (e.g. year or year & iboot)
		grp_cols = names(sp_nm)[names(sp_nm) %in% c("YEAR","iboot")]
	# Now calculate prop overlap & max nm for each combination of grouping columns
		prop_over = aggregate(sp_nm[,"NM",drop=FALSE], sp_nm[,grp_cols,drop=FALSE], FUN = nm_overlap, cnt_nm = cnt_nm)
	# Find combinations where values don't meet criteria
		fail_inds = which(prop_over$NM[,"prop_over"] < min_overlap | prop_over$NM[,"max_nm"] > max_nm | is.nan(prop_over$NM[,"prop_over"]) | is.nan(prop_over$NM[,"max_nm"]) )
	# Find which rows in sp_nm correspond to these dodgy year/boot combinations
		rm_inds = which(apply(sp_nm[,grp_cols,drop=FALSE],1,paste,collapse="_") %in% apply(prop_over[fail_inds,grp_cols,drop = FALSE],1,paste,collapse="_"))
	# If any combinations fail then exclude then from sp_nm before returning remaining
	if(length(rm_inds) > 0){
		sp_nm = sp_nm[-rm_inds,]
	}
	return(sp_nm)
}

check_nm_old = function(sp_data,sp_nm, max_nm = 0.7, min_overlap = 0.3, tp_col = "WEEKNO", season_lim = c(1,26)){
	# Create small wrapper function used later to check nm across all year/bootstrap combinations
	nm_overlap = function(cur_nm, cnt_nm){
		prop_over = (2 - sum(abs(cur_nm - cnt_nm)))/2
		max_nm = max(cur_nm)
		ret_obj = c(prop_over = prop_over, max_nm = max_nm)
		return(ret_obj)
	}
	# Get indices of sp_data that are withing flight period
		inds = which(sp_data$WEEKNO >= season_lim[1] & sp_data$WEEKNO <= season_lim[2])
	# Estimate a normalised flight period directly from counts by estimating mean weekly/daily count across all years
		wk_cnt = tapply(sp_data$COUNT[inds], sp_data[inds,tp_col], mean, na.rm = TRUE)
	# Setup object to hold normalised count flight period that spans entire season (may be weeks/days without counts)
		cnt_nm = setNames(rep(0,length(season_lim[1]:season_lim[2])),season_lim[1]:season_lim[2])
	# Normalise weekly mean count
		cnt_nm[names(wk_cnt)] = wk_cnt / sum(wk_cnt)
	# figure out which grouping columns to work over when applying nm checks (e.g. year or year & iboot)
		grp_cols = names(sp_nm)[names(sp_nm) %in% c("YEAR","iboot")]
	# Now calculate prop overlap & max nm for each combination of grouping columns
		prop_over = aggregate(sp_nm[,"NM",drop=FALSE], sp_nm[,grp_cols,drop=FALSE], FUN = nm_overlap, cnt_nm = cnt_nm)
	# Find combinations where values don't meet criteria
		fail_inds = which(prop_over$NM[,"prop_over"] < min_overlap | prop_over$NM[,"max_nm"] > max_nm)
	# Find which rows in sp_nm correspond to these dodgy year/boot combinations
		rm_inds = which(apply(sp_nm[,grp_cols,drop=FALSE],1,paste,collapse="_") %in% apply(prop_over[fail_inds,grp_cols,drop = FALSE],1,paste,collapse="_"))
	# If any combinations fail then exclude then from sp_nm before returning remaining
	if(length(rm_inds) > 0){
		sp_nm = sp_nm[-rm_inds,]
	}
	return(sp_nm)
}
