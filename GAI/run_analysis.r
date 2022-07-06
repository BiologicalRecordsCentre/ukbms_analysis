# Load required packages
library(mgcv)
library(speedglm)

# Load in required functions from scripts
	# NOTE will need to set wd to location of these files or alter the paths below to point to their location
source("estimate_ci.r")
source("shared_funcs.r")
source("create_data_files.r")
source("estimate_fp.r")
source("bootstrapping.r")


# NOW YOU NEED TO LOAD IN YOUR COUNT DATA SO THAT THE FUNCTION CAN CREATE THE REQUIRED DATA.FILES

# Load in count data as data.frame the data frame should have columns SPECIES, SITE, YEAR, WEEKNO (or alternatively a column called DAYNO), COUNT. At present all of these fields should be numeric fields (e.g. numeric species and site codes) as this is the structure of the UKBMS data

# If you have any site index only type data (e.g. a single abundance value/metric for the species/year/site comb rather than a series of visit abundance values) and not the raw data. Examples of this type of data in the UKBMS are egg counts/timed counts that are done once a year at the peak of abundance. The site index data should be in a data.frame called si_data that has the same columns as above but where the COUNT column is replaced by a column called SINDEX that contains the site index value/metric

# RUN ANALYSIS

# Simple version for a single species using lower level functions and no parallel processing)
#-----------------------------------------------------------------------

# Note for this to work species data should be in data.frame sp_data and site index only data in si_data (if present)

# Stage 1 - Determine normalised flight periods for each year can alter defaults below to change way function works
	sp_nm = est_fp_spline_all(sp_data, verbose = FALSE, season_lim = c(1,26), add_anchor = TRUE, tp_col = "WEEKNO", min_visit = NULL, min_pos = NULL)
	
# Stage 2 - estimate collated index
	sp_ci = gai_coll_ind(sp_data = sp_data, nm_data = sp_nm, si_data = si_data, methodCI = "GLM_N", season_lim = c(1,26), tp_col = "WEEKNO")
	
	# Alternative version of stage 2 where the site indices calculated during stage 2 are outputted (save as an rds file)
	sp_ci = gai_coll_ind(sp_data = sp_data, nm_data = sp_nm, si_data = si_data, methodCI = "GLM_N", season_lim = c(1,53), tp_col = "WEEKNO", sindex_fname = "GAI_sindex_data.rds")
		# Now read the site index data back in from the RDS file produced above
			pred_si = readRDS("GAI_sindex_data.rds")
	
# END OF SIMPLE VERSION


	
# Version using higher levels functions to run analysis for all species (no bootstrapping)
#------------------------------------------------------------------------------------
	# Directory into which data files will be created
		data_dpath = "~/BMS data"
		# Create directory if it doesn't already exist
		if(!file.exists(data_dpath)){
			dir.create(data_dpath)
		}
	# Directory into which results of current analysis will saved
		out_dpath = "~/GAI Results"

	# Build species data files (as required by code)
		# Vector containing all species codes
			spp_codes = c(12,48) # replace with your species codes
		# Need to loop through species
		for(sp in spp_codes){
			# Code to read in species abundance data for current species (note section above detailing format that data needs to have SPECIES, SITE, YEAR, WEEKNO, COUNT. Needs to save this to an object called sp_data
				# INSERT YOUR CODE HERE
					sp_data = extract_visit_totals(spp = sp, years = 1976:2016, channel = channel)
				
			# Code to read in species site index only data (if any exists) into object called si_data otherwise set si_data to null
				# INSERT YOUR CODE HERE
					si_data = extract_sindex_only(spp = sp, years = 1976:2016, channel = channel)
					
			# Species count data (will create a single file for each species rather than year files which old script used to)
				create_data_files(sp_data = sp_data, dpath = data_dpath, yr_files = FALSE)

			
			# Site index only data
			if(!is.null(si_data)){
				if(nrow(si_data) > 0){
					create_data_files(sp_data = si_data, dpath = data_dpath, yr_files = FALSE)
				}
			}
		}
		
	# Stage 1 - Normalised flight period estimation
	all_nm = fit_fp_splines(dpath_data = data_dpath, dpath_out = out_dpath, spp = spp_codes, yrs = NULL, parallel = TRUE, parallel_pkg = "parallel", n_cpus = NULL, nboot = NULL, season_lim = c(1,26), add_anchor = TRUE, tp_col="WEEKNO")
	
	# Stage 2 - Estimate collated index
		# NOTE code below will attempt to run in parallel using all but one of your cpus which can
	all_ci = fit_gai(dpath_data = data_dpath, dpath_fp_res = out_dpath, spp = spp_codes, yrs = NULL, parallel = TRUE, parallel_pkg = "parallel", n_cpus = NULL, nboot = NULL, boot_only = FALSE, tp_col = "WEEKNO", season_lim = c(1,26), save_sindex = TRUE)


# Version using higher levels functions to run analysis for all species but this time including bootstrapping (will need species data files as created in above example and the code here only replaced the last few couple of lines in the previous example where stage 1 or stage 2 are performed

		# Stage 1 - Normalised flight period estimation
	all_nm = fit_fp_splines(dpath_data = data_dpath, dpath_out = out_dpath, spp = spp_codes, yrs = NULL, parallel = TRUE, parallel_pkg = "parallel", n_cpus = 2, nboot = 10, batch_nboot = 10, season_lim = c(1,26), add_anchor = TRUE, tp_col="WEEKNO")
	
	# Stage 2 - Estimate collated index
		# NOTE code below will attempt to run in parallel using all but one of your cpus which can
	all_ci = fit_gai(dpath_data = data_dpath, dpath_fp_res = out_dpath, spp = spp_codes, yrs = NULL, parallel = TRUE, parallel_pkg = "parallel", n_cpus = 2, nboot = 10, batch_nboot = 10, boot_only = FALSE, tp_col = "WEEKNO", season_lim = c(1,26), save_sindex = TRUE)
#------------------------------------------------------------------------------------