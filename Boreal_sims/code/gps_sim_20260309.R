library(tidyverse)
library(terra)
library(glmmTMB)
library(reproducible)
library(foreach)
library(doParallel)
library(data.table)

source(file.path("code", "sim_paths.R"))

# Read stuff in ####
data_path = file.path('inputs')
if (!dir.exists(data_path)) {
  dir.create(data_path)
  
  ntland = prepInputs(url = 'https://drive.google.com/file/d/1Xp2x3kTX0HteuBmio-tG1mOHoQvml39C/view?usp=share_link',
                      fun = 'terra::rast',
                      destinationPath = data_path,
                      targetFile = 'ntland.tif')
  
  # # Peter note 2026/03/09 04:39 PM PDT: These might need to be terra::wrap()'d for the sharing as RDS file to work. Make a note about this for Julie tomorrow
  # NTyearly = prepInputs(url = 'https://drive.google.com/file/d/1eprxug3JCNf3c2q7pItUszzjaGu1_fC9/view?usp=drive_link',
  #                       destinationPath = data_path,
  #                       fun = 'readRDS')
  # 
  # # Peter note 2026/03/09 04:40 PM PDT: Similar predicament here except maybe even worse because I can't even run str() on it. Difference may be because this one is a SpatVectorCollection and the other is a SpatRaster?
  # NT5yearly = prepInputs(url = 'https://drive.google.com/file/d/19Srjt6yKTM-lJAYfzG5hWnh6f9AQ0bT7/view?usp=drive_link',
  #                        destinationPath = data_path,
  #                        fun = 'readRDS')
  
  # Peter note 2026/03/09 04:42 PM PDT: This is weirdly a list of length 1028, but the values (kappa, mu, shape, scale) appear to just be repeated. The list names have a bunch of random gibberish on there. Not sure what's up but can certainly at least just reduce to the first 4 elements?
  NTdistparams = prepInputs(url = 'https://drive.google.com/file/d/1QGTdkrx_FKgmtsMwKly53uZXBzYcL44i/view?usp=drive_link',
                            destinationPath = data_path,
                            fun = 'readRDS')
  
  # Peter note 2026/03/09 04:49 PM PDT: Welp, this one just completely fails. It gives the following error message: "error reading from connection". Did some Googling and got very little help. Questions to ask Julie for troubleshooting: 1) what packages did she have loaded when making this? Is there anything other than a model data table in there? Maybe if I load some other packages it'll work? 2) Does it work on her computer (if she downloads from GDrive)? What about if she just reads the file that's presumably on her computer? Could it be some sort of corruption from the transfer?
  ntmodel = prepInputs(url = 'https://drive.google.com/file/d/11td8Wbn2m_TbHVJJppgpfpuiE2cmQ1nP/view?usp=share_link',
                       destinationPath = data_path,
                       fun = 'load',
                       targetFile = 'ntmodel.Rdata')
  ntmod = ntmodel$ntmod
  
  # Peter note 2026/03/09 04:52 PM PDT: The only one that actually worked! Miraculous
  ntstudyarea = prepInputs(url = 'https://drive.google.com/file/d/1YOsRhBImlNuoAU4Jkdkz9tMeTfPaF_jq/view?usp=drive_link',
                           destinationPath = data_path)
} else {
  
  ntland = rast(file.path(data_path, "ntland.tif"))
  NTdistparams = readRDS(file.path(data_path, "NTdistparams.rds"))
  ntstudyarea = vect(file.path(data_path, "ntstudyarea.shp"))
  load(file.path(data_path, "ntmodel.RData"))
  
}


# If the directory already exists, assume we're running this and have already read in the files.
NTdistparams = NTdistparams[1:4] # for some reason it adds extra parametrs on but they appear to be all the same, so I'm just going to keep the first four

# Prep the arguments for the function ####

## Step length and turning angle distributions ####
# Helper function to extract values from NTdistparams, which has some weird list names but not to worry
get_value_by_name = function(obj, name) as.numeric(obj[str_detect(names(obj), name)])

# named and ordered according to documentation in sim_paths.R
this_move_pars = c(gamma_shape = get_value_by_name(NTdistparams, "shape"), 
                   gamma_rate = 1 / get_value_by_name(NTdistparams, "scale"), 
                   vm_mu = get_value_by_name(NTdistparams, "mu"), 
                   vm_kappa = get_value_by_name(NTdistparams, "kappa"))

## Covariates: the function I've written takes the covariate rasters pre-processed (i.e., pre-multiplied by the beta's to make a series of SSF "maps"), separated by interactions with steplength, turning angle, etc. This saves a bit of time when extracting the values.
# Helper function (for now) to deal with NA values in data
fix_na = function(x, replace_val = 0) ifelse(is.na(x), replace_val, x)

# get all beta coefficients from model
ntmod_sum = summary(ntmod)
beta_all = ntmod_sum$coefficients$cond[, 1] %>% fix_na

lsl_var = str_detect(names(beta_all), "log\\(sl")
cta_var = str_detect(names(beta_all), "cos\\(ta")
lsl_cta_var = lsl_var & cta_var # interaction between log(sl) and cos(ta)
base_var = !(lsl_var | cta_var) # anything without sl or ta is just a regular old beta
# remove interaction from lsl and cta variables
lsl_var = lsl_var & !lsl_cta_var
cta_var = cta_var & !lsl_cta_var

# Helper function that takes model outputs and multiplies them onto rasters to get desired iSSF outputs
get_metric_raster = function(include_vars,
                             start_or_end = c("start", "end"),
                             metric_regex = c("I\\(log\\(sl_ \\+ 1\\)\\)(:?)", "I\\(cos\\(ta_\\)\\)(:?)")) {
  
  start_or_end = paste0("_", match.arg(start_or_end))
  metric_regex = match.arg(metric_regex)
  
  beta_this = beta_all[include_vars]
  # remove "_start" because they all have it
  beta_this_names = str_replace_all(names(beta_this), start_or_end, "") %>% 
    # also remove the step length interaction part so we can focus on what the covariate actually is
    str_replace_all(metric_regex, "")
  beta_this_names = ifelse(nchar(beta_this_names) == 0, "intercept", beta_this_names)
  names(beta_this) = beta_this_names
  
  ntland_manip = do.call(c, lapply(beta_this_names, function(nm) {
    
    message("beginning ", nm) # just for debugging because this takes a few seconds - I know not the most efficient way to do this with an apply but it works. Could revise once we expand to larger extent
    if (nm == "intercept") {
      # get any raster with the same extent and resolution as ntland, fill with 0, add the coefficient, and return!
      return(ntland[[1]] * 0 + beta_this[nm])
    }
    
    if (str_detect(nm, "^I\\(log")) {
      # if this covariate was log-transformed, we can tell by the name, and we do it in the raster
      ntland_nm = str_split_i(nm, "(\\(| )", 3) # remove the "I(log(... + 1))" stuff from the name so we can access the correct raster layer
      return(beta_this[nm] * log1p(ntland[[ntland_nm]]))
    }
    
    # if we've gotten here, it presumably means we just return the base raster value because there are no other I() functions defined
    beta_this[nm] * ntland[[nm]]
    
  }))
  
  sum(ntland_manip)
  
}

W_issf_base = get_metric_raster(base_var, start_or_end = "end")
W_issf_sl = get_metric_raster(lsl_var)
W_issf_ta = get_metric_raster(cta_var, metric_regex = "I\\(cos\\(ta_\\)\\)(:?)")

# combine
W_issf_all = c(W_issf_base, W_issf_sl, W_issf_ta)

# Run the function to simulate paths ####
ncores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")) * as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))
if (!is.finite(ncores)) ncores = 1
message("Number of cores: ", ncores)
registerDoParallel(cores = ncores)

NP = 50000 # hopefully this is enough for a map!
np_per_core = ceiling(NP / ncores)

out = sim_paths_issf(n_lists = ncores,
                     n_paths_per_list = np_per_core,
                     n_steps_per_path = 675, 
                     move_pars = this_move_pars,
                     R_env = W_issf_all,
                     R_mask = !is.na(W_issf_base), # so the caribou don't leave the domain where covariates are defined
                     coef_interact = c("identity", "steplog", "angle"), # see function documentation for more info on this
                     coef_point_vals = c("end", "start", "end"), # so the code knows whether to look at the start or end point for each portion of the model
                     log_step_cos_angle_par = beta_all[lsl_cta_var],
                     log_sl_offset = 1, # because all the logs in this model are log(... + 1)
                     n_print = 10)

# FOR CLUSTER ONLY
out_dir = "/scratch/pt1/borealcaribou/outputs"

if (dir.exists(out_dir)) {
  
  DATE_OUT = str_sub(str_replace_all(Sys.time(), "-", ""), 1, 8)
  saveRDS(out, file.path(out_dir, paste0("paths_NWT_", DATE_OUT, ".rds")))
  out_bind = do.call(rbind, out)
  
  v = vect(out_bind, geom = c("x", "y"))
  rv = rasterize(v, W_issf_base, fun = sum, background = 0)
  
  writeRaster(rv, file.path(out_dir, paste0("TUD_NWT_", DATE_OUT, ".tif")), overwrite = TRUE)

}