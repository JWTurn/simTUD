# Simulate random paths based on SSF / iSSF input
#
# n_lists: integer > 0; number of independent lists of paths to make. Only relevant if computing in parallel (this will break the work up into "n_lists" threads / cores / processes) otherwise default value of 1 should be fine
# n_paths_per_list: integer > 0; number of paths per list (these will all be included in the same data.frame with a path_id column)
# n_steps_per_path: integer > 0; how many steps are included in each path? Note that the time units here are simply "steps", and you specify how long these steps are based on the inputs to the "move_pars" argument.
# n_rand: integer > 0; at each step, the algorithm chooses the animal's next step from a set of random points; this argument chooses how many to pick from. Higher values will have more computational cost but will be more accurate. Default is probably fine but can adjust as needed.
# move_pars: named numeric vector; contains parameters necessary for simulating random step lengths and turning angles, names include "gamma_shape", "gamma_rate", "vm_mu", "vm_kappa". If you're using an exponential distribution for step lengths, set gamma_shape to 1 and gamma_rate to 1 / mean step length. If local_gibbs is TRUE (it shouldn't be, if you're using iSSF outputs), only looks at the first parameter which is the radius size
# R_env: a SpatRasterDataset (from the terra package) containing all environmental covariates. SpatRasterDatasets basically allow you to generate a 2-D array of rasters and this is useful here because the rows can represent time / seasons (if habitat selection does not vary across seasons, use only one row) and the columns can represent a covariate / interaction type (see "coef_interact"). Everything should be on log scale (e.g., linear combination of iSSF parameter outputs) because we will take the exp() inside this function. Typically, it is sufficient to have one row (i.e., no seasonal variation) and two columns (one representing habitat quality from non-step-length terms, and one representing habitat quality from step length terms). It's most efficient to perform these raster operations before calling this function so you only have to do it once.
# R_mask: SpatRaster object; should be 0 or 1 (or NA) everywhere; determines whether certain areas are "in bounds" for steps. This allows for boundaries that are not simply a convex polygon around the study area (e.g., barriers that can't be crossed).
# R_x0: SpatRaster object; determines sampling probabilities for initial location. Should be in log form. If not supplied, will use first component of R_env.
# coef_values: 3-D array of coefficient values for each simulation (in case we want path outputs to vary by path). First dimension is path # (length equal to n_paths_per_list; this allows for individual variation within the simulations), second dimension is equal to season (length equal to n_seasons / # of rows of R_env), third dimension is equal to # of columns of R_env / # of different covariates. If you've pre-calculated the habitat quality values for inclusion in R_env, these can stay as the default value (1). For covariates interacted with step length, turning angle, etc, see "coef_interact".
# coef_interact: character vector with length equal to the # of columns of "coef_values", also the # of layers in each element of "R_env". Indicates whether coefficient values interact with nothing ("identity"; for traditional selection coefficients), "step", "steplog", "angle", or "daylight". This tells the algorithm which components of R_env need to be multiplied by step length, turn angle, etc.
# coef_point_vals: character vector with length equal to the # of columns of "coef_values", also the # of layers in each element of "R_env". Each value is either "end" (the default; for covariates we want to extract at the end point), or "start" (covariates we want to extract from the beginning of the step).
# step_cos_angle_par: numeric matrix with (# of paths per list) rows and (# of unique seasons) columns; interaction parameter between step length ("rate" parameter of gamma distribution) and cosine of turning angle ("kappa" parameter of von Mises distribution). If only one value is supplied the function will apply it to all paths / seasons as necessary.
# log_step_cos_angle_par: numeric matrix with (# of paths per list) rows and (# of unique seasons) columns; interaction parameter between log(step length) ("shape" parameter of gamma distribution) and cosine of turning angle ("kappa" parameter of von Mises distribution). If only one value is supplied the function will apply it to all paths / seasons as necessary.
# step_day_par: numeric matrix with (# of paths per list) rows and (# of unique seasons) columns; interaction parameter between step length ("rate" parameter of gamma distribution) and daylight_var. If only one value is supplied the function will apply it to all paths / seasons as necessary.
# log_step_day_par: numeric matrix with (# of paths per list) rows and (# of unique seasons) columns; interaction parameter between log(step length) ("shape" parameter of gamma distribution) and daylight_var. If only one value is supplied the function will apply it to all paths / seasons as necessary.
# cos_angle_day_par: numeric matrix with (# of paths per list) rows and (# of unique seasons) columns; interaction parameter between cos(turn angle) ("kappa" parameter of von Mises distribution) and daylight_var. If only one value is supplied the function will apply it to all paths / seasons as necessary.
# memory_par: numeric; parameter for selection of previously visited areas; can also be a vector if parameter value is to be different for each path. Usuaully 0 is sufficient but sometimes can help make simulations more realistic. If only one value is supplied the function will apply it to all paths / seasons as necessary.
# mem_fun: character; what method are we using to calculate memory? See the code itself (specifically, line 283 and below) where memory is calculated for more details
# scale_mem_value: numeric > 0; what do we multiply the memory values by? Depends on how they were scaled to fit the original model that estimated memory_par. If mem_fun == "d_initial_norm_int", then this is a vector of length 2 with scaling constants for normalized and un-normalized distances
# local_gibbs: logical; do we simulate using the SSF or using the local Gibbs algorithm? Usually OK to leave this as false.
# x0: either NULL (default; will sample random initial locations) or a list of data.frames / numeric matrices with two columns depicting initial locations
# x0_bounds: SpatVector polygon object; space from which to sample initial locations. If R_env or R_mask are provided and this is NULL, will use extent of those instead.
# daylight_var: numeric vector with length equal to n_steps_per_path; contains information on some quantity (e.g., time of day) at each time step. If there are no "time of day" effects in the model, this can stay as the default (0) value.
# season_var: integer vector; all values should be between 1 and the # of layers in the covariate rasters ("R_..."). Tells us which index to grab from at each point in time. If there is no seasonality this will always be 1 (the default)
# n_check_mask: integer > 0; how many points along a step do we check for being out of bounds? This can be useful in, for example, areas with complex topography and long timesteps where animals may cross a steep cliff or gulley to go from one good habitat patch to another. If the midpoint, or any of a series of points along the path between these two locations, is in an area deemed by R_mask to be out of bounds, we will replace this step with a newly re-drawn one.
# step_factor: numeric > 0; how much do we multiply randomly generated steps by to get them to their true values? Can be relevant if move_pars are on a different scale from the units of the rasters (e.g., metres versus kilometres)
# min_sl: numeric >= 0; minimum step length to be simulated at each step in m. FOR NOW, only relevant if local_gibbs is FALSE
# log_sl_offset: numeric; if log step lengths were shifted during model fitting, how should we adjust them accordingly here?
# n_cores: integer > 0; how many cores do we divide the work into? Should be at most n_lists. If it's less than or equal to 1 we don't use parallelization
# n_print: integer > 0; how often do we print an update on the progress of the function? If it's greater than n_steps_per_path then we never print anything
# ret_list: logical; only really relevant when n_lists == 1, and in that case, tells the function whether to return a list (with 1 element) or a data.frame
# ...: additional arguments to memory function (whatever it is). Currently not implemented
#
# Returns a data.frame (or list of data.frames) for paths
sim_paths_issf = function(n_lists = 1, 
                          n_paths_per_list = 1, 
                          n_steps_per_path = 100, 
                          n_rand = 10,
                          move_pars = c(gamma_shape = 1, gamma_rate = 1, vm_mu = 0, vm_kappa = 0),
                          R_env = NULL,
                          R_mask = NULL,
                          R_x0 = NULL,
                          coef_values = array(1, c(n_paths_per_list, 1, 1)),
                          coef_interact = rep("identity", dim(coef_values)[2]),
                          coef_point_vals = rep("end", dim(coef_values)[2]),
                          step_cos_angle_par = 0,
                          log_step_cos_angle_par = 0,
                          step_day_par = 0,
                          log_step_day_par = 0,
                          cos_angle_day_par = 0,
                          memory_par = 0,
                          mem_fun = c("KDE", "d_initial", "d_initial_int", "d_initial_rel", "d_initial_norm", "d_initial_norm_int", "cos_dir_hrc"),
                          scale_mem_value = 1,
                          local_gibbs = FALSE,
                          x0 = NULL,
                          x0_bounds = NULL, # TO DO: Could make an option for this to be a list?
                          daylight_var = numeric(n_steps_per_path),
                          season_var = rep(1, n_steps_per_path),
                          n_check_mask = 1,
                          step_factor = 1,
                          min_sl = 0,
                          log_sl_offset = 0,
                          n_cores = n_lists, 
                          n_print = n_steps_per_path + 1,
                          ret_list = (n_lists > 1),
                          ...) {
  
  # TO DO: Change the for loop so it doesn't rbind every time
  
  require(doParallel)
  require(terra)
  require(circular)
  require(tidyverse)
  require(data.table)
  require(abind)
  
  mem_fun = match.arg(mem_fun)
  if (mem_fun == "d_initial_norm_int" && length(scale_mem_value) == 1) stop("Incorrect format for scale_mem_value")
  if (n_cores > 1) registerDoParallel(cores = min(n_cores, n_lists))
  if (local_gibbs && min_sl > 0) {
    warning("local_gibbs is set to TRUE and min_sl is set > 0. These are incompatible for now and thus ignoring min_sl value.")
    min_sl = 0
  }
  
  if (is.null(R_env) & is.null(R_mask)) {
    # Make a blank raster with the default resolution
    R_env = sds(rast(vals = 0, crs = "+init=epsg:4326"))
    # vals = 1 here because we want the animal to be able to go everywhere
    R_mask = rast(vals = 1, crs = "+init=epsg:4326")
  } else if (is.null(R_env)) {
    # Make a new raster using the template from R_mask
    R_env = sds(rast(vals = 0, extent = ext(R_mask), resolution = res(R_mask), crs = crs(R_mask)))
  } else if (is.null(R_mask)) {
    # Make a mask that's 1 everywhere (i.e., everywhere's good)
    R_mask = rast(vals = 1, extent = ext(R_env), resolution = res(R_env), crs = crs(R_env))
  }
  if (!is(R_env, "SpatRasterDataset")) R_env = sds(R_env)
  
  n_seasons = max(season_var)
  # Repeat R_env if necessary
  if (n_seasons > 1 && length(R_env) == 1) R_env = do.call(sds, rep(list(R_env), n_seasons))
  
  # Get total # of layers in each season. They should be all the same!
  n_lyr_r_env = nlyr(R_env)
  if (any(n_lyr_r_env != n_lyr_r_env[1])) stop("All elements of 'R_env' must have the same number of layers. Please supply a value for this argument that meets this criterion.")
  n_coefs = n_lyr_r_env[1]
  
  # Fix coef_values if necessary
  if (dim(coef_values)[2] != n_seasons | dim(coef_values)[3] != n_coefs) {
    correct_cv_dim = c(n_paths_per_list, n_seasons, n_lyr_r_env)
    warning("coef_values has the wrong dimensions. Replacing with correct dimensions but removing all existing values and replacing with 1. Please re-run with proper dimensions ([", paste0(correct_cv_dim, collapse = ", "), "]) if you want to apply user-selected values.")
    coef_values = array(1, correct_cv_dim)
  }
  
  # Repeat all of these things as necessary
  coef_interact = rep_ifnecessary(coef_interact, n_coefs)
  coef_point_vals = rep_ifnecessary(coef_point_vals, n_coefs)
  step_cos_angle_par = rep_ifnec_2d(step_cos_angle_par, n_seasons, n_paths_per_list, n_rand)
  log_step_cos_angle_par = rep_ifnec_2d(log_step_cos_angle_par, n_seasons, n_paths_per_list, n_rand)
  step_day_par = rep_ifnec_2d(step_day_par, n_seasons, n_paths_per_list, n_rand)
  log_step_day_par = rep_ifnec_2d(log_step_day_par, n_seasons, n_paths_per_list, n_rand)
  cos_angle_day_par = rep_ifnec_2d(cos_angle_day_par, n_seasons, n_paths_per_list, n_rand)
  memory_par = rep_ifnec_2d(memory_par, n_seasons, n_paths_per_list, n_rand)
  
  start_vals = coef_point_vals == "start"
  end_vals = coef_point_vals == "end"
  
  if (is.null(x0_bounds) & is.null(x0)) x0_bounds = vect(ext(R_mask), crs = crs(R_mask))
  
  i_dist_prop = seq(1/n_check_mask, 1, length.out = n_check_mask) # get proportions to check mask
  
  if (n_cores > 1) {
    `%dothis%` = `%dopar%`
  } else {
    `%dothis%` = `%do%`
  }
  
  nr = n_rand * n_paths_per_list
  
  list_of_all_paths = foreach(yyy = 1:n_lists, .errorhandling = "stop") %dothis% {
    
    if (is.null(x0)) {
      
      # Random start locations within boundary
      xy = spatSample(x0_bounds, size = nr)
      xy = crds(xy)
      
      # Extract values from the raw RSF raster to pick the random initial locations in a way that represents habitat quality
      if (is.null(R_x0)) {
        this_coef_vals = apply(adrop(coef_values[, 1, , drop = FALSE], drop = 2), 2, rep, each = n_rand)
        this_coef_matrix = coef_interact_matrix(n_rows = nr,
                                                n_cols = n_coefs,
                                                coef_interact = coef_interact,
                                                st_vals = numeric(nr) + ("steplog" %in% coef_interact))
        t_ssf = apply(terra::extract(R_env[season_var[1]], xy) * this_coef_vals * this_coef_matrix, 1, sum)
      } else {
        t_ssf = terra::extract(R_x0, xy)[, 1]
      }
      t_zero = terra::extract(R_mask, xy)[ ,1]
      
      ud <- exp(t_ssf) * t_zero
      # convert NA values to 0
      ud[!is.finite(ud)] = 0
      keep = sample(1:nrow(xy), size = n_paths_per_list, replace = FALSE, prob = ud) 
      xy_st = xy[keep, , drop = FALSE]
      
    } else {
      xy_st = x0[[yyy]]
    }
    
    # template data.frame
    tmp0 = data.table(path_id = 1:n_paths_per_list, 
                      i_strata = 0,
                      x_first = as.numeric(xy_st[, 1, drop = TRUE]), 
                      y_first = as.numeric(xy_st[, 2, drop = TRUE]),
                      x_prev = as.numeric(xy_st[, 1, drop = TRUE]),
                      y_prev = as.numeric(xy_st[, 2, drop = TRUE]), 
                      x = NA, 
                      y = NA, 
                      t_var = daylight_var[1],
                      direction_prev = round(runif(n_paths_per_list, -pi, pi), 3))
    
    path_all = data.table()
    
    for (i in 1:n_steps_per_path) {
      
      tmp1 = tmp0[, c("i_strata", "t_var") := .(i, daylight_var[i])]
      
      path_nums = rep(tmp1$path_id, each = n_rand)
      
      if (local_gibbs) {
        # Generate random intermediate location from disk around current location. Note that for each path the same intermediate location is used. This is my best interpretation of Michelot et al. (2019) who describe selecting one intermediate point per timestep for simulation, and then selecting the final location using a sampling algorithm (in our case, this means picking n_rand points).
        r_circle_r = runif(n_paths_per_list, 0, move_pars[1] * step_factor)
        r_circle_theta = runif(n_paths_per_list, -pi, pi)
        
        r_circle_x = r_circle_r * cos(r_circle_theta)
        r_circle_y = r_circle_r * sin(r_circle_theta)
        
        # Sample randomly from the truncated target distribution i.e., by approximating with a bunch of random locations
        r_location_r = runif(nr, 0, move_pars[1] * step_factor)
        r_location_theta = runif(nr, -pi, pi)
        
        r_location_x = rep(r_circle_x + tmp1$x_prev, each = n_rand) + r_location_r * cos(r_location_theta)
        r_location_y = rep(r_circle_y + tmp1$y_prev, each = n_rand) + r_location_r * sin(r_location_theta)
        r_location_he = atan2(r_location_y - rep(tmp1$y_prev, each = n_rand), r_location_x - rep(tmp1$x_prev, each = n_rand))
        
        r_steplengths = sqrt(rep(r_circle_r, each = n_rand)^2 + r_location_r^2)
        r_turn_angles = r_location_he - rep(tmp1$direction_prev, each = n_rand)
        move_metrics_df = data.table(step = r_steplengths, angle = r_turn_angles, path_id = path_nums)
      } else {
        # Here is the more traditional method for generating random steps based on step length & turning angle
        r_steplengths = rgamma(length(path_nums), move_pars["gamma_shape"], move_pars["gamma_rate"])
        r_turn_angles = as.numeric(rvonmises(length(path_nums), circular(move_pars["vm_mu"]), move_pars["vm_kappa"]))
        move_metrics_df = data.table(step = r_steplengths * step_factor, angle = r_turn_angles, path_id = path_nums)
      }
      
      tmp2a = tmp1[move_metrics_df, on = "path_id", nomatch = 0, allow.cartesian = TRUE][
        ,direction := (direction_prev + angle) %% (2*pi)][
          ,c("x", "y") := .(x_prev + step * cos(direction), y_prev + step * sin(direction))]
      
      # Filter steps that cross inhospitible habitat - check midway points for fast steps as they are longer
      zero_mx = sapply(i_dist_prop, function(prop) {
        tt_x <- prop * tmp2a$step * cos(tmp2a$direction) + tmp2a$x_prev
        tt_y <- prop * tmp2a$step * sin(tmp2a$direction) + tmp2a$y_prev
        terra::extract(R_mask, cbind(tt_x, tt_y))[,1]
      })
      zero_mx[!is.finite(zero_mx)] = 0
      
      tmp2b = tmp2a[, step_good := ifelse(rowSums(zero_mx) == length(i_dist_prop) & r_steplengths >= min_sl, 1, 0)][
        , step_good_all := sum(step_good), by = path_id] # do all path_id's have at least one good step?
      # Identify whether we need to "re-roll" any paths because none of the steps ended up inside the mask
      # TO DO: This code currently seems as if "all" should be "any" but a) that would take longer and b) seems too much. It still works as intended as it is now but could maybe clean this up a bit in the future.
      while(any(tmp2b$step_good_all == 0)) {
        tmp2b_bad = tmp2b[step_good == 0]
        tmp2b_good = tmp2b[step_good == 1] # we keep these
        
        n_bad_locs = nrow(tmp2b_bad)
        
        if (local_gibbs) {
          # Generate random intermediate location from disk around current location
          r_circle_r = runif(n_paths_per_list, 0, move_pars[1] * step_factor)
          r_circle_theta = runif(n_paths_per_list, -pi, pi)
          
          r_circle_x = r_circle_r * cos(r_circle_theta)
          r_circle_y = r_circle_r * sin(r_circle_theta)
          
          # Sample randomly from the truncated target distribution i.e., by approximating with a bunch of random locations
          r_location_r = runif(n_bad_locs, 0, move_pars[1] * step_factor)
          r_location_theta = runif(n_bad_locs, -pi, pi)
          
          bad_locs_table = table(factor(paste0("p", tmp2b_bad$path_id), levels = paste0("p", 1:n_paths_per_list)))
          
          r_location_x = rep(r_circle_x + tmp1$x_prev, bad_locs_table) + r_location_r * cos(r_location_theta)
          r_location_y = rep(r_circle_y + tmp1$y_prev, bad_locs_table) + r_location_r * sin(r_location_theta)
          r_location_he = atan2(r_location_y - rep(tmp1$y_prev, bad_locs_table), r_location_x - rep(tmp1$x_prev, bad_locs_table))
          
          r_steplengths_bad = sqrt(rep(r_circle_r, bad_locs_table)^2 + r_location_r^2)
          r_turn_angles_bad = r_location_he - tmp2b_bad$direction_prev
          move_metrics_df_bad = data.table(step = r_steplengths_bad, angle = r_turn_angles_bad, path_id = tmp2b_bad$path_id)
        } else {
          r_steplengths_bad = rgamma(n_bad_locs, move_pars["gamma_shape"], move_pars["gamma_rate"])
          r_turn_angles_bad = as.numeric(rvonmises(n_bad_locs, circular(move_pars["vm_mu"]), move_pars["vm_kappa"]))
          move_metrics_df_bad = data.table(step = r_steplengths_bad * step_factor, angle = r_turn_angles_bad, path_id = tmp2b_bad$path_id)
        }
        
        tmp2a_bad = tmp1[move_metrics_df_bad, on = "path_id", nomatch = 0, allow.cartesian = TRUE][
          ,direction := (direction_prev + angle) %% (2*pi)][
            ,c("x", "y") := .(x_prev + step * cos(direction), y_prev + step * sin(direction))]
        
        zero_mx_bad = sapply(i_dist_prop, function(prop) {
          tt_x <- prop * tmp2a_bad$step * cos(tmp2a_bad$direction) + tmp2a_bad$x_prev
          tt_y <- prop * tmp2a_bad$step * sin(tmp2a_bad$direction) + tmp2a_bad$y_prev
          terra::extract(R_mask, cbind(tt_x, tt_y))[,1]
        })
        zero_mx_bad[!is.finite(zero_mx_bad)] = 0
        # Sometimes if there is only one "bad" value the sapply will return a vector and not a matrix. I couldn't find a solution for this; maybe there is one but for now this check will suffice. So far I am not including the check outside of the while loop as there will basically never be only one location (would require n_paths == 1 and n_rand == 1)
        if (!is.matrix(zero_mx_bad)) zero_mx_bad = t(zero_mx_bad)
        
        tmp2b_bad_new = tmp2a_bad[, step_good := ifelse(rowSums(zero_mx_bad) == length(i_dist_prop) & r_steplengths_bad >= min_sl, 1, 0)][
          , step_good_all := sum(step_good), by = path_id]
        tmp2b = rbind(tmp2b_good, tmp2b_bad_new)[order(path_id)]
      }
      
      # tmp2b = tmp2b[, step := ifelse(step == 0, 1e-5, step)]
      
      # Get covariate values at all proposed steps
      xy_vals_extract_end = cbind(tmp2b$x, tmp2b$y)
      xy_vals_extract_sta = cbind(tmp2b$x_prev, tmp2b$y_prev)
      
      # extract all values, then multiply by necessary covariates, then sum together
      coef_interact_matrix_i = coef_interact_matrix(nr, 
                                                    coef_interact = coef_interact,
                                                    st_vals = tmp2b$step,
                                                    ta_vals = cos(tmp2b$angle),
                                                    dl_vals = rep(daylight_var[i], nr),
                                                    lsl_offset = log_sl_offset)
      this_coef_vals = apply(adrop(coef_values[, season_var[i], , drop = FALSE], drop = 2), 2, rep, each = n_rand)
      this_coef_final = coef_interact_matrix_i * this_coef_vals
      this_r_env = R_env[season_var[i]]
      if (any(start_vals)) {
        log_p_start = apply((terra::extract(this_r_env, xy_vals_extract_sta) * this_coef_final)[, start_vals, drop = FALSE], 1, sum)
      } else {
        log_p_start = numeric(nr)
      }
      if (any(end_vals)) {
        log_p_end = apply((terra::extract(this_r_env, xy_vals_extract_end) * this_coef_final)[, end_vals, drop = FALSE], 1, sum)
      } else {
        log_p_end = numeric(nr)
      }
      log_p = log_p_start + log_p_end

      # TO DO: See if it's worth rearranging this so these extra parameters come as a list or something.
      tmp2b = tmp2b[, log_step_adjusted := log(step) + log_sl_offset][
        , logp := log_p + daylight_var[i] * (step_day_par[, season_var[i]] * step + log_step_day_par[, season_var[i]] * log_step_adjusted + cos_angle_day_par[, season_var[i]] * cos(angle)) + log_step_adjusted * cos(angle) * log_step_cos_angle_par[, season_var[i]] + step * cos(angle) * step_cos_angle_par[, season_var[i]]]
      
      # add in memory. The i > 1 condition is necessary because at that point in the loop, path_all doesn't have anything in it yet so it will cause an error
      if (i > 1 & any(memory_par[[season_var[i]]] != 0)) {
        if (mem_fun == "KDE") {
          # Include as a model covariate a kernel density estimate (KDE) of previous locations. 
          tmp2b$memory = unlist(sapply(unique(tmp2b$path_id), function(id) {
            this_path = tmp2b[path_id == id]
            prev_pts = path_all[path_id == id] # gets all previous points
            prev_pts_trk = amt::make_track(prev_pts, .x = x, .y = y, crs = crs(R_rsf))
            get_KDE_value(this_path, prev_pts_trk)
          })) * scale_mem_value * memory_par[[season_var[i]]]
        } else if (mem_fun == "d_initial") {
          # distance from initial location
          tmp2b$memory = unlist(sapply(unique(tmp2b$path_id), function(id) {
            this_path = tmp2b[path_id == id]
            first_loc = path_all[path_id == id & i_strata == 1]
            this_path_sp = vect(this_path, geom = c("x", "y"), crs = crs(R_rsf))
            first_loc_sp = vect(first_loc, geom = c("x_prev", "y_prev"), crs = crs(R_rsf))
            as.numeric(terra::distance(this_path_sp, first_loc_sp))
          })) * scale_mem_value * memory_par[[season_var[i]]]
        } else if (mem_fun == "d_intial_q") {
          # Distance from the initial location but with a linear and quadratic term
          tmp2b = tmp2b[, d_first := sqrt((x - x_first) ^ 2 + (y - y_first) ^ 2)][
            , memory := d_first * memory_par[[season_var[i]]][1] + d_first ^ 2 * memory_par[[season_var[i]]][2]]
        } else if (mem_fun == "d_initial_int") {
          # distance from initial location, interacted with the previous distance from initial location (i.e., if you're not far away from your initial location the effect of "distance from initial location" may not mean much)
          tmp2b = tmp2b[, d_hrc_raw := sqrt((x - x_first)^2 + (y - y_first)^2) * scale_mem_value][
            , d_hrc_prev_raw := sqrt((x_prev - x_first)^2 + (y_prev - y_first)^2) * scale_mem_value][
              , memory := d_hrc_raw * memory_par[[season_var[i]]][1] + 
                d_hrc_raw * d_hrc_prev_raw * memory_par[[season_var[i]]][2]]
        } else if (mem_fun == "d_initial_rel") {
          # distance between the range centre distances for two most recent locations; i.e., is the animal moving further away from the centre?
          tmp2b$memory = unlist(sapply(unique(tmp2b$path_id), function(id) {
            this_path = tmp2b[path_id == id]
            first_loc = path_all[path_id == id & i_strata == 1]
            this_path_sp = vect(this_path, geom = c("x", "y"), crs = crs(R_rsf))
            this_path_prev_sp = vect(this_path, geom = c("x_prev", "y_prev"), crs = crs(R_rsf))
            first_loc_sp = vect(first_loc, geom = c("x_prev", "y_prev"), crs = crs(R_rsf))
            as.numeric(terra::distance(this_path_sp, first_loc_sp) - terra::distance(this_path_prev_sp, first_loc_sp))
          })) * scale_mem_value * memory_par[[season_var[i]]]
        } else if (mem_fun == "d_initial_norm") {
          # squared distance between home range centre and current location
          tmp2b = tmp2b[, memory := exp(-((sqrt((x - x_first)^2 + (y - y_first)^2)) * scale_mem_value)^2) * memory_par[[season_var[i]]]]
        } else if (mem_fun == "d_initial_norm_int") {
          tmp2b = tmp2b[, d_hrc_raw := sqrt((x - x_first)^2 + (y - y_first)^2) * scale_mem_value[2]][
            , d_hrc_prev_raw := sqrt((x_prev - x_first)^2 + (y_prev - y_first)^2) * scale_mem_value[2]][
              , d_hrc_norm_raw := exp(-(d_hrc_raw * scale_mem_value[1])^2)][
                , memory := d_hrc_norm_raw * memory_par[[season_var[i]]][1] + 
                  d_hrc_norm_raw * d_hrc_prev_raw * memory_par[[season_var[i]]][2] + 
                  d_hrc_norm_raw * d_hrc_prev_raw^2 * memory_par[[season_var[i]]][3]]
        } else if (mem_fun == "cos_dir_hrc") {
          # instead measures direction from home range centre (i.e. is it moving away or towards it)
          tmp2b$memory = as.numeric(sapply(unique(tmp2b$path_id), function(id) {
            this_path = tmp2b[path_id == id]
            first_loc = path_all[path_id == id & i_strata == 1]
            cos(atan2(first_loc$y - this_path$y_prev, first_loc$x - this_path$x_prev) - this_path$direction)
          })) * scale_mem_value * memory_par[[season_var[i]]]
        } else {
          # nothing else has been implemented yet
          tmp2b$memory = 0
        }
      } else {
        tmp2b$memory = 0
      }
      
      # Add in step effects, correct NA values for memory if necessary
      tmp2c = tmp2b[, step_good := ifelse(is.finite(memory), step_good, 0)][
        , exp_lp := exp(logp + memory) * step_good] 
      # do we need the " * step_good "? Shouldn't it be 1 always?
      tmp2c$exp_lp[!is.finite(tmp2c$exp_lp)] = 0
      tmp2c = tmp2c[, sum_exp_lp := sum(exp_lp), by = path_id][
        ,exp_lp := ifelse(sum_exp_lp == 0, step_good, exp_lp)]
      # In the event that the probablities all round off and become 0 (seems quite rare) we just pick from them equally
      
      tmp3b = tmp2c[, .SD[sample(.N, 1, prob = exp_lp)], by = path_id]
      
      tmp_save = tmp3b[, list(path_id, i_strata, x, y, x_prev, y_prev, t_var, step, direction)]
      path_all = rbind(path_all, tmp_save)
      
      # for use in the next iteration
      tmp0_names = names(tmp0)
      tmp0 = tmp3b[, c("x_prev", "y_prev", "x", "y", "direction_prev") := .(x, y, NA, NA, direction)][
        , ..tmp0_names]
      
      if (i %% n_print == 0) message("Completed step ", i, " of ", n_steps_per_path)
      
    }  # End of i (individual steps)
    
    path_all[, i_rep := yyy]
    
  }
  
  if (n_lists == 1 & !ret_list) return(list_of_all_paths[[1]]) # just give a data.frame
  list_of_all_paths
  
}

# Helper function for sim_paths
#
# n_row: number of rows for desired matrix
# n_col: number of columns for desired matrix. Only matters if coef_interact is not specified so can often be left equal to 1
# coef_interact: character vector with length equal to the # of columns of "coef_values", also the # of layers in each element of "R_env". Indicates whether coefficient values interact with nothing ("identity"; for traditional selection coefficients), "step", "steplog", "angle", or "daylight".
# st_vals: values of step lengths
# ta_vals: values of turning angles (or cosines of turning angles; whatever's fed in)
# dl_vals: values of daylight (or cosines of daylight; whatever's fed in)
# lsl_offset: offset value for log step length
coef_interact_matrix = function(n_rows,
                                n_cols = 1,
                                coef_interact = rep("identity", n_cols),
                                st_vals = numeric(n_rows) + 1e-4, # so we don't have log(0)
                                ta_vals = numeric(n_rows),
                                dl_vals = numeric(n_rows),
                                lsl_offset = 0) {
  
  if (n_cols != length(coef_interact)) n_cols = length(coef_interact)
  
  raw_vals = Vectorize(function(vv, stv, slv, tav, dlv) switch(vv, identity = 1, step = stv, steplog = slv, angle = tav, daylight = dlv),
                       vectorize.args = c("vv", "stv", "slv", "tav", "dlv"))(rep(coef_interact, each = n_rows), 
                                                                             rep(st_vals, n_cols), 
                                                                             rep(log(st_vals), n_cols) + lsl_offset, 
                                                                             rep(ta_vals, n_cols), 
                                                                             rep(dl_vals, n_cols))
  
  matrix(raw_vals, n_rows, n_cols)
  
}

# Helper function for sim_paths that repeats v by length.out but only if a) length.out is greater than 1 and b) v is length 1. Throws an error if v does not have length 1 or equal to length.out.
#
# v: vector
# length.out: integer > 0
#
# Returns a vector of length 'length.out'
rep_ifnecessary = function(v, length.out) {
  if (!(length(v) %in% c(1, length.out))) stop("Supplied value 'v' should either have length 1 or length equal to 'length.out'.")
  if (length.out > 1 && length(v) == 1) return(rep(v, length.out))
  v
}

# Helper function for sim_paths that repeats v in a 2-D matrix
#
# v: numeric value, vector, or matrix
# ncol_out: desired number of columns of matrix
# nrow_out: desired number of rows of matrix before being repeated
# desired number of repetitions for each row
#
# Returns a matrix with (nrow_out * nrep_out) rows by ncol_out columns
rep_ifnec_2d = function(v, ncol_out, nrow_out, nrep_out) {
  
  if (!is.matrix(v)) v = as.matrix(v)
  
  # repeat across columns
  if (!(ncol(v) %in% c(1, ncol_out))) stop("v must either be a vector / 1-column matrix or a matrix with ncol_out columns")
  if (ncol(v) != ncol_out) v_col = do.call(cbind, rep(list(v), ncol_out)) else v_col = v
  
  # repeat across rows
  if (!(nrow(v_col) %in% c(1, nrow_out))) stop("v_col must either be a vector of length (matrix with # of rows) 1 or nrow_out")
  if (nrow(v_col) != nrow_out) v_col_row = apply(v_col, 2, rep, nrow_out) else v_col_row = v_col
  
  # repeat finally
  apply(v_col_row, 2, rep, each = nrep_out)
  
}