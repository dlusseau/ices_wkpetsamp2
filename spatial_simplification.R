# Generate fake-spatial pbycatch or fishing effort by area
# Change p_bycatch or effort based on areas. I.e., if a species has a bycatch hotspot in a specific location, or if there is a hotspot where fishing effort is concentrated.
# For pbycatch: This value can be multiplied by a baseline probability for a given species, or can provide a vector of probabilities alone (but make sure they're realistic values!)
# For effort: This can be a multiplier for the mean effort by boat. You can further concentrate effort in one area by changing the number of boats, but that is not covered here.
area_multiplier <- function(nareas, hotspot_alpha = NA, alpha_vec = NA) {
  if (is.na(alpha_vec) & is.na(hotspot_alpha)) {
    stop("Specify an alpha value, either in the form of a single hotspot with hotspot_alpha, or a vector of values in alpha_vec")
  }
  if (!is.na(alpha_vec) & length(alpha_vec != nareas)) {
    stop("alpha_vec must be of length nareas")
  }

  # unless a vector of alpha values is specified, there is a single hotspot
  if (is.na(alpha_vec)) {
    alpha_vec <- c(hotspot_alpha, rep(1, times = nareas - 1))
  }
  prop_by_area <- extraDistr::rdirichlet(
    n = 1,
    alpha = alpha_vec
  ) # /rdirichlet
  return(prop_by_area)
}

# test
# area_multiplier(nareas = 6,hotspot_alpha = 5)
