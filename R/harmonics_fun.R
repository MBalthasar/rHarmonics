#' Harmonic modelling
#'
#' This function enables the user to model different number of cycles per year
#'
#' @param user_vals A vector with numeric values.
#'
#' @param user_dates A Vector with Date objects. See \code{\link[base]{as.Date}}.
#'
#' @param harmonic_deg Numeric. The number of cycles per year (harmonic degree)
#'                     that should be modelled.
#'
#' @param ref_date (optional) A Date object. Default is 1970-01-01.
#'
#' @return A numeric vector with the fitted values.
#'
#' @details To calculate the harmonic fitted curve of a periodic signal,
#'          ordinary least squares regressions are computed using coupled
#'          sine and cosine curves on time-series data. The underlying algorithm
#'          is based on Shumway & Stoffer (2017) equations 4.1 – 4.2.
#'
#' @references Shumway, R. H., & Stoffer, D. S. (2017). Time series analysis and its
#'             applications: with r examples. Springer.
#'
#' @examples
#'
#' library(ggplot2)
#'
#' # Load sample NDVI time-series data.frame
#' ndvi_df <- base::readRDS(system.file(package = "rHarmonics",
#'                                      "extdata", "MODIS_NDVI_TimeSeries.rda"))
#'
#' # Apply harmonic function using 3 cycles per year
#' fitted_3rd_deg <- harmonics_fun(user_vals = ndvi_df$ndvi,
#'                                 user_dates = ndvi_df$dates,
#'                                 harmonic_deg = 3)
#'
#' # Combine fitted values with original df
#' combined_df <- cbind(ndvi_df, fitted_3rd_deg)
#' names(combined_df) <- c("dates", "ndvi", "fitted")
#'
#' # Plot original data with fitted values
#' ggplot2::ggplot() +
#'   geom_point(data=combined_df, aes(x = dates, y = ndvi), color='black', size=1) +
#'   geom_line(data=combined_df, aes(x = dates, y = fitted), colour = "red", size=1)+
#'   labs(title="MODIS NDVI TimeSeries with fitted values", x="Date", y="NDVI") +
#'   theme(plot.title = element_text(hjust = 0.5, face="bold", size=14),
#'         axis.text=element_text(size=10),
#'         axis.title=element_text(size=12),
#'         legend.title=element_text(size=13, face="bold"),
#'         legend.text=element_text(size=13))+
#'   ylim(0,1)
#'
#'
#' # Apply harmonic function on multi-layer raster file
#'
#' library(raster)
#'
#' # Load sample cloud- and snow-free MODIS TimeSeries
#' sample_raster <- raster::brick(system.file(package = "rHarmonics",
#'                                            "extdata",
#'                                            "MODIS_NDVI_stack.tif"))
#'
#' fitted_raster <- raster::calc(sample_raster,
#'                               function(x){
#'                                  harmonics_fun(x, user_dates = ndvi_df$dates,
#'                                                harmonic_deg = 3)
#'                               })
#'
#' @export
harmonics_fun <- function(user_vals, user_dates, harmonic_deg, ref_date){

  ### For every missing output parameter set the default ###

  if (missing(user_vals)){
    stop("Values must be provided.")
  }
  if (missing(user_dates)){
    stop("Dates must be provided.")
  }
  # Check if dates are in "Date"-format
  if (class(user_dates) != "Date"){
    stop("Dates must be provided as 'Date' objects.")
  }
  if (missing(harmonic_deg)){
    stop("Harmonic degree must be provided.")
  }
  if (missing(ref_date)){
    ref_date <- as.Date("1970-01-01")
  } else if (class(ref_date) != "Date"){
    stop("Reference date must be provided as a 'Date' object.")
  }

  # If user vals are only NA, output is same as input

  if (all(is.na(user_vals))) {
    print("User values consist of NA values only. Output is same as input")
    return(user_vals)

    # Otherwise apply harmonic analysis
  } else {

    ### Calculate difference to ref_date in radians ###

    # Start for loop
    for (i in 1:length(user_dates)){
      current_date <- user_dates[i]
      # Calculate the difference in days
      current_diff_days <- as.numeric(difftime(current_date, ref_date), units="days")
      # Convert to years
      current_diff_years <- current_diff_days/365.25
      # Convert to radians
      current_diff_radians <- current_diff_years * 2 * pi
      if (i == 1){
        my_radians <- current_diff_radians
      } else {
        my_radians <- c(my_radians, current_diff_radians)
      }
      rm(i, current_date, current_diff_days, current_diff_years, current_diff_radians)
    }

    ### Caculate sines and cosines ###

    # Define constant
    my_cons <- rep(1, times = length(my_radians))
    # Create sines and cosines data frames with one column for each harmonic
    my_sin <- data.frame(matrix(nrow = length(my_radians), ncol = harmonic_deg))
    my_cos <- data.frame(matrix(nrow = length(my_radians), ncol = harmonic_deg))
    # Create names for df
    sin_names <- rep("sin_", times=harmonic_deg)
    cos_names <- rep("cos_", times=harmonic_deg)
    my_seq <- as.character(seq(1,harmonic_deg, by=1))
    sin_names <- paste(sin_names, my_seq, sep="")
    cos_names <- paste(cos_names, my_seq, sep="")
    # Add names to df
    names(my_sin) <- sin_names
    names(my_cos) <- cos_names
    # Calculate sines and cosines for each harmonic
    for (j in 1:harmonic_deg){
      # Calculate current sines and cosines by multiplying the radians with the current
      # degree of harmonic and then apply the sine/cosine function
      current_sines <- sin(my_radians * j)
      current_cosines <- cos(my_radians * j)
      # Fill data frames with values
      my_sin[,j] <- current_sines
      my_cos[,j] <- current_cosines
      # remove redundant variables
      rm(j, current_cosines, current_sines)
    }
    # Create df from the dependent and all independent variables
    df_for_reg <- cbind(user_vals, my_cons, my_radians, my_sin, my_cos)

    ### Apply Ordinary Least Squares Regression ###

    # Ordinary Least Squares Regression
    my_reg <- stats::lm(formula = user_vals ~ ., data = df_for_reg)
    # Get coefficients (constant is intercept value)
    # For coefficient and radians it's easy ...
    cons_coef <- my_reg$coefficients[1]
    t_coef <- my_reg$coefficients[3]
    # ... for the sines and cosines selecting the right columns is a little more complicated
    # Start with 4 because the first three values are the ndvi, my_cons and my_radians
    # First define the start and stop column for the sines and cosines ...
    sin_start <- 4
    sin_end <- 4 + harmonic_deg - 1
    cos_start <- 4 + harmonic_deg
    cos_end <- 4 + harmonic_deg + harmonic_deg -1
    # ... and then subset the data accordingly
    sin_coef <- my_reg$coefficients[c(sin_start:sin_end)]
    cos_coef <- my_reg$coefficients[c(cos_start:cos_end)]

    ### Calculate fitted values ###

    # multiply independent variables with the coefficients
    df_for_reg[,2] <- df_for_reg[,2] * cons_coef
    df_for_reg[,3] <- df_for_reg[,3] * t_coef
    # for loop multiplying the factor for each harmonic degree
    for (k in 1:harmonic_deg){
      # for the sines define i + 3 because the first sine column is at the 4th position
      df_for_reg[,k + 3] <- df_for_reg[,k + 3] * sin_coef[k]
      # for the cosines define i + 3 + harmonic_deg to get to the first cosine column
      df_for_reg[,k + 3 + harmonic_deg] <- df_for_reg[,k + 3 + harmonic_deg] * cos_coef[k]
      # remove reduntant variables
      rm(k)
    }
    # calculate sum (fitted value) of the multplied independent variables
    fitted <- rowSums(df_for_reg[,c(2:ncol(df_for_reg))], na.rm = TRUE)
    # return fitted values
    return(fitted)
  }
}
