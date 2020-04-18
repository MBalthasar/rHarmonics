# rHarmonics
R package for harmonic modelling of time-series data.

## Installation

To install the current version, use `devtools`.

```R
devtools::install_github("MBalthasar/rHarmonics")
```

## Available Functions

The following functions are currently available and tested on Windows 10.

* `harmonics_fun()` This function enables the user to model different number of cycles per year for a given time series data set.

## Example

In this example a fitted curve using three cylces per year based on a NDVI MODIS timeseries is computed:

```R
library(ggplot2)

# Load sample NDVI time-series data.frame
ndvi_df <- base::readRDS(system.file(package = "rHarmonics",
                                     "extdata", "MODIS_NDVI_TimeSeries.rda"))

# Apply harmonic function using 3 cycles per year
fitted_3rd_deg <- harmonics_fun(user_vals <- ndvi_df$ndvi,
                                user_dates = ndvi_df$dates,
                                harmonic_deg <- 3)

# Combine fitted values with original df
combined_df <- cbind(ndvi_df, fitted_3rd_deg)
names(combined_df) <- c("dates", "ndvi", "fitted")

# Plot original data with fitted values
ggplot2::ggplot() +
  geom_point(data=combined_df, aes(x = dates, y = ndvi), color='black', size=1) +
  geom_line(data=combined_df, aes(x = dates, y = fitted), colour = "red", size=1)+
  labs(title="MODIS NDVI TimeSeries with fitted values", x="Date", y="NDVI") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=14),
        axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        legend.title=element_text(size=13, face="bold"),
        legend.text=element_text(size=13))+
  ylim(0,1)
```

<img src="images/Rplot01.png" width=1000>

The function can also be applied on a multi-layer raster stack to create a cloud-interpolated psuedo times-series data set:

```R
# Apply harmonic function on multi-layer raster file

library(raster)

# Load sample cloud- and snow-free MODIS TimeSeries
sample_raster <- raster::brick(system.file(package = "rHarmonics",
                                           "extdata",
                                           "MODIS_NDVI_stack.tif"))

# Convert raster stack to matrix
sample_matrix <- as.matrix(sample_raster)

# Flip rows and columns
# -> Each column represents the NDVI time-series in one pixel
sample_matrix_T <- t(sample_matrix)

# Each column represents one cell
mat_fitted <- apply(sample_matrix_T, 2, harmonics_fun,
                    user_dates = ndvi_df$dates,
                    harmonic_deg <- 3)

# Define new raster
new_raster <- sample_raster

# Fill values on new raster layer based on data.frame
for (i in 1:ncol(mat_fitted)) {
  # Fill each cell of the new raster with the fitted values
  # stored in each column of the matrix
  new_raster[i] <- mat_fitted[,i]
}
```
