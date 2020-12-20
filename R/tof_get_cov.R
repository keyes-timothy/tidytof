# Function to calculate, adjust, and invert the covariance matrix for the reference populations
tof_get_cov <- function(data) {
    cov(data) %>% 
    solve()
}
