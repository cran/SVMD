
#' @title Spearman Variational Mode Decomposition
#' @description  Optimal number of modes of Variational Mode Decomposition (VMD) using Spearman's rank correlation coefficient
#' @param series The input time series signal to be decomposed.
#' @param alpha The balancing parameter of the data-fidelity constraint. Default is 2000.
#' @param tau Time-step of the dual ascent (pick 0 for noise-slack). Default is 0.
#' @param DC If TRUE, the first mode is put and kept at DC (0 frequency). Default is FALSE.
#' @param init Mode initialization (1 = all omegas start at 0). Default is 1.
#' @param tol Convergence tolerance criterion. Default is 1e-7.
#' @param threshold The correlation coefficient threshold to determine the optimal number of modes. Default is 0.997.
#' @param max_modes The maximum number of modes to consider. Default is 10.
#' @param verbose Logical, if TRUE, prints detailed messages about the decomposition process.
#'
#' @return Returns a list containing the optimal number of modes, reconstructed signal, and additional outputs from the VMD process:
#' \itemize{
#'   \item \code{optimal_K}: The optimal number of modes.
#'   \item \code{reconstructed_signal}: The reconstructed signal from the selected modes.
#'   \item \code{imfs}: Intrinsic Mode Functions (IMFs) obtained from SVMD.
#'   \item \code{u_hat}: Estimated envelopes of the modes.
#'   \item \code{omega}: Frequencies of the modes.
#' }
#'
#' @examples{
#' # Example data generation:
#' # Set the number of observations
#' N <- 300
#' # Set a random seed for reproducibility
#' set.seed(123)
#' # Generate random uniform values
#' rand_unif <- runif(n = N, min = 0, max = 1.0)
#' # Create the components of the time series
#' sig1 <- 6 * rand_unif
#' sig2 <- sin(8 * pi * rand_unif)  # Using sine function
#' sig3 <- 0.5 * sin(40 * pi * rand_unif)  # Using sine function
#' # Combine the components to form the final signal
#' signal <- sig1 + sig2 + sig3
#' # Apply the sVMD function to the signal
#' result <- sVMD(signal)
#'}
#' @export
#'
#' @references
#' Yang, H., Cheng, Y., and Li, G. (2021). A denoising method for ship radiated noise based on Spearman variational mode decomposition, spatial-dependence recurrence sample entropy, improved wavelet threshold denoising, and Savitzky-Golay filter. Alexandria Engineering Journal, 60(3), 3379-3400
sVMD <- function(series, alpha = 2000, tau = 0, DC = FALSE, init = 1, tol = 1e-7, threshold = 0.997, max_modes = 10, verbose = FALSE) {
  # Initialize K
  K <- 1

  # Function to reconstruct series from IMFs
  reconstruct_series <- function(imfs) {
    return(rowSums(imfs))
  }

  # Store IMFs and additional outputs for the optimal number of modes
  optimal_imfs <- NULL
  optimal_u_hat <- NULL
  optimal_omega <- NULL

  while (TRUE) {
    # VMD decomposition
    result <- vmd(series, alpha = alpha, tau = tau, K = K, DC = DC, init = init, tol = tol)

    if (!is.null(result$u)) {
      if (verbose) {
        message("IMFs dimensions: ", paste(dim(result$u), collapse = " x "))
        if (!is.null(result$u_hat)) {
          message("U_hat dimensions: ", paste(dim(result$u_hat), collapse = " x "))
        }
        if (!is.null(result$omega)) {
          message("Omega dimensions: ", length(result$omega))
        }
      }
    } else {
      if (verbose) {
        message("IMFs are NULL")
      }
    }

    # Extract IMFs, estimated envelopes, and frequencies
    u <- result$u
    u_hat <- result$u_hat
    omega <- result$omega

    # Reconstruct series from IMFs
    reconstructed_series <- reconstruct_series(u)

    # Calculate Spearman correlation
    corr <- cor(series, reconstructed_series, method = "spearman")

    # Check the correlation coefficient threshold
    if (corr >= threshold || K >= max_modes) {
      optimal_imfs <- u  # Store IMFs for the optimal number of modes
      optimal_u_hat <- u_hat  # Store estimated envelopes
      optimal_omega <- omega  # Store frequencies
      break
    }

    # If threshold is not met, increase K and repeat
    K <- K + 1
  }

  # Return the optimal K, reconstructed series, and additional outputs
  Result <- list(optimal_K = K, reconstructed_series = reconstructed_series, imfs = optimal_imfs, u_hat = optimal_u_hat, omega = optimal_omega)
  return(Result)
}


