# create a new smoothing function

# LOGIC
# Baseline Correction: Moving minimum find the "floor."
# SG Smoothing: One pass with a wide window (n=41) creates the curvy line.
# Shape Detection: The 2nd derivative identifies all candidate "humps."
# Targeting: The algorithm calculates the absolute difference |Apex_{RT} - Target_{RT}|.
# Selection: The hump with the smallest difference is chosen.
# Integration: Boundaries are set at X\% of that specific hump's height.

#' @author Eewon Tai


library(signal)

# --- 1. THE ALGORITHM FUNCTIONS ---

# Step 0: Moving Minimum Baseline
simple_baseline <- function(y, window = 60) {
  L <- length(y)
  base <- sapply(1:L, function(i) {
    start <- max(1, i - window/2)
    end <- min(L, i + window/2)
    min(y[start:end])
  })
  return(base)
}

# The Processor: Baseline -> SG Smooth -> 2nd Deriv -> RT-Based Selection
peak_processor_rt_targeted <- function(x, y,
                                       expected_rt,     # The library RT for your compound
                                       base_window = 60,
                                       sg_window = 41,
                                       stringency = 0.4,
                                       edge_percent = 0.05,
                                       p = 3) {  # p: 1,2,3 (degree of polynomial curve fitted)

  # A. Baseline Correction
  base <- simple_baseline(y, window = base_window)
  y_corr <- y - base

  # B. SG Smoothing (The Curvy Line)
  if (sg_window %% 2 == 0) sg_window <- sg_window + 1
  y_smooth <- as.vector(sgolayfilt(y_corr, p = p, n = sg_window))

  # C. 2nd Derivative (Inflection Detection)
  dy2 <- as.vector(sgolayfilt(y_corr, p = p, n = sg_window, m = 2))
  threshold <- stringency * min(dy2, na.rm = TRUE)

  # D. Find all "Peaky" Candidates
  potential_indices <- which(dy2 < threshold & y_smooth > (max(y_smooth) * 0.05))

  if (length(potential_indices) == 0) return(list(y_smooth=y_smooth, base=base, peak=NULL))

  # Group contiguous indices into individual peak candidates
  peak_groups <- split(potential_indices, cumsum(c(1, diff(potential_indices) > 5)))

  peaks_list <- lapply(peak_groups, function(group) {
    apex_idx <- group[which.max(y_smooth[group])]
    peak_height <- y_smooth[apex_idx]
    cutoff <- peak_height * edge_percent

    # Boundary Expansion (Step 3)
    l_idx <- apex_idx; while(l_idx > 1 && y_smooth[l_idx] > cutoff) l_idx <- l_idx - 1
    r_idx <- apex_idx; while(r_idx < length(y_smooth) && y_smooth[r_idx] > cutoff) r_idx <- r_idx + 1

    # Area - trapezoid
    p_range <- l_idx:r_idx
    area <- sum(diff(x[p_range]) * (y_smooth[p_range][-1] + y_smooth[p_range][-length(p_range)]) / 2)

    return(data.frame(apex_rt=x[apex_idx], left_rt=x[l_idx], right_rt=x[r_idx],
                      height=peak_height, area=area))
  })

  all_candidates <- do.call(rbind, peaks_list)

  # --- E. NEW SELECTION LOGIC: Closest to Expected RT ---
  # Calculate distance from each candidate's apex to our target RT
  all_candidates$rt_diff <- abs(all_candidates$apex_rt - expected_rt)

  # Sort by distance and pick the one with the smallest difference
  winner <- all_candidates[order(all_candidates$rt_diff), ][1, ]

  return(list(y_smooth = y_smooth, base = base, peak = winner, all_found = all_candidates))
}

# --- 2. IMPROVED RANDOMIZED TEST EXECUTION ---
# The second half of the script is a Simulation to prove the algorithm works.
# Creates a Fake Metabolite: A random peak between RT 45 and 55.
# Adds 'Distractors': It adds up to 4 other peaks that might be much taller than your target to try and 'trick' the computer.
# Adds Noise & Tilt: It adds random static and a tilting baseline.
# The Visualization: * Grey dots: Raw messy data.
# Green dashed line: The identified baseline.
# Blue line: The smoothed 'clean' signal.
# Red shaded area: The specific peak the algorithm chose because it was closest to your target (RT 50).


# Reset plotting layout to 2x2
# par(mfrow=c(2,2))
#
# for(trial in 1:4) {
#   t_vec <- seq(0, 100, length.out = 1000)
#
#   # 1. Randomize our "Metabolite of Interest" (Target)
#   # We expect it at RT 50, but let's say it shifts slightly between samples
#   true_rt <- runif(1, 45, 55)
#   target_h <- runif(1, 0.5, 1.5)
#   target_w <- runif(1, 1.5, 3.5)
#   target_peak <- target_h * exp(-(t_vec - true_rt)^2 / (2 * target_w^2))
#
#   # 2. Add "Distractor" Peaks (Interference)
#   # These are peaks that are NOT our compound but appear in the same range
#   n_distractors <- sample(0:4, 1)
#   distractor_signal <- rep(0, 1000)
#   if(n_distractors > 0) {
#     for(d in 1:n_distractors) {
#       # Distractors can be anywhere from RT 10 to 90
#       d_rt <- runif(1, 10, 90)
#       # Ensure distractors aren't exactly on top of our target for the test
#       if(abs(d_rt - 50) < 5) d_rt <- d_rt + 10
#
#       d_h <- runif(1, 0.2, 3.0) # Can be much bigger than our target!
#       d_w <- runif(1, 1, 5)
#       distractor_signal <- distractor_signal + d_h * exp(-(t_vec - d_rt)^2 / (2 * d_w^2))
#     }
#   }
#
#   # 3. Add Random Baseline Lift and Noise
#   # Baseline could be high, low, or tilted
#   base_lift <- (t_vec * runif(1, -0.01, 0.01)) + runif(1, 0.1, 1.0)
#   noise_level <- runif(1, 0.01, 0.05)
#   raw_y <- target_peak + distractor_signal + base_lift + rnorm(1000, 0, noise_level)
#
#   # 4. RUN YOUR NEW ALGORITHM
#   # We tell it we expect our compound at RT 50
#   res <- peak_processor_rt_targeted(t_vec, raw_y, expected_rt = 50,
#                                     base_window = 60, sg_window = 41, stringency = 0.4,
#                                     edge_percent = 0.05, p = 2)  # change params
#
#   # 5. PLOTTING
#   plot(t_vec, raw_y, col="grey80", pch=16, cex=0.3,
#        main=paste("Trial", trial, "| Found RT:", round(res$peak$apex_rt, 1)),
#        xlab="Time", ylab="Intensity")
#
#   # The "Curvy Line" (Blue) and "Floor" (Green)
#   lines(t_vec, res$base, col="darkgreen", lty=2, lwd=1)
#   lines(t_vec, res$y_smooth + res$base, col="blue", lwd=1.5)
#
#   # Mark the "Target RT" we were looking for
#   abline(v = 50, col="black", lty=3)
#
#   # Shade the Peak that the algorithm selected as "Closest to 50"
#   if(!is.null(res$peak)) {
#     p_rt <- t_vec[t_vec >= res$peak$left_rt & t_vec <= res$peak$right_rt]
#     p_y  <- (res$y_smooth + res$base)[t_vec >= res$peak$left_rt & t_vec <= res$peak$right_rt]
#     polygon(c(min(p_rt), p_rt, max(p_rt)), c(min(res$base), p_y, min(res$base)),
#             col=rgb(1, 0, 0, 0.3), border="red")
#
#     # Label the height and area
#     text(res$peak$apex_rt, res$peak$height + res$base[1],
#          labels=paste("A:", round(res$peak$area,1)), pos=3, cex=0.8, col="red")
#   }
# }
#
# # Reset layout
# par(mfrow=c(1,1))
