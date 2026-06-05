# Example usage
H <- 28  # Hits
Fa <- 72  # False alarms
M <- 23  # Misses
CR <- 2680 # Correct rejections

# Define ETS function
ets <- function(H, Fa, M, CR) {
  # Calculate expected hits due to chance
  H_ar = (H + Fa) * (H + M) / (H + Fa + M + CR)

  # Calculate ETS
  ETS = (H - H_ar) / (H + Fa + M - H_ar)

  return(ETS)
}


ets(H, Fa, M, CR)

ets_alt <- function(H, Fa, M, CR) {
  # Numerator
  numerator <- H * CR - Fa * M

  # Denominator
  denominator <- H * Fa + H * M + H * CR + Fa^2 + M^2 + Fa * M + Fa * CR + M * CR

  # Calculate ETS
  ETS <- numerator / denominator

  return(ETS)
}
ets_alt(H, Fa, M, CR)

heidke <- function(H, Fa, M, CR) {
  # Calculate Heidke Skill Score (HSS)
  HSS = 2*(H * CR - Fa * M) / ((H + Fa) * (M + CR) + (H + M) * (Fa + CR))

  return(HSS)
}
heidke(H, Fa, M, CR)

pierce <- function(H, Fa, M, CR) {
  # Calculate Pierce Skill Score (PSS)
  PSS = (H / (H + M)) - (Fa / (CR + Fa))

  return(PSS)
}
pierce(H, Fa, M, CR)

