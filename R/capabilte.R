#' Indices de capabilité du procédé (Cp, Cpk, niveau sigma)
#'
#' Évalue l'aptitude d'un procédé à produire dans les tolérances, en
#' supposant que les mesures suivent une loi normale (Laplace-Gauss).
#'
#' \itemize{
#'   \item \eqn{C_p = (USL - LSL) / (6\sigma)} : capabilité potentielle
#'     (dispersion seule).
#'   \item \eqn{C_{pk} = \min((USL - \mu)/(3\sigma), (\mu - LSL)/(3\sigma))} :
#'     capabilité réelle (dispersion et centrage).
#'   \item Niveau sigma = \eqn{3 C_{pk}}.
#' }
#' Une seule limite peut être fournie (spécification unilatérale) : Cp vaut
#' alors NA et Cpk est calculé sur la limite disponible.
#'
#' @param x Vecteur numérique des mesures.
#' @param LSL Limite de spécification inférieure (NA si absente).
#' @param USL Limite de spécification supérieure (NA si absente).
#' @return Une liste de classe "capabilite" avec mu, sigma, n, Cp, Cpk,
#'   Cpl, Cpu, niveau_sigma, ppm et verdict.
#' @examples
#' set.seed(1)
#' x <- rnorm(200, mean = 10, sd = 0.03)
#' capabilite(x, LSL = 9.9, USL = 10.1)
#' @export
#' @importFrom stats sd pnorm
capabilite <- function(x, LSL = NA, USL = NA) {

  # --- Vérifications ---
  if (!is.numeric(x)) stop("x doit être un vecteur numérique.")
  x <- x[!is.na(x)]
  if (length(x) < 2) stop("Il faut au moins 2 mesures valides.")
  if (is.na(LSL) && is.na(USL))
    stop("Fournir au moins une limite de spécification (LSL ou USL).")
  if (!is.na(LSL) && !is.na(USL) && LSL >= USL)
    stop("LSL doit être strictement inférieure à USL.")

  mu    <- mean(x)
  sigma <- sd(x)
  if (sigma == 0) stop("Écart-type nul : indices non définis.")

  # --- Indices unilatéraux ---
  Cpu <- if (is.na(USL)) NA_real_ else (USL - mu) / (3 * sigma)
  Cpl <- if (is.na(LSL)) NA_real_ else (mu - LSL) / (3 * sigma)

  # --- Indices principaux ---
  Cp  <- if (is.na(LSL) || is.na(USL)) NA_real_ else (USL - LSL) / (6 * sigma)
  Cpk <- min(Cpu, Cpl, na.rm = TRUE)

  # --- Taux de non-conformité (loi normale) ---
  p_bas  <- if (is.na(LSL)) 0 else pnorm(LSL, mu, sigma)
  p_haut <- if (is.na(USL)) 0 else 1 - pnorm(USL, mu, sigma)

  verdict <- if (Cpk >= 1.67) "Excellent (Cpk >= 1.67)"
  else if (Cpk >= 1.33) "Capable (1.33 <= Cpk < 1.67)"
  else if (Cpk >= 1.00) "Juste capable (1.00 <= Cpk < 1.33)"
  else                  "Non capable (Cpk < 1.00)"

  res <- list(
    mu = mu, sigma = sigma, n = length(x),
    LSL = LSL, USL = USL,
    Cp = Cp, Cpk = Cpk, Cpl = Cpl, Cpu = Cpu,
    niveau_sigma = 3 * Cpk,
    ppm = (p_bas + p_haut) * 1e6,
    verdict = verdict
  )
  class(res) <- "capabilite"
  res
}

#' @export
print.capabilite <- function(x, ...) {
  cat("Capabilité du procédé (loi normale)\n")
  cat(sprintf("  n = %d | moyenne = %.4f | écart-type = %.4f\n",
              x$n, x$mu, x$sigma))
  if (!is.na(x$Cp))
    cat(sprintf("  Cp  = %.3f\n", x$Cp))
  cat(sprintf("  Cpk = %.3f\n", x$Cpk))
  cat(sprintf("  Niveau sigma = %.2f | non-conformité = %.0f ppm\n",
              x$niveau_sigma, x$ppm))
  cat("  Verdict : ", x$verdict, "\n", sep = "")
  if (!is.na(x$Cp) && x$Cp - x$Cpk > 0.2)
    cat("  Remarque : Cp nettement supérieur à Cpk, le procédé est décentré.\n")
  invisible(x)
}

