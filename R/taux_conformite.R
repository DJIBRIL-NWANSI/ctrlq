#' Taux de conformité sous la loi de Laplace-Gauss
#'
#' Estime la moyenne et l'écart-type d'une série de mesures, puis calcule
#' la proportion de pièces conformes et non conformes par rapport aux
#' limites de spécification, en supposant X ~ N(mu, sigma^2).
#'
#' @param x Vecteur numérique des mesures (ex. diamètre en mm).
#' @param LSL Limite de spécification inférieure.
#' @param USL Limite de spécification supérieure.
#' @param graphique Logique. Si TRUE, trace la courbe de Gauss et les zones.
#' @return Une liste de classe "taux_conformite" avec mu, sigma, n,
#'   p_conforme, p_sous_LSL, p_sur_USL, p_non_conforme et ppm.
#' @examples
#' set.seed(1)
#' x <- rnorm(200, mean = 10, sd = 0.05)
#' taux_conformite(x, LSL = 9.9, USL = 10.1)
#' @importFrom stats sd pnorm dnorm
#' @importFrom graphics plot polygon abline
#' @export
taux_conformite <- function(x, LSL, USL, graphique = TRUE) {

  # --- Vérifications ---
  if (!is.numeric(x)) stop("x doit être un vecteur numérique.")
  x <- x[!is.na(x)]
  if (length(x) < 2) stop("Il faut au moins 2 mesures valides.")
  if (LSL >= USL) stop("LSL doit être strictement inférieure à USL.")

  # --- Estimation des paramètres de la loi normale ---
  mu    <- mean(x)
  sigma <- sd(x)

  # --- Probabilités via la fonction de répartition ---
  p_sous_LSL <- pnorm(LSL, mean = mu, sd = sigma)
  p_sur_USL  <- 1 - pnorm(USL, mean = mu, sd = sigma)
  p_conforme <- 1 - p_sous_LSL - p_sur_USL

  res <- list(
    mu = mu, sigma = sigma, n = length(x),
    LSL = LSL, USL = USL,
    p_conforme     = p_conforme,
    p_sous_LSL     = p_sous_LSL,
    p_sur_USL      = p_sur_USL,
    p_non_conforme = 1 - p_conforme,
    ppm            = (1 - p_conforme) * 1e6
  )
  class(res) <- "taux_conformite"

  # --- Graphique : courbe de Gauss + zones hors tolérance ---
  if (graphique) {
    t <- seq(mu - 4 * sigma, mu + 4 * sigma, length.out = 500)
    d <- dnorm(t, mu, sigma)
    plot(t, d, type = "l", lwd = 2,
         xlab = "Mesure", ylab = "Densité",
         main = "Loi normale et limites de spécification")
    g <- t <= LSL; polygon(c(t[g], max(t[g]), min(t[g])),
                           c(d[g], 0, 0), col = "red", border = NA)
    h <- t >= USL; polygon(c(t[h], max(t[h]), min(t[h])),
                           c(d[h], 0, 0), col = "red", border = NA)
    abline(v = c(LSL, USL), lty = 2, col = "red")
    abline(v = mu, lty = 3)
  }

  res
}

#' @export
print.taux_conformite <- function(x, ...) {
  cat("Taux de conformité (loi normale)\n")
  cat(sprintf("  n = %d | moyenne = %.4f | écart-type = %.4f\n",
              x$n, x$mu, x$sigma))
  cat(sprintf("  Spécifications : [%.4f ; %.4f]\n", x$LSL, x$USL))
  cat(sprintf("  Conforme       : %.2f %%\n", 100 * x$p_conforme))
  cat(sprintf("  Non conforme   : %.2f %% (%.0f ppm)\n",
              100 * x$p_non_conforme, x$ppm))
  invisible(x)
}

