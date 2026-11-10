#' Test de normalité des mesures
#'
#' Vérifie l'hypothèse sur laquelle reposent les autres fonctions du
#' package : les mesures suivent-elles une loi de Laplace-Gauss ?
#' Combine un test formel (Shapiro-Wilk), des indicateurs de forme
#' (asymétrie et aplatissement) et deux graphiques (histogramme avec
#' courbe de Gauss, diagramme quantile-quantile).
#'
#' @param x Vecteur numérique des mesures (3 à 5000 valeurs).
#' @param alpha Seuil de risque du test (par défaut 0.05).
#' @param graphique Logique. Si TRUE, trace les deux graphiques.
#' @return Une liste de classe "test_normalite" avec n, mu, sigma,
#'   W, p_value, asymetrie, aplatissement, normal et verdict.
#' @examples
#' set.seed(1)
#' test_normalite(rnorm(150, 10, 0.05))
#' test_normalite(rexp(150))
#' @export
#' @importFrom stats shapiro.test sd dnorm qqnorm qqline
#' @importFrom graphics par hist lines
test_normalite <- function(x, alpha = 0.05, graphique = TRUE) {

  # --- Vérifications ---
  if (!is.numeric(x)) stop("x doit être un vecteur numérique.")
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 3)    stop("Il faut au moins 3 mesures valides.")
  if (n > 5000) stop("Shapiro-Wilk est limité à 5000 mesures.")
  if (sd(x) == 0) stop("Toutes les mesures sont identiques.")
  if (alpha <= 0 || alpha >= 1) stop("alpha doit être dans ]0 ; 1[.")

  mu    <- mean(x)
  sigma <- sd(x)

  # --- Test de Shapiro-Wilk ---
  sw <- shapiro.test(x)

  # --- Indicateurs de forme (0 pour une loi normale) ---
  z <- (x - mu) / sigma
  asymetrie     <- mean(z^3)
  aplatissement <- mean(z^4) - 3

  normal <- sw$p.value >= alpha
  verdict <- if (normal) {
    "Hypothèse de normalité non rejetée : les outils du package sont applicables."
  } else {
    "Normalité rejetée : prudence, les taux et indices calculés peuvent être faux."
  }

  res <- list(
    n = n, mu = mu, sigma = sigma,
    W = unname(sw$statistic), p_value = sw$p.value,
    asymetrie = asymetrie, aplatissement = aplatissement,
    alpha = alpha, normal = normal, verdict = verdict
  )
  class(res) <- "test_normalite"

  # --- Graphiques : histogramme + Gauss, et QQ-plot ---
  if (graphique) {
    ancien <- par(mfrow = c(1, 2))
    on.exit(par(ancien), add = TRUE)

    hist(x, freq = FALSE, col = "lightblue", border = "white",
         main = "Histogramme et loi normale ajustée",
         xlab = "Mesure", ylab = "Densité")
    t <- seq(min(x), max(x), length.out = 300)
    lines(t, dnorm(t, mu, sigma), col = "red", lwd = 2)

    qqnorm(x, main = "Diagramme quantile-quantile", pch = 16,
           xlab = "Quantiles théoriques", ylab = "Quantiles observés")
    qqline(x, col = "red", lwd = 2)
  }

  res
}

#' @export
print.test_normalite <- function(x, ...) {
  cat("Test de normalité (Shapiro-Wilk)\n")
  cat(sprintf("  n = %d | moyenne = %.4f | écart-type = %.4f\n",
              x$n, x$mu, x$sigma))
  cat(sprintf("  W = %.4f | p-value = %.4f (seuil %.2f)\n",
              x$W, x$p_value, x$alpha))
  cat(sprintf("  Asymétrie = %.3f | Aplatissement = %.3f\n",
              x$asymetrie, x$aplatissement))
  cat("  ", x$verdict, "\n", sep = "")
  invisible(x)
}
