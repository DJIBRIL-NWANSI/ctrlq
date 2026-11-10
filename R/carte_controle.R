#' Carte de contrôle de Shewhart (moyennes ou valeurs individuelles)
#'
#' Trace une carte de contrôle fondée sur la loi normale : la moyenne d'un
#' sous-groupe suit N(mu, sigma^2 / n), donc les limites naturelles sont
#' centre +/- 3 sigma (environ 99,73 % des points en régime maîtrisé).
#'
#' \itemize{
#'   \item n > 1 : carte des moyennes (X-barre), sigma estimé par
#'     \eqn{\bar{R}/d_2}.
#'   \item n = 1 : carte des valeurs individuelles, sigma estimé par
#'     \eqn{\overline{MR}/1.128}.
#' }
#'
#' Deux règles de détection sont appliquées :
#' \enumerate{
#'   \item un point au-delà des limites de contrôle (+/- 3 sigma) ;
#'   \item 8 points consécutifs du même côté de la ligne centrale.
#' }
#'
#' @param x Vecteur numérique des mesures, dans l'ordre de production.
#' @param taille_echantillon Taille des sous-groupes (entier de 1 à 10).
#' @param graphique Logique. Si TRUE, trace la carte.
#' @return Une liste de classe "carte_controle" avec centre, LCL, UCL,
#'   sigma_estime, points, hors_limites et serie_8.
#' @examples
#' set.seed(1)
#' x <- c(rnorm(100, 10, 0.05), rnorm(25, 10.12, 0.05))
#' carte_controle(x, taille_echantillon = 5)
#' @export
#' @importFrom graphics plot abline points text
carte_controle <- function(x, taille_echantillon = 5, graphique = TRUE) {

  # --- Vérifications ---
  if (!is.numeric(x)) stop("x doit être un vecteur numérique.")
  x <- x[!is.na(x)]
  n <- taille_echantillon
  if (n != round(n) || n < 1 || n > 10)
    stop("taille_echantillon doit être un entier entre 1 et 10.")

  k <- floor(length(x) / n)
  if (k < 2) stop("Il faut au moins 2 sous-groupes complets.")
  if (length(x) > k * n) {
    warning("Mesures excédentaires ignorées (sous-groupe incomplet).")
    x <- x[seq_len(k * n)]
  }

  # --- Statistiques par sous-groupe ---
  d2 <- c(1.128, 1.128, 1.693, 2.059, 2.326, 2.534,
          2.704, 2.847, 2.970, 3.078)   # constantes d2 pour n = 1..10

  if (n == 1) {
    pts       <- x
    sigma_hat <- mean(abs(diff(x))) / d2[1]
  } else {
    M         <- matrix(x, nrow = k, ncol = n, byrow = TRUE)
    pts       <- rowMeans(M)
    R         <- apply(M, 1, function(v) max(v) - min(v))
    sigma_hat <- mean(R) / d2[n]
  }

  # --- Limites de contrôle (loi normale de la moyenne) ---
  centre <- mean(pts)
  sigma_pts <- sigma_hat / sqrt(n)
  UCL <- centre + 3 * sigma_pts
  LCL <- centre - 3 * sigma_pts

  # --- Règle 1 : points hors limites ---
  hors_limites <- which(pts > UCL | pts < LCL)

  # --- Règle 2 : 8 points consécutifs du même côté ---
  r     <- rle(sign(pts - centre))
  fin   <- cumsum(r$lengths)
  debut <- fin - r$lengths + 1
  serie_8 <- integer(0)
  for (i in seq_along(r$lengths)) {
    if (r$values[i] != 0 && r$lengths[i] >= 8)
      serie_8 <- c(serie_8, (debut[i] + 7):fin[i])
  }

  res <- list(
    centre = centre, LCL = LCL, UCL = UCL,
    sigma_estime = sigma_hat,
    taille_echantillon = n, k = k,
    points = pts,
    hors_limites = hors_limites,
    serie_8 = serie_8
  )
  class(res) <- "carte_controle"

  # --- Graphique ---
  if (graphique) {
    ylim <- range(c(pts, LCL, UCL))
    plot(pts, type = "b", pch = 16, ylim = ylim,
         xlab = "Sous-groupe", ylab = "Valeur",
         main = "Carte de contrôle de Shewhart")
    abline(h = centre, lwd = 2)
    abline(h = c(LCL, UCL), lty = 2, col = "red", lwd = 2)
    text(1, UCL, "LSC", pos = 3, col = "red")
    text(1, LCL, "LIC", pos = 1, col = "red")
    alertes <- sort(unique(c(hors_limites, serie_8)))
    if (length(alertes) > 0)
      points(alertes, pts[alertes], pch = 19, col = "red", cex = 1.5)
  }

  res
}

#' @export
print.carte_controle <- function(x, ...) {
  cat("Carte de contrôle de Shewhart\n")
  cat(sprintf("  %d sous-groupes de taille %d\n", x$k, x$taille_echantillon))
  cat(sprintf("  Centre = %.4f | LIC = %.4f | LSC = %.4f\n",
              x$centre, x$LCL, x$UCL))
  if (length(x$hors_limites) == 0 && length(x$serie_8) == 0) {
    cat("  Procédé sous contrôle statistique.\n")
  } else {
    if (length(x$hors_limites) > 0)
      cat("  Points hors limites : ",
          paste(x$hors_limites, collapse = ", "), "\n", sep = "")
    if (length(x$serie_8) > 0)
      cat("  Série de 8 points du même côté : sous-groupes ",
          paste(x$serie_8, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

