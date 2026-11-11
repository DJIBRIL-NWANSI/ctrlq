#' Animation de l'effet du décentrage sur la capabilité
#'
#' Fait glisser la moyenne d'un procédé N(mu, sigma^2) vers l'une des
#' limites de tolérance et affiche, image par image, la courbe de Gauss,
#' la zone hors tolérance (en rouge), Cpk et le taux de non-conformité.
#' L'écart-type reste constant : seul le centrage change, ce qui montre
#' qu'un procédé précis (Cp élevé) peut quand même produire des rebuts.
#'
#' Sans \code{fichier}, l'animation se joue dans la fenêtre graphique.
#' Avec \code{fichier}, elle est enregistrée en GIF (package gifski requis).
#'
#' @param LSL Limite de spécification inférieure.
#' @param USL Limite de spécification supérieure.
#' @param sigma Écart-type du procédé (constant).
#' @param mu_depart Moyenne initiale (par défaut, le centre de l'intervalle).
#' @param decalage_max Décalage total de la moyenne (par défaut un quart de
#'   l'intervalle de tolérance). Peut être négatif.
#' @param n_images Nombre d'images de l'aller (au moins 2).
#' @param aller_retour Si TRUE, la moyenne revient à sa position de départ.
#' @param fichier Chemin d'un GIF à créer, ou NULL pour jouer à l'écran.
#' @param delai Durée d'une image, en secondes.
#' @param largeur,hauteur Dimensions du GIF en pixels.
#' @return Invisiblement, une liste avec les moyennes (mu), les Cpk, les
#'   taux de non-conformité en ppm et le chemin du fichier.
#' @examples
#' animation_decentrage(LSL = 9.9, USL = 10.1, sigma = 0.03, delai = 0)
#' \dontrun{
#' animation_decentrage(9.9, 10.1, 0.03, fichier = "decentrage.gif")
#' }
#' @export
#' @importFrom graphics plot polygon abline text
#' @importFrom stats dnorm pnorm
animation_decentrage <- function(LSL, USL, sigma,
                                 mu_depart = (LSL + USL) / 2,
                                 decalage_max = (USL - LSL) / 4,
                                 n_images = 40, aller_retour = TRUE,
                                 fichier = NULL, delai = 0.1,
                                 largeur = 800, hauteur = 500) {

  # --- Vérifications ---
  if (!is.numeric(c(LSL, USL, sigma, mu_depart, decalage_max)))
    stop("LSL, USL, sigma, mu_depart et decalage_max doivent \u00eatre num\u00e9riques.")
  if (LSL >= USL) stop("LSL doit \u00eatre strictement inf\u00e9rieure \u00e0 USL.")
  if (sigma <= 0) stop("sigma doit \u00eatre strictement positif.")
  if (n_images < 2) stop("n_images doit valoir au moins 2.")

  # --- Trajectoire de la moyenne ---
  mus <- seq(mu_depart, mu_depart + decalage_max, length.out = n_images)
  if (aller_retour) mus <- c(mus, rev(mus[-c(1, n_images)]))

  # --- Calculs analytiques (aucune donnée simulée) ---
  cpk_f <- function(mu) min(USL - mu, mu - LSL) / (3 * sigma)
  ppm_f <- function(mu) (pnorm(LSL, mu, sigma) + 1 - pnorm(USL, mu, sigma)) * 1e6

  # --- Axes fixes pour toutes les images ---
  xlim <- range(c(LSL - sigma, USL + sigma, mus - 4 * sigma, mus + 4 * sigma))
  ymax <- dnorm(0, 0, sigma) * 1.3

  dessiner <- function(mu) {
    t   <- seq(xlim[1], xlim[2], length.out = 600)
    d   <- dnorm(t, mu, sigma)
    cpk <- cpk_f(mu)
    couleur <- if (cpk >= 1.33) "darkgreen" else if (cpk >= 1) "darkorange" else "red"

    plot(t, d, type = "l", lwd = 3, xlim = xlim, ylim = c(0, ymax),
         xlab = "Mesure", ylab = "Densit\u00e9",
         main = "Effet du d\u00e9centrage sur la capabilit\u00e9")
    bas <- t <= LSL
    if (sum(bas) > 1)
      polygon(c(t[bas], max(t[bas]), min(t[bas])), c(d[bas], 0, 0),
              col = "red", border = NA)
    haut <- t >= USL
    if (sum(haut) > 1)
      polygon(c(t[haut], max(t[haut]), min(t[haut])), c(d[haut], 0, 0),
              col = "red", border = NA)
    abline(v = c(LSL, USL), lty = 2, col = "red", lwd = 2)
    abline(v = mu, lty = 3, lwd = 2)
    text(mean(xlim), ymax * 0.97,
         sprintf("mu = %.3f   |   Cpk = %.2f   |   %.0f ppm hors tol\u00e9rance",
                 mu, cpk, ppm_f(mu)),
         col = couleur, font = 2, cex = 1.2)
  }

  # --- Lecture à l'écran ou export GIF ---
  if (is.null(fichier)) {
    for (mu in mus) {
      dessiner(mu)
      Sys.sleep(delai)
    }
  } else {
    if (!requireNamespace("gifski", quietly = TRUE))
      stop("Le package gifski est n\u00e9cessaire : install.packages(\"gifski\")")
    gifski::save_gif({ for (mu in mus) dessiner(mu) },
                     gif_file = fichier, width = largeur, height = hauteur,
                     delay = delai, progress = FALSE)
    message("Animation enregistr\u00e9e : ", fichier)
  }

  invisible(list(mu = mus,
                 Cpk = vapply(mus, cpk_f, numeric(1)),
                 ppm = vapply(mus, ppm_f, numeric(1)),
                 fichier = fichier))
}
