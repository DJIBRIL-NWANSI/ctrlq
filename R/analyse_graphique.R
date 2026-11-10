#' Analyse graphique complète d'un procédé
#'
#' Produit un tableau de bord de 4 panneaux à partir de mesures
#' supposées suivre une loi normale (Laplace-Gauss) :
#' \enumerate{
#'   \item courbe de Gauss avec les zones hors tolérance
#'     (\code{\link{taux_conformite}}) ;
#'   \item diagramme quantile-quantile (vérification de la normalité) ;
#'   \item carte de contrôle de Shewhart (\code{\link{carte_controle}}) ;
#'   \item synthèse chiffrée : conformité, Cp, Cpk, niveau sigma,
#'     test de Shapiro-Wilk.
#' }
#'
#' @param x Vecteur numérique des mesures, dans l'ordre de production.
#' @param LSL Limite de spécification inférieure.
#' @param USL Limite de spécification supérieure.
#' @param taille_echantillon Taille des sous-groupes de la carte de contrôle.
#' @return Invisiblement, une liste avec les résultats des quatre analyses
#'   (conformite, capabilite, carte, normalite).
#' @examples
#' set.seed(1)
#' x <- c(rnorm(100, 10, 0.04), rnorm(25, 10.06, 0.04))
#' analyse_graphique(x, LSL = 9.9, USL = 10.1)
#' @export
#' @importFrom graphics par plot.new text
#' @importFrom stats qqnorm qqline
analyse_graphique <- function(x, LSL, USL, taille_echantillon = 5) {

  if (!is.numeric(x)) stop("x doit être un vecteur numérique.")
  x <- x[!is.na(x)]
  if (missing(LSL) || missing(USL))
    stop("LSL et USL sont nécessaires pour l'analyse graphique.")

  # --- Calculs (sans graphique) ---
  conf <- taux_conformite(x, LSL, USL, graphique = FALSE)
  cap  <- capabilite(x, LSL, USL)
  norm <- test_normalite(x, graphique = FALSE)

  # --- Mise en page 2 x 2 ---
  ancien <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  on.exit(par(ancien), add = TRUE)

  # Panneau 1 : courbe de Gauss et zones hors tolérance
  taux_conformite(x, LSL, USL, graphique = TRUE)

  # Panneau 2 : QQ-plot
  qqnorm(x, main = "Vérification de la normalité", pch = 16,
         xlab = "Quantiles théoriques", ylab = "Quantiles observés")
  qqline(x, col = "red", lwd = 2)

  # Panneau 3 : carte de contrôle
  carte <- carte_controle(x, taille_echantillon, graphique = TRUE)

  # Panneau 4 : synthèse chiffrée
  plot.new()
  titre <- "SYNTHÈSE"
  lignes <- c(
    sprintf("n = %d | moyenne = %.4f | écart-type = %.4f",
            conf$n, conf$mu, conf$sigma),
    sprintf("Conformes : %.2f %% (%.0f ppm hors tolérance)",
            100 * conf$p_conforme, conf$ppm),
    sprintf("Cp = %.2f | Cpk = %.2f | niveau sigma = %.2f",
            cap$Cp, cap$Cpk, cap$niveau_sigma),
    cap$verdict,
    sprintf("Shapiro-Wilk : p = %.3f (%s)", norm$p_value,
            if (norm$normal) "normalité OK" else "normalité rejetée"),
    if (length(carte$hors_limites) + length(carte$serie_8) == 0)
      "Carte : procédé sous contrôle"
    else
      sprintf("Carte : %d alerte(s), procédé instable",
              length(unique(c(carte$hors_limites, carte$serie_8))))
  )
  text(0.5, 0.95, titre, font = 2, cex = 1.3)
  text(0.02, seq(0.8, 0.1, length.out = length(lignes)),
       lignes, adj = 0, cex = 0.95)

  invisible(list(conformite = conf, capabilite = cap,
                 carte = carte, normalite = norm))
}
