#' Cognitive risk sample correlation matrix
#'
#' Reported sample correlation matrix from a cross-sectional study on
#' cognitive risk and resilience factors in remitted depression patients, from
#' Hoorelebeke, Marchetti, DE Schryver and Koster (2016).
#' The study was conducted with 69 participants, and the correlation matrix
#' consists of six variables. The variables are as follows:
#'
#' * `BRIEF_WM`: working memory complaints, a self-report measure of perceived cognitive control
#' * `PASAT_ACC`: PASAT accuracy, performance on behavioural measure of congitive control
#' * `Adapt ER`: self-report adaptive emotion regulation strategies
#' * `Maladapt ER`: self-report maladaptive emotion regulation strategies
#' * `Resilience`: self-report resilience
#' * `Resid Depress`: self-report residual depressive symptoms
#'
#' @docType data
#'
#' @usage data(riskcor)
#'
#' @format A 6 by 6 correlation matrix
#'
#' @keywords datasets
#'
#' @source <https://ars.els-cdn.com/content/image/1-s2.0-S0165032715313252-mmc1.pdf>
#'
#' @importFrom Rdpack reprompt
#' @references
#'     \insertRef{hoorelbeke2016interplay}{SEset}
#' @examples
#' data(riskcor)
#' print(rownames(riskcor))
#' print(riskcor)
"riskcor"
