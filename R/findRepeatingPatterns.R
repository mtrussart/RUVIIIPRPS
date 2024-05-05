#' Find repeating factors or characters in a vector.

#' @author Ramyar Molania

#' @param vec A vector of factors or characters.
#' @param n.repeat Numeric. Indicates the minimum repeat of individual factors ro characters in the vector.

#' @return A vector of factors or characters that are repeated at least "n.repeat" times.

findRepeatingPatterns <- function(vec, n.repeat) {
    char.counts <- table(vec)
    repeated.chars <- subset(char.counts, char.counts >= n.repeat)
    return(names(repeated.chars))
}
