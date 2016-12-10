exploreRob <- function(df_d, lbl, f_vect) {
    if (lbl %in% df_d[, 1] == FALSE) {
        stop("Invalid label")
    }
    for (f in f_vect) {
        if (f %in% df_d[, 2] == FALSE) {
            stop("Invalid f set")
        }
    }
    df_d_sub <- df_d[df_d[, 1] == lbl, ]
    df_d_sub <- df_d_sub[df_d_sub[, 2] %in% f_vect, 2:3]
    message("ANOVA test for label '", lbl, "' and f values ", paste0(as.character(f_vect), 
        collapse = ", "))
    anv <- aov(df_d_sub[, 2] ~ df_d_sub[, 1], data = df_d_sub)
    invisible(anv)
}

exploreRep <- function(df_d, lbl_vect, f) {
    if (f %in% df_d[, 2] == FALSE) {
        stop("Invalid f")
    }
    for (lbl in lbl_vect) {
        if (lbl %in% df_d[, 1] == FALSE) {
            stop("Invalid lbl set")
        }
    }
    df_d_sub <- df_d[df_d[, 2] == f, ]
    df_d_sub <- df_d_sub[df_d_sub[, 1] %in% lbl_vect, c(1, 3)]
    message("ANOVA test for labels '", paste0(as.character(lbl_vect), collapse = "', '"), 
        "' and f value ", f)
    anv <- aov(df_d_sub[, 2] ~ df_d_sub[, 1], data = df_d_sub)
    invisible(anv)
}

