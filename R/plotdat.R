.gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

.bracket <- function(x, width, y, height) {
    data.frame(x = (c(0, 1, 4, 5, 6, 9, 10)/10 - 0.5) * (width) + x, y = c(0, 1, 
        1, 2, 1, 1, 0)/2 * (height) + y)
}
plotsamExplorer <- function(dat, save = FALSE, filename = NULL, p.depth = 0.9, font.size = 3.5, 
    anova = TRUE, x.lab = NULL, y.lab = NULL, leg.lab = NULL) {
    pgg <- NULL
    addlines <- function(xc, yc, cl, lab = NULL, font.size) {
        annotate("text", label = lab, fontface = "bold", x = xc, y = yc, size = font.size, 
            colour = cl)
    }
    ytmp <- max(dat$Value)
    xtmp <- min(dat$Variable)
    
    L <- length(dat$Value)
    grp <- names(table(dat$Label))
    dpth <- as.numeric(names(table(dat$Variable)))
    temp1 <- which((dpth >= p.depth & dpth < 1) == TRUE)
    
    dpth <- dpth[temp1]
    nmfactor <- names(table(dat$Label))
    lenfactor <- length(table(dat$Label))
    xm1 <- min(dpth)
    xm1 <- which(xm1 == as.numeric(names(table(dat$Variable))))
    xm2 <- length(table(dat$Variable)) - 1
    cols = .gg_color_hue(lenfactor)
    new.df <- data.frame(fa = rep("# DE genes", L), f = as.factor(dat$Variable), 
        value = dat$Value, group = as.factor(dat$Label))
    colar <- "coral"
    colar2 <- "seagreen3"
    # xm1 <- 11.5 xm2 <- 14.5
    pgg <- ggplot(new.df, aes(x = f, y = value, fill = group, color = group)) + 
	geom_boxplot(notch = FALSE, , position = position_dodge(width = 0), outlier.size = 1) + 
        scale_fill_manual(name = "Annotation", values = cols) + scale_color_manual(guide = "none", 
        values = rep("gray40", lenfactor)) + ylab("")  + theme(legend.position = c(0.1, 
        0.85)) 
    #+ scale_fill_discrete(name = "Annotation")
    if (!is.null(leg.lab)) {
        pgg <- pgg + scale_fill_discrete(name = "Annotation", labels = leg.lab)
    }
    if (!is.null(x.lab)) {
        pgg <- pgg + xlab(x.lab)
    }
    if (!is.null(y.lab)) {
        pgg <- pgg + ylab(y.lab)
    }
    # scale_x_discrete(name='xxx') scale_y_discrete(name='yyy')
    
    # print(xm1)
    if (anova) {
        pv <- c()
        crd <- c()
        for (i in seq_along(grp)) {
            p <- summary(exploreRob(dat, grp[i], dpth))
            pv <- c(pv, round(p[[1]][["Pr(>F)"]][1], digits = 7))
            a <- which(dat$Label == grp[i])
            xtmp <- max(dat$Variable[a])
            ytmp <- max(dat$Value[a])
            crd <- rbind(crd, c(xtmp, ytmp))
        }
        
        cnt = 0
        for (i in seq_len(lenfactor)) {
            pgg <- pgg + addlines(length(table(dat$Variable)) - 3 + crd[i, 1], crd[i, 
                2] + 5, cols[i], lab = paste0("p-value", pv[i]), font.size)
            cnt = cnt - 50
        }
        # print(dpth)
        for (i in seq_len(lenfactor)) {
            tt <- which(dat$Label == nmfactor[i])
            # print(tt) for(j in seq_along(dpth)){
            tt1 <- which(dat[tt, ]$Variable >= dpth[1])
            ymn <- min(dat[tt, ]$Value[tt1])
            # print(ymn)
            ymx <- max(dat[tt, ]$Value[tt1])
            # }
            
            pgg <- pgg + annotate("rect", xmin = xm1 - 0.5, xmax = xm2 + 0.5, ymin = ymn, 
                ymax = ymx, alpha = 0.2)
        }
    }
    
    # annotate('text', label = 'r = 3', fontface='bold', x = 13.0, y = 90,
    # face='bold', size = font.size, colour = 'black')
    
    plot(pgg)
    # grid.brackets(425, 190, 340, 190, lwd=1, h=0.02)
    
    
    if (save) {
        fn <- NULL
        if (!is.null(filename)) {
            fn <- paste0(getwd(), "/", filename, ".pdf", sep = "")
        } else {
            fn <- paste0("plot", runif(1, min = 1000, max = 5000), ".pdf", sep = "")
            fn <- paste0(getwd(), "/", fn, sep = "")
        }
        ggsave(filename = fn, plot = pgg)
        message("file is saved as: ", fn)
    }
  invisible(pgg)
}
