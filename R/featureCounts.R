samExplore <- function(..., subsample_d = 1, N_boot = 1, countboot = c("all", "Assigned", 
    "Unassigned_Ambiguity", "Unassigned_MultiMapping", "Unassigned_NoFeatures", "Unassigned_Unmapped", 
    "Unassigned_MappingQuality", "Unassigned_FragmentLength", "Unassigned_Chimera", 
    "Unassigned_Secondary", "Unassigned_Nonjunction", "Unassigned_Duplicate")) {
    
    .printParam(subsample_d, N_boot)
    
    fc <- NULL
    fc <- featureCounts(...)
    if (is.null(fc)) {
        stop("no output returned")
    }
    ts <- vector("character", length(fc$targets))
    resultlist <- vector("list", length(fc$targets))
    for (i in seq_along(fc$targets)) {
        fctmp <- list()
        fctmp[["counts"]] <- as.matrix(fc$counts[, i])
        fctmp[["stat"]] <- fc$stat[, c(1, (i + 1))]
        message("Sub sampling: ",fc$targets[i])
        message("Sub sample N: ", i)
        restmp <- .featureCountsBoot(fctmp, subsample_d = subsample_d, N_boot = N_boot, 
            countboot = countboot)
        
        resultlist[[i]] <- restmp[[1]]
        ts[i] <- restmp[[2]]
    }
    names(resultlist) <- fc$targets
    names(ts) <- fc$targets
    list(bootres = resultlist, target.size = ts, featuremain = fc)
}


.featureCountsBoot <- function(flcount, subsample_d = 1, N_boot = 10, countboot = "all") {
    ts <- seq <- NULL
    totalcount <- as.numeric(flcount$stat[, 2])
    names(totalcount) <- as.vector(flcount$stat[, 1])
    if (countboot[1] == "all") {
        ts <- round(sum(totalcount) * subsample_d)
        seqtmp <- as.vector(flcount$counts)
        names(seqtmp) <- rownames(flcount$counts)
        seq <- c(seqtmp, totalcount[-which(names(totalcount) == "Assigned")])
    } else {
        countboot <- unique(c("Assigned", countboot))
        ts <- round(sum(totalcount[countboot]) * subsample_d)
        seqtmp <- as.vector(flcount$counts)
        names(seqtmp) <- rownames(flcount$counts)
        tmpx <- totalcount[countboot]
        tmpx <- tmpx[-which(names(tmpx) == "Assigned")]
        seq <- c(seqtmp, tmpx)
    }
    bootmat <- c()
    
    for (i in seq_len(N_boot)) {
        tmp <- thinCounts(seq, target.size = ts)
        bootmat <- cbind(bootmat, tmp[rownames(flcount$counts), ])
    }
    list(bootmat, target.size = unique(ts, sum(tmp)))
    
}
.printParam <- function(depth = 1, boot = 1) {
    message("Simulated Depth: ", depth * 100, "%")
    message("Number of sub samples: ", boot)
}
