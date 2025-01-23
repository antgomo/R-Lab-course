CellDMC
function (beta.m, pheno.v, frac.m, adjPMethod = "fdr", 
          adjPThresh = 0.05, cov.mod = NULL, sort = FALSE, mc.cores = 1) 

# CellDMC(betas.fis, pD.w0$R, cc.df,mc.cores = 6,cov.mod = Phen) 
pheno.v <- pD.w0$R
beta.m <- betas.fis; rm(betas.fis)
frac.m <- cc.df
cov.mod <- Phen
  
# start function 
# checkings 
  if (sum(is.na(pheno.v)) > 0) stop("No NA allowed in pheno.v!")
  if (ncol(beta.m) != length(pheno.v))  stop("Number of columns of beta.m should equal to length of pheno.v!")
  if (ncol(beta.m) != nrow(frac.m)) stop("Number of columns of beta.m should equal to number of rows of frac.m!")
  if (length(colnames(frac.m)) != ncol(frac.m)) stop("Pls assign correct names of cell-type to frac.m")
# look if betas are 0 to 1
  is.beta <- ((min(beta.m) >= 0) & (max(beta.m) <= 1))
# check if binary pheno
  if (nlevels(factor(pheno.v)) == 2) {
    message("Binary phenotype detected. Predicted change will be 1 - 0.")
    pheno.v <- factor(pheno.v)
  }
# Now pheno is "0" (NR) and "1" (R)

  if (!is.factor(pheno.v) & !is.character(pheno.v)) message("pheno.v is not factor or character. Treating as continuous variables.")
  # False
  
  # complete model matrix  
  design <- model.matrix(~frac.m + pheno.v:frac.m)[, -1]
  
  # look for rank
  if (Matrix::rankMatrix(design) < ncol(design)) {stop("The design matrix is not full ranked.\nThis means that you coundn't make inference for all cell-types in your fraction matrix.\n         This is usally casued by fractions of a cell-type of one pheno type are all 0 or some fractions in one pheno type are paralle to that of another cell-type.\n         You might use which(colSums(model.matrix(~ frac.m + pheno.v:frac.m)[, -1]) == 0) to find the cell type.")}
  
  # add covariates to model, if they exist
if (!is.null(cov.mod)) 
    design <- cbind(design, cov.mod[, -1])
  
  library(stringr)
  # change the names of the interacting coefficients
  IntNames.v <- str_c(colnames(frac.m), "Pheno") # cell type names
  colnames(design)[(1 + ncol(frac.m)):(2 * ncol(frac.m))] <- IntNames.v
  
  # apply the linear model
  detectCores() # we have 8
  mc.cores <- 6 # apply 
  allCoe.m <- do.call(rbind, mclapply(seq_len(nrow(beta.m)), 
                                      function(i) {
                                        beta.v <- beta.m[i, ]
                                        Int.o <- lm(beta.v ~ ., data = data.frame(design))
                                        IntCoe.m <- summary(Int.o)$coe[IntNames.v, ]
                                        IntCoe.v <- unlist(apply(IntCoe.m, 1, function(x) list(x)))
                                        names(IntCoe.v) <- NULL
                                        return(IntCoe.v)
                                      }, mc.preschedule = TRUE, mc.cores = mc.cores, mc.allow.recursive = TRUE))
  
  
  coe.ld <- lapply(seq_len(ncol(frac.m)), function(j) {
    idx <- ((j - 1) * 4 + 1):((j - 1) * 4 + 4)
    tmp.m <- allCoe.m[, idx]
    tmp.m <- cbind(tmp.m, p.adjust(tmp.m[, 4], method = adjPMethod))
    if (is.beta) {
      tmp.m[which(tmp.m[, 1] > 1), 1] <- 1
      tmp.m[which(tmp.m[, 1] < -1), 1] <- -1
    }
    colnames(tmp.m) <- c("Estimate", "SE", "t", 
                         "p", "adjP")
    rownames(tmp.m) <- rownames(beta.m)
    return(data.frame(tmp.m))
  })
  names(coe.ld) <- colnames(frac.m)
  dmct.m <- matrix(rep(0, ncol(frac.m) * nrow(beta.m)), ncol = ncol(frac.m))
  dmct.idx <- which(sapply(coe.ld, "[[", "adjP") < 
                      adjPThresh)
  dmct.m[dmct.idx] <- sign(sapply(coe.ld, "[[", "Estimate")[dmct.idx])
  dmc.v <- ifelse(rowAlls(dmct.m == 0), 0, 1)
  dmct.m <- cbind(dmc.v, dmct.m)
  colnames(dmct.m) <- c("DMC", colnames(frac.m))
  rownames(dmct.m) <- rownames(beta.m)
  if (sort) 
    coe.ld <- lapply(coe.ld, function(x) x[order(x$p), ])
  return(list(dmct = dmct.m, coe = coe.ld))
}
