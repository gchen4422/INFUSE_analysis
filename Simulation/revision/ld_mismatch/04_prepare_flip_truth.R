#!/usr/bin/env Rscript
.infuse_project <- Sys.getenv("INFUSE_PROJECT", "/path/to/baseline_simulation")
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)
base <- .infuse_project
module <- file.path(.infuse_work, "LD_mismatch_sensitivity")
out <- file.path(module,"inputs")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
truth <- list(); sources <- character()
for(c in 1:3) for(locus in 1:100) for(h in 1:2) {
  zpath <- function(rho) file.path(base,"More_annotes",paste0("causal_num_",c),
    "summary_data",paste0("flipped_0.01_rho_",rho),
    sprintf("CAUSAL_%d_LOCI_%d_h2_%d_flipped_0.01_rho_%d.txt",c,locus,h,rho))
  z <- fread(zpath(4)); positive <- fread(zpath(1))
  sources <- c(sources,zpath(4),zpath(1))
  stopifnot(identical(z$RSID,positive$RSID),identical(z$Signal,positive$Signal))
  for(k in 1:2) {
    col <- paste0("zscore_",k)
    original <- sign(positive[[col]])*abs(z[[col]])
    stopifnot(all(pmin(abs(positive[[col]]-original),
      abs(positive[[col]]-.5*original))<1e-10))
    idx <- which(abs(z[[col]]-original)>1e-10)
    stopifnot(length(idx)==ceiling(.01*nrow(z)))
    truth[[length(truth)+1L]] <- data.table(causal_index=c,locus=locus,
      h2_index=h,ancestry=k,row_index=idx,RSID=z$RSID[idx])
  }
}
fwrite(rbindlist(truth),file.path(out,"flip_truth.tsv"),sep="\t")
fwrite(data.table(path=sources,md5=unname(tools::md5sum(sources))),
  file.path(out,"flip_truth_sources.tsv"),sep="\t")
