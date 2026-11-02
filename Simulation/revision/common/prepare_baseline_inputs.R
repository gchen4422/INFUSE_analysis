.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
library(data.table)
setDTthreads(1)
alignment <- file.path(.infuse_work, "baseline_figure2_alignment_20260907")
source(file.path(.infuse_repo, "Simulation/revision/common/baseline_sources.R"))
setting <- commandArgs(TRUE)[1]
spec <- figure2_baseline_spec(setting)
e <- new.env(); load(spec$legacy_cache, e)
d <- as.data.table(e$data_all)
stopifnot(nrow(d) == 1314114L, !anyDuplicated(d[, .(causal_num,h2,locus,RSID)]))
n_cells <- 0L
for(k in 1:3) for(h in 1:2) for(l in 1:100) {
  stem <- sprintf("CAUSAL_%d_LOCI_%d_h2_%d",k,l,h)
  result <- file.path(spec$result_root,paste0("causal_num_",k),spec$result_subdir,paste0("Res_",stem,".RData"))
  truth <- file.path(spec$data_root,paste0("causal_num_",k),spec$data_subdir,stem)
  z <- fread(truth); fenv <- new.env(); load(result,fenv)
  idx <- which(d$causal_num == k & d$h2 == h & d$locus == l)
  idx <- idx[match(z$RSID, d$RSID[idx])]
  stopifnot(length(idx) == nrow(z), !anyNA(idx), all(d$Signal[idx] == z$Signal))
  candidate <- list(
    MAFPC_Either=fenv$MESuSiE_res_mlk_annot_subset$pip,
    MAFPC_WB=fenv$MESuSiE_res_mlk_annot_subset$pip_config[,1],
    MAFPC_BB=fenv$MESuSiE_res_mlk_annot_subset$pip_config[,2],
    MAFPC_Shared=fenv$MESuSiE_res_mlk_annot_subset$pip_config[,3],
    MESuSiE_Either=fenv$MESuSiE_res$pip,
    MESuSiE_WB=fenv$MESuSiE_res$pip_config[,1],
    MESuSiE_BB=fenv$MESuSiE_res$pip_config[,2],
    MESuSiE_Shared=fenv$MESuSiE_res$pip_config[,3],
    SuSiE_Either=pmax(fenv$susie_EU$pip,fenv$susie_BB$pip),
    SuSiE_Shared=pmin(fenv$susie_EU$pip,fenv$susie_BB$pip),
    SuSiE_WB=fenv$susie_EU$pip,SuSiE_BB=fenv$susie_BB$pip)
  for(suffix in c(".mcmc.paintor",".fun.mcmc.paintor")) {
    p <- fread(file.path(dirname(result),paste0(stem,suffix)))
    stopifnot(!anyDuplicated(p$RSID),setequal(p$RSID,z$RSID))
    nm <- if(suffix == ".mcmc.paintor") "Paintor_PIP" else "Paintor_annot_pip"
    candidate[[nm]] <- p$Posterior_Prob[match(z$RSID,p$RSID)]
  }
  for(nm in names(candidate)) set(d, i=idx, j=nm, value=as.numeric(candidate[[nm]]))
  n_cells <- n_cells + 1L
}
dir.create(dirname(spec$verified_cache),recursive=TRUE,showWarnings=FALSE)
saveRDS(list(data_all=as.data.frame(d),baseline=spec,verified_cells=n_cells),
        spec$verified_cache,compress="gzip")
cat("Prepared baseline PIPs:",spec$verified_cache,"\n")
