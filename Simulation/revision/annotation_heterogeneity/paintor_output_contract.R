# Shared PAINTOR output parsing and optimizer-status helpers.
paintor_valid_pip <- function(pip) {
  is.numeric(pip) && length(pip)>0L && all(is.finite(pip)) &&
    all(pip >= -1e-8 & pip <= 1+1e-8)
}

paintor_pip_status <- function(pip) {
  if(!paintor_valid_pip(pip)) return("invalid_pip")
  if(all(pmin(pmax(pip,0),1)==0)) "all_zero_no_selection" else "positive_pip"
}

paintor_log_diagnostics <- function(path) {
  txt <- if(file.exists(path)) readLines(path,warn=FALSE) else character()
  it <- grep("Enrichment estimates at iteration [0-9]+",txt,value=TRUE)
  iter <- if(length(it)) max(as.integer(sub(".*iteration ([0-9]+).*","\\1",it))) else NA_integer_
  list(log_available=file.exists(path),
       optimizer_error_logged=any(grepl("Optimization error encountered",txt,fixed=TRUE)),
       last_logged_iteration=iter,
       at_default_iteration_limit=if(is.na(iter)) NA else iter>=10L,
       convergence_assessment="not_certified_by_output_validation")
}

# Reference compute_set_pip returns an empty selection for a zero vector:
# 0/0 yields NaN, and which(cumsum(...) < coverage) has no indices.
# Make that same outcome explicit, preserving all cells in the denominator.
legacy_paintor_set <- function(pip, truth, coverage=0.95) {
  if(!paintor_valid_pip(pip)) stop("Invalid PAINTOR PIPs.")
  if(length(truth)!=length(pip) || anyNA(truth) ||
     length(coverage)!=1L || !is.finite(coverage) || coverage<=0 || coverage>1)
    stop("Invalid PAINTOR truth or coverage.")
  pip <- pmin(pmax(as.numeric(pip),0),1)
  selected <- integer()
  if(sum(pip)>0) {
    normalized <- pip/sum(pip)
    ordering <- order(normalized,decreasing=TRUE)
    selected <- ordering[cumsum(normalized[ordering])<coverage]
  }
  n <- length(selected);tp <- sum(as.logical(truth)[selected]);nc <- sum(as.logical(truth))
  data.table::data.table(Power=if(nc)tp/nc else 0,
                        FDR=if(n)(n-tp)/n else 0,
                        True_selected=tp,Size=n)
}
