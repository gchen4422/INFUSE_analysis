.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
# Use allele-matched baseline186 annotations from the same genome build.
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
root<-normalizePath(Sys.getenv("AOU_UNION_ROOT", file.path(.infuse_work, "INFUSE_variant_availability/real_data_AoU_union")),mustWork=TRUE)
path<-file.path(root,"summary/evaluation_annotation_rows.tsv.gz")
a<-fread(path)
m<-fread(file.path(root,"inputs/frozen_manifests/region_manifest.tsv"))
features<-c("Conserved_LindbladToh_common","H3K4me3_Trynka_common")
stopifnot(all(features %in% readLines(file.path(root,"inputs/baseline_186/annotation_names.txt"))),
 !anyDuplicated(m[,.(trait,locus)]),!anyDuplicated(a[,.(trait,locus,SNP)]))
oldcols<-setdiff(names(a),features)
old<-copy(a[,..oldcols])
chr<-m$CHR[match(paste(a$trait,a$locus),paste(m$trait,m$locus))]
stopifnot(!anyNA(chr))
values<-matrix(NA_real_,nrow(a),length(features),dimnames=list(NULL,features))
provenance<-list()
for(k in sort(unique(chr))){
 f<-file.path(root,"inputs/baseline_186",paste0("chr",k,".tsv"))
 b<-fread(f,select=c("SNP",features))
 stopifnot(!anyDuplicated(b$SNP))
 ii<-which(chr==k);hit<-match(a$SNP[ii],b$SNP)
 stopifnot(!anyNA(hit))
 values[ii,]<-as.matrix(b[hit,..features])
 provenance[[length(provenance)+1L]]<-data.table(CHR=k,records=length(ii),source=f,md5=unname(tools::md5sum(f)))
}
stopifnot(all(is.finite(values)),all(values %in% c(0,1)))
for(f in features){
 if(f %in% names(a))stopifnot(identical(as.numeric(a[[f]]),values[,f]))
 set(a,j=f,value=values[,f])
}
stopifnot(identical(old,a[,..oldcols]))
fwrite(a,path,sep="\t",na="NA")
fwrite(rbindlist(provenance),file.path(root,"summary/reference_enrichment_annotation_provenance.tsv"),sep="\t")
