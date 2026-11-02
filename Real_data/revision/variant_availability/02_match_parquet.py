#!/usr/bin/env python3
"""Resolve baseline annotations by rsID, chromosome and unordered allele pair.

Writes coverage and matched records; never chooses a missing-data policy.
"""
import csv
import json
import os
from collections import defaultdict
from pathlib import Path
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

root = Path(os.environ.get("AOU_UNION_ROOT", str(Path(os.environ.get("INFUSE_REVISION", "/path/to/revision_output")) / "INFUSE_variant_availability/real_data_AoU_union"))).resolve()
src = Path(os.environ.get("INFUSE_AOU_INPUTS", "/path/to/aou_inputs"))
baseline = Path(os.environ.get("INFUSE_BASELINE_PARQUET", "/path/to/baseline_annotation_parquet"))
with (root/'qc/locus_manifest.tsv').open() as f:
    manifest = list(csv.DictReader(f,delimiter='\t'))
bychr = defaultdict(dict)
per_region = {}
for row in manifest:
    union = {}
    for anc in ('eur','afr'):
        path = src/f"results_{row['trait']}"/'summary_data'/f"{anc}_loci_{row['locus']}_summ"
        with path.open() as f:
            for ss in csv.DictReader(f):
                key = (int(ss['CHR']),ss['SNP'],tuple(sorted((ss['ALT'],ss['REF']))))
                bychr[key[0]][key[1]] = key[2]
                union.setdefault(key,set()).add(anc)
    per_region[(row['trait'],row['locus'])] = union
outdir = root/'inputs/baseline_186'
outdir.mkdir(parents=True,exist_ok=True)
status = {}
allqc = []
feature_names = None
for chrom in sorted(bychr):
    source = baseline/f'baselineLF2.2.UKB.{chrom}.annot.parquet'
    pf = pq.ParquetFile(source)
    features = [x for x in pf.schema_arrow.names if x not in ('CHR','SNP','BP','A1','A2')]
    assert len(features)==187 and 'GTEx_eQTL_MaxCPP_common' in features
    fitnames = [x for x in features if x!='GTEx_eQTL_MaxCPP_common']
    if feature_names is None:
        feature_names=fitnames
    assert feature_names==fitnames
    wanted = pa.array(list(bychr[chrom]),type=pa.string())
    records = defaultdict(list)
    seen = set()
    for batch in pf.iter_batches(batch_size=65536):
        table = pa.Table.from_batches([batch])
        table = table.filter(pc.is_in(table['SNP'],value_set=wanted))
        for r in table.to_pylist():
            snp = r['SNP']; seen.add(snp)
            if int(r['CHR'])!=chrom or tuple(sorted((r['A1'],r['A2'])))!=bychr[chrom][snp]:
                continue
            vals = tuple(r[n] for n in fitnames)
            records[snp].append((r['BP'],r['A1'],r['A2'],vals))
    with (outdir/f'chr{chrom}.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n')
        w.writerow(['SNP','CHR','annotation_BP','A1','A2']+fitnames)
        for snp,alleles in bychr[chrom].items():
            key = (chrom,snp,alleles)
            matches=records[snp]
            if not matches:
                status[key]='allele_mismatch' if snp in seen else 'absent_rsid'
            elif len({x[3] for x in matches})>1:
                status[key]='ambiguous_annotation'
            else:
                status[key]='matched'
                bp,a1,a2,vals=matches[0]
                w.writerow([snp,chrom,bp,a1,a2]+list(vals))
    allqc.append(dict(CHR=chrom,input_rsids=len(bychr[chrom]),matched=sum(status[(chrom,s,a)]=='matched' for s,a in bychr[chrom].items()),
      allele_matched_duplicate_rsids=sum(len(v)>1 for v in records.values()),source=str(source)))
    print(allqc[-1],flush=True)
with (outdir/'annotation_names.txt').open('w') as f:
    f.write('\n'.join(feature_names)+'\n')
with (root/'qc/annotation_match_status.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['trait','locus','SNP','CHR','alleles','availability','annotation_status'])
    for (trait,locus),union in per_region.items():
        for key,ancs in union.items():
            w.writerow([trait,locus,key[1],key[0],'/'.join(key[2]),'+'.join(sorted(ancs)),status[key]])
with (root/'qc/annotation_coverage_by_locus.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['trait','locus','union_n','matched','absent_rsid','allele_mismatch','ambiguous_annotation'])
    for (trait,locus),union in per_region.items():
        w.writerow([trait,locus,len(union)]+[sum(status[k]==label for k in union) for label in ('matched','absent_rsid','allele_mismatch','ambiguous_annotation')])
with (root/'qc/parquet_annotation_provenance.json').open('w') as f:
    json.dump(dict(match_key='chromosome+rsid+unordered_alleles; no BP matching across builds',
      excluded_feature='GTEx_eQTL_MaxCPP_common',fitting_features=feature_names,chromosomes=allqc),f,indent=2)
