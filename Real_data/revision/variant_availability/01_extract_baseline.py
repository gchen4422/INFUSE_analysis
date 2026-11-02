#!/usr/bin/env python3
"""Stream existing baselineLF2.2 records by rsID; do not compare hg19/hg38 BP.

This creates coverage diagnostics only. Missing/conflicting records are not
silently imputed or dropped from the intended fine-mapping universe.
"""
import csv
import json
import os
from pathlib import Path

root = Path(os.environ.get("AOU_UNION_ROOT", str(Path(os.environ.get("INFUSE_REVISION", "/path/to/revision_output")) / "INFUSE_variant_availability/real_data_AoU_union"))).resolve()
src = Path(os.environ.get("INFUSE_AOU_INPUTS", "/path/to/aou_inputs"))
baseline = Path(os.environ.get("INFUSE_BASELINE_ANNOTATIONS", "/path/to/baselineLF2.2.UKB")) / "annot_combine/combined_baselineLF2.2.UKB.annot"
with (root/'qc/locus_manifest.tsv').open() as f:
    manifest = list(csv.DictReader(f,delimiter='\t'))
variants = {}
for row in manifest:
    for anc in ('eur','afr'):
        path = src/f"results_{row['trait']}"/'summary_data'/f"{anc}_loci_{row['locus']}_summ"
        with path.open() as f:
            for ss in csv.DictReader(f):
                v = (ss['CHR'], ss['POS'], ss['ALT'], ss['REF'])
                variants.setdefault(ss['SNP'],set()).add(v)
conflict = {k:sorted(v) for k,v in variants.items() if len(v)!=1}
with (root/'qc/rsid_input_identity_conflicts.json').open('w') as f:
    json.dump(conflict,f,indent=2)
if conflict:
    raise ValueError(f'{len(conflict)} rsIDs map to inconsistent input identities; review before matching')
found, duplicated, conflicting = {}, {}, set()
outdir = root/'inputs/baseline_records'
outdir.mkdir(exist_ok=True,parents=True)
with baseline.open() as f, (outdir/'matched_records.tsv').open('w') as out:
    reader = csv.reader(f,delimiter='\t')
    header = next(reader)
    assert len(header)==191 and header[:4]==['CHR','BP','SNP','CM']
    exclude = header.index('GTEx_eQTL_MaxCPP_common')
    fitting = [x for i,x in enumerate(header) if i>=4 and i!=exclude]
    assert len(fitting)==186
    writer = csv.writer(out,delimiter='\t',lineterminator='\n')
    writer.writerow(header)
    for record in reader:
        snp = record[2]
        if snp not in variants:
            continue
        assert len(record)==191
        input_chr = next(iter(variants[snp]))[0]
        if record[0] != input_chr:
            conflicting.add(snp)
        vals = tuple(record[4:])
        if snp in found:
            duplicated[snp] = duplicated.get(snp,1)+1
            if found[snp]!=vals:
                conflicting.add(snp)
        else:
            found[snp]=vals
        writer.writerow(record)
with (outdir/'fitting_annotation_names.txt').open('w') as f:
    f.write('\n'.join(fitting)+'\n')
with (root/'qc/baseline_coverage.json').open('w') as f:
    json.dump(dict(source=str(baseline),match_key='rsID; chromosome consistency checked; no cross-build BP matching',
        input_unique_rsids=len(variants),matched_unique_rsids=len(found),missing_unique_rsids=len(set(variants)-set(found)),
        duplicated_matched_rsids=len(duplicated),conflicting_matched_rsids=len(conflicting),
        fitting_annotation_count=186,excluded_annotation='GTEx_eQTL_MaxCPP_common'),f,indent=2)
with (root/'qc/baseline_missing_rsids.txt').open('w') as f:
    f.write('\n'.join(sorted(set(variants)-set(found)))+'\n')
with (root/'qc/baseline_conflicting_rsids.txt').open('w') as f:
    f.write('\n'.join(sorted(conflicting))+'\n')
with (root/'qc/baseline_coverage.json').open() as f:
    print(f.read())
