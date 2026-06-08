import sys
import re

if len(sys.argv) != 3:
    sys.exit("Usage: python make_primer3_input.py candidate_6_indel.tsv candidate_6_flank_1001bp.fa > primer3_input.txt")

tsv = sys.argv[1]
fasta = sys.argv[2]

candidates = []
with open(tsv) as f:
    header = f.readline().strip().split("\t")
    for line in f:
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        row = dict(zip(header, parts))
        row["Pos"] = int(row["Pos"])
        row["Delta_length"] = int(row["Delta_length"])
        candidates.append(row)

def read_fasta(path):
    name = None
    seqs = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seqs).upper()
                name = line[1:].split()[0]
                seqs = []
            else:
                seqs.append(line)
        if name is not None:
            yield name, "".join(seqs).upper()

for name, seq in read_fasta(fasta):
    # samtools faidx header is usually like Chr08:34568896-34569896
    m = re.match(r"([^:]+):(\d+)-(\d+)", name)
    if not m:
        sys.stderr.write(f"Skip header without region format: {name}\n")
        continue

    chrom = m.group(1)
    start = int(m.group(2))
    end = int(m.group(3))

    matched = [x for x in candidates if x["Chr"] == chrom and start <= x["Pos"] <= end]
    if len(matched) != 1:
        sys.stderr.write(f"Cannot uniquely match {name}, matched {len(matched)} candidates\n")
        continue

    c = matched[0]
    target_start = c["Pos"] - start   # Primer3 uses 0-based position
    target_len = 1

    # Avoid placing primers directly on the InDel site or very close to it
    exclude_start = max(0, target_start - 20)
    exclude_end = min(len(seq), target_start + 21)
    exclude_len = exclude_end - exclude_start

    seq_id = f'{c["Marker"]}_{chrom}_{c["Pos"]}_delta{c["Delta_length"]}'

    print(f"SEQUENCE_ID={seq_id}")
    print(f"SEQUENCE_TEMPLATE={seq}")
    print(f"SEQUENCE_TARGET={target_start},{target_len}")
    print(f"SEQUENCE_EXCLUDED_REGION={exclude_start},{exclude_len}")
    print("PRIMER_TASK=generic")
    print("PRIMER_PICK_LEFT_PRIMER=1")
    print("PRIMER_PICK_RIGHT_PRIMER=1")
    print("PRIMER_PICK_INTERNAL_OLIGO=0")
    print("PRIMER_NUM_RETURN=10")
    print("PRIMER_PRODUCT_SIZE_RANGE=180-350")
    print("PRIMER_OPT_SIZE=22")
    print("PRIMER_MIN_SIZE=18")
    print("PRIMER_MAX_SIZE=25")
    print("PRIMER_OPT_TM=60.0")
    print("PRIMER_MIN_TM=58.0")
    print("PRIMER_MAX_TM=62.0")
    print("PRIMER_MIN_GC=40.0")
    print("PRIMER_MAX_GC=60.0")
    print("PRIMER_MAX_POLY_X=4")
    print("PRIMER_MAX_NS_ACCEPTED=0")
    print("PRIMER_EXPLAIN_FLAG=1")
    print("=")
