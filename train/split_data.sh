#!/bin/bash

if [ $# -lt 1 ]; then
    echo "usage: $0 features"
    exit 1
fi
input="$1"
prefix=`basename $input`
prefix=${prefix%.gz}
prefix=${prefix%.tsv}

for group in protein-other protein-snp nonpro-other  nonpro-snp; do
    echo "##$group" > ./features/${prefix}.${group}.tsv
    zcat ${input} | grep "^#" >> ./features/${prefix}.${group}.tsv
done

zcat ${input}|grep -v "^#" | awk '$(NF-1) != "NaN"' | awk '$18 != 0' >> ./features/${prefix}.protein-other.tsv
zcat ${input}|grep -v "^#" | awk '$(NF-1) != "NaN"' | awk '$18 == 0' >> ./features/${prefix}.protein-snp.tsv
zcat ${input}|grep -v "^#" | awk '$(NF-1) == "NaN"' | awk '$18 != 0' >> ./features/${prefix}.nonpro-other.tsv
zcat ${input}|grep -v "^#" | awk '$(NF-1) == "NaN"' | awk '$18 == 0' >> ./features/${prefix}.nonpro-snp.tsv