# Create directory
refDir=/nfs/tubro/umms-eicarpen/references
mkdir -p $ref

# Download reference file
curl --output ${ref}/refdata-gex-GRCh38-2020-A.tar.gz \
https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz

# Uncompress zipped file
tar -xvf ${ref}/refdata-gex-GRCh38-2020-A.tar.gz \
-C ${ref}/
rm ${ref}/refdata-gex-GRCh38-2020-A.tar.gz