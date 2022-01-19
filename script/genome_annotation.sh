# !/bin/bash
set -eu

# Help message.
echo ""
echo "Script to run all annotations on a genome. Output will be in the current"
echo "directory."
echo ""
echo "Usage:"
echo "gene_annotation.sh genus species prefix fasta matrix threads"
echo ""


# Parse variables.
genus=$1
species=$2
prefix=$3
fasta=$4
threads=$5
annotation_pipeline=prokka

# Create output directory in the current directory.

mkdir -p annotation
mkdir -p $annotation_pipeline
mkdir -p checkv
mkdir -p virsorter2
mkdir -p vibrant
mkdir -p tmp

# Compute GC content and GC skew
echo "Compute GC content and GC skew with dnaglider..."
dnaglider -fasta $fasta -stride 1000 -window 5000 -threads $threads \
    -out annotation/"$prefix"_gc.bed
dnaglider -fasta $fasta -stride 1000 -window 5000 -threads $threads \
    -out annotation/"$prefix"_gcskew.bed -fields GCSKEW

# Run annotation pipeline.
if [ $annotation_pipeline = prokka ]
then
    # Run prokka annotation pipeline
    echo "Run prokka anotation pipeline..."
    prokka --outdir ./prokka --prefix $prefix --addgenes --genus $genus \
        --species $species --usegenus --force --cpus $threads $fasta

# TODO: Correct it to make it work.
elif [ $annotation_pipeline = pgap ]
then
    # Run pgap annotation pipeline
    echo "Run pgap anotation pipeline..."
    # TODO: write input.yaml and others files.
    wdir=$(pwd)
    cd /data/pgap/
    /home/rsg/repo/pgap/pgap.py -o $wdir/pgap -n --cpus $threads \
        --ignore-all-errors $wdir/tmp/input.yaml
    cd $wdir
fi

# Extract features from gff
echo "Extract features from gff..."

# Extract CDS from gff to bed
grep -i "CDS" $annotation_pipeline/*gff | sed 's/product=/\t/' | \
    sed 's/;protein_id=/\t/' | \
    awk -v FS="\t" -v OFS="\t" '{print $1,$4,$5,$10,$6,$7}' \
    > annotation/"$prefix"_cds.bed

# Extract tRNA
awk -v OFS="\t" '$3=="tRNA" {print $1,$4,$5,$3}' $annotation_pipeline/*.gff \
    > annotation/"$prefix"_tRNA.bed

# Extract rRNA
awk -v OFS="\t" '$3=="rRNA" {print $1,$4,$5,$3}' $annotation_pipeline/*.gff \
    > annotation/"$prefix"_rRNA.bed

# Find structural sequences sites
script_path=$(dirname "$0")
echo "Find parS sites..."
python3 $script_path/pars_seqfinder.py $fasta annotation/"$prefix"_pars.bed
echo "Find matS sites..."
python3 $script_path/mats_seqfinder.py $fasta annotation/"$prefix"_mats.bed

# Detect provirus sequences with checkV.
echo "Detect provirus sequences with checkV"
# Split fasta in 250kb fragments as checkV detect only one provirus per 
# fragment.
seqkit sliding -c -s 250000 -W 250000 $fasta -o tmp/"$prefix"_split.fa
# Run checkV
checkv end_to_end -t $threads tmp/"$prefix"_split.fa checkv/
# Extract provirus positions from the output.
grep "Yes" checkv/contamination.tsv | sed 's/[:,-]/\t/g' | \
    awk '{if ($11 == "host" && $12 == "viral" && $13 == "host") \
            {print $1"\t"$2+$19"\t"$2+$20} \
        else {if ($11 == "host" && $12 == "viral") \
            {print $1"\t"$2+$17"\t"$2+$18} \
            else {if ($11 == "viral" && $12 == "host") \
                {print $1"\t"$2+$15"\t"$2+$16}}}}' | \
    sed 's/_sliding//g' > annotation/"$prefix"_checkv.bed

# Run virsorter2
echo "Launch virsorter2"
virsorter run -w virsorter2 -i $fasta --min-length 5000 \
    --min-score 0.5 -j $threads --rm-tmpdir all --keep-original-seq
cut -f1,4,5,27,28 virsorter2/final-viral-boundary.tsv |
    awk 'NR != 1 {print $0}' > annotation/"$prefix"_virsorter2.bed

# Run VIBRANT
echo "Launch VIBRANT"
VIBRANT_run.py -i $fasta -t $threads -folder vibrant
name=$(basename $fasta .fa)
cut -f1,6,7 vibrant/VIBRANT_$name/VIBRANT_results_$name/VIBRANT_integrated_prophage_coordinates_$name.tsv | \
    awk 'NR != 1 {print $0"\tvibrant"}' > annotation/"$prefix"_vibrant.bed

# Delete temporary files
rm -r tmp