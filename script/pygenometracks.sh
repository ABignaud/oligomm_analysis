# !/bin/bash
set -eu

# Defined variables
prefix=$1
matrix_file=$2
ori_value=$3
ter_value=$4
max_cmap=$5
label_bac=$6
script_path=$(dirname "$0")

# Create a pyGenome Tracks initiator file from the template
echo "Display plot with pyGenometracks"
chr=$(tail -n 1 annotation/"$prefix"_gc.bed | cut -f1)
end=$(tail -n 1 annotation/"$prefix"_gc.bed | cut -f3)
sed 's/prefix/'$prefix'/' $script_path/../utils/tracks_template_wo_vline.ini | \
    sed 's/matrix_file/..\/cool\/'$matrix_file'/' | \
    sed 's/chr_length/'$end'/'  | \
    sed 's/ori_value/'$ori_value'/' | \
    sed 's/ter_value/'$ter_value'/' | \
    sed 's/max_cmap/'$max_cmap'/' | \
    sed 's/label_bac/'$label_bac'/' \
    > annotation/"$prefix"_tracks.ini
pyGenomeTracks --region $chr:1-$end --tracks annotation/"$prefix"_tracks.ini \
    --out "$prefix"_wo_line.pdf
    