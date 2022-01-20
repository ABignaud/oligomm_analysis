# OligoMM12 Analysis

Analysis notebook for the characterization of the oligoMM 12 gut metapopulation.

This repository described the analysis from this paper:
[In vitro and in vivo characterization of bacterial and prophages 3D organizations in the OMM12 consortium.](DOI) Q. Lamy-Besnier, A. Bignaud, J. Garneau, M. Titecat, A. von Strempel, B. Stecher-Letsch, R. Koszul, L. Debarbieux, M. Marbouty, Journal, 2022.

The analysis from the paper are described in the notebook repository. Useful functions and a command line interface are available in the script repository.

## Installation

A quick installation to launch the comman line tools with most of the dependencies installed can be done with the conda environment file. The python package needs to be installed as editable to access correctly the bash scripts.

```bash
git clone https://github.com/ABignaud/oligomm_analysis.git
cd oligomm_analysis
conda env create --file oligomm.yml
conda activate oligomm
pip install -e .
```

If you want to run the annotation workflow, you have to make sure that the following software have been correctly installed (especially the database for the annotation pipelines), otherwise you can launch it with the `--skip-annotation` parameter to avoid ro install them:

- dnaglider v0.0.5 (<https://github.com/cmdoret/dnaglider>)
- prokka v1.14.6 (<https://github.com/tseemann/prokka>)
- seqkit v2.1.0 (<https://bioinf.shenwei.me/seqkit/>)
- checkV v0.8.6 (<https://bitbucket.org/berkeleylab/checkv/src/master/>)
- virsorter v2.2.3 (<https://github.com/jiarong/VirSorter2>)
- VIBRANT v3.6 (<https://github.com/AnantharamanLab/VIBRANT>)

## Quick start

To launch the command line function you need to be in a directory named `genus_species` with the folders of the HiC matrix in `graal` format computed by hicstuff inside the repository. Then, the oligomm pipeline can be launch with the following command:

```bash
oligomm -o 0,1,,,2,3,,4,5 -c 16 -b 5000 -B 10000 -i 10 \
    -l in_vitro_R1,in_vitro_R2,in_vivo_2019_R1,in_vivo_2019_R2,in_vivo_2020_R1,in_vivo_2020_R2 \
    MM18 MM81 MM1 M11 MM10 MM4
```

For more detailed explanation of the analysis done by the pipeline you can look at the following [notebook](notebook/OMM12_01_Analysis_of_one bacteria.ypnb).
