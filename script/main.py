#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import bacchus.antidiagonal as bca
import bacchus.hic as bch
import bacchus.io as bcio
import bacchus.oriter as bco
import bacchus.plot as bcp
import click
import cooler
import hicstuff.distance_law as hcdl
import hicstuff.hicstuff as hcs
import hicstuff.io as hio
import matplotlib.pyplot as plt
import script.plot_matrices as cma
import os
import numpy as np
import pandas as pd
import subprocess
from tinycov.tinycov import covplot
from os.path import isfile, join


@click.command()
@click.argument("hicstuff_folder", nargs=-1, type=click.Path())
@click.option(
    "-A",
    "--skip_annotation",
    is_flag=True,
    type=bool,
    help="If enable, skip annotation workflow.",
)
@click.option(
    "-b", "--bin_size", default=2000, type=int, help="Binning size of the map."
)
@click.option(
    "-B",
    "--bin_size_ratio",
    default=0,
    type=int,
    help="Binning size of the log ratio.",
)
@click.option("-c", "--cpus", default=16, type=int, help="Threads to used.")
@click.option(
    "-C",
    "--skip_coverage",
    is_flag=True,
    type=bool,
    help="If enable, skip coverage workflow.",
)
@click.option(
    "-D",
    "--skip_distance_law",
    is_flag=True,
    type=bool,
    help="If enable, skip distance law workflow.",
)
@click.option(
    "-H",
    "--skip_hicrep",
    is_flag=True,
    type=bool,
    help="If enable, skip hicrep workflow.",
)
@click.option(
    "-i",
    "--iteration",
    default=10,
    type=int,
    help="Number of iteration for serpentine.",
)
@click.option(
    "-l",
    "--labels",
    default="Blue",
    type=str,
    help="Label to use to save the maps separated by commas if more than one.",
)
@click.option(
    "-M",
    "--skip_contact_map",
    is_flag=True,
    type=bool,
    help="If enable, skip contact map workflow.",
)
@click.option(
    "-o",
    "--order",
    default="0",
    type=str,
    help="Order of matrix for ratio. 0 based values.",
)
@click.option(
    "-P",
    "--skip_pygenometrack",
    is_flag=True,
    type=bool,
    help="If enable, skip pygenometrack workflow.",
)
@click.option(
    "-S",
    "--skip_antidiagonal",
    is_flag=True,
    type=bool,
    help="If enable, skip antidiagonal workflow.",
)
@click.option(
    "-t",
    "--threshold",
    default="auto",
    type=str,
    help="Threshold for serpentine (integer). If auto is given it will guess a threshold.",
)
def main(
    hicstuff_folder,
    bin_size,
    bin_size_ratio,
    cpus,
    labels,
    order,
    iteration,
    skip_annotation,
    skip_antidiagonal,
    skip_contact_map,
    skip_coverage,
    skip_distance_law,
    skip_hicrep,
    skip_pygenometrack,
    threshold,
):
    # Defined matrices files and fragment_file and distance_law files
    fragment_file = join(hicstuff_folder[0], "fragments_list.txt")
    matrice_files = []
    distance_law_files = []
    bam_files = []
    for folder in hicstuff_folder:
        matrice_files.append(
            join(folder, "abs_fragments_contacts_weighted.txt")
        )
        distance_law_files.append(join(folder, "distance_law.txt"))
        if os.path.isfile(join(folder, "tmp", "for.sorted.bam")):
            bam_files.append(join(folder, "tmp", "for.sorted.bam"))
        else:
            bam_files.append(join(folder, "tmp", "for.bam"))

    # Build a list with the order of matrices to compare.
    order = order.split(",,,")
    for i, v in enumerate(order):
        order[i] = v.split(",,")
        for j, v in enumerate(order[i]):
            order[i][j] = list(map(int, v.split(",")))

    labels = labels.split(",")
    genus, species = os.path.basename(os.getcwd()).split("_")
    prefix = genus[0] + "_" + species

    # Plot the contact maps.
    os.makedirs("contact_map", exist_ok=True)

    # Defined which part of the pipeline will be done.
    contact_map = not skip_contact_map
    coverage = not skip_coverage
    distance_law = not skip_distance_law
    hicrep = not skip_hicrep
    annotation = not skip_annotation
    pygenometrack = not skip_pygenometrack
    antidiagonal = not skip_antidiagonal

    if coverage:
        os.makedirs("cov", exist_ok=True)
        # Compute coverage with tinycov for all the maps.
        cov_files = []
        for i, bam_file in enumerate(bam_files):
            label = labels[i]
            output_fig = join("cov", f"{prefix}_{label}_cov.pdf")
            output_file = join("cov", f"{prefix}_{label}_cov.bed")
            title = prefix + "_" + label
            covplot(
                bam=bam_file,
                out=output_fig,
                res=5000,
                skip=500,
                name=title,
                text=output_file,
            )
            plt.close()
            cov_files.append(output_file)

    if annotation:
        # Annotate the genome and create the pygenometrack file.
        fasta = f"{genus[0].upper()}_{species}.fa"
        script_path = join(
            os.path.dirname(os.path.abspath(__file__)), "genome_annotation.sh"
        )
        cmd_annotation = (
            f"{script_path} {genus} {species} {prefix} {fasta} {cpus}"
        )
        annotation_process = subprocess.Popen(cmd_annotation, shell=True)
        out, err = annotation_process.communicate()

    if pygenometrack:
        k = 1
        oriter_file = join("annotation", f"{prefix}_oriter.bed")
        ratio_file = join("cov", f"{prefix}_ratio_ori_ter.txt")
        # Compute the ori ter position if not already available.
        if not coverage:
            cov_files = []
            for i, bam_file in enumerate(bam_files):
                label = labels[i]
                output_file = join("cov", f"{prefix}_{label}_cov.bed")
                cov_files.append(output_file)

        # Extract data:
        gc_skew_data = pd.read_csv(
            join("annotation", f"{prefix}_gcskew.bed"),
            sep="\t",
            header=0,
            names=["chr", "start", "end", "val"],
        )
        pars_data = pd.read_csv(
            join("annotation", f"{prefix}_pars.bed"),
            sep="\t",
            header=None,
            names=["chr", "start", "end", "type", "seq", "err"],
        )

        ori, ter = bco.detect_ori_ter(gc_skew_data, pars_data)

        with open(oriter_file, "w") as out:
            chrom = f"{genus[0].upper()}_{species}"
            out.write("%s\t%s\t%s\t%s\n" % (chrom, ori - 1, ori + 1, "ori"))
            out.write("%s\t%s\t%s\t%s\n" % (chrom, ter - 1, ter + 1, "ter"))

        with open(ratio_file, "w") as out:
            for i, cov_file in enumerate(cov_files):
                label = labels[i]
                cov_data = pd.read_csv(
                    cov_file,
                    sep="\t",
                    header=None,
                    names=["chr", "start", "end", "val"],
                )
                ratio, ori_cov, ter_cov = bco.compute_oriter_ratio(
                    cov_data, ori, ter
                )
                out.write("%s\t%s\t%s\t%s\n" % (label, ori, ter, ratio))
                if i == k:
                    ori_plot = ori_cov * 1.5
                    ter_plot = ter_cov * 0.5

        # Compute the matrix in cool format and balanced it.
        label_bac = labels[k]
        cool_out = join("cool", f"{label_bac}_{bin_size // 1000}kb.cool")
        # Check if cool file already exist.
        M, frags = bcio.binned_map(matrice_files[k], fragment_file, bin_size)
        if not isfile(cool_out):
            os.makedirs("cool", exist_ok=True)
            hio.save_cool(
                cool_out,
                M,
                frags,
                metadata={"hicstuff": "3.1.0", "bin-type": "fixed"},
            )
        # Balance matrix
        cool_matrix = cooler.Cooler(cool_out)
        cooler.balance.balance_cooler(
            cool_matrix, mad_max=50, min_nnz=0, store=True
        )

        # Look for vmax
        M = hcs.normalize_sparse(M, norm="ICE", n_mad=10)
        M = bch.get_symmetric(M)
        M = M.toarray()
        max_cmap = np.nanpercentile(M, 99)
        script_path = join(
            os.path.dirname(os.path.abspath(__file__)), "pygenometracks.sh"
        )
        cmd_pgt = f"{script_path} {prefix} {os.path.basename(cool_out)} {ori_plot} {ter_plot} {max_cmap} {label_bac}"
        pgt_process = subprocess.Popen(cmd_pgt, shell=True)
        out, err = pgt_process.communicate()

    if contact_map:
        # Plot the matrices and their log ratio according to the order.
        if len(order) == 2:
            os.makedirs("log_ratio", exist_ok=True)
            # flatten the list
            l1 = sum(order[0], [])
            l2 = sum(order[1], [])
            matrice_files_1 = [matrice_files[a] for a in l1]
            matrice_files_2 = [matrice_files[a] for a in l2]
            label = [
                "-".join([labels[a] for a in l1]),
                "-".join([labels[a] for a in l2]),
            ]
            print(
                "Plot the matrices and their log ratio: {0} - {1}".format(
                    label[0], label[1]
                )
            )
            cma.plot_matrices(
                matrice_files_1,
                matrice_files_2,
                bin_size,
                bin_size_ratio,
                cpus,
                fragment_file,
                label,
                iteration,
                threshold,
            )
            if len(order[1]) >= 2:
                matrice_files_1 = [matrice_files[a] for a in order[1][0]]
                matrice_files_2 = [matrice_files[a] for a in order[1][1]]
                label = [
                    "-".join([labels[a] for a in order[1][0]]),
                    "-".join([labels[a] for a in order[1][1]]),
                ]
                print(
                    "Plot the matrices and their log ratio: {0} - {1}".format(
                        label[0], label[1]
                    )
                )
                cma.plot_matrices(
                    matrice_files_1,
                    matrice_files_2,
                    bin_size,
                    bin_size_ratio,
                    cpus,
                    fragment_file,
                    label,
                    iteration,
                    threshold,
                )
                if len(order[1][0]) == 2:
                    matrice_files_1 = [matrice_files[order[1][0][0]]]
                    matrice_files_2 = [matrice_files[order[1][0][1]]]
                    label = [
                        "-".join([labels[order[1][0][0]]]),
                        "-".join([labels[order[1][0][1]]]),
                    ]
                    print(
                        "Plot the matrices and their log ratio: {0} - {1}".format(
                            label[0], label[1]
                        )
                    )
                    cma.plot_matrices(
                        matrice_files_1,
                        matrice_files_2,
                        bin_size,
                        bin_size_ratio,
                        cpus,
                        fragment_file,
                        label,
                        iteration,
                        threshold,
                    )
                if len(order[1][1]) == 2:
                    matrice_files_1 = [matrice_files[order[1][1][0]]]
                    matrice_files_2 = [matrice_files[order[1][1][1]]]
                    label = [
                        "-".join([labels[order[1][1][0]]]),
                        "-".join([labels[order[1][1][1]]]),
                    ]
                    print(
                        "Plot the matrices and their log ratio: {0} - {1}".format(
                            label[0], label[1]
                        )
                    )
                    cma.plot_matrices(
                        matrice_files_1,
                        matrice_files_2,
                        bin_size,
                        bin_size_ratio,
                        cpus,
                        fragment_file,
                        label,
                        iteration,
                        threshold,
                    )
                if len(order[1]) == 3:
                    matrice_files_1 = [matrice_files[a] for a in order[1][0]]
                    matrice_files_2 = [matrice_files[a] for a in order[1][2]]
                    label = [
                        "-".join([labels[a] for a in order[1][0]]),
                        "-".join([labels[a] for a in order[1][2]]),
                    ]
                    print(
                        "Plot the matrices and their log ratio: {0} - {1}".format(
                            label[0], label[1]
                        )
                    )
                    cma.plot_matrices(
                        matrice_files_1,
                        matrice_files_2,
                        bin_size,
                        bin_size_ratio,
                        cpus,
                        fragment_file,
                        label,
                        iteration,
                        threshold,
                    )
                    matrice_files_1 = [matrice_files[a] for a in order[1][1]]
                    matrice_files_2 = [matrice_files[a] for a in order[1][2]]
                    label = [
                        "-".join([labels[a] for a in order[1][1]]),
                        "-".join([labels[a] for a in order[1][2]]),
                    ]
                    print(
                        "Plot the matrices and their log ratio: {0} - {1}".format(
                            label[0], label[1]
                        )
                    )
                    cma.plot_matrices(
                        matrice_files_1,
                        matrice_files_2,
                        bin_size,
                        bin_size_ratio,
                        cpus,
                        fragment_file,
                        label,
                        iteration,
                        threshold,
                    )
                    if len(order[1][2]) == 2:
                        matrice_files_1 = [matrice_files[order[1][2][0]]]
                        matrice_files_2 = [matrice_files[order[1][2][1]]]
                        label = [
                            "-".join([labels[order[1][2][0]]]),
                            "-".join([labels[order[1][2][1]]]),
                        ]
                        print(
                            "Plot the matrices and their log ratio: {0} - {1}".format(
                                label[0], label[1]
                            )
                        )
                        cma.plot_matrices(
                            matrice_files_1,
                            matrice_files_2,
                            bin_size,
                            bin_size_ratio,
                            cpus,
                            fragment_file,
                            label,
                            iteration,
                            threshold,
                        )
        if len(order[0]) == 2:
            os.makedirs("log_ratio", exist_ok=True)
            matrice_files_1 = [matrice_files[a] for a in order[0][0]]
            matrice_files_2 = [matrice_files[a] for a in order[0][1]]
            label = [
                "-".join([labels[a] for a in order[0][0]]),
                "-".join([labels[a] for a in order[0][1]]),
            ]
            print(
                "Plot the matrices and their log ratio: {0} - {1}".format(
                    label[0], label[1]
                )
            )
            cma.plot_matrices(
                matrice_files_1,
                matrice_files_2,
                bin_size,
                bin_size_ratio,
                cpus,
                fragment_file,
                label,
                iteration,
                threshold,
            )
            if len(order[0][1]) == 2:
                matrice_files_1 = [matrice_files[order[0][1][0]]]
                matrice_files_2 = [matrice_files[order[0][1][1]]]
                label = [
                    "-".join([labels[order[0][1][0]]]),
                    "-".join([labels[order[0][1][1]]]),
                ]
                print(
                    "Plot the matrices and their log ratio: {0} - {1}".format(
                        label[0], label[1]
                    )
                )
                cma.plot_matrices(
                    matrice_files_1,
                    matrice_files_2,
                    bin_size,
                    bin_size_ratio,
                    cpus,
                    fragment_file,
                    label,
                    iteration,
                    threshold,
                )
        if len(order[0][0]) == 2:
            os.makedirs("log_ratio", exist_ok=True)
            matrice_files_1 = [matrice_files[order[0][0][0]]]
            matrice_files_2 = [matrice_files[order[0][0][1]]]
            label = [
                "-".join([labels[order[0][0][0]]]),
                "-".join([labels[order[0][0][1]]]),
            ]
            print(
                "Plot the matrices and their log ratio: {0} - {1}".format(
                    label[0], label[1]
                )
            )
            cma.plot_matrices(
                matrice_files_1,
                matrice_files_2,
                bin_size,
                bin_size_ratio,
                cpus,
                fragment_file,
                label,
                iteration,
                threshold,
            )
        # Case where there is only one map.
        if len(order) == 1 and len(order[0]) == 1 and len(order[0][0]) == 1:
            print("Plot the matrice: {0}".format(labels[0]))
            M1 = bcio.build_map(matrice_files, fragment_file, bin_size)
            bcp.contact_map(
                M1,
                binning=bin_size,
                title=labels[0],
                out_file=join(
                    "contact_map",
                    f"{labels[0]}_{bin_size / 1000}kb.pdf",
                ),
            )
            distance_law = False
            hicrep = False

    if distance_law:
        # Compute the distance law.
        print("Compute distance law...")
        # Defined output:
        output_file_img = "distance_law.pdf"
        output_file_table = "distance_law.txt"

        # Make new lists for the modified distance law.
        length_files = len(distance_law_files)
        xs = [None] * length_files
        ps = [None] * length_files
        names = [None] * length_files
        for i in range(length_files):
            xs[i], ps[i], names[i] = hcdl.import_distance_law(
                distance_law_files[i]
            )
        names = [name[0] for name in names]
        inf = 3000
        sup = max(max(xs[0], key=len))
        # Iterate on the different file given by the user.
        for i in range(length_files):
            # Make the average
            xs[i], ps[i] = hcdl.average_distance_law(xs[i], ps[i], sup, False)
        # Normalize and make the derivative
        ps = hcdl.normalize_distance_law(xs, ps, inf, sup)
        hcdl.plot_ps_slope(xs, ps, labels, output_file_img, inf, sup)
        hcdl.export_distance_law(xs, ps, labels, output_file_table)

    if hicrep:
        # Compute HicRep
        print("Compute HiCRep...")

        # Save the matrices as cool format:
        os.makedirs("cool", exist_ok=True)
        cool_files = []
        binning = 10_000  # Large binning for HiCreppy
        # Import matrix
        for i, matrix_file in enumerate(matrice_files):
            cool_out = join("cool", f"{labels[i]}_{binning // 1000}kb.cool")
            # Create cool file if doesn't exist.
            if not isfile(cool_out):
                M, frags = binned_map(matrix_file, fragment_file, 10000)
                matrix.append(M.tocoo())
                hio.save_cool(
                    cool_out,
                    M,
                    frags,
                    metadata={"hicstuff": "3.1.0", "bin-type": "fixed"},
                )
            cool_files.append(cool_out)

        corr_matrix = bch.get_hicreppy(cool_files, subsample=0, h=None)

        out_file = "hicreppy.pdf"
        bcp.hicreppy_plot(corr_matrix, labels, out_file)

    if antidiagonal:
        # Compute antidiagonal
        print("Compute antidiagonal...")

        # Search parS
        pars_file = join("annotation", f"{prefix}_pars.bed")
        if isfile(pars_file):
            pars = np.array(
                pd.read_csv(pars_file, sep="\t", header=None).iloc[:, 1]
            )
        else:
            pars = None

        # Search for ori
        ori_file = join("cov", f"{prefix}_ratio_ori_ter.txt")
        if isfile(ori_file):
            ori = pd.read_csv(ori_file, sep="\t", header=None).iloc[1, 1]
        else:
            ori = None

        # Run antidiagonal analysis.
        all_values = []
        for i, matrice_file in enumerate(matrice_files):
            M = bcio.build_map([matrice_file], fragment_file, bin_size)
            all_values.append(bca.compute_antidiagonal(M, full=False))

        # Display the antidiagonal plot
        out_file = f"{prefix}_antidiagonal.pdf"
        bcp.antidiagonal_scalogram(
            all_values,
            binning=bin_size,
            title=f"{genus[0].upper()} {species}",
            ori_pos=ori,
            pars_pos=pars,
            labels=labels,
            out_file=out_file,
        )


if __name__ == "__main__":
    main()
