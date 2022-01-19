#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import bacchus.hic as bch
import bacchus.io as bcio
import bacchus.plot as bcp
import hicstuff.hicstuff as hcs
import numpy as np
import serpentine as srp
from os.path import join


def log_ratio_multiple(
    matrix_files_1,
    matrix_files_2,
    fragment_file,
    bin_size,
    iterations=10,
    threshold=10,
    cpus=16,
):
    """Function to compute the log ratio of two matrices and use serpentine to
    return a smoother log ratio.
    """
    for i, matrix_file_1 in enumerate(matrix_files_1):
        M, frags = bcio.binned_map(matrix_file_1, fragment_file, bin_size)
        if i == 0:
            M1 = M
        else:
            M1 += M
    for i, matrix_file_2 in enumerate(matrix_files_2):
        M, _ = bcio.binned_map(matrix_file_2, fragment_file, bin_size)
        if i == 0:
            M2 = M
        else:
            M2 += M
    M1 = M1.tocoo()
    M2 = M2.tocoo()
    m1 = M1.sum()
    m2 = M2.sum()
    if m1 < m2:
        M2 = hcs.subsample_contacts(M2, int(m1))
    else:
        M1 = hcs.subsample_contacts(M1, int(m2))
    subsample = min(m1, m2)
    print("Matrices have been subsampled to {0} contacts.".format(subsample))
    # M1 = M1.tocsr() - M1.tocsr().multiply(
    #     get_win_density(M1, win_size=3, sym_upper=True) < 2 / 9
    # )
    # M2 = M2.tocsr() - M2.tocsr().multiply(
    #     get_win_density(M2, win_size=3, sym_upper=True) < 2 / 9
    # )
    M1 = bch.get_symmetric(M1).toarray()
    M2 = bch.get_symmetric(M2).toarray()
    if threshold == "auto":
        trend, threshold = srp.MDbefore(M1, M2, show=False)
    else:
        threshold = int(threshold)
    M1_serp, M2_serp, ratio_srp = srp.serpentin_binning(
        M1,
        M2,
        parallel=cpus,
        triangular=False,
        threshold=threshold,
        minthreshold=threshold / 5,
        iterations=iterations,
        verbose=False,
    )
    ratio_log = np.log2(M2) - np.log2(M1)
    return M1_serp, M2_serp, ratio_log, ratio_srp, subsample


def plot_matrices(
    matrice_files_1,
    matrice_files_2,
    bin_size,
    bin_size_ratio,
    cpus,
    fragment_file,
    label,
    iteration,
    threshold,
):
    # Defined bin_size for the log ratio if necessary.
    if bin_size_ratio == 0:
        bin_size_ratio = bin_size

    # Compute log ratio and smoothen it with serpentine.
    M1_serp, M2_serp, ratio_log, ratio_srp, subsample = log_ratio_multiple(
        matrice_files_1,
        matrice_files_2,
        fragment_file,
        bin_size_ratio,
        iteration,
        threshold,
        cpus,
    )

    # Plot the log_ratio
    bcp.contact_map_ratio(
        ratio_srp,
        ratio_log,
        binning=bin_size_ratio,
        ratio=True,
        out_file=join(
            "log_ratio",
            "{0}_{1}_{2}kb.pdf".format(
                label[0], label[1], bin_size_ratio / 1000
            ),
        ),
        title=f"{label[0]} - {label[1]}",
    )

    # Plot the maps with subsampling according to the smallest matrix.
    M1 = bcio.build_map(
        matrice_files_1, fragment_file, bin_size, subsample=subsample
    )
    bcp.contact_map(
        M1,
        binning=bin_size,
        title=label[0],
        out_file=join(
            "contact_map",
            "{0}_{1}kb_subsample.pdf".format(label[0], bin_size / 1000),
        ),
    )
    M2 = bcio.build_map(
        matrice_files_2, fragment_file, bin_size, subsample=subsample
    )
    bcp.contact_map(
        M2,
        binning=bin_size,
        title=label[1],
        out_file=join(
            "contact_map",
            "{0}_{1}kb_subsample.pdf".format(label[1], bin_size / 1000),
        ),
    )
    M1 = bcio.build_map(matrice_files_1, fragment_file, bin_size)
    bcp.contact_map(
        M1,
        binning=bin_size,
        title=label[0],
        out_file=join(
            "contact_map", "{0}_{1}kb.pdf".format(label[0], bin_size / 1000)
        ),
    )
    M2 = bcio.build_map(matrice_files_2, fragment_file, bin_size)
    bcp.contact_map(
        M2,
        binning=bin_size,
        title=label[1],
        out_file=join(
            "contact_map", "{0}_{1}kb.pdf".format(label[1], bin_size / 1000)
        ),
    )
