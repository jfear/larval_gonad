"""Common elements used for the scRNA-Seq data.

This is a place to store functions and variables that are used over and over
again with the scRNA-seq data.

"""

CLUSTER_ANNOT = {
    0: 'Late 1º Spermatocytes (0)',
    1: 'Late Cyst Cells (1)',
    2: 'Spermatogonia (2)',
    3: 'Early 1º Spermatocytes (3)',
    4: 'Late Cyst Cells (4)',
    5: 'Late Cyst Cells (5)',
    6: 'Early Cyst Cells (6)',
    7: 'Terminal Epithelium(7)',
    8: 'Pigment (8)',
    9: 'Unknown (9)',
}

CLUSTER_ORDER = [
    'Spermatogonia (2)',
    'Early 1º Spermatocytes (3)',
    'Late 1º Spermatocytes (0)',
    'Early Cyst Cells (6)',
    'Late Cyst Cells (1)',
    'Late Cyst Cells (4)',
    'Late Cyst Cells (5)',
    'Terminal Epithelium(7)',
    'Pigment (8)',
    'Unknown (9)',
]


# CLUSTER_ANNOT = {
#     0: 'Late Primary Spermatocytes (0)',
#     1: 'Early Somatic Cyst Cells (1)',
#     2: 'Late Somatic Cyst Cells (2)',
#     3: 'Late Somatic Cyst Cells (3)',
#     4: 'Spermatogonia (4)',
#     5: 'Terminal Epithelium (5)',
#     6: 'Mid Primary Spermatocytes (6)',
#     7: 'Late Somatic Cyst Cells (7)',
#     8: 'Early Primary Spermatocytes (8)',
#     9: 'Pigment Cells (9)',
#     10: 'Early Somatic Cyst Cells (10)',
#     11: 'Unknown (11)',
# }
#
# CLUSTER_ORDER = [
#     'Spermatogonia (4)',
#     'Early Primary Spermatocytes (8)',
#     'Mid Primary Spermatocytes (6)',
#     'Late Primary Spermatocytes (0)',
#     'Early Somatic Cyst Cells (1)',
#     'Early Somatic Cyst Cells (10)',
#     'Late Somatic Cyst Cells (2)',
#     'Late Somatic Cyst Cells (3)',
#     'Late Somatic Cyst Cells (7)',
#     'Terminal Epithelium (5)',
#     'Pigment Cells (9)',
#     'Unknown (11)',
# ]
