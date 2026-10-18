# WDBC raw data

Files `wdbc.data` and `wdbc.names` are the original UCI Wisconsin Diagnostic Breast Cancer files, bundled for offline MATLAB execution.

Source: https://archive.ics.uci.edu/dataset/17/breast+cancer+wisconsin+diagnostic
Dataset DOI: https://doi.org/10.24432/C5DW2B
License listed by UCI: CC BY 4.0, https://creativecommons.org/licenses/by/4.0/

569 rows, 212 malignant (M) and 357 benign (B). Each row contains an identifier, diagnosis and 30 numeric features. Identifier is excluded; M maps to +1 and B to −1.

Features 1–10 are the means, 11–20 the standard errors, and 21–30 the worst values of these ten measurements, in this order:

1. Radius
2. Texture
3. Perimeter
4. Area
5. Smoothness
6. Compactness
7. Concavity
8. Concave points
9. Symmetry
10. Fractal dimension

`np_wdbc.m` reads these files, generates the stratified splits and records original row indices, sample identifiers, training means and training sample standard deviations.
