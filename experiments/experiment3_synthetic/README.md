# Experiment 3: phân loại tổng hợp

Chạy `RUN_EXPERIMENT3_SYNTHETIC`; `RUN_EXPERIMENT3` chạy cả tổng hợp và WDBC.
Cấu hình: `CONFIG_EXPERIMENT3.m`.

- n=100,200; prevalence 0.3,0.5; 400 train,200 validation,2000 test.
- Gaussian features; 10 reference coefficients bằng 1/sqrt(10); 5% label flips.
- Free và affine (5 independent homogeneous equations + 1 duplicated row).
- CRN,CRN0,SSQP: 6 positive weights mỗi split.
- SVM l1/l2: 6 matched coefficients và 34 independent tuning candidates mỗi split.
- 30 paired realizations; tổng 23520 runs; bảng 29,30 và hai hình Experiment 3.
