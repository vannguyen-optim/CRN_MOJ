# Experiment 3: WDBC

Chạy `RUN_WDBC` từ thư mục chính. Cấu hình: `CONFIG_WDBC.m`.
Dữ liệu gốc có sẵn ở `data/wdbc.data`, `data/wdbc.names`.

- 569 observations,30 features; malignant=+1,benign=-1; identifier bị loại.
- 30 stratified paired splits: 341 train,114 validation,114 test.
- Class counts: train(127,214),validation(42,72),test(43,71).
- Seed 20271005+s, s=0,...,29; offsets 10000 và 20000.
- Mean/sample standard deviation lấy từ train; bias không phạt.
- Mỗi split: 18 fractional,12 matched convex,68 tuned convex runs.
- Tổng 2940 runs; bảng 35 và 34; tuning và matched-comparison CSV riêng.
