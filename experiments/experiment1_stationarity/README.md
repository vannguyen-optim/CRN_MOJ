# Experiment 1: stationarity và nghiệm tham chiếu

Chạy `RUN_EXPERIMENT1` từ thư mục chính. Cấu hình: `CONFIG_EXPERIMENT1.m`.

- Miền [-2,2]^2; 171 trọng số dương trên lưới mẫu số 20.
- Hai điểm đầu (1.5,-1.5), (-1.5,1.5) và zero-start diagnostic riêng.
- CRN, CRN0, SSQP: tổng 1539 runs.
- Nghiệm tham chiếu: xét zero, endpoints và tất cả positive real roots
  của đa thức theo từng tọa độ; lưu root/value tolerances.
- Xuất success, iterations, time và reference gaps; bảng 27, hai hình Experiment 1.
