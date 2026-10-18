# Experiment 2: hồi quy thưa

Chạy `RUN_EXPERIMENT2` từ thư mục chính. Cấu hình: `CONFIG_EXPERIMENT2.m`.

- d=100,200,500; N1=200,N2=100; 2000 evaluation observations mỗi task.
- Support mỗi vector có 10 entries; giao 0,5,10; noise standard deviation 0.1.
- Free, box [-5,5]^d, affine Cz=0 (6 rows, rank 5, một duplicated row).
- Gamma=0.001,0.01,0.05; 30 paired realizations mỗi cấu hình.
- Free: CRN,CRN0,IReNA,HpgSRN; box/affine: CRN,CRN0,SSQP.
- Tổng 8100 runs; bảng 28 và hai hình Experiment 2.
