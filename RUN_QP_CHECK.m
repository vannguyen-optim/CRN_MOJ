function result=RUN_QP_CHECK()
% Fast native-backend check, including the exact check that failed previously.
root=CRN_SETUP();cfg=CRN_CONFIG_PUBLICATION();o=cfg.options;
backend=np_native_optim();
fprintf('quadprog: %s\n',backend.quadprog);
fprintf('linprog: %s\n',backend.linprog);
np_work_budget('start',Inf);
[y,solved,meta]=np_qp(eye(2),[-2;-3],[],[],[1 1],1,zeros(2,1),ones(2,1),[.5;.5],o);
point_error=Inf;if isequal(size(y),[2 1]),point_error=norm(y-[0;1]);end
result=struct('backend',backend,'y',y,'solved',solved,'point_error',point_error,'meta',meta);
folder=fullfile(root,'results_checks');if ~exist(folder,'dir'),mkdir(folder);end
save(fullfile(folder,'qp_check.mat'),'result','-v7');
fprintf('Exit flag: %g; y=%s; point error=%.3e.\n',meta.flag,mat2str(y,16),point_error);
if ~solved||point_error>=1e-7
 error('CRN:QPCheck','Native QP check failed. %s\n%s',meta.identifier,meta.message);
end
fprintf('PASS: Native MATLAB QP and bounds.\n');
end
