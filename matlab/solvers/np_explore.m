function [y,u,ok,s]=np_explore(P,x,o,center)
% Finite-budget stabilized IRL1 proposals. Direct acceptance is checked outside.
if nargin<4,center=x;end
v=x;T=P.L+o.eta+1;
s=struct('qp',0,'certificate_qp',0,'steps',0,'certificate_calls',0, ...
 'model_certificate_calls',0,'polishing',{{}});ok=false;u=[];y=x;
for r=1:o.explore_steps
 [~,g]=np_eval(P,v);g=g+o.eta*(v-center);
 assert(strcmp(o.explore_policy,'stable_irl1'),'Use the publication exploration policy');
 w=P.alpha.*P.p.*max(abs(v),o.rw_eps).^(P.p-1);
 [yn,~,modelok,mw,mcw,polish]=np_l1_model(P,v,g,T,w,[],o);
 s.qp=s.qp+mw;s.certificate_qp=s.certificate_qp+mcw;s.steps=r;
 s.model_certificate_calls=s.model_certificate_calls+double(polish.model_certificate_called);
 if ~isempty(polish.raw),s.polishing{end+1}=polish;end
 if ~modelok,return;end
 v=yn;
end
y=v;[u,~,ok,cw]=np_certificate(P,y,[],o);
s.certificate_qp=s.certificate_qp+cw;s.certificate_calls=1;
end
