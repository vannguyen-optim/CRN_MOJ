function value=np_work_budget(action,limit)
% Per-worker QP reservation counter; every solver call starts its own run.
% Standalone mathematical probes default to an unbounded counter.
persistent maximum used;
if isempty(maximum),maximum=Inf;used=0;end
switch action
 case 'start'
  assert(limit>=0 && (isinf(limit)||limit==floor(limit)),'Invalid QP solve budget');
  maximum=limit;used=0;value=0;
 case 'reserve'
  value=used<maximum;if value,used=used+1;end
 case 'used',value=used;
 case 'remaining',value=max(0,maximum-used);
 otherwise,error('Unknown work-budget operation: %s',action);
end
end
