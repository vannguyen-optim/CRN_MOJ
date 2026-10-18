function U=np_null(C,d)
if isempty(C),U=eye(d);else,U=null(full(C));end
end
