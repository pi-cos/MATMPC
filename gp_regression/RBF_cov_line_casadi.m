function cov_line = RBF_cov_line_casadi(X_tr, x, l)
% shape:
% X_tr: [N_data_tr, D]
% x: [1,D]
% l: [1,D]
cov_line = exp(-sum(((X_tr-repmat(x,size(X_tr,1),1)).^2)./repmat((l.^2),size(X_tr,1),1),2))';
% X_tr = X_tr./l;
% x = x./l;
% cov_line = exp(-sum(X_tr.^2 + x.^2 -2*X_tr.*x, 2))';
end