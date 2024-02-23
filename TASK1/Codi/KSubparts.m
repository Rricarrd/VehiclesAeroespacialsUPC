function [KNN,KND, KDN, KDD] = KSubparts(K, NeumannDOF, DirichlettDOF)

KNN = sparse(K(NeumannDOF,NeumannDOF));
KND = sparse(K(NeumannDOF,DirichlettDOF));
KDN = sparse(K(DirichlettDOF,NeumannDOF));
KDD = sparse(K(DirichlettDOF,DirichlettDOF));
end

