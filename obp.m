%oblique projection  A/BC 
% obp_ABCd is calculated by using the definition
% obp_ABClq is calculated by using LQ decomposition
function [obp_ABCd,obp_ABClq]=obp(A,B,C)
obp_ABCd=A*[C' B']*pinv([C;B]*[C' B'])*[C;0*B];
[L11,L21,L22,L31,L32,L33,Q1T,Q2T,Q3T]=lqBCA(B,C,A);
obp_ABClq=L32*inv(L22)*[L21 L22]*[Q1T;Q2T];
