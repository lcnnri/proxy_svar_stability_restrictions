function D = dupl (n)
% DUPL Duplication matrix
% returns duplication matrix D_n. Output is sparse.
%
% see Lutkepohl (1991) New Introduction To Multiple Time Series Analysis
%   appendix A.12.2
%=========================================================================%
m   = n * (n + 1) / 2;
nsq = n^2;
r   = 1;
a   = 1;
v   = zeros(1, nsq);
for i = 1:n
   v(r:r + i - 2) = i - n + cumsum(n - (0:i-2));
   r = r + i - 1;
   
   v(r:r + n - i) = a:a + n - i;
   r = r + n - i + 1;
   a = a + n - i + 1;
end
D = full(sparse(1:nsq, v, 1, nsq, m));
end % eof
