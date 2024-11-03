function [u,bi] = HaraseF2(m,N)

%   Reference:
%
%   S. Harase, A table of short-period Tausworthe generators 
%   for Markov chain quasi-Monte Carlo, 
%   J. Comput. Appl. Math., 384 (2021), p. 113136.
% 

sigma = [70,179,146,139,5192,1028,12749,20984,72349,92609];

g = sigma(m-9);
switch m-9
    case 1
        a = [1,0,0,0,0,0,1,1,0,1];
    case 2
        a = [1,1,0,0,1,0,0,1,1,0,1];
    case 3
        a = [1,1,1,1,1,0,0,1,0,0,1,1];
    case 4
        a = [1,1,1,0,1,0,0,0,1,0,1,1,1];
    case 5
        a = [1,0,1,0,1,1,0,1,1,1,1,0,1,1];
    case 6
        a = [1,1,0,1,1,0,0,1,1,1,0,1,0,1,1];
    case 7
        a = [1,1,0,1,0,1,1,1,1,1,0,0,1,0,0,1];
    case 8
        a = [1,0,1,1,1,0,0,0,0,1,0,1,1,0,0,0,1];
    case 9
        a = [1,1,0,1,0,1,1,0,1,0,1,0,0,0,1,1,0,1];
    case 10
        a = [1,0,1,1,0,1,1,1,1,0,0,0,1,1,0,0,1,0,0];
    case 11
        a = [1,1,1,0,1,0,1,0,1,1,1,0,0,1,1,1,0,0,1,0];
end
a = logical(a);
a = find(a==1)-1;
a = flip(m-a);
mm = N-1;
w = 32;

b = false(mm,1);
bi = false(w,mm);

b(1) = 1;
for i = (m+1):mm
    b(i) = mod(sum(b(i-a)), 2);
end
b = [b;b(1:w-1)];

for i = 1:mm
    temp1 = mod((i-1)*g+1,mm);
    temp1 = temp1+(temp1==0)*mm;
    temp2 = temp1+w-1;
    bi(:,i) = b(temp1:temp2);
end
u = (2.^(-(1:w)))*bi;

% Points before randomization

% u = u(:);
% s = gcd(d,N-1);
% if s == 1
%     U = (reshape(repmat(u,d,1),d,N-1))';
%     Ub = reshape(repmat(bi,1,d),[w,d,N-1]);
% else
%     U = zeros(N-1,d);
%     Ub = false(w,d,N-1);
%     bb = (N-1)/s;
%     ss = d/s;
%     Ua = [repmat(u,ss,1);u(1:s)];
%     Uab = [repmat(bi,1,ss),bi(:,1:s)];
%     for i = 1:s
%         U((i-1)*bb+1:i*bb,:) = (reshape(Ua(i:ss*(N-1)+i-1),d,bb))';
%         Ub(:,:,(i-1)*bb+1:i*bb) = reshape(Uab(:,i:ss*(N-1)+i-1),[w,d,bb]);
%     end
% end