S=randn(4096);
M=ones(size(S));
K=gaussian(-3:3, 0, 2);
K1=ones(2)/4;
%K1=gaussian(-12:12, 0, 4);
K=K/sum(K);
K1=K1/sum(K1);
rtot=zeros(size(S));
m0=M;
s0=S;
tic
for k=1:6
    s1=conv_corrected(s0, K, true); 
    m1=conv2_separable(m0, K );
    r1=conv_corrected((s1-s0).^2, K1, true);
    rtot=r1(2:2:end, 2:2:end)+rtot(2:2:end, 2:2:end);
    s0=s1(2:2:end, 2:2:end);
    m0=m1(2:2:end, 2:2:end);
    mean(rtot(m0>0.75)./m0(m0>0.75))
end
toc

tic
N=3*2^(k-1);
N0=2.^k;
K2=gaussian(-N:N, 0, 2^(k-1));
%K2=ones(N0)/(N0^2);
SS=conv_corrected(S, K2, true);
RR=conv_corrected((SS-S).^2, ones(N0), true);
MM=conv2_separable(M, K2);
SS1=SS(N0:N0:end, N0:N0:end);
RR1=RR(N0:N0:end, N0:N0:end);
MM1=MM(N0:N0:end, N0:N0:end);
mean(RR1(MM1>0.75)./MM1(MM1>.75))
toc