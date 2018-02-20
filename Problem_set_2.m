N=100;
S=zeros(N,N);
I=zeros(N,N);
a=0.1;
b=1;
c=0.5;
K=30;
D=1;

for i=1:N
    for j=1:N-1
    
    S(i,j+1)=S(i,j)+deltat*(b*(I(i,j)+S(i,j))-c*S(i,j)-S(i,j)*(I(i,j)+S(i,j))/K-a*S(i,j)*I(i,j)+D(S(i+1,j)+S(i-1,j)-2*S(i,j)));
    %Tillsätt boundary conditions
    end
end