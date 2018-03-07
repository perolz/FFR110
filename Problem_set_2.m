clc, clear all
N=100;
iterations=700;
S=zeros(N,iterations);
I=zeros(N,iterations);
a=0.1;
b=1;
c=0.5;
K=30;
D=1;
S(1,1)=b/a;
I(1,1)=(-b+a*b*K-a*c*K)/a;
deltat=0.1;
deltax=1;
for j=1:iterations-1
    for i=1:N
        if(i==1)
            S(i,j+1)=S(i,j)+deltat*(b*(I(i,j)+...
                S(i,j))-c*S(i,j)-S(i,j)*(I(i,j)+S(i,j))/K-a*S(i,j)*I(i,j)+D*(S(i+1,j)+S(N,j)-2*S(i,j))/deltax);
            I(i,j+1)=I(i,j)+deltat*(-c*I(i,j)-I(i,j)*(I(i,j)+S(i,j))/K+a*S(i,j)*I(i,j)...
                +D*(I(i+1,j)+I(N,j)-2*I(i,j))/deltax);
        elseif i==N
            S(i,j+1)=S(i,j)+deltat*(b*(I(i,j)+...
                S(i,j))-c*S(i,j)-S(i,j)*(I(i,j)+S(i,j))/K-a*S(i,j)*I(i,j)+D*(S(1,j)+S(i-1,j)-2*S(i,j))/deltax);
            I(i,j+1)=I(i,j)+deltat*(-c*I(i,j)-I(i,j)*(I(i,j)+S(i,j))/K+a*S(i,j)*I(i,j)...
                +D*(I(1,j)+I(i-1,j)-2*I(i,j))/deltax);
        else
            I(i,j+1)=I(i,j)+deltat*(-c*I(i,j)-I(i,j)*(I(i,j)+S(i,j))/K+a*S(i,j)*I(i,j)...
                +D*(I(i+1,j)+I(i-1,j)-2*I(i,j))/deltax);
            S(i,j+1)=S(i,j)+deltat*(b*(I(i,j)+...
                S(i,j))-c*S(i,j)-S(i,j)*(I(i,j)+S(i,j))/K-a*S(i,j)*I(i,j)+D*(S(i+1,j)+S(i-1,j)-2*S(i,j))/deltax);
        end
    end
end
X=0:deltat:iterations*deltat-deltat;
Y=0:deltax:N*deltax-deltax;
fig=figure(1);
contourf(X,Y,S,20)
colorbar
title('Susceptibles as a function of time and space')
xlabel('Time (t)')
ylabel('Space (x)')
tmp1=find(S(50,:)>2,1,'first');
tmp2=find(S(50,:)>4,1,'first');
tmp3=find(S(50,:)>6,1,'first');
wavespeed=(S(50,tmp1+1)-S(50,tmp1))/deltat
saveas(fig,'Susceptibles.png')
fig2=figure(2);
contourf(X,Y,I,20)
colorbar
title('Infectives as a function of time and space')
xlabel('Time (t)')
ylabel('Space (x)')
tmp1=find(S(50,:)>2,1,'first');
wavespeedInfectives=(S(50,tmp1+1)-S(50,tmp1))/deltat;
saveas(fig2,'Infectives.png')