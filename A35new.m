%%%%35
clear all
th=zeros(2,2);
r=zeros(2,2);
[Ix,Iy,Iz]=prodop(1/2,1);
a=[1;0];b=[0;1];p=(a+b)/sqrt(2);m=(a-b)/sqrt(2);
% tqb(:,1)= kron(kron(a,a),a);tqb(:,2)= kron(kron(a,a),b);tqb(:,3)= kron(kron(a,b),a);tqb(:,4)= kron(kron(a,b),b);
% tqb(:,5)= kron(kron(b,a),a);tqb(:,6)= kron(kron(b,a),b);tqb(:,7)= kron(kron(b,b),a);tqb(:,8)= kron(kron(b,b),b);
% hi= (332.5*(10)^6)*2*Ix(:,:,1);
% hf= (332.5*(10)^6)*eye(2);
col = ['c' 'm' 'k' 'b'];
hi= ((10)^6)*2*Iz(:,:,1);
hf= ((10)^6)*eye(2);
T=10e-6;
% T=1;
n=9;%
M=8;% number of steps=8
m=0:8;
t= linspace(0,T,n);
s=t/T;
t=t(2)-t(1);
for i=1:9
 h(:,:,i)= (1-s(i))*hi+s(i)*hf;
 [v,e]=eig(h(:,:,i));
 e(:,:,i)=e;
end
u=zeros(2,2,8);
for i=1:M
h1(:,:,i)= (1-m(i+1)/M)*hi+(m(i+1)/M)*hf;
u(:,:,i)= expm(-1i*(1-m(i+1)/M)*hi*t)*expm(-1i*(m(i+1)/M)*hf*t);
u1(:,:,i)= expm(-1i*h1(:,:,i)*t);
u2 = u(:,:,i)
[th,r]=cart2pol(real(u2),imag(u2))
%  det(u2);
end  
return

%  for k=1:2
% for l=1:2  
%     
% [th(k,l),r(k,l)]=cart2pol(real(u2(k,l)),imag(u2(k,l)))
% 
% end
%  end
% for i=1:n
% th(:,:)
% r(:,:)
% 
%  [v,d]=eig(h(:,:,i));
%  v1(:,:,i)= v;d1(:,:,i)= d;
%  v=v1(:,:,i);
%  
% %  [a1 k]=min(diag(d1(:,:,i)));
% %  grst(:,i)= v(:,k); 
% %  st(:,1)=v1(:,:,1);
%  st(:,i)=v(:,1);
%  st(:,i);
% end
% % 
% % op1= a*a';op2= b*b';
% % for j=1:n-1




    
%     grst(:,j+1)= u(:,:,j+1)*grst(:,j);
%     grst(:,j+1)= grst(:,j+1)/0.7071;
 %     st(:,j+1)= u(:,:,j+1)*st(:,j);
%end
%    st= st/-0.7071;
% here I have calculated phases manually by compairing state at every time
% point with cos(phi1)+isin(phi1). thij corresponds to the phase associated
% with ith state. values of j corresponds to the trignometric functions cosine and sine functions. 
% th111 = 0.9375;th112 = 0.9375;th121 = 0.9375;th122 = 0.9375;the(2)=0.9375;
% th211 = 1.57;th212 = 1.57;th221 = 1.57;th222 = 1.57;the(3)=1.57;
% th1=pi-1.266627772030418;th2=pi-1.266627772030418; the(4)=pi-1.266627772030418;
% th1=pi-1.266627772030418;th2=pi-1.266627772030418;the(5)=pi-1.266627772030418;
% th1 = 1.57;th2 = 1.57;the(6)=pi/2;
% th1=0.937492959136589;th2=0.937492959136589;the(7)=0.937492959136589;
% th1 = 2*pi;th2 = 2*pi; the(8)=2*pi;
% th1 = 2*pi-1.249996813450404;th2 = 2*pi-1.249996813450404;the(9)=2*pi-1.249996813450404;
 
% for i=1:n
%     st(:,i)= st(:,i)*exp(-1i*the(i));
%     rho(:,:,i)= st(:,i)*st(:,i)';
% end
% st=st*0.7071;
%    rho(:,:,8)= st(:,8)*st(:,8)';
%   p(1,8) = trace(rho(:,:,8)*op1);p(2,8) = trace(rho(:,:,8)*op2);
%   bar3(p)

% r21=real(st(1,2));r22=real(st(2,2));i21=imag(st(1,2));i22=imag(st(2,2));
% r11=rho(1,1,1);r12=rho(2,2,1); 
% r21=rho(1,1,1);r22=rho(2,2,1); 
% r31=rho(1,1,1);r32=rho(2,2,1); 
% r41=rho(1,1,1);r42=rho(2,2,1); 
% r51=rho(1,1,1);r52=rho(2,2,1); 
% r61=rho(1,1,1);r62=rho(2,2,1); 
% r71=rho(1,1,1);r72=rho(2,2,1); 
% r81=rho(1,1,1);r82=rho(2,2,1); 

% h1=[r21 r22];h2=[i21 i22];
% h1=[r11 r12];%%% dm elements at point 2
% h2=[r21 r22];%%% dm elements at point 2
% h3=[r31 r32];%%% dm elements at point 2
% h4=[r41 r42];%%% dm elements at point 2
% h5=[r51 r52];%%% dm elements at point 2
% h6=[r61 r62];%%% dm elements at point 2
% h7=[r71 r72];%%% dm elements at point 2
% h8=[r81 r82];%%% dm elements at point 2

% r41=real(st(1,4));r42=real(st(2,4));i41=imag(st(1,4));i42=imag(st(2,4));
% h3=[r41 r42];h4=[i41 i42];
% r21=rho(1,1,4);r22=rho(2,2,4); 
% hr4=[r41 r42];hi4=[i41 i42];
% 
% r61=real(st(1,6));r62=real(st(2,6));i61=imag(st(1,6));i62=imag(st(2,6));
% % h5=[r41 r22];h6=[i21 i22];
% hr6=[r41 r22];hi6=[i21 i22];
% 
% r81=real(st(1,8));r82=real(st(2,8));i81=imag(st(1,8));i82=imag(st(2,8));
% % h7=[r81 r82];h8=[i81 i82];
% hr8=[r81 r82];hi8=[i81 i82];


% % subplot(4,2,1);hist(h1);subplot(4,2,2);hist(h2);
% % subplot(4,2,3);bar3(real(rho(:,:,2)));subplot(4,2,4);bar3(imag(rho(:,:,2)));
% % subplot(4,2,5);bar3(real(st(:,6)));subplot(4,2,6);bar3(imag(st(:,6)));
% % subplot(4,2,7);bar3(real(st(:,8)));subplot(4,2,8);bar3(imag(st(:,8)));

% figure(1)
% subplot(4,2,1);hist(hr1);subplot(4,2,3);hist(hr3);subplot(4,2,5);hist(hr5);subplot(4,2,7);hist(hr7);
% subplot(4,2,2);hist(hr2);subplot(4,2,4);hist(hr4);subplot(4,2,6);hist(hr6);subplot(4,2,8);hist(hr8);

for i=1:2
plot(s,e(i,:),col(i))
hold on
end
% subplot(4,2,1);bar3(imag(rho(:,:,2)));subplot(4,2,2);bar3(imag(rho(:,:,2)));
% subplot(4,2,3);bar3(real(rho(:,:,2)));subplot(4,2,4);bar3(imag(rho(:,:,2)));
% subplot(4,2,5);bar3(real(st(:,6)));subplot(4,2,6);bar3(imag(st(:,6)));
% subplot(4,2,7);bar3(real(st(:,8)));subplot(4,2,8);bar3(imag(st(:,8))*10^15);
% h=subplot(4,2,1)
% p=get(h,'daspectratio')
% p(1)=p(1)+0.25;p(2)=p(2)+0.25;
% % p(3)=p(3)+0.05;p(1)=p(1)+0.05;p(2)=p(2)+0.05;
% set(h,'pos',p)