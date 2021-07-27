% V=(0:10:400);
% r=(0:10:400);
% [VV,rr]=meshgrid(V,r);
% I= log2(1+(VV+1).*(exp(2.*rr)))/2;
% figure
% surf(VV,rr,I);
% grid on

% rr=100;

 I1=zeros(100,100,4);
 I2=zeros(100,100,4);
 kr= linspace(0,100,100);
 ks= linspace(1,4,1);
 tt= linspace(0.01,1,100);
 s=tt;
 Z= [1 10^1];
%  k=1;
% Y= 1;


 % defining 'μ' function for 3-d plot (k=10 parties)
%  for k=2
 for rr= 1:2
     for i=1:100
        for sig1=1:100
I1(sig1,i,rr)= (0.5)*log(1+(sig1/((.5).*exp(-2.*Z(rr))+(1/s(i))-1))); %even state
I2(sig1,i,rr)=(0.5)*log(1+(sig/((.5)*(exp(-2.*Z(rr)))+((4.5)*((1/s(i))-1))))); %odd state
        end
     end
 end
 
 [X,Y]=meshgrid(tt,kr);
 for j=1:2
  surf(X,Y,I1(:,:,j));
  end
%  figure('DefaultAxesFontSize',14)

subplot(2,2,1);surf(X,Y,I1(:,:,1));
title('r=1','FontSize',18)
 ylabel('Σ','FontSize', 18)
 xlabel('\eta','FontSize', 18)
zlabel('I (A,S)','FontSize', 18)
 xlim([0,1])
 ylim ([0,20])
%zlim ([0,5])
subplot(2,2,2);surf(X,Y,I1(:,:,2));
title('r=10','FontSize',18)
ylabel('Σ','FontSize', 18)
xlabel('\eta','FontSize', 18)
zlabel('I (A,S)','FontSize', 18)
xlim([0,1])
ylim ([0,20])


subplot(2,2,3);surf(X,Y,I2(:,:,1));
title('r=1','FontSize',18)
 ylabel('Σ','FontSize', 18)
 xlabel('\eta','FontSize', 18)
zlabel('I (A,S)','FontSize', 18)
 xlim([0,1])
 ylim ([0,20])
%zlim ([0,5])

subplot(2,2,4);surf(X,Y,I2(:,:,2));
title('r=10','FontSize',18)
ylabel('Σ','FontSize', 18)
xlabel('\eta','FontSize', 18)
zlabel('I (A,S)','FontSize', 18)
xlim([0,1])
ylim ([0,20])
% zlim ([0,1])

% subplot(2,2,3);surf(X,Y,I1(:,:,3));
% title('k=1000','FontSize',18)
% zlim([0,5])
% subplot(2,2,4);surf(X,Y,I1(:,:,4));
% title('k=10000','FontSize',18)
% zlim([0,5])