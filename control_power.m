% V=(0:10:400);
% r=(0:10:400);
% [VV,rr]=meshgrid(V,r);
% I= log2(1+(VV+1).*(exp(2.*rr)))/2;
% figure
% surf(VV,rr,I);
% grid on

% rr=100;

 I1=zeros(100,100,2);
 I2=zeros(100,100,2);
   
 kr= linspace(0,30,100);
 ks= linspace(1,4,1);
 tt= linspace(0.01,100,100);
 s=tt/100;
 Z= [1 10];
 SZ= [10 100];
 k=1;
Y= ceil(k/2);


 % defining 'μ' function for 3-d plot (k=10 parties)
 %after revealing initial information
 % A general representation of control
 for rr= 1:2
     for i=1:100
        for sig1=1:100
            
I1(sig1,i,rr)= (0.5)*log(1+(SZ(1)/((.5).*exp(-2.*Z(rr))+kr(sig1)+((1/s(i))-1))));
I2(sig1,i,rr)= (0.5)*log(1+(SZ(1)/((.5).*exp(-2.*Z(rr))+kr(sig1)+(4.5)*((1/s(i))-1))));
      
     end
 end
 end
 [X,Y]=meshgrid(s,kr);
 for j=1:2
    
  surf(X,Y,I1(:,:,j));
  surf(X,Y,I2(:,:,j));  
 end
%  figure('DefaultAxesFontSize',14)

subplot(2,2,1);surf(X,Y,I1(:,:,1));
title('r=1,Even','FontSize',18)
ylabel('Σ_{control}','FontSize', 18)
xlabel('\eta','FontSize', 18)
zlabel('I (A,S)','FontSize', 18)
% xlim([0,1])
ylim ([0,30])
%zlim ([0,5])
subplot(2,2,2);surf(X,Y,I2(:,:,1));
title('r=1,Odd','FontSize',18)
ylabel('Σ_{control}','FontSize', 18)
xlabel('\eta','FontSize', 18)
zlabel('I (A,S)','FontSize', 18)
% xlim([0,1])
 ylim ([0, 30])
%zlim ([0,5])
subplot(2,2,3);surf(X,Y,I1(:,:,2));
title('r=10,Even','FontSize',18)
ylabel('Σ_{control}','FontSize', 18)
xlabel('\eta','FontSize', 18)
zlabel('I (A,S)','FontSize', 18)
% xlim([0,1])
 ylim ([0,30])
zlim ([0,5])
subplot(2,2,4);surf(X,Y,I2(:,:,2));
title('r=10,Odd','FontSize',18)
ylabel('Σ_{control}','FontSize', 18)
xlabel('\eta','FontSize', 18)
zlabel('I (A,S)','FontSize', 18)
% xlim([0,1])
ylim ([0,30])
zlim ([0,5])