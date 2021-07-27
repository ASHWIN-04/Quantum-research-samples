sig = 10;

r=1;

 I=zeros(7,100); %initiating Matrix for mutual information 4-mode GHZ
 I1=zeros(7,100); %initiating Matrix for mutual information after revealing control values
 I11=zeros(7,100); %initiating Matrix for mutual information after revealing control values
 
 I2=zeros(7,100); % for 2 mode squeezed state
 I3=zeros(7,100);%initiating Matrix for mutual information after revealing control values
 d=zeros(7,100);
 d1=zeros(7,100);
 kr= linspace(2,8,7); % no. of parties
%  ks= linspace(1,4,1);
 tt= linspace(.0001,1,100);
 s=tt;
 %J=randperm(10);

 J=[6     1     7     4    11     9     5     8     3    10     2]; % Control input by Charlie- position quadrature
 K=[11     1     9     7     8     4    10     6     5     2     3];% Control input by Charlie- momentum quadrature
 L= [0 19 16 30 13 29 8 17]; % Simplification of combined quadrature value
 M= [0,19,18,27,30,32,28,29,22];
 
 markers = ['o' '+' '*' '.' 'x' 'd' 'p' 's'];
 

Z=[0 1 0 1 0 1 0 1 0 1 0 1 0 1 0]; % reducing to odd and even cases

% function for  mutual info after result (k- parties)

 for k=kr
%      Y(k)=ceil(k/2);
% k
% Z(k)
% Z(k-1)
for i=1:100 
I(k,i)= Z(k)*(0.5)*log(1+(sig/((.5)*exp(-2*r)+L(k)+((1/s(i))-1)))) +Z(k-1)*(0.5)*log(1+(sig/((.5)*(exp(-2*r))+L(k)+((4.5)*((1/s(i))-1)))));
I1(k,i)= Z(k)*(0.5)*log(1+(sig/((.5)*exp(-2*r)+((1/s(i))-1)))) +Z(k-1)*(0.5)*log(1+(sig/((.5)*(exp(-2*r))+((4.5)*((1/s(i))-1)))));
I11(k,i)= Z(k)*(0.5)*log(1+(sig/((.5)*exp(-2*r)))) + Z(k-1)*(0.5)*log(1+(sig/((.5)*(exp(-2*r)))))

I2(k,i)= (0.5)*log(1+(sig/((.5)*(exp(-2*r)+exp(2*r))+M(k)+((2.5)*((1/s(i))-1)))));% two mode squeezed
I3(k,i)= (0.5)*log(1+(sig/((.5)*(exp(-2*r)+exp(2*r))+(2.5)*((1/s(i))-1))));% two mode squeezed 

end
 end 
 
 
 for k=kr
for i=1:100
d(k,i)=I(k,i)-I2(k,i);
d1(k,i)=I1(k,i)-I3(k,i);
end
 end 
 % function plot for noise and mutual info after result (k- parties)

  h2=figure(2);
  set(gca, 'FontSize', 16, 'LineWidth', 1)
  
  P2=plot(tt,d1(2,:),markers(2),'DisplayName','Party2','LineWidth',1,'MarkerSize',5);
  hold on;
  P3=plot(tt,d1(3,:),markers(3),'DisplayName','Party3','LineWidth',1,'MarkerSize',5);
  hold on;
  P4=plot(tt,d1(4,:),markers(4),'DisplayName','Party4','LineWidth',1,'MarkerSize',5);
  hold on;
  P5=plot(tt,d1(5,:),markers(5),'DisplayName','Party5','LineWidth',1,'MarkerSize',5);
  hold on;
  P6=plot(tt,d1(6,:),markers(6),'DisplayName','Party6','LineWidth',1,'MarkerSize',5);
  hold on;
  P7=plot(tt,d1(7,:),markers(7),'DisplayName','Party7','LineWidth',1,'MarkerSize',5);
  hold on;
  P8=plot(tt,d1(8,:),markers(8),'DisplayName','Party8','LineWidth',1,'MarkerSize',5);
  hold on;
  
 
  
%   P2=plot(tt,I(2,:),col(2),'DisplayName','Even','LineWidth',1.5,'MarkerSize',8);
%   hold on;
%   P3=plot(tt,I(3,:),col(3),'DisplayName','Part3','LineWidth',1.5,'MarkerSize',8);
%   hold on;
%   P4=plot(tt,I(4,:),col(4),'DisplayName','Party4','LineWidth',1.5,'MarkerSize',8);
%   hold on;
%   P5=plot(tt,I(5,:),col(5),'DisplayName','Party5','LineWidth',1.5,'MarkerSize',8);
%   hold on;
%   P6=plot(tt,I(6,:),col(6),'DisplayName','Party6','LineWidth',1.5,'MarkerSize',8);
%   hold on;
%   P7=plot(tt,I(7,:),col(7),'DisplayName','Odd','LineWidth',1.5,'MarkerSize',8);
%   hold on;
%   P8=plot(tt,I(8,:),col(8),'DisplayName','Party8','LineWidth',1.5,'MarkerSize',8);
%   hold on;
%  
  
  
  
  
  % function plot for noise and mutual info after result and control(k- parties)
  
%  h2=figure(2);
%   set(gca, 'FontSize', 16, 'LineWidth', 1)
%   
%   
%   P2=plot(tt,I1(2,:),col(2),'DisplayName','Party2','LineWidth',1.5,'MarkerSize',8);
%   hold on;
%   P3=plot(tt,I1(3,:),col(3),'DisplayName','Part3','LineWidth',1.5,'MarkerSize',8);
%   hold on;
%   P4=plot(tt,I1(4,:),col(4),'DisplayName','Party4','LineWidth',1.5,'MarkerSize',8);
%   hold on;
%   P5=plot(tt,I1(5,:),col(5),'DisplayName','Party5','LineWidth',1.5,'MarkerSize',8);
%   hold on;
%   P6=plot(tt,I1(6,:),col(6),'DisplayName','Party6','LineWidth',1.5,'MarkerSize',8);
%   hold on;
%   P7=plot(tt,I1(7,:),col(7),'DisplayName','Party7','LineWidth',1.5,'MarkerSize',8);
%   hold on;
%   P8=plot(tt,I1(8,:),col(8),'DisplayName','Party8','LineWidth',1.5,'MarkerSize',8);
%   hold on;



 xlabel('\eta','FontSize', 18)
 ylabel('I(A,S)','FontSize', 18)
%  legend([P2 P3 P4 P5 P6 P7 P8], 'Second','Third','Fourth','Fifth','Sixth','Seventh','Eighth')
 legend([ P4 P7], 'Even','Odd')

% for k=kr
%  p(k)=plot(tt,I(k,:),col(k),'DisplayName','Party','LineWidth',1.5,'MarkerSize',8)
%   hold on
% end
%  % Random encoding by Charlie
%  M=J+K;
 
 %  Sum of encoding
%  for n=2:11  
%     d=1;
%     G(n)= M(n);
% while d < n
%   G(n) = G(n)+ M(n-d);
%   d=d+1;
% end
%  end
