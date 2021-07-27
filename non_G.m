clear r;
N=100;   % matrix dimension
r = 0.5; % squeezing parameter
ra = .5; % displacement parameter
% l=1;   % number photon subtracted
% m=1;   % number photon added
% n=N-1; % Fock state photon number
% p=2;



%col = ['c' 'm' 'k' 'b'];
% w=1;
% T=1;

t=1001;
% v=4;
eta=linspace(-5*pi,5*pi,t);
%ea=linspace(-pi,pi,v);
% ea=(-pi:pi/30:pi); % phase parameter in cat state
%ea1=(-pi:pi/4:pi) % phase parameter

theta = 0;             % parameter in  displacement operator
phi=0;                 % parameter in  squeezing operator
a = diag(sqrt(1:N-1),1); % annihilation  Operator


%Ttt=(0:1:50);


VC=zeros(2,2);

% scat1=zeros(N,t);
% scat12=zeros(N,5);






%Squeezing
zta = r* exp(1i*phi);       %complex parameter in  squeezing operator
S = expm((1/2)*(conj(zta)*(a)^2-(zta)*(a')^2)); %Squeezing Operator



%gen = normalized_state(a^l*(a')^m*S*D*Fock(0,N));



%Coherent States

alpha  = ra* exp(1i*theta);              %parameter for +alpha coherent state 
alpham = -ra* exp(1i*theta);              %parameter for -alpha coherent state 
alphapc  = 1i*ra* exp(1i*theta);              %parameter for +alpha coherent state 
alphamc = -1i*ra* exp(1i*theta);              %parameter for -alpha coherent state 

D = expm(alpha*a'- conj(alpha)*a);     %Diplacement operator  for alpha
Dm = expm(alpham*a'- conj(alpham)*a);    %Diplacement operator for -alpha
Dpc = expm(alphapc*a'- conj(alphapc)*a); %Diplacement operator  for ialpha
Dmc = expm(alphamc*a'- conj(alphamc)*a); %Diplacement operator  for -ialpha



coh1 = normalized_state(D*Fock(0,N)); %coherent state alpha
coh2 = normalized_state(Dm*Fock(0,N));  %coherent state -alpha
coh3 = normalized_state(Dpc*Fock(0,N)); %coherent state ialpha
coh4 = normalized_state(Dmc*Fock(0,N)); %coherent state -ialpha

% In12=coh1'*coh2;%inner product
% In13=coh1'*coh3;
% In14=coh3'*coh4;
% In23=coh2'*coh3;
% In24=coh2'*coh4;
% In34=coh3'*coh4;
% 
% fIn12=(abs(coh1'*coh2))^2%Fidelity
% fIn13=(abs(coh1'*coh3))^2
% fIn14=(abs(coh3'*coh4))^2
% fIn23=(abs(coh2'*coh3))^2
% fIn24=(abs(coh2'*coh4))^2
% fIn34=(abs(coh3'*coh4))^2
% 
% 
%PASCS

Pas1 = normalized_state(a*a'*D*Fock(0,N)); %coherent state alpha
Pas2 = normalized_state(a*a'*Dm*Fock(0,N));  %coherent state -alpha
Pas3 = normalized_state(a*a'*Dpc*Fock(0,N)); %coherent state ialpha
Pas4 = normalized_state(a*a'*Dmc*Fock(0,N)); %coherent state -ialpha

PIn12= Pas1'*Pas2;
PIn13= Pas1'*Pas3;
PIn14= Pas3'*Pas4;
PIn23= Pas2'*Pas3;
PIn24= Pas2'*Pas4;
PIn34= Pas3'*Pas4;

%PASCS






F=[];
F1=[];

for tt=eta  
%for tt= -pi/4
%CAT STATE

% scat1= normalized_state(coh1+(expm(1i*tt))*coh2);

% for j=1:v
% scat12(:,j)= normalized_state(coh1+(expm(1i*ea(j)))*coh2);
% SIn1=scat12(:,j)'*scat12(:,j);
% end





% scat= normalized_state((cos(eta(i))*coh1+ sin(eta(i))*coh2)/sqrt(1+sin(2*eta(i))*exp(-2*(alpha)^2))); %general cat state
scat= ((cos(tt)*coh1 + sin(tt)*coh2)/sqrt(1+sin(2*tt)*exp(-2*(alpha)^2)));

SDM = scat*(scat)';               % cat state density matrix

nc= trace(a'*a*SDM); % average number of photons in general cat state


% 
%two mode entangled 

% Tsc1=normalized_state( kron(coh1,coh2) +kron(coh2,coh1));
% Tsc2=normalized_state( kron(coh1,coh2) -kron(coh2,coh1));
% Tsc3=normalized_state( kron(coh1,coh1) +kron(coh2,coh2));
% Tsc4=normalized_state( kron(coh1,coh1) -kron(coh2,coh2));

% DTsc=Tsc1*Tsc1';

%two mode entangled

%CAT STATE 






% %THERMAL STATE
% 
% nb=nc/100;                     %average number of photons
% x=nb/(nb+1);
% rth=zeros(N,N);  
% l=(0:1:N-1);
% 
% for p=l
%     rth = rth + (1-x)*(x^(p))*Fock(p,N)*Fock(p,N)'; %thermal state
% end
% 
% ntb = trace(a'*a*rth);
% 
% 
% dsrth= D*S*rth*S'*D';  % displaced squeezed thermal state
% 
% nt= trace(a'*a*dsrth); % average number of photons in displaced squeezed thermal state
% 
%THERMAL STATE
   
    



% %CoVariance Matrix 
% 

Var(SDM,N);
% 
% %CoVariance Matrix 


% IVC=inv(Var,N);


% nG Measures

% rho=SDM;   %nG state
% tao=dsrth; %reference Gaussian
% 
% Mu1=trace(rho^2);
% Mu2=trace(tao^2);
% K=trace(rho*tao);
% 
% D = (Mu1 + Mu2 - 2*K)/2;
% del_a= D/Mu1;      %hilbert-schmidt distance
% del_b=trace(rho*logm(rho))-trace(rho*logm(tao)); %Relative entropy


sx=sqrt(det(Var,N));
del_b = ((sx/2)+ 0.5)*log((sx/2)+0.5)-((sx/2)- 0.5)*log((sx/2)- 0.5);% - S(ρ),(S(ρ)=0,pure state) nG measure pure state(Relative entropy)
%del_b= (x+1)*log2(x+1)-x*log2(x);
% del_b = -x*log2(x)-(1-x)*log2(1-x);
%kurt = trace(den*(q^4))-3*((trace(den*(q^2))^2));
% nG Measures



% F=[F,del_a];
F1=[F1,del_b];

end
   


% vac = [1; zeros(N-1,1)];
% P = S*vac;                           %Squeezed Vacuum State
% ss = P/sqrt(sum(conj(P).*P));        %Normalized Squeezed Vacuum State

% cs1 = D*vac;
% cs = cs1/sqrt(sum(conj(cs1).*cs1));  %coherent state
% rho_c=cs*cs';                         %density matrix for coherent state
% rho_s=ss*ss';                         %density matrix for squeezed state
% gen1 = a^l*a'^m*S*D*vac
% gen = a^l*a'^m*S*D*vac/sqrt(sum(conj(gen1).*gen1));

%gen = normalized_state(a^l*(a')^m*S*D*Fock(0,N));

% ben = normalized_state(S*D*Fock(0,N));
% rho_t = expm((-w*(a')*a)/(T))/(trace(expm((-w*(a')*a)/(T))));%thermal-statedm

% rho=gen*gen';
% tao=ben*ben';
% Pn=Fock(n,N)'*rho*Fock(n,N);

% Reproduction of Mandel Q parameter
% avg=gen'*(a'*a)*gen
% sq_avg=gen'*((a')^2*a^2)*gen
% Q=(sq_avg/avg)-avg
% 
% 
% % dl = real(gen'*(a')^p*(a)^p*gen - (gen'*(a')*a*gen)^p);
%F=[F,Pn];
% F1=[F1,S];
%end

%%%FIGURE%%%%%%
 clf;
 h1=figure(1);
 set(gca, 'FontSize', 16, 'LineWidth', 1)
 %hold on
% plot(ea,F,'r','LineWidth',1.5,'MarkerSize',8)
%  xlabel('\eta','FontSize', 16)
%  ylabel('δ_a','FontSize', 18)
% 
 plot(eta,F1,'LineWidth',1.5,'MarkerSize',8)
 xlabel('\eta','FontSize', 16)
 ylabel('δ_b','FontSize', 18)
% %legend('numerical','numerical','analytic','analytic')
% bar(F)
 ylim([-5,5])
  xlim([-pi,pi])
print(h1,'-depsc','pnda10.eps')


%plot(eta,del_a)

function VC = Var(den,N)
a = diag(sqrt(1:N-1),1); % annihilation  Operator

q= (a + a')/sqrt(2);
m= (a-a')/1i*sqrt(2);

aq= trace(q*den);
ap= trace(m*den);


%r= (q*m + m*q)/2;% use when aq,ap are not 0

aqs= trace(q*q*den);

aps= trace(m*m*den);

apq= trace(m*q*den);

aqp= trace(q*m*den);

V11 = aqs-((aq)^2);
V12 = (apq + aqp)/2 - ap*aq;
V21 = (apq + aqp)/2 - ap*aq;
V22 = aps-((ap)^2);


VC=[V11 V12 ; V21 V22];% Variance Matrix Cat state
end


function state=Fock(n,N)
%Fock(n,N) = Fock state |n> in N-dim space

 a=diag(sqrt(1:N-1),1);
 
 vac=[1;zeros(N-1,1)];
 state=normalized_state(a'^n*vac);
end