
% 
% alpha=0.5; % displacement parameter
% r=0.5 ;     % squeezing parameter
% 
% 
M=50; 
Hyp=zeros(M,1);% Hilbert Space dimension
Fock(0,M);

t=101;
eta=linspace(0.3,0.7,t);
% 
% theta = 0;             % parameter in  displacement operator
% phi= 0;                   % squeezing parameter
% 
a = diag(sqrt(1:M-1),1); % annihilation  Operator
% 
 VC=zeros(2,2);
% 
% m=0;
% n=0;



%Defining parameters for hypergeometric state


% et=0.5; 



%  %Coherent States
% 
%  
% D = expm(alpha*a'- conj(alpha)*a);       %Diplacement operator  for alpha
% Dm = expm((-alpha)*a'- conj(-alpha)*a);    %Diplacement operator for -alpha
% 
% coh1 = normalized_state(D*Fock(0,N));  %coherent state alpha
% coh2 = normalized_state(Dm*Fock(0,N)); %coherent state -alpha





F=[];
% F=zeros(t,t);
% del_b=zeros(t,t,1);
% D1=zeros(t);
for et= eta  %prob parameter btw 0 and 1

k=M*et^(-1);
l=M*(1-et)^(-1);

W=[k l];
L=max(W)+1;

Ln= L*et;
Ln1=L*(1-et);
 

%Hypergeometric State
for n1=0:M-1
    
    Mn=M-n1;

% C1=nchoosek(Ln,n1);
% C2=nchoosek(Ln1,Mn);
% C3=nchoosek(L,M);
C1=gamma(Ln+1)/(gamma(Ln-n1+1)*gamma(n1+1));
C2=gamma(Ln+1)/(gamma(Ln-Mn+1)*gamma(Mn+1));
C3=gamma(L+1)/(gamma(L-M+1)*gamma(M+1));

%Hypergeomeric coeffecient;
HetL=(sqrt(C1*C2))*((sqrt(C3))^(-1));
 Hyp= Hyp + HetL*Fock(n1,M);
end
Hyp1=normalized_state(Hyp);
dHyp=Hyp1*Hyp1';

%Photon added Hypergeometric State
PAHyp= normalized_state(a'*Hyp1);
dPAHyp=PAHyp*PAHyp';

%Photon subtracted Hypergeometric State
PSHyp= normalized_state(a*Hyp1);
dPSHyp=PSHyp*PSHyp';


%Photon added then subtracted Hypergeometric State
PASHyp= normalized_state(a*a'*Hyp1);
dPASHyp=PASHyp*PASHyp';

% %COHERENT STATE  
%  D1 = expm(tt*a'- conj(tt)*a);
% coh11 = normalized_state(D1*Fock(0,N));  %coherent state alpha
% % COHERENT STATE 
% 
% %PASCS
% PScoh11=normalized_state((a^n)*(a'^m)*(D1*Fock(0,N)));
% Den_pscoh = PScoh11*(PScoh11)';
% %PASCS
% 
% 
% % %SQueezing
% % for r= eta
% % D1 = expm(tt*a'- conj(tt)*a);
% zta = r* exp(1i*phi);       %complex parameter in  squeezing operator
%  S = expm((1/2)*(conj(zta)*(a)^2-(zta)*(a')^2)); %Squeezing Operator
% % SQ = normalized_state((a^n)*((a')^m)*S*D1*Fock(0,N));
% % DSQ=SQ*(SQ)';
% % 
% % %SQueezing
% 
% 
% 
% % General state
% genS=normalized_state((a^n)*(a'^m)*(S*D1*S*D1*Fock(0,N)));
% DG=genS*(genS)';
% 
% % squeezed displaced and superposition
% genDS = normalized_state((a^n)*(a'^m)*(S*D1*Fock(0,N)));
% DG0=genDS*(genDS)';
% 
% %photon added squeezed
% genPS=normalized_state((S*D1*Fock(0,N)));
% DG1=genPS*(genPS)';
% 
% % CAT STATE
% 
% scat= ((cos(tt)*coh1 + sin(tt)*coh2)/sqrt(1+sin(2*tt)*exp(-2*(alpha)^2)));  %general cat state
% den_scat = scat*(scat)';    % cat state density matrix
% 
% %photon_added_cat_state
% pacs=normalized_state((a^n)*((a')^m)*(cos(tt)*coh1 + sin(tt)*coh2)/sqrt(1+sin(2*tt)*exp(-2*(alpha)^2)));
% den_pacs=pacs*(pacs)';

%CoVariance Matrix 
den=dHyp; 
q= (a + a')/sqrt(2);    % position Operator
p= (a-a')/1i*sqrt(2);   % momentum operator


aq= trace(q*den);       % avg position 
ap= trace(p*den);       % avg momentum



aqs= trace(q*q*den);

aps= trace(p*p*den);

apq= trace(p*q*den);

aqp= trace(q*p*den);

V11 = aqs-((aq)^2);
V12 = (apq + aqp )/2 - ap*aq;
V21 = (apq + aqp)/2 - ap*aq;
V22 = aps-((ap)^2);

VC=[V11 V12 ; V21 V22]; % Variance Matrix Cat state/Thermal state




% nG Measures
x=sqrt(det(VC));
del_b = (x/2+ 0.5)*log(x/2+0.5)-(x/2- 0.5)*log(x/2- 0.5); % As (S(ρ)=0 for pure state) nG measure pure state(Relative entropy)

% nG Measures

%F1(:,:)=del_b;

F=[F,del_b];
end



 clf;
 h1=figure(1);
 set(gca, 'FontSize', 16, 'LineWidth', 1)
%  plotfunc3d(del_b,tt,r)
% [X,Y]=ndgrid(tt,r);
% V=F(X,Y);
% figure
% surf(X,Y,V)
%  title('del_b','tt','r')
plot(eta,F,'r','LineWidth',1.5,'MarkerSize',8); grid on;
 xlabel('\eta','FontSize', 16)
%  zlabel('\eta2','FontSize', 18)
 ylabel('\delta_b','FontSize', 18)


 xlim([0.3,0.7])
  ylim([0,0.3])
 %zlim([-2*pi,2*pi])

 
 %%%%%%%%%%%%%%%%functions%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function result=normalized_state(state)
%normalized_state (state) - transforms to normalized column state (miran)

result=column_vector(state)/sqrt(sum(conj(state).*state));
end
 
 function result=column_vector(state)
% |state> writes as a column vector
%
% row vector -> column vector
% column vector -> the same

dim=size(state);
if dim(1)==1
state=state';
end;

result=state;
end



function state=Fock(n,N)
%Fock(n,N) = Fock state |n> in N-dim space

 a=diag(sqrt(1:N-1),1);
 
 vac=[1;zeros(N-1,1)];
 state=normalized_state(a'^n*vac);
end


%%%%%%%%%%%%%%%%functions%%%%%%%%%%%%%%%%%%%%%%%%%%%