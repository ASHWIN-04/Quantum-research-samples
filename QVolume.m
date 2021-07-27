ket0=[1 ;0];
ket1=[0 ;1];
bra0= ctranspose(ket0);
r00=ket0*bra0;%1-d
r11=ket1*ctranspose(ket1);
sigx =  [0 1 ; 1 0];
sigy=   [0 -1i ; 1i 0];
sigz =  [1 0 ; 0 -1];
H = [1 1; 1 -1]/sqrt(2);
I = eye(2);
If = eye(8);
Ir=eye(4);
sigxx =  kron(sigx,I);
sigxxx = kron(I,sigx);
sigzz =  kron(sigz,I);
sigzzz = kron(I,sigz);
sigyy  = kron(sigy,I);
sigyyy = kron(I,sigy);

sigx1 = kron(sigxx,I);
sigx2 = kron(sigxxx,I);
sigx3 = kron(I,sigxxx);

sigz1=kron(sigzz,I);
sigz2=kron(sigzzz,I);
sigz3=kron(I,sigzzz);

xy= kron(sigx,sigy);
xz= kron(sigx,sigz);

%%%%%%%%%%%%%%%%%%%%% Quantum Volume %%%%%%%%%%%%%%%%% %Q1- part-1, h2=1

Rx1=[cos(0.1*pi/2)  -1i*sin(0.1*pi/2); -1i*sin(0.1*pi/2) cos(0.1*pi/2)];
Rx2=[cos(0.2*pi/2)  -1i*sin(0.2*pi/2); -1i*sin(0.2*pi/2) cos(0.2*pi/2)];
Rx3=[cos(-0.6*pi/2)  -1i*sin(-0.6*pi/2); -1i*sin(-0.6*pi/2) cos(-0.6*pi/2)];
Rx4=[cos(0.2*pi/2)  -1i*sin(0.2*pi/2); -1i*sin(0.2*pi/2) cos(0.2*pi/2)];

Rz1=[1  0; 0 exp(1i*0.4*pi)];
Rz2=[1  0; 0 exp(1i*0.4*pi)];

Ry1=[cos(pi/2)  -sin(pi/2);sin(pi/2) cos(pi/2)];
Ry2=[cos(0.5*pi/2) -sin(0.5*pi/2);sin(0.5*pi/2) cos(0.5*pi/2)];


Rxy1= kron(Rx1,Ry1);
Rxz1= kron(Rx2,Rz1);
Rxy2= kron(Rx3,Ry2);
Rxz2= kron(Rx4,Rz2);

Cn = [1 0 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0];%control qubit above
Cn1 = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0];%control qubit below
HI = kron(H,I);
SW = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];

ket00=kron(ket0,ket0);
ket01=kron(ket0,ket1);
bra00=ctranspose(ket00);

rhoi= ket00*bra00;

U=Cn1*Rxz2*Cn1*HI*Cn1*Rxy2*SW*Cn1*Rxz1*Cn1*HI*Cn1*Rxy1;
kket=U*ket00;

rhof=U*rhoi*ctranspose(U);

[v,d]=eig(rhof);% this outputs the state and probability distribution

 %heavy output : h2=1, form above eigenvalue of density matrix

%%%%%%%%%%%%%%%%%%%%% toffoli %%%%%%%%%%%%%%%%% %Q2 part-1
To=[1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 ;0 0 1 0 0 0 0 0; 0 0 0 0 0 0 0 1 ; 0 0 0 0 1 0 0 0; 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 1 0; 0 0 0 1 0 0 0 0];%control qubit above
To1=[1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 ;0 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0 ; 0 0 0 0 1 0 0 0; 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 0 1; 0 0 0 0 0 0 1 0];%control qubit below
CCZ= [1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 ;0 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0 ; 0 0 0 0 1 0 0 0; 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 -1];
ket000=kron(ket00,ket0);
% IIH=  kron(Ir,H);
% HII=  kron(H,Ir);
% 
% X=IIH*CCZ*IIH;
% Y=IIH*To1*IIH;

Tk=To*ket000;
rtk=Tk*ctranspose(Tk);
rti= [.403 0 0 0 0 0 0 0; 0 .263 0 0 0 0 0 0 ;0 0 .01 0 0 0 0 0; 0 0 0 .008 0 0 0 0 ; 0 0 0 0 .191 0 0 0; 0 0 0 0 0 .110 0 0; 0 0 0 0 0 0 .012 0; 0 0 0 0 0 0 0 .003];     %obtained from experiment


F= (trace(sqrtm(rtk)*sqrtm(rti)))^2; %Fidelity


%%%%%%%%%%%%%%%%%%%%% K-State %%%%%%%%%%%%%%%%% %Q3
% https://github.com/qiskit-community/qiskit-community-tutorials/blob/master/awards/teach_me_qiskit_2018/w_state/W%20State%201%20-%20Multi-Qubit%20Systems.ipynb
ket100=kron(ket1,ket00);
ket001=kron(ket00,ket1);
ket010=kron(ket01,ket0);

th=pi/2;
x=pi/2;
y=pi/2;
ph=1/sqrt(2) ;         %0.785398163;             %1.91063323622;
ph1=   1/sqrt(3);         % 0.955316618 ;
uu=[cos(th/2)  -exp(1i*x)*sin(th/2); exp(1i*y)*sin(th/2) exp(1i*(x+y))*cos(th/2)];
RY=[cos(ph/2)  -sin(ph/2); sin(ph/2) cos(ph/2)];
RY1=[cos(-ph/2)  -sin(-ph/2); sin(-ph/2) cos(-ph/2)];
RY2=[cos(ph1/2)  -sin(ph1/2); sin(ph1/2) cos(ph1/2)];
RY3=[cos(-ph1/2)  -sin(-ph1/2); sin(-ph1/2) cos(-ph1/2)];



CU=kron(I,r00)+kron(uu,r11);
CZ=kron(I,r00)+kron(sigz,r11);
CX=kron(I,r00)+kron(sigx,r11)
CY=kron(I,r00)+kron(sigy,r11);

XC=kron(r00,I)+kron(r11,sigx)

U1= sigx1*kron(Cn,I)*kron(I,Cn)*kron(CU,I)*kron(RY,Ir);
U1*ket000;
% rJ=U1*ket000*ctranspose(U1*ket000);
U2= To*kron(I,CX)*To*kron(I,Cn1)*kron(Cn1,I)*To*kron(Ir,RY)*kron(I,CZ)*kron(Ir,RY1)*kron(kron(I,RY2),I)*kron(CZ,I)*kron(kron(I,RY3),I);
jj=U2*ket100;% w-state kron(CZ,Rz1)


p0=pi/4;
p1=-pi/3;
p2=0;
p3=-pi/4;
Rz0=[1  0; 0 exp(1i*p0)];
Rz=[1  0; 0 exp(1i*p1)];
Rzz=[1  0; 0 exp(1i*p2)];
Rzzz=[1  0; 0 exp(1i*p3)];
Rzz1=sigx*Rzz*sigx;
Rzz2=sigx*Rzzz*sigx;
R2z2=Rz*Rzz1;
R2z3=Rz0*Rzz2;
RX=kron(CZ,R2z2)*kron(R2z3,CZ);
jj1=RX*jj;


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















