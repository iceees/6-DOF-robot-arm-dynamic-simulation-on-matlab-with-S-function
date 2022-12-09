function [sys,x0,str,ts]=plant(t,x,u,flag)
% S-function for continuous state equation?
switch flag,  
%Initialization
  case 0,     
    [sys,x0,str,ts]=mdlInitializeSizes;
case 1,       
    sys=mdlDerivatives(t,x,u);
%Outputs
  case 3,
    sys=mdlOutputs(t,x,u);
%Unhandled flags
  case {2, 4, 9 }
    sys = [];
%Unexpected flags
  otherwise
    error(['Unhandled flag = ',num2str(flag)]);
end

%-------mdlInitializeSizes -------------------%
function [sys,x0,str,ts]=mdlInitializeSizes
global  g                 
sizes = simsizes;
sizes.NumContStates  = 12;   
sizes.NumDiscStates  = 0;   
sizes.NumOutputs     = 12;   
sizes.NumInputs      =6;    
sizes.DirFeedthrough = 0;   
sizes.NumSampleTimes = 0;   
sys=simsizes(sizes);
x0=[0 0 0 0 0 0 0 0 0 0 0 0];
str=[];
ts=[];
%% 
q=sym('q',[1 6]);
dq=sym('dq',[1 6]);
ddq=sym('ddq',[1 6]);
ddq(1)=0;
 ddq(2)=0;
 ddq(3)=0;
 ddq(4)=0;
 ddq(5)=0;
 ddq(6)=0;
 S=zeros(6,1);
 KK=zeros(6,15);
 G=zeros(6,1);
 M=zeros(6,6);
 CC=zeros(6,6);


%-------mdlDerivatives -------------------%
function sys=mdlDerivatives(t,x,u)
global  g  q dq ddq
%% 
g=9.8;      %  9.8
alpha=[0 pi/2 0 0 -pi/2 pi/2];
a=[0 0 0.264 0.237 0 0];
d=[0.144 0 -0.0075 0.114 0.114 0.067];
thet=[0 pi/2 0 -pi/2 0 0];
dh=[alpha;a;d;thet]';
Pc(:,1) = [0.0316 -3.1464 -13.8983]*10^-3;
Pc(:,2) = [131.5620 -0.0210 112.1840]*10^-3;
Pc(:,3) = [190.3840 0.0410 17.1800]*10^-3;
Pc(:,4) = [0.0886 21.0083 -2.5014]*10^-3;
Pc(:,5) = [-0.0886 -21.0083 -2.5014]*10^-3;
Pc(:,6) = [0 0 8.0000]*10^-3;
m = [2.920; 6.787; 2.450; 1.707; 1.707; 0.176]; 
Ic = zeros(3,3,6);
Ic(:,:,1)  = [42.614 0.046 0.062; 0.046 41.164 -1.386; 0.062 -1.386 31.883]*10^-4;
Ic(:,:,2)  = [100.7 -1.8 1.6; -1.8 1100.8 0; 1.6 0 1087.1]*10^-4;
Ic(:,:,3)  = [31.45 0.48 7.23; 0.48 172.41 -0.15; 7.23 -0.15 166.82]*10^-4;
Ic(:,:,4)  = [20.92 -0.061 0.078; -0.061 16.808 0.992; 0.078 0.992 19.75]*10^-4;
Ic(:,:,5)  = [20.92 -0.061 -0.078; -0.061 16.808 -0.992; -0.078 -0.992 19.75]*10^-4;
Ic(:,:,6)  = [0.9296 0 0; 0 0.9485 0; 0 0 1.5925]*10^-4; 
 tau=u(1:6);        
 q(1)=x(1);
 q(2)=x(3);
 q(3)=x(5);
 q(4)=x(7);
 q(5)=x(9);
 q(6)=x(11);
 dq(1)=x(2);
 dq(2)=x(4);
 dq(3)=x(6);
 dq(4)=x(8);
 dq(5)=x(10);
 dq(6)=x(12);
  [R,P] = compute_frame_transform(dh,q);
%   [F,N] = outsideEquation(dh,dq,ddq,Pc,m,Ic,g,R,P);
 [G,M]=calculteGandM(dh,Pc,m,Ic,g,R,P);
 [CC]=calculteCCandKK(dh,Pc,m,Ic,g,R,P,G);
 for i=1:6
   Cc(:,i)=CC(:,i).*dq(i);  
 end
 [KK]=calculteKK(dh,G,Pc,m,Ic,g,R,P,CC);
K1 = [0,0,0,0,0,0]';
K2=KK*[dq(1),0,0,0,0,dq(3),dq(4),dq(5),dq(6),0,0,0,0,0,0]';
K3 = KK*[0,dq(1),0,0,0,0,0,0,0,dq(4),dq(5),dq(6),0,0,0]';
K4 = KK*[0,0,dq(1),0,0,0,0,0,0,0,0,0,dq(5),dq(6),0]';
K5 = KK*[0,0,0,dq(1),0,0,0,0,0,0,0,0,0,0,dq(6)]';
K6 = KK*[0,0,0,0,dq(1),0,0,0,0,0,0,0,0,0,0]';
k=[K1,K2,K3,K4,K5,K6];
C = k + Cc;
%% 
 S=inv(M)*(tau-C*dq'-G);  
 ddq(1)=S(1);
 ddq(2)=S(2);
 ddq(3)=S(3);
 ddq(4)=S(4);
 ddq(5)=S(5);
 ddq(6)=S(6);
sys(1)=x(2);
sys(2)=S(1);
sys(3)=x(4);
sys(4)=S(2);
sys(5)=x(6);
sys(6)=S(3);
sys(7)=x(8);
sys(8)=S(4);
sys(9)=x(10);
sys(10)=S(5);
sys(11)=x(12);
sys(12)=S(6);
%-------mdlOutputs -------------------%
function sys=mdlOutputs(t,x,u)
sys(1)=x(1);  
sys(2)=x(2);  
sys(3)=x(3);  
sys(4)=x(4);  
sys(5)=x(5);  
sys(6)=x(6);  
sys(7)=x(7);  
sys(8)=x(8);  
sys(9)=x(9);  
sys(10)=x(10);  
sys(11)=x(11);  
sys(12)=x(12); 
