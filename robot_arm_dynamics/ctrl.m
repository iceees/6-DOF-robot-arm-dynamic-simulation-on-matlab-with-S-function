function [sys,x0,str,ts] = ctrl(t,x,u,flag)
switch flag,  
case 0,       
    [sys,x0,str,ts]=mdlInitializeSizes;
case 3,       
    sys=mdlOutputs(t,x,u);
case {2,4,9}
    sys=[];  
otherwise  
    error(['Unhandled flag = ',num2str(flag)]);
end

%-------mdlInitializeSizes -------------------%
function [sys,x0,str,ts]=mdlInitializeSizes
sizes = simsizes;         
sizes.NumOutputs     = 6; 
sizes.NumInputs      = 18; 
sizes.DirFeedthrough = 1; 
sizes.NumSampleTimes = 1; 
sys = simsizes(sizes);    

x0  = [];                
str = [];                 
ts  = [0 0];              
q=sym('q',[1 6]);
dq=sym('dq',[1 6]);
ddq=sym('ddq',[1 6]);


%-------mdlOutputs(t,x,u) ------------------%
function sys=mdlOutputs(t,x,u)  
g=9.8;
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

R1=u(1);     
dr1=0;        
R2=u(2);      
dr2=0;        
R3=u(3);      
dr3=0;        
R4=u(4);      
dr4=0;       
R5=u(5);      
dr5=0;        
R6=u(6);      
dr6=0;        

x(1)=u(7);    
x(2)=u(8);    
x(3)=u(9);    
x(4)=u(10);    
x(5)=u(11);    
x(6)=u(12);    
x(7)=u(13);    
x(8)=u(14);    
x(9)=u(15);    
x(10)=u(16);    
x(11)=u(17);    
x(12)=u(18);    

e1=R1-x(1);   
e2=R2-x(3);   
e3=R3-x(5);   
e4=R4-x(7);   
e5=R5-x(9);   
e6=R6-x(11);  
e=[e1;e2;e3;e4;e5;e6];    
 
de1=dr1-x(2); 
de2=dr2-x(4); 
de3=dr3-x(6); 
de4=dr4-x(8); 
de5=dr5-x(10); 
de6=dr6-x(12); 
de=[de1;de2;de3;de4;de5;de6]; 

Kp=[50 0 0 0 0 0; 0 50 0 0 0 0 ; 0 0 50 0 0 0 ; 0 0 0 50 0 0 ; 0 0 0 0 50 0; 0 0 0 0 0 50]; 
Kd=[50 0 0 0 0 0; 0 50 0 0 0 0 ; 0 0 50 0 0 0 ; 0 0 0 50 0 0 ; 0 0 0 0 50 0; 0 0 0 0 0 50]; 

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
[G,M]=calculteGandM(dh,Pc,m,Ic,g,R,P);
tol=Kp*e+Kd*de+G; 

sys(1)=tol(1);  
sys(2)=tol(2);  
sys(3)=tol(3);  
sys(4)=tol(4);  
sys(5)=tol(5);  
sys(6)=tol(6);  
