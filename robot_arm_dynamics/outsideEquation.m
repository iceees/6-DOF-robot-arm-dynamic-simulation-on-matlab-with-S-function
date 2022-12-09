function [F,N] = outsideEquation(dh,dtheta,ddtheta,Pc,m,Ic,g,R,P)
%%外推公式，
number = size(dh,1);
Z = [0,0,1]';
w0 = zeros(3,1); dw0 = zeros(3,1);dv0 = [0;0;g];%基座的角速度和速度
%为什么没有速度的推导   因为直接写加速度的表达式，速度表达式用不到，但是角速度表达式可以用到
%% 1-n外推公式  速度推到好理解
w(:,1) = R(:,:,1)' * w0 + dtheta(1) * Z;
dw(:,1) = R(:,:,1)' * dw0 + cross(R(:,:,1)' * w0, dtheta(1) * Z) + ddtheta(1) * Z;
dv(:,1) = R(:,:,1)' * (cross(dw0,P(:,1)) + cross(w0,cross(w0, P(:,1))) + dv0);
dvc(:,1) = cross(dw(:,1), Pc(:,1))+cross(w(:,1), cross(w(:,1), Pc(:,1))) + dv(:,1);

for i = 1:number-1
   w(:,i+1) = R(:,:,i+1)' * w(:,i) + dtheta(i+1) * Z;
   dw(:,i+1) = R(:,:,i+1)' * dw(:,i) + cross(R(:,:,i+1)' * w(:,i), dtheta(i+1) * Z)+ ddtheta(i+1) * Z;
   dv(:,i+1) = R(:,:,i+1)' * (cross(dw(:,i), P(:,i+1)) + cross(w(:,i), cross(w(:,i), P(:,i+1))) + dv(:,i));
   dvc(:,i+1) = cross(dw(:,i+1), Pc(:,i+1)) + cross(w(:,i+1), cross(w(:,i+1), Pc(:,i+1))) + dv(:,i+1);
end

%力推导
for i = 1:number
   F(:,i)=m(i)*dvc(:,i);  %i关节质心出合力
   N(:,i)=Ic(:,:,i) * dw(:,i) + cross(w(:,i), Ic(:,:,i) * w(:,i));%为什么力矩会有两部分
end

end