%% RUN THIS FUNCTION NEW VERSION...
%function [pt,y1]=CHW4_main()
clear all;
CHW4;
syms p34;
N=9;
pts3d1=[pts3d(:,1:8) pts3d(:,12)];
pts2d1=[pts2d(:,1:8) pts2d(:,12)];
save('2d_used_points.mat','pts2d1');
%create A matrix
for i=1:N
A(i,:)=[pts3d1(:,i)' 0 0 0 0 -pts3d1(1:3,i)'.*pts2d1(1,i)];
end
for i=N+1:2*N
A(i,:)=[0 0 0 0 pts3d1(:,i-N)' -pts3d1(1:3,i-N)'.*pts2d1(2,i-N)];
end
%--------------------------------------------------------------------------
%create y matrix in respect of p34
y(1:N,1)=pts2d1(1,1:N)'.*p34;
y(N+1:2*N,1)=pts2d1(2,1:N)'.*p34;
%--------------------------------------------------------------------------
Z=(A.'*A)\A.'*y; %Z
P34=double(solve(Z(11)^2+Z(10)^2+Z(9)^2-1,p34)); %solve contraint ||q3||=1 to evaluate p34


y1(1:N,1)=pts2d1(1,1:N)'.*P34(2); %determine y matrix by p34
y1(N+1:2*N,1)=pts2d1(2,1:N)'.*P34(2); 
Z1=double((A.'*A)\A.'*y1);  %evaluate Z matrix contain all pij
%-------------------------------------------------------------------------
%evaluate again specefic dots on 2d image by project 3d points by p~ matrix 
%(reprojected points)
for i=1:N
    xi(:,i)=[Z1(1:4)';Z1(5:8)';Z1(9:11)' P34(2)]*pts3d1(:,i);
    xi(:,i)=xi(:,i)./xi(3,i);
end
%plot dots again on image by red circles
figure(2); plot(xi(1,:),xi(2,:),'O','Color','R');
title('"O" labled are projected points using P~ and "*" labled are initial points in 2d');
save('2d_reprojct_pts.mat','xi');
%--------------------------------------------------------------------------
% specify pij values in p~ matrix 
p1=Z1(1:3);
p14=Z1(4);
p2=Z1(5:7);
p24=Z1(8);
p3=Z1(9:11);
p34=P34(2);

%evaluate intrinsic and extrinsic camera parameters 

%EVALUATE TETA
teta=acosd(-dot(cross(p1,p3),cross(p2,p3))/(norm(cross(p2,p3))*norm(cross(p1,p3))));

r3=p3;
tz=p34;
u0=p1.'*p3;
v0=p2.'*p3;
alphau(1)=norm(cross(p1,p3),2)*sind(teta);
alphau(2)=-norm(cross(p1,p3),2)*sind(teta);
alphav(1)=norm(cross(p2,p3),2)*sind(teta);
alphav(2)=-norm(cross(p2,p3),2)*sind(teta);
r1(:,1)=(1/sind(teta))*(cross(p3,cross(p1,p3)/norm(cross(p1,p3),2))+cross(p3,cross(p2,p3)/norm(cross(p2,p3),2))*cosd(teta));
r1(:,2)=-(1/sind(teta))*(cross(p3,cross(p1,p3)/norm(cross(p1,p3),2))+cross(p3,cross(p2,p3)/norm(cross(p2,p3),2))*cosd(teta));
r2(:,1)=cross(p3,cross(p2,p3)/norm(cross(p2,p3),2));
r2(:,2)=-cross(p3,cross(p2,p3)/norm(cross(p2,p3),2));
tx=(p14+(p24-v0*p34)*alphau(1)*cosd(teta)/alphav(1)-u0*p34)/alphau(1);
ty=(sind(teta)/alphav(1))*(p24-v0*p34);
%--------------------------------------------------------------------------
%evaluate optical center coordinate
R=[r1(:,1).';r2(:,1).';r3.'];
T=[tx;ty;tz];
co=-R\T; %optical center coordinate in old coordinate system=[0,0,0] and we know:
            % Mold=T+R*Mnew therefore Mnew=-inv(R)*T=new-optical-center;
save('intrc_OptCent.mat','co','alphav','alphau','u0','v0');

%BUILD NEW P~ MATRIX WITH NEW INTRINSIC AND EXTRINSIC PARAM.
pt2=[alphau(1)*r1(:,1)'+u0*r3',alphau(1)*tx+u0*tz;alphav(1)*r2(:,1)'+v0*r3',alphav(1)*ty+v0*tz;r3',tz];

f11=pt2(1,1:3);
f12=pt2(2,1:3);
f13=pt2(3,1:3);
%CHECK CONSTRAINT 2 FROM 3.23 BOOK
check=dot(cross(f11,f13),cross(f12,f13));

% for say that has did pass ball from goal line
% select four point on ball on 2d picture and then use aech point in
% m~=p~*M~ and we know
% v=alphav*y/z+v0 and z>=0 
% so (y/z)=(v-v0)/alphav ----
% solve 3 equation with three unknown parameters and get results

%ball_coor=[1661,1661,1661,1661;981.8,981.8,981.8,981.8];%coordinates of four selected point on ball
%v=[981.8,981.8,981.8,981.8];
ball_coor=[1591,1595,1595,1604;675,690,660,679];%coordinates of four selected point on ball
v=[675,690,660,679];
%ball_coor=[967,967,967,967;1213,1213,1213,1213];
pt=pt2;  % P~ matrix
save('Ptild.mat','pt2');
syms y1 x1 z1
Mt=[x1;y1;z1;1]; 
mt=pt*Mt;
mt1=mt(1:2)./mt(3);
opt_c=-(pt(1:3,1:3)\pt(1:3,4));
save('opt_c.mat','opt_c');
for i=1:4
    yz(i)=(v(i)-v0)/alphav(2);
    out(i)=solve(mt1(1)-ball_coor(1,i),mt1(2)-ball_coor(2,i),y1/z1-yz(i),x1,y1,z1);
end
x1=double(out(1).x1+out(2).x1+out(3).x1+out(4).x1)/4; %determine mean of four point to increase accuracy

y1=double(out(1).y1+out(2).y1+out(3).y1+out(4).y1)/4; % y coordinate if negative then ball pass from goal line  
                                                                %else is
                                                                %not goal
z1=double(out(1).z1+out(2).z1+out(3).z1+out(4).z1)/4;
save('coords_ball.mat','x1','y1','z1');

%end