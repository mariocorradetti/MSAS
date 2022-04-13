function []=sat(ang1,ang2)
% Function that creates the shape of LUMIO
% Body
A=[0.5 0.5 0.5];
B=[-0.5 0.5 0.5];
C=[-0.5 -0.5 0.5];
D=[0.5 -0.5 0.5];
E=[0.5 0.5 -0.5];
F=[-0.5 0.5 -0.5];
G=[-0.5 -0.5 -0.5];
H=[0.5 -0.5 -0.5];

X=[A(1) B(1) C(1) D(1); D(1) C(1) G(1) H(1); G(1) H(1) E(1) F(1); E(1) F(1) B(1) A(1); B(1) C(1) G(1) F(1); A(1) D(1) H(1) E(1)]';
Y=[A(2) B(2) C(2) D(2); D(2) C(2) G(2) H(2); G(2) H(2) E(2) F(2); E(2) F(2) B(2) A(2); B(2) C(2) G(2) F(2); A(2) D(2) H(2) E(2)]';
Z=[A(3) B(3) C(3) D(3); D(3) C(3) G(3) H(3); G(3) H(3) E(3) F(3); E(3) F(3) B(3) A(3); B(3) C(3) G(3) F(3); A(3) D(3) H(3) E(3)]';
Co=[99/255,99/255,99/255]; 

xcg_b=0;
ycg_b=0;
zcg_b=0;
W=226;
H=340;
L=226;
Xbody= xcg_b+X*W;
Ybody= ycg_b+Y*L;
Zbody= zcg_b+Z*H;

fill3(Xbody,Ybody,Zbody,Co,'FaceAlpha',1)
hold on
axis equal
xlabel('x [cm]')
ylabel('y [cm]')
zlabel('z [cm]')


%% Yoke
Ly=198;
Wy=0.2;
Hy=100;

y=[0 1 0];

xy= xcg_b;
yy1= ycg_b+L/2+Ly/2;
yy2=-yy1;
zy= zcg_b;

scal=[Wy Ly Hy];
A=rotatevect([0.5 0.5 0.5].*scal,ang1,y);
B=rotatevect([-0.5 0.5 0.5].*scal,ang1,y);
C=rotatevect([-0.5 -0.5 0.5].*scal,ang1,y);
D=rotatevect([0.5 -0.5 0.5].*scal,ang1,y);
E=rotatevect([0.5 0.5 -0.5].*scal,ang1,y);
F=rotatevect([-0.5 0.5 -0.5].*scal,ang1,y);
G=rotatevect([-0.5 -0.5 -0.5].*scal,ang1,y);
H=rotatevect([0.5 -0.5 -0.5].*scal,ang1,y);


X=[A(1) B(1) C(1) D(1); D(1) C(1) G(1) H(1); G(1) H(1) E(1) F(1); E(1) F(1) B(1) A(1); B(1) C(1) G(1) F(1); A(1) D(1) H(1) E(1)]';
Y=[A(2) B(2) C(2) D(2); D(2) C(2) G(2) H(2); G(2) H(2) E(2) F(2); E(2) F(2) B(2) A(2); B(2) C(2) G(2) F(2); A(2) D(2) H(2) E(2)]';
Z=[A(3) B(3) C(3) D(3); D(3) C(3) G(3) H(3); G(3) H(3) E(3) F(3); E(3) F(3) B(3) A(3); B(3) C(3) G(3) F(3); A(3) D(3) H(3) E(3)]';
Co=[163/255,163/255,163/255];

Xy1=xy+X;
Yy1=yy1+Y;
Zy1=zy+Z;


A=rotatevect([0.5 0.5 0.5].*scal,ang2,y);
B=rotatevect([-0.5 0.5 0.5].*scal,ang2,y);
C=rotatevect([-0.5 -0.5 0.5].*scal,ang2,y);
D=rotatevect([0.5 -0.5 0.5].*scal,ang2,y);
E=rotatevect([0.5 0.5 -0.5].*scal,ang2,y);
F=rotatevect([-0.5 0.5 -0.5].*scal,ang2,y);
G=rotatevect([-0.5 -0.5 -0.5].*scal,ang2,y);
H=rotatevect([0.5 -0.5 -0.5].*scal,ang2,y);

X=[A(1) B(1) C(1) D(1); D(1) C(1) G(1) H(1); G(1) H(1) E(1) F(1); E(1) F(1) B(1) A(1); B(1) C(1) G(1) F(1); A(1) D(1) H(1) E(1)]';
Y=[A(2) B(2) C(2) D(2); D(2) C(2) G(2) H(2); G(2) H(2) E(2) F(2); E(2) F(2) B(2) A(2); B(2) C(2) G(2) F(2); A(2) D(2) H(2) E(2)]';
Z=[A(3) B(3) C(3) D(3); D(3) C(3) G(3) H(3); G(3) H(3) E(3) F(3); E(3) F(3) B(3) A(3); B(3) C(3) G(3) F(3); A(3) D(3) H(3) E(3)]';

Xy2=xy+X;
Yy2=yy2-Y;
Zy2=zy+Z;

fill3(Xy2,Yy2,Zy2,Co,'FaceAlpha',1)
fill3(Xy1,Yy1,Zy1,Co,'FaceAlpha',1)
hold on

%% Solar array
Lp=650;
Wp=5;
Hp=205;

xp= xcg_b;
yp1= ycg_b+L/2+Ly+Lp/2;
yp2=-yp1;
zp= zcg_b;


y=[0 1 0];

scal=[Wp Lp Hp];
A=rotatevect([0.5 0.5 0.5].*scal,ang1,y);
B=rotatevect([(-0.5) 0.5 0.5].*scal,ang1,y);
C=rotatevect([-0.5 -0.5 0.5].*scal,ang1,y);
D=rotatevect([0.5 -0.5 0.5].*scal,ang1,y);
E=rotatevect([0.5 0.5 -0.5].*scal,ang1,y);
F=rotatevect([-0.5 0.5 -0.5].*scal,ang1,y);
G=rotatevect([-0.5 -0.5 -0.5].*scal,ang1,y);
H=rotatevect([0.5 -0.5 -0.5].*scal,ang1,y);

X=[A(1) B(1) C(1) D(1); D(1) C(1) G(1) H(1); G(1) H(1) E(1) F(1); E(1) F(1) B(1) A(1); B(1) C(1) G(1) F(1); A(1) D(1) H(1) E(1)]';
Y=[A(2) B(2) C(2) D(2); D(2) C(2) G(2) H(2); G(2) H(2) E(2) F(2); E(2) F(2) B(2) A(2); B(2) C(2) G(2) F(2); A(2) D(2) H(2) E(2)]';
Z=[A(3) B(3) C(3) D(3); D(3) C(3) G(3) H(3); G(3) H(3) E(3) F(3); E(3) F(3) B(3) A(3); B(3) C(3) G(3) F(3); A(3) D(3) H(3) E(3)]';
Co='b';


Xp1=xp+X;
Yp1=yp1+Y;
Zp1=zp+Z;


A=rotatevect([0.5 0.5 0.5].*scal,ang2,y);
B=rotatevect([(-0.5) 0.5 0.5].*scal,ang2,y);
C=rotatevect([-0.5 -0.5 0.5].*scal,ang2,y);
D=rotatevect([0.5 -0.5 0.5].*scal,ang2,y);
E=rotatevect([0.5 0.5 -0.5].*scal,ang2,y);
F=rotatevect([-0.5 0.5 -0.5].*scal,ang2,y);
G=rotatevect([-0.5 -0.5 -0.5].*scal,ang2,y);
H=rotatevect([0.5 -0.5 -0.5].*scal,ang2,y);

X=[A(1) B(1) C(1) D(1); D(1) C(1) G(1) H(1); G(1) H(1) E(1) F(1); E(1) F(1) B(1) A(1); B(1) C(1) G(1) F(1); A(1) D(1) H(1) E(1)]';
Y=[A(2) B(2) C(2) D(2); D(2) C(2) G(2) H(2); G(2) H(2) E(2) F(2); E(2) F(2) B(2) A(2); B(2) C(2) G(2) F(2); A(2) D(2) H(2) E(2)]';
Z=[A(3) B(3) C(3) D(3); D(3) C(3) G(3) H(3); G(3) H(3) E(3) F(3); E(3) F(3) B(3) A(3); B(3) C(3) G(3) F(3); A(3) D(3) H(3) E(3)]';

Xp2=xp+X;
Yp2=yp2-Y;
Zp2=zp+Z;

fill3(Xp2,Yp2,Zp2,Co,'FaceAlpha',0.75)
fill3(Xp1,Yp1,Zp1,Co,'FaceAlpha',0.75)
hold on


%% FUNCTIONS 
function [vect2]=rotatevect(vect1,theta,k)
%%RODRIGUES'S FORMULA
vect2=vect1*cosd(theta)+cross(k,vect1)*sind(theta)+k*dot(k,vect1)*(1-cosd(theta));
end
end