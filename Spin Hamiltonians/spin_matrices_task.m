run Physical_Constants.m
Sx=0.5.*[0,1;1,0];
Sy=-0.5i.*[0,1;-1,0];
Sz=0.5.*[1,0;0,-1];
Splus_single=Sx+1i.*Sy;
Sminus_single=Sx-1i*Sy;
function t=Commutator(mat1,mat2)
    t=mat1*mat2-mat2*mat1;
end
%Commutator(Sx,Sy)==1i*Sz
%Commutator(Sy,Sz)==1i*Sx
%Commutator(Sz,Sx)==1i*Sy
%Commutator(Sz,(Sx*Sx+Sy*Sy+Sz*Sz))
%Commutator(Sz,Splus_single)
%Commutator(Sz,Sminus_single)
%Commutator(Splus_single,Sminus_single)
E=[1,0;0,1];
S1x=kron(Sx,E);
S2x=kron(E,Sx);
S1y=kron(Sy,E);
S2y=kron(E,Sy);
S1z=kron(Sz,E);
S2z=kron(E,Sz);
S1plus=[0,0,1,0;0,0,0,1;0,0,0,0;0,0,0,0];
S2plus=[0,1,0,0;0,0,0,0;0,0,0,1;0,0,0,0];
S1minus=[0,0,0,0;0,0,0,0;1,0,0,0;0,1,0,0];
S2minus=[0,0,0,0;1,0,0,0;0,0,0,0;0,0,1,0];
Splus=S1plus+S2plus;
Sminus=S1minus+S2minus;
S1xS2x=S1x*S2x;
S1yS2y=S1y*S2y;
S1zS2z=S1z*S2z;
S1S2=(S1x+S1y+S1z)*(S2x+S2y+S2z);
J=5E6;
B=100E-3;
g1=2;
g2=2.2;
field_independent=J*hbar.*(S1xS2x+S1yS2y+S1zS2z)/planck;
field_dependent=bmagn.*B.*(g1.*S1z+g2.*S2z)/planck;
H=field_independent+field_dependent;
[vectors,values]=eig(H);
vectors;
values;
for i=1:4
    energy(i)=values(i,i);
end
energy;
energy(4)-energy(2);%E42
energy(4)-energy(3);%E43
energy(3)-energy(1);%E31
energy(2)-energy(1);%E21
E42=zeros(1,100);
E43=zeros(1,100);
E31=zeros(1,100);
E21=zeros(1,100);
B=zeros(1,100);
for b=1:100
    field_dependent=bmagn.*10*b*1E-3.*(g1.*S1z+g2.*S2z)/planck;
    H=field_independent+field_dependent;
    [vectors,values]=eig(H);
    for i=1:4
        energy(i)=values(i,i);
    end
    B(b)=10*b;
    E42(b)=energy(4)-energy(2);
    E43(b)=energy(4)-energy(3);
    E31(b)=energy(3)-energy(1);
    E21(b)=energy(2)-energy(1);
end

%plot(B,E42,'r');%same as E31 as expected
hold on;
%plot(B,E43,'g');%same as E21 as expected
hold on;
%plot(B,E31,'y');
hold on;
%plot(B,E21,'b');

B42=interp1(E42,B,9E9);
B43=interp1(E43,B,9E9);
B31=interp1(E31,B,9E9);
B21=interp1(E21,B,9E9);%in mT
%Field=B42
field_independent=J*hbar.*(S1xS2x+S1yS2y+S1zS2z)/planck;
field_dependent=bmagn.*Field*1E-3.*(g1.*S1z+g2.*S2z)/planck;
H=field_independent+field_dependent;
[vectors,values]=eig(H);
vectors;
values;
H;
vectors(:,1)%column
vectors(1,:)%row
%vectors(n,:)*H*vectors(:,m) gives nHm where n,m=1,4
%a42
%compute eigenvectors with B42
Field=B42;
field_independent=J*hbar.*(S1xS2x+S1yS2y+S1zS2z)/planck;
field_dependent=bmagn.*Field*1E-3.*(g1.*S1z+g2.*S2z)/planck;
H=field_independent+field_dependent;
[vectors,values]=eig(H);
a42=vectors(4,:)*(S1x+S2x)*vectors(:,2);

Field=B43;
field_independent=J*hbar.*(S1xS2x+S1yS2y+S1zS2z)/planck;
field_dependent=bmagn.*Field*1E-3.*(g1.*S1z+g2.*S2z)/planck;
H=field_independent+field_dependent;
[vectors,values]=eig(H);
a43=vectors(4,:)*(S1x+S2x)*vectors(:,3);

Field=B31;
field_independent=J*hbar.*(S1xS2x+S1yS2y+S1zS2z)/planck;
field_dependent=bmagn.*Field*1E-3.*(g1.*S1z+g2.*S2z)/planck;
H=field_independent+field_dependent;
[vectors,values]=eig(H);
a31=vectors(3,:)*(S1x+S2x)*vectors(:,1);

Field=B21;
field_independent=J*hbar.*(S1xS2x+S1yS2y+S1zS2z)/planck;
field_dependent=bmagn.*Field.*(g1.*S1z+g2.*S2z)/planck;
H=field_independent+field_dependent;
[vectors,values]=eig(H);
a21=vectors(2,:)*(S1x+S2x)*vectors(:,1);

a43^2+a42^2+a31^2+a21^2;
new_B=0.1:1:1000;
%plot(new_B,a43^2*LS_lorentzian(new_B,B43,5,1));
hold on;
%plot(new_B,a42^2*LS_lorentzian(new_B,B42,5,1));
hold on;
plot(new_B,a31^2*LS_lorentzian(new_B,B31,5,1)+a21^2*LS_lorentzian(new_B,B21,5,1),'b');
hold on;
%plot(new_B,a21^2*LS_lorentzian(new_B,B21,5,1),'b');