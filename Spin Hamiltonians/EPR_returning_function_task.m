function t=EPR(coupling)
    J=coupling*1E6
    Sx=0.5.*[0,1;1,0];
    Sy=-0.5i.*[0,1;-1,0];
    Sz=0.5.*[1,0;0,-1];
    Splus_single=Sx+1i.*Sy;
    Sminus_single=Sx-1i*Sy;
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
    g1=2;
    g2=2.2;
    bmagn=9.274009995000001e-24;
    hbar=1.054571817000000e-34;
    planck=6.626070150000000e-34;
    %calculate Bnm
    E42=zeros(1,100);
    E43=zeros(1,100);
    E31=zeros(1,100);
    E21=zeros(1,100);
    B=zeros(1,100);
    for b=1:100
        field_independent=J*hbar.*(S1xS2x+S1yS2y+S1zS2z)/planck;
        field_dependent=bmagn.*10*b*1E-3.*(g1.*S1z+g2.*S2z)/planck;
        H=field_independent+field_dependent;
        [vectors,values]=eig(H);
        energy=zeros(1,4);
        for i=1:4
            energy(i)=values(i,i);
        end
        B(b)=10*b;
        E42(b)=energy(4)-energy(2);
        E43(b)=energy(4)-energy(3);
        E31(b)=energy(3)-energy(1);
        E21(b)=energy(2)-energy(1);
    end
    B42=interp1(E42,B,9E9);
    B43=interp1(E43,B,9E9);
    B31=interp1(E31,B,9E9);
    B21=interp1(E21,B,9E9);
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
    t=[B42,B43,B31,B21;a42,a43,a31,a21];
end

M=EPR(2500);

B=100:1:600;
%plot(B,M(2,1)^2*LS_lorentzian(B,M(1,1),5,1));%42
hold on;
%plot(B,M(2,2)^2*LS_lorentzian(B,M(1,2),5,1));%43
hold on;
%plot(B,M(2,3)^2*LS_lorentzian(B,M(1,3),5,1));%31
hold on;
%plot(B,M(2,4)^2*LS_lorentzian(B,M(1,4),5,1));%21

plot(B,(M(2,1)^2*LS_lorentzian(B,M(1,1),5,1)+M(2,2)^2*LS_lorentzian(B,M(1,2),5,1)+M(2,3)^2*LS_lorentzian(B,M(1,3),5,1)+M(2,4)^2*LS_lorentzian(B,M(1,4),5,1)))