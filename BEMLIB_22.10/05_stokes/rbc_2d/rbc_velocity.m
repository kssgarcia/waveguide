function [U,V] = rbc_velocity ...
 ...
  (NSG,X,Y,Xint ...
  ,Axint,Bxint,Cxint ...
  ,Ayint,Byint,Cyint ...
  ,Dfx,Dfy ...
  ,vnx,vny)

%==================================
% solve the integral equation for the
% interfacial velocity
%==================================

%---
% prepare
%---

global NGL
global Iflow
global wall
global shrt
global mu1
global mu2

%---
% equal viscosities?
%---

Iflag = 1;

if(abs(mu1-mu2)>0.00001*mu1)
  Iflag = 2;
end

if(Iflag==2)
 MAT=zeros(2*NSG,2*NSG);
end
  
%====================================
% compute the potentials at the nodes
%====================================

%---
for node=1:NSG
%---

 U(node) = 0.0;
 V(node) = 0.0;

 Dfx0 = Dfx(node);
 Dfy0 = Dfy(node);

%---
 for j=1:NSG  % run over segments
%---

 j1 = j+1;

 Ising = 0;
 if(node==j) Ising = 1; end           % first end-point
 if(node==j+1) Ising = 2; end         % second end-point
 if(node==1 & j==NSG) Ising = 2; end   % second end-point

 %~~~
 if(Iflag==1) % equal viscosities
 %~~~

 [Wx,Wy] = rbc_slp_spline  ...
 ...
  (X(node),Y(node) ...
  ,X(j),Y(j)...
  ,X(j1),Y(j1)...
  ,Xint(j),Xint(j1) ...
  ,Axint(j),Bxint(j),Cxint(j) ...
  ,Ayint(j),Byint(j),Cyint(j) ...
  ,Dfx0,Dfy0 ...
  ,Dfx(j),Dfx(j1) ...
  ,Dfy(j),Dfy(j1) ...
  ,Ising);

 %~~~
 elseif(Iflag==2)
 %~~~

  [Wx,Wy ...
  ,DLPxx1, DLPxy1 ...
  ,DLPyx1, DLPyy1 ...
  ,DLPxx2, DLPxy2 ...
  ,DLPyx2, DLPyy2 ...
  ...
  ] = rbc_sdlp_spline  ...
  ...
  (X(node),Y(node) ...
  ,X(j),Y(j)...
  ,X(j1),Y(j1)...
  ,Xint(j),Xint(j1) ...
  ,Axint(j),Bxint(j),Cxint(j) ...
  ,Ayint(j),Byint(j),Cyint(j) ...
  ,Dfx0,Dfy0 ...
  ,Dfx(j),Dfx(j1) ...
  ,Dfy(j),Dfy(j1) ...
  ,vnx(j),vnx(j1) ...
  ,vny(j),vny(j1) ...
  ,Ising);

  nNSG = node+NSG;

  MAT(node,j    ) = MAT(node,j    ) + DLPxx1;
  MAT(node,j+NSG) = MAT(node,j+NSG) + DLPyx1;
  MAT(nNSG,j    ) = MAT(nNSG,j    ) + DLPxy1;
  MAT(nNSG,j+NSG) = MAT(nNSG,j+NSG) + DLPyy1;

  k = j1;

  if(j==NSG) % point NSG+1 is mapped to point 1
   k = 1;
  end

   MAT(node,k    ) = MAT(node,k    ) + DLPxx2;
   MAT(node,k+NSG) = MAT(node,k+NSG) + DLPyx2;
   MAT(nNSG,k    ) = MAT(nNSG,k    ) + DLPxy2;
   MAT(nNSG,k+NSG) = MAT(nNSG,k+NSG) + DLPyy2;

 %~~~
  end
 %~~~

  U(node) = U(node)+Wx;
  V(node) = V(node)+Wy;

%---
  end  % over segments
%---

  U(node) = -U(node)/(4*pi*mu1);
  V(node) = -V(node)/(4*pi*mu1);

%---
end  % over nodes
%---

%---
% add the incident flow
%---

 for i=1:NSG

  if(Iflow==1)
   U(i) = U(i)+shrt*Y(i);
  elseif(Iflow==2)
   U(i) = U(i)+shrt*(Y(i)-wall);
  elseif(Iflow==3)
   U(i) = U(i)+shrt*(Y(i)-wall);
  end
 end

 %~~~
 if(Iflag==2)
 %~~~

 lambda = mu2/mu1;
 MAT = -(1.0-lambda)*MAT/(4.0*pi);

 for i=1:2*NSG
  MAT(i,i) = MAT(i,i)+0.5*(1.0+lambda);
 end

 for i=1:NSG
  RHS(i) = U(i);
  RHS(i+NSG) = V(i);
 end

%MAT=eye(2*NSG,2*NSG);

 SOL = RHS/MAT';

 for i=1:NSG
  U(i) = SOL(i);
  V(i) = SOL(i+NSG);
 end

 %~~~
 end
 %~~~

%---
% wrap
%---

U(NSG+1)=U(1);
V(NSG+1)=V(1);

%---
% done
%---

return
