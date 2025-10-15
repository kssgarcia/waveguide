function [f,s,q] = flsphere_ode(ndiv,smax,capls,xcl,scl,slope);

%===============================
% solve ODEs for a semi-infinite
% axisymmetric meniscus
% using RK4
%
% scl: sigma of contact line
%
% will integrate from scl to smax
%================================

dsg = (smax-scl)/ndiv;

dsgh = 0.5*dsg;

s(1) = scl;  % starting point
f(1) = xcl;
q(1) = slope;

%---
% integrate
%---

for i=1:ndiv

   fp = q(i);
   tmp = 1.0+q(i)*q(i);
   tmg = sqrt(tmp)/capls;
   qp = tmp*( -fp/s(i) + f(i)*tmg );
   fp1 = fp;
   qp1 = qp;

   s(i+1) = s(i)+   dsgh;
   f(i+1) = f(i)+fp*dsgh;
   q(i+1) = q(i)+qp*dsgh;

   fp = q(i+1);
   tmp = 1.0+q(i+1)*q(i+1);
   tmg = sqrt(tmp)/capls;
   qp= tmp*(-fp/s(i+1) + f(i+1)*tmg);
   fp2 = fp;
   qp2 = qp;

   f(i+1) = f(i)+fp*dsgh;
   q(i+1) = q(i)+qp*dsgh;

   fp = q(i+1);
   tmp = 1.0+q(i+1)*q(i+1);
   tmg = sqrt(tmp)/capls;
   qp = tmp*(-fp/s(i+1) + f(i+1)*tmg);
   fp3 = fp;
   qp3 = qp;

   s(i+1) = s(i)+   dsg;
   f(i+1) = f(i)+fp*dsg;
   q(i+1) = q(i)+qp*dsg;

   fp = q(i+1);
   tmp = 1.0+q(i+1)*q(i+1);
   tmg = sqrt(tmp)/capls;
   qp = tmp*(-fp/s(i+1) + f(i+1)*tmg);
   fp4 = fp;
   qp4 = qp;

   f(i+1) = f(i) + (fp1+2*fp2+2*fp3+fp4)*dsg/6.0;
   q(i+1) = q(i) + (qp1+2*qp2+2*qp3+qp4)*dsg/6.0;

end

%---
% done
%---

return
