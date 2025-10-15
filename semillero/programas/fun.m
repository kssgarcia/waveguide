function  u = fun (t,x)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

u =0;
nk=50;
for k=0:nk 
    ter = (4/pi^2) * ((-1)^k/(2*k+1)^2) * cos((2*k+1)*pi*t) * sin((2*k+1)*pi*x); 
    u=u+ter ;
end

end