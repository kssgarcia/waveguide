    function erf= erfun(x)

%=============================
% error function
% computed by approximation
%=============================

  t = 1.0D0/(1.0D0+0.3275911*abs(x));

%------------------------------
%  complementary error function
%------------------------------

    erfunc = exp(-x*x)*t*( 0.254829592 ...
                      +t*(-0.284496736 ...
                      +t*( 1.421413741 ...
                      +t*(-1.453152027 ...
                      +t*  1.061405429))));

   if(x< 0.0) erfunc = 2.0-erfunc; end

%----------------
%  error function
%----------------

   erf = 1.0D0-erfunc;

%-----
% done
%-----

return
