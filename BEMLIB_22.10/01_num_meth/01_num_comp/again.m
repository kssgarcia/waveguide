function again(level)

%=========
% demonstration of recursive function calling
%========

if(level>0)
 level = level-1;
 level
 again(level);
 level
end

%---
% done
%---

return
