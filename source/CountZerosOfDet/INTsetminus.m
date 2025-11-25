function out = INTsetminus(A,B)
%do A\B where B is sorted intervals contained in A

out = [];

b = B(1); 
if length(B) ~= 1
    c = B(2);
else
    c = B(1);
end

if A.inf <= b.inf
    out = infsup(A.inf,b.inf);
else
    out = infsup(b.sup,c.inf);
end

for i = 2:length(B)-1

    b = B(i); c = B(i+1);
    out = [out;infsup(b.sup,c.inf)];
    
end
   
if A.sup >= c.sup
    out = [out;infsup(c.sup,A.sup)];
end


end

