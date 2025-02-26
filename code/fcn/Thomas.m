function [q] = thomas(aa,bb,cc,dd)
bet(1)=bb(1);
gam(1)=dd(1)/bb(1);
n=length(bb);
for i=2:n
bet(i)=bb(i)-(aa(i)*cc(i-1)/bet(i-1));
gam(i)=(dd(i)-aa(i)*gam(i-1))/bet(i);
end 
q(n)=gam(n);
for i=n-1:-1:1
   q(i)=gam(i)-(cc(i)*q(i+1)/bet(i));   
end 
