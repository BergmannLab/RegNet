function f=fkjellstrom2(x)
  [n,m] = size(x);
  f=zeros(m,1);
  h = 0.01.*((cos(x+1.982) + cos(2*x+5.720) + cos(3*x + 1.621) + cos(4*x + 0.823) + cos(5*x + 3.222))); 
  for i=1:m
    f(i) = prod(1+h(:,i));
  end
  
  