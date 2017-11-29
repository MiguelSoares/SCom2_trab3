function y=mary(levels,m,n)
level_values=-(levels-1):2:(levels-1);
x=randsrc(1,m,level_values);
y=x(1)*ones(1,n);
for i=2:m
    y=[y x(i)*ones(1,n)];
end

end