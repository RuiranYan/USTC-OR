function y = f(x)
%Rosenbrock 函数
    n = length(x);
    y = 0;
    for i=1:n-1
        y = y + (100*(x(i+1)-x(i)^2)^2+(x(i)-1)^2);
    end
end

