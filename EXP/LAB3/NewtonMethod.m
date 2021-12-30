%可调参数
eps = 1e-12;
alpha_=1;
rho=0.25;
sigma=0.5;
x=[-1,0,0,0,0];
%NewtonMethod
itr = 0;
while itr<1000
    t = sym('t',size(x));
    hess = hessian(f(t),t);
    grd = gradient(f(t));
    if norm(double(subs(grd,t,x)))<eps
        break
    end
    d = -double(subs(hess,t,x))\double(subs(grd,t,x));
    d = d';
    % WolfePowell非精确一维搜索
    alpha = WolfePowell(alpha_,rho,sigma,@f,x,d);
    x = x+alpha.*d;
    itr = itr+1;
    if mod(itr,100)==0
        itr
    end
end
itr
x
f(x)