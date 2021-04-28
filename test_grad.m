Z = Power_temp;
m = J' * mm;

x  = inv(conj(Z)) *conj(m);
f(x,Z,m)
df(x,Z,m)

xt = rand(30,1)+2*rand(30,1)*1j;
f(xt,Z,m)


function res = f(x,Z,m)
    res = -conj(x)'*Z*conj(x)+2*real(m'*conj(x));
end

function res = df(x,Z,m)
    res = norm(-2*conj(Z)*x+2*conj(m));
end


