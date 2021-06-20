using LinearAlgebra

function NewtonMethodSystem(x,N)
    for j=1:N
        x=x-J(x)\F(x)
        println("j=",j,"   x=",x,"   |F(x)|=",norm(F(x)))
    end
end

function F(x)
    return [x[1]^2+x[2]^2-2*x[1]-2*x[2]+1; x[1]+x[2]-2*x[1]*x[2]]
end

function J(x)
    return [2*x[1]-2 2*x[2]-2; 1-2*x[2] 1-2*x[1]]
end

function Broyden(x,A,N)
    for j=1:N
        xNew=x-A\F(x)
        s=xNew-x
        y=F(xNew)-F(x)
        A=A+((y-A*s)/norm(s)^2)*s'
        x=xNew
        println("j=",j,"   x=",x,"   |F(x)|=",norm(F(x)))
    end
end

function BFGS(x,AInv,N)
    for j=1:N
        xNew=x-AInv*F(x)
        s=xNew-x
        y=F(xNew)-F(x)
        AInv=AInv+((s'*y+y'*AInv*y)/(s'*y)^2)*(s*s')-(AInv*y*s'+s*y'*AInv)/(s'*y)
        x=xNew
        println("j=",j,"   x=",x,"   |F(x)|=",norm(F(x)))
    end
end

NewtonMethodSystem([1.9;0.6],6)
println("****************************")
Broyden([1.0,0.6],[1.0 0.0;0.0 1.0],13)
println("****************************")
BFGS([1.9;0.6],[1.0 0.0;0.0 1.0],14)
