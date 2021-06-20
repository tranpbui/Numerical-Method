function BisectionMethod(a,b,N)
    for j=1:N
        c=0.5*(b+a)
        if f(c)*f(a)>0
            a=c
        else
            b=c
        end
        println("j=",j,"   c=",c,"   f(c)=",f(c))
    end
end

function NewtonMethod(x,N)
    for j=1:N
        x=x-f(x)/df(x)
        println("j=",j,"   x=",x,"   f(x)=",f(x))
    end
end

function ModifiedNewtonMethod(x,alpha,N)
    for j=1:N
        x=x-alpha*f(x)/df(x)
        println("j=",j,"   x=",x,"   f(x)=",f(x))
    end
end

function f(x)
    return x^5-2.0
end

function df(x)
    return 5.0*x^4
end

NewtonMethod(2^(1/5)*(0.5+1.0im),10)
