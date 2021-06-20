function RER(a,b,N)
    deltax=(b-a)/N
    x=a+deltax
    sum=0.0
    for j=1:N
        sum+=f(x)
        x+=deltax
    end
    return sum*deltax
end

function TrapR(a,b,N)
    deltax=(b-a)/N
    x=a+deltax
    sum=0.0
    for j=1:N-1
        sum+=f(x)
        x+=deltax
    end
    return (sum+0.5*(f(a)+f(b)))*deltax
end

function TrapR2D(a,b,c,d,Nx,Ny)
    deltax=(b-a)/Nx
    deltay=(d-c)/Ny
    sum=0.0
    x=a
    for j=1:Nx-1
        x+=deltax
        sum+=0.5*(f2D(x,c)+f2D(x,d))
        y=c
        for k=1:Ny-1
            y+=deltay
            sum+=f2D(x,y)
        end
    end
    y=c
    for k=1:Ny-1
        y+=deltay
        sum+=0.5*(f2D(a,y)+f2D(b,y))
    end
    return(0.25*(f2D(a,c)+f2D(b,c)+f2D(a,d)+f2D(b,d))+sum)*deltax*deltay
end

function SimpR(a,b,N)
    deltax=(b-a)/N
    x=a+deltax
    simp2=0.0
    simp4=0.0
    for j=1:N/2-1
        simp4+=f(x)
        x+=deltax
        simp2+=f(x)
        x+=deltax
    end
    simp4+=f(x)
    return (2.0*simp2+4.0*simp4+f(a)+f(b))*deltax/3.0
end

function MCInt(a,b,N)
    bminusa=b-a
    sum=0.0
    for j=1:N
        sum+=f(a+bminusa*rand())
    end
    return sum*bminusa/N
end

function MCInt2D(a,b,c,d,Nx,Ny)
    bminusa=b-a
    dminusc=d-c
    sum=0.0
    for j=1:Nx
        for k=1:Ny
            sum+=f2D(a+bminusa*rand(),c+dminusc*rand())
        end
    end
    return sum*bminusa*dminusc/Nx/Ny
end

function f2D(x,y)
    return sin(x^2+y^3)
end

#### Code starts here
#### Problem 2a
function f(x)
    return tan(x)
end
a=0.0
b=1.0
exactans=-log(cos(1.0))

println("  ")
println("Problem 2a: Value Table")
for N=[10,100,1000,10000,100000,1000000]
    RERans=RER(a,b,N)
    TrapRans=TrapR(a,b,N)
    SimpRans=SimpR(a,b,N)
    MCIntans=MCInt(a,b,N)    
    println(N,"  ",[RERans,TrapRans,SimpRans,MCIntans])
end
    
println("***********************************************")
println("Problem 2a: Error Table")
for N=[10,100,1000,10000,100000,1000000]
    RERans=RER(a,b,N)
    TrapRans=TrapR(a,b,N)
    SimpRans=SimpR(a,b,N)
    MCIntans=MCInt(a,b,N)
    println(N,"  ",abs.([RERans,TrapRans,SimpRans,MCIntans].-exactans))
end

#### Problem 2b
function f(x)
    return (abs(x-sqrt(2)))^(1/3)
end
a=-5.0
b=5.0
exactans=3/4*((5-sqrt(2))^(4/3)+(5+sqrt(2))^(4/3))

println("  ")
println("Problem 2b: Value Table")
for N=[10,100,1000,10000,100000,1000000]
    RERans=RER(a,b,N)
    TrapRans=TrapR(a,b,N)
    SimpRans=SimpR(a,b,N)
    MCIntans=MCInt(a,b,N)    
    println(N,"  ",[RERans,TrapRans,SimpRans,MCIntans])
end
    
println("***********************************************")
println("Problem 2b: Error Table")
for N=[10,100,1000,10000,100000,1000000]
    RERans=RER(a,b,N)
    TrapRans=TrapR(a,b,N)
    SimpRans=SimpR(a,b,N)
    MCIntans=MCInt(a,b,N)
    println(N,"  ",abs.([RERans,TrapRans,SimpRans,MCIntans].-exactans))
end

#### Problem 2c
function f(x)
    return sin(2*pi*cos(x))
end
a=-pi
b=pi
exactans=0.0

println("  ")
println("Problem 2c: Value Table")
for N=[10,100,1000,10000,100000,1000000]
    RERans=RER(a,b,N)
    TrapRans=TrapR(a,b,N)
    SimpRans=SimpR(a,b,N)
    MCIntans=MCInt(a,b,N)    
    println(N,"  ",[RERans,TrapRans,SimpRans,MCIntans])
end
    
println("***********************************************")
println("Problem 2c: Error Table")
for N=[10,100,1000,10000,100000,1000000]
    RERans=RER(a,b,N)
    TrapRans=TrapR(a,b,N)
    SimpRans=SimpR(a,b,N)
    MCIntans=MCInt(a,b,N)
    println(N,"  ",abs.([RERans,TrapRans,SimpRans,MCIntans].-exactans))
end

#### Problem 3
a=2.0
b=3.0
c=0.0
d=1.0
exactans=0.027368544838490438156

println("  ")
println("Problem 3: Value Table")
for N=[10,100,1000,10000]
    TrapR2Dans=TrapR2D(a,b,c,d,N,N)
    MCInt2Dans=MCInt2D(a,b,c,d,N,N)
    println(N,"  ",[TrapR2Dans,MCInt2Dans])
end

println("***********************************************")
println("Problem 3: Error Table")
for N=[10,100,1000,10000]
    TrapR2Dans=TrapR2D(a,b,c,d,N,N)
    MCInt2Dans=MCInt2D(a,b,c,d,N,N)
    println(N,"  ",abs.([TrapR2Dans,MCInt2Dans].-exactans))
end
