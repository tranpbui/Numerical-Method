#initialize
    N=10
    A=rand(N,N);
    B=rand(N,N);
    a=rand(N);
    b=rand(N);

    AInv=inv(A)
    LU=lu(A)
#end initialize

N=1000
A=rand(N,N);
B=rand(N,N);
a=rand(N);
b=rand(N);

println("**************************")
println("Time to add two vectors of length N=",N)
@time a+b
println("Time to multiply a matrix by a vector of length N=",N)
@time A*a
println("Time to multiply two square matrices with size N=",N)
@time A*B

