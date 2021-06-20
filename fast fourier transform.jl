using FFTW
using DelimitedFiles

### Problem 1

N=128
L=2.0
x=collect(-div(N,2):div(N,2)-1)*(L/N)
f=x.^2
println(fft(f)./N)

ssss

### Problem 2
println(" ")
println("***Let's now start problem 2.***")
N=16
L=2.0*pi
x=collect(-div(N,2):div(N,2)-1)*(L/N)
g=cos.(2.0*x)
h=cos.(14.0*x)
println(fft(g)./N)
println(fft(h)./N)

### Problem 3 create data file

N=1024
L=20.0
x=collect(-div(N,2):div(N,2)-1)*(L/N)
y=zeros(N)
for k=1:2:9
    y.+=sqrt(k)*sin.(k*pi*x)
end
for k=2:2:10
    y.+=sqrt(k)*cos.(k*pi*x)
end
y.+=2021.0
println(y)
yHat=fft(y)./N
for j=1:N
    if abs(real.(yHat[j]))<1.0E-12
        yHat[j]=1.0im*imag(yHat[j])
    end
    if abs(imag.(yHat[j]))<1.0E-12
        yHat[j]=1.0*real(yHat[j])
    end
end
myFile=open("FourierData.out","w")
writedlm(myFile,yHat)
close(myFile)

### Problem 3 Solution

println(" ")
println("***Let's now start problem 3.***")

myFile=open("FourierData.out","r")
mHat=readdlm(myFile,'\t',Complex{Float64},'\n')
close(myFile)
m=ifft(mHat)*length(mHat)
println(real(m))
