# Problem 1

using LinearAlgebra
A=[-2 0; 2 -2; 0 2]
b=[1;2;3]

M,N=size(A)
if N>M
    ATA=transpose(A)*A
    SingVals=reverse(sqrt.(eigvals(ATA)))
    V=reverse(eigvecs(ATA),dims=2)
    Sigma=zeros(Float64,size(A))
    SigmaCross=zeros(Float64,size(transpose(A)))
    U=zeros(Float64,M,M)
    for k=1:M
        U[:,k]=A*V[:,k]/SingVals[k]
        Sigma[k,k]=SingVals[k]
        SigmaCross[k,k]=1.0/SingVals[k]
    end
else
    AAT=A*transpose(A)
    SingVals=reverse(sqrt.(eigvals(AAT)))
    U=reverse(eigvecs(AAT),dims=2)
    Sigma=zeros(Float64,size(A))
    SigmaCross=zeros(Float64,size(transpose(A)))
    V=zeros(Float64,N,N)
    for k=1:N
        V[:,k]=transpose(U[:,k])*A/SingVals[k]
        Sigma[k,k]=SingVals[k]
        SigmaCross[k,k]=1.0/SingVals[k]
    end
end

println("1(a) The system does not have a solution because it is overdetermined.")
println("1(b) The singular value decomposition has the following components")
println("U=",U)
println("Sigma=",Sigma)
println("V=",V)
println(" ")
println("1(c) The least squares solution is: ",V*SigmaCross*transpose(U)*b)
println("1(d) The minimal value of the residual is: ",norm(A*V*SigmaCross*transpose(U)*b-b)^2)


# Problem 2

function MySVD(A)
    ATA=transpose(A)*A
# Julia orders the eigenvalues from smallest to largest, so their order needs to be switched
# when they are put into the matrix Sigma.
# In order to make sure the signs of the eigenvectors line up we define U in terms of V
    SingVals=reverse(sqrt.(abs.(eigvals(ATA))))
    V=reverse(eigvecs(ATA),dims=2)
    Sigma=zeros(Float64,size(A))
    U=zeros(Float64,size(A)[1],size(A)[1])
    for k=1:size(A)[2]
        U[:,k]=A*V[:,k]/SingVals[k]
        Sigma[k,k]=SingVals[k]
    end
    return U,Sigma,V
end

# Problem 3

println(" ")
println(" ")

function MyTruncatedSVD(A,s)
    ATA=transpose(A)*A
    SingVals=reverse(sqrt.(abs.(eigvals(ATA))))
    V=reverse(eigvecs(ATA),dims=2)
    Sigma=zeros(Float64,size(A))
    U=zeros(Float64,size(A)[1],size(A)[1])
    counter=0
    for k=1:size(A)[2]
        U[:,k]=A*V[:,k]/SingVals[k]
        Sigma[k,k]=SingVals[k]
        if SingVals[k]>s
            counter+=1
        end
    end
    UTrunc=U[:,1:counter]
    SigmaTrunc=Sigma[1:counter,1:counter]
    VTrunc=V[:,1:counter]
    return UTrunc,SigmaTrunc,VTrunc
end

MyMatrixA=[1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
s=1.0
TruncU,TruncSigma,TruncV=MyTruncatedSVD(MyMatrixA,s)
println("Problem 3")
println("The original matrix (as an example) is: ",MyMatrixA)
println("The singular values are: ",reverse(sqrt.(abs.(eigvals(transpose(MyMatrixA)*MyMatrixA)))))
println("The SVDtruncated matrix with s=",s," is: ",TruncU*TruncSigma*transpose(TruncV))
s=1.5
TruncU,TruncSigma,TruncV=MyTruncatedSVD(MyMatrixA,s)
println(TruncU)
ppp
println("The SVDtruncated matrix with s=",s," is: ",TruncU*TruncSigma*transpose(TruncV))
println("Note that using s=1 doesn't set any singular values to zero, so the SVD truncated matrix is the same as the original")
println("When s=1.5, one singular value gets set to zero and that causes the SVD truncated matrix to be different than the original")

println(" ")
println(" ")

# Problem 4

using FileIO
using Images

function MyImageTruncatedSVD(A,NumberDeleted)
    ATA=transpose(A)*A
    SingVals=reverse(sqrt.(abs.(eigvals(ATA))))
    V=reverse(eigvecs(ATA),dims=2)
    Sigma=zeros(Float64,size(A))
    U=zeros(Float64,size(A)[1],size(A)[1])
    for k=1:rank(A)
        U[:,k]=A*V[:,k]/SingVals[k]
        Sigma[k,k]=SingVals[k]
    end
    UTrunc=U[:,1:end-NumberDeleted]
    SigmaTrunc=Sigma[1:end-NumberDeleted,1:end-NumberDeleted]
    VTrunc=V[:,1:end-NumberDeleted]
    SingValueRatio=sum(SingVals[1:end-NumberDeleted])/sum(SingVals[1:end])
    return UTrunc,SigmaTrunc,VTrunc,SingValueRatio
end

MyJPG=load("TOG.JPG",view=false)
MyRGB=convert(Array{Float64},channelview(MyJPG))
R=MyRGB[1,:,:];
G=MyRGB[2,:,:];
B=MyRGB[3,:,:];

NumberDeleted=900

RTruncU,RTruncSigma,RTruncV,RSingValueRatio=MyImageTruncatedSVD(R,NumberDeleted)
GTruncU,GTruncSigma,GTruncV,GSingValueRatio=MyImageTruncatedSVD(G,NumberDeleted)
BTruncU,BTruncSigma,BTruncV,BSingValueRatio=MyImageTruncatedSVD(G,NumberDeleted)

RTrunc=RTruncU*RTruncSigma*transpose(RTruncV)
GTrunc=GTruncU*GTruncSigma*transpose(GTruncV)
BTrunc=BTruncU*BTruncSigma*transpose(BTruncV)

TruncImage=zeros(Float64,size(MyRGB));
TruncImage[1,:,:]=RTrunc
TruncImage[2,:,:]=GTrunc
TruncImage[3,:,:]=BTrunc

# Thanks to Gabriel for showing us the following command!
save("TruncatedTOG.jpg", map(clamp01nan, colorview(RGB,TruncImage)))

println("I discared 900 of 976 singular values and the image still looks decent.")
println("The percent of singular valuesness I kept in the red channel was: ",RSingValueRatio*100)
println("The percent of singular valuesness I kept in the green channel was: ",GSingValueRatio*100)
println("The percent of singular valuesness I kept in the blue channel was: ",BSingValueRatio*100)

