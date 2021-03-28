# 变量赋值以及简单运算（教材p3）
A=52.738
a=3.2477e2
B=A+a
a=3.2477e-2
B=A+a
S="Hello"       # matlab使用单引号，此处是区别。

# 复数计算(教材p4)
2*(300-2*im)^2+1/sqrt(5+2*im)

A=[1,3,-1,0,7,2]    #输入向量
B=.~(A.>2)          #A中不大于2的元素的位置，相比matlab多了一个点。
D=A.*B              #A中不大于2的元素
C=(A.>0).& (A.<3)   #A中大于0且小于3的元素位置
E=A.*C

# 例8 （教材8页）
using GR
x=0:0.1:6     
A=[x 4x]
B=sin.(A)
plot(x,B[:,1],x,B[:,2])
plot(x,B[:,1],"r:+")


π                   #圆周率
eps()               #数值最小间隔

#单缝衍射的模拟
using Makie
b = 0.1 #缝宽，毫米
n = 1000 #分成n份
λ = 5.6e-4
θm = π/180
θ=-θm:0.001*θm:θm
I=zeros(length(θ))
for i=1:length(I)
    ϕ = 2π*b*sin(θ[i])/λ
    dϕ = ϕ/n
    A=1/(n+1)
    for j=1:n
        A += cos(j*dϕ)/n
    end
    I[i]=A^2
end
lines(θ,I)
II=repeat(I,1,10)
heatmap(II,colormap=:greys)