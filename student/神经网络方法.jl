
using GR
include("参数以及定义.jl")


# 神经网络方法
using Flux
m=Dense(2,1)                #神经网络模型，这里仅用到一个神经元，有两个输入，一个输出。未设置激活函数。
loss(x,y)=Flux.mse(m(x),y)  #损失函数即l=(y-ŷ)^2,ŷ=m(x)是神经网络的输出。
opt = ADAM(0.0001)          #神经网络训练算法
# 从已知的波动方程数据中学习y=cos(ωt-kx+ϕ)
tspan=0:τ:2000τ
Loss=Float64[]
# 训练神经网络
for n=1:100000
    ϕ=rand()*π                  #随机设置初始相位
    t=rand(tspan)               #随机选择时间
    x=rand(X)                   #随机选择x坐标
    y1=cos(ω*t-k*x+ϕ)           #t时刻质点x的位移
    y2=cos(ω*t-k*(x+h)+ϕ)       #t时刻质点x+h的位移
    y3=cos(ω*(t+τ)-k*(x+2h)+ϕ)  #t+τ时刻质点x+2h的位移，监督学习中，这相当于答案
    data=[([y1,y2],y3)]         #神经网络输入数据是y1和y2，目标数据是y3。
    Flux.train!(loss,params(m),data,opt)
    if n%100==0 
        l=loss([y1,y2],y3).data[]
        println("step=",n," loss=",l)
        push!(Loss,l) #当n等于100的整数倍时，输出神经网络的损失函数值
    end
end
plot(Loss) #训练过程中损失函数的变化情况
#训练完毕后神经网络m就可以用于预测了

E=Float64[]
function absorb_nn()
    y0=zeros(N)                             
    y1=zeros(N)                           
    y =zeros(N)
    l=y1[N-2]
    r=y1[N-1]
    for t=0:τ:32T
        y[2:N-1]=c*(y1[3:N]+y1[1:N-2])+2*(1-c)*y1[2:N-1]-y0[2:N-1]
        y[1]=gsp(t)
        y[N]=m([l,r]).data[] #利用神经网络预测边界值
        y0 .= y1
        y1 .= y
        l=y1[N-2]
        r=y1[N-1]
        plot(X,y)
        xlim([0,L])
        ylim([-2,2])

        v=y1 .- y0
        Ek=sum(v.^2)            # 波的动能计算 
        Ep=c*sum(diff(y).^2)    # 势能计算  

        push!(E,Ek+Ep) #波的无量纲能量
    end
end
absorb_nn()
ylim([0,maximum(E)])
plot(E)
入射波能量=E[length(E)÷2]
反射波能量=E[end]
反射率=反射波能量/入射波能量
println("反射率=",反射率)

