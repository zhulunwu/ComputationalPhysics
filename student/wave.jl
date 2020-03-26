# 波动微分方程求解，即Y_tt-u^2*Y_xx=0,u是波速
using GR
# 波的参数
const λ=0.2     #波长
const f=10      #频率
const u=λ*f     #波速
const ω=2π*f    #角频率
const k=2π/λ    #波数
const T=1/f     #周期
const w=2T      #脉冲宽度

const h=λ/25
const α=0.8
const c=α^2
const τ=α*h/u #如此设置是为了方便后面的计算，可以看出，一个步长需要两个时间步来完成。
const L=16λ
const X=0:h:L
const N=length(X)

#脉冲振动曲线
gsp(t)=sin(ω*t)exp(-(t-6T)^2/w^2)
t=0:0.02T:12T
y=gsp.(t)
plot(t,y)   

#模拟1：行波传播，固定边界，反射.
function wave()      
    y0=zeros(N)                             
    y1=zeros(N)                           
    y =zeros(N)
    for t=0:τ:32T 
        y[2:N-1]=c*(y1[3:N]+y1[1:N-2])+2*(1-c)*y1[2:N-1]-y0[2:N-1]
        y[1]=gsp(t)
        y0 .= y1
        y1 .= y
        plot(X,y)
        xlim([0,L])
        ylim([-2,2]) 
    end
end
wave()

#吸收边界方法一：根据波动规律直接获得边界位移值
#根据波的传播规律，右端点的振动位移是左端点传播L/u的结果。
#理论丰满，现实骨感。由于计算存在误差，吸收并不很理想，只是说明原理。
#况且需要知道波的表达式，现实中没有用处。
function absorb_demo()
    y0=zeros(N)                             
    y1=zeros(N)                           
    y =zeros(N)
    for t=0:τ:32T 
        y[2:N-1]=c*(y1[3:N]+y1[1:N-2])+2*(1-c)*y1[2:N-1]-y0[2:N-1]
        y[1]=gsp(t)
        y[N]=gsp(t-L/u)
        y0 .= y1
        y1 .= y
        plot(X,y)
        xlim([0,L])
        ylim([-2,2]) 
    end
end
absorb_demo()

# 外插方法
function absorb_interplot()
    y0=zeros(N)                             
    y1=zeros(N)                           
    y =zeros(N)
    l=y1[N-2]
    r=y1[N-1]
    for t=0:τ:32T
        y[2:N-1]=c*(y1[3:N]+y1[1:N-2])+2*(1-c)*y1[2:N-1]-y0[2:N-1]
        y[1]=gsp(t)
        y[N]=r*(2-α)-l*(1-α)
        y0 .= y1
        y1 .= y
        l=y1[N-2]
        r=y1[N-1]
        plot(X,y)
        xlim([0,L])
        ylim([-2,2])
    end
end
absorb_interplot()

# 神经网络方法
using Flux
m=Dense(2,1)                #神经网络模型，这里仅用到一个神经元，有两个输入，一个输出。未设置激活函数。
loss(x,y)=Flux.mse(m(x),y)  #损失函数即l=(y-ŷ)^2,ŷ=m(x)是神经网络的输出。
opt = ADAM(0.0001)          #神经网络训练算法
# 从已知的波动方程数据中学习y=cos(ωt-kx+ϕ)
tspan=0:τ:2000τ
# 训练神经网络
for n=1:20000
    ϕ=rand()*π                  #随机设置初始相位
    t=rand(tspan)               #随机选择时间
    x=rand(X)                   #随机选择x坐标
    y1=cos(ω*t-k*x+ϕ)           #t时刻质点x的位移
    y2=cos(ω*t-k*(x+h)+ϕ)       #t时刻质点x+h的位移
    y3=cos(ω*(t+τ)-k*(x+2h)+ϕ)  #t+τ时刻质点x+2h的位移，监督学习中，这相当于答案
    data=[([y1,y2],y3)]         #神经网络输入数据是y1和y2，目标数据是y3。
    Flux.train!(loss,params(m),data,opt)
    n%100==0 && @info(loss([y1,y2],y3)) #当n等于100的整数倍时，输出神经网络的损失函数值
end
#训练完毕后神经网络m就可以用于预测了

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
    end
end
absorb_nn()
