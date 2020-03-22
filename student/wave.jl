# 波动微分方程求解，即Y_tt-u^2*Y_xx=0,u是波速
using GR
# 波的参数
const λ=0.2     #波长
const f=10      #频率
const u=λ*f     #波速
const ω=2π*f    #角频率
const T=1/f     #周期
const w=2T      #脉冲宽度

const h=λ/25
const α=0.5
const c=α^2
const τ=α*h/u
const L=16λ
const x=0:h:L
const N=length(x)

#脉冲振动曲线
gsp(t)=sin(ω*t)exp(-(t-6T)^2/w^2)
t=0:0.02T:12T
y=gsp.(t)
plot(t,y)   

#模拟1：行波传播，固定边界，反射和驻波
function wave()      
    y0=zeros(N)                             
    y1=zeros(N)                           
    y =zeros(N)
    for t=0:τ:32T #增加时间到1000τ可以观察反射现象
        y[2:N-1]=c*(y1[3:N]+y1[1:N-2])+2*(1-c)*y1[2:N-1]-y0[2:N-1]
        y[1]=gsp(t)
        y0 .= y1
        y1 .= y
        plot(x,y)
        xlim([0,L])
        ylim([-2,2]) 
    end
end
wave()

#吸收边界方法一：根据波动规律直接获得边界位移值
#根据波的传播规律，右端点的振动位移是左端点传播L/u的结果
function absorb1()
    y0=zeros(N)                             
    y1=zeros(N)                           
    y =zeros(N)
    for t=0:τ:32T 
        y[2:N-1]=c*(y1[3:N]+y1[1:N-2])+2*(1-c)*y1[2:N-1]-y0[2:N-1]
        y[1]=gsp(t)
        y[N]=gsp(t-16T-τ)
        y0 .= y1
        y1 .= y
        plot(x,y)
        xlim([0,L])
        ylim([-2,2]) 
    end
end
absorb1()

# 上述方法理论上正确，但是由于计算误差，效果并不好。
# 可以直接利用计算结果，然后根据振动方程预测
function absorb2()
    y0=zeros(N)                             
    y1=zeros(N)                           
    y =zeros(N)
    d =zeros(3) # 延迟2个时间步长，刚好对应1个空间步长。
    for t=0:τ:32T
        y[2:N-1]=c*(y1[3:N]+y1[1:N-2])+2*(1-c)*y1[2:N-1]-y0[2:N-1]
        y[1]=gsp(t)
        d=circshift(d,1)
        d[1]=y[N-1]
        y[N]=d[3]
        y0 .= y1
        y1 .= y
        plot(x,y)
        xlim([0,L])
        ylim([-2,2])
    end
end
absorb2()