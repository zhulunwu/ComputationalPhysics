# 波动微分方程求解，即Y_tt-u^2*Y_xx=0,u是波速
# 波的参数
const λ=0.2   #波长
const f=10    #频率

#模拟1：行波传播，固定边界，反射和驻波
function wave()      
    u=λ*f       #波速
    ω=2*pi*f    #角频率
    h=λ/25
    α=0.8;c=α^2
    τ=α*h/u
    N=8*25      #8个波长，每个波长25个h
    L=N*h
    x=range(0,L,length=N)
    y0=zeros(N)                             
    y1=zeros(N)                           
    y =zeros(N)
    for t=0:τ:400τ
        y[2:N-1]=c*(y1[3:N]+y1[1:N-2])+2*(1-c)*y1[2:N-1]-y0[2:N-1]
        y[1]=sin(ω*t)
        y0 .= y1
        y1 .= y
        plot(x,y)
        ylim([-2,2]) 
        sleep(0.1)
    end
end
wave()

#吸收边界方法一：根据波动方程解获得边界位移值
#根据波的传播规律，右端点的振动位移是左端点传播L/u的结果
function absorb()
    u=λ*f       #波速
    ω=2*pi*f    #角频率
    h=λ/25
    α=0.8;c=α^2
    τ=α*h/u
    N=8*25      #8个波长，每个波长25个h
    L=N*h
    x=range(0,L,length=N)
    y0=zeros(N)            
    y1=zeros(N)                           
    y =zeros(N)
    for t=0:τ:400τ
        if t>L/u
            y[N]=sin(ω*(t-L/u))
        end
        y[2:N-1]=c*(y1[3:N]+y1[1:N-2])+2*(1-c)*y1[2:N-1]-y0[2:N-1]
        y[1]=sin(ω*t)
        y0 .= y1
        y1 .= y
        plot(x,y)
        ylim([-2,2]) 
        sleep(0.01)
    end
end
absorb()