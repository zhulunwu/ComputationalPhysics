
using GR
include("参数以及定义.jl")

#吸收边界方法一：根据波动规律直接获得边界位移值
#根据波的传播规律，右端点的振动位移是左端点传播L/u的结果。
#理论丰满，现实骨感。由于计算存在误差，吸收并不很理想，只是说明原理。
#况且需要知道波的表达式，现实中没有用处。
# 波动微分方程求解，即Y_tt-u^2*Y_xx=0,u是波速

E=Float64[]
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
       
        v=y1 .- y0
        Ek=sum(v.^2)    
        Ep=c*sum(diff(y).^2) #势能计算  

        push!(E,Ek+Ep) #波的无量纲能量
    end
end
absorb_demo()
ylim([0,maximum(E)])
plot(E)
入射波能量=E[length(E)÷2]
反射波能量=E[end]
反射率=反射波能量/入射波能量
println("反射率=",反射率)
