#半波损失的动力学模拟
#基本模型是质点和无质量弹簧连接。节点用质点的质量差异来模拟。
#质点之间的受力和弹簧的长度相关。

using DifferentialEquations
using WGLMakie

#一维数组表示质点的平衡坐标
Δx=1    #无量纲
L=2Δx
x0=0:Δx:L

#向右为受力的正方向
# function wave(ds,s,p,t)
    # ds[1]=s[2]
    # ds[3]=s[4]
    # ds[5]=s[6]
    # ds[2]=p*(s[3]-s[1]-1)   #大于零则表明弹簧处于拉长状态，左边质点受力为正
    # ds[4]=p*(s[1]-s[3]+1)+p*(s[5]-s[3]-1)
    # ds[6]=p*(s[3]-s[5]+1)
    # ds=[s[2],p*(s[3]-s[1]-1),s[4],p*(s[1]-s[3]+1)+p*(s[5]-s[3]-1),s[5],p*(s[3]-s[5]+1)]
    # ds[1]=s[2]
    # ds[2]=-s[1]
    # ds=[s[2],-s[1]] 
# end
halfwave(s,p,t)=[s[2],p*(s[3]-s[1]-1),s[4],p*(s[1]-s[3]+1)+p*(s[5]-s[3]-1),s[5],p*(s[3]-s[5]+1)]
tspan = [0.0,50.0]
s0=[0,0,1,0,2,0]
p=0.01
prob = ODEProblem(halfwave,s0,tspan,p)
sol = solve(prob)

ps=hcat(sol.u...)
lines(sol.t,ps[1,:])