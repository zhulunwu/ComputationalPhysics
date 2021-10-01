using DifferentialEquations
using Makie

β=0.5
function surfaceTensor!(du,u,p,t)
    du[1] = cos(u[3])
    du[2] = sin(u[3])
    du[3] = 2+β*u[2]-sin(u[3])/u[1]
end
u0 = [0.1;0.0;0.0]
tspan = (0.0,1.5)
prob = ODEProblem(surfaceTensor!,u0,tspan)
sol = solve(prob)

ps=hcat(sol.u...)
lines(ps[1:2,:])

# 给定三点，输出由这三个点确定园的曲率半径.数据格式是三行两列
# 任意园的方程为(x-a)^2+(y-b)^2=r^2，三个参数，三个方程可以定出。
function radius(points)
    x1,y1=points[1,:]
    x2,y2=points[2,:]
    x3,y3=points[3,:]
    A1=x1^2+y1^2
    A2=x2^2+y2^2
    A3=x3^2+y3^2
    B1=0.5(A1-A2)
    B2=0.5(A1-A3)
    x12=x1-x2
    y12=y1-y2
    x13=x1-x3
    y13=y1-y3 
    c1=x12*y13-x13*y12
    c2=y12*x13-y13*x12
    if abs(c1)>eps() && abs(c2)>eps() 
        a=(B1*y13-B2*y12)/c1  
        b=(B1*x13-B2*x12)/c2
        r=sqrt((x1-a)^2+(y1-b)^2)
        return r
    else
        return "∞"
    end
end

# 测试上述函数。数据由下述函数产生。由于只关心曲率半径，所以圆心坐标随机。
function circle_points(r)
    a=rand()
    b=rand()
    x1=rand()
    y1=sqrt(r^2-(x1-a)^2)
    x2=rand()
    y2=sqrt(r^2-(x2-a)^2)
    x3=rand()
    y3=sqrt(r^2-(x3-a)^2)
    return [x1 y1;x2 y2;x3 y3]   
end
# 测试结果正确
ps=[[332.5 261.25];[426.25 282.5];[481.25 326.25]]
