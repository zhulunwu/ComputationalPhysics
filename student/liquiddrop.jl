#读入数据
using DelimitedFiles
pixdata=readdlm("pixdata.txt")

#数据转化成物理单位
k=7.0e-6
data=pixdata*k

#坐标转换
xs=data[:,1]
offset=0.5*(maximum(xs)-minimum(xs))
data[:,1] = data[:,1] .- offset
ys=data[:,2]
data[:,2]=data[:,2] .- minimum(ys)

# 简单画个图看看
# using Makie
# lines(data[:,1],data[:,2])

#截取右半数据来和ode进行比较
rdata=data[data[:,1] .> 0,:]
# lines(rdata[:,1],rdata[:,2])

#构造s向量
sl=first(size(rdata))
se=zeros(sl) #s of experiment
for i=2:sl
    x1=rdata[i-1,1]
    y1=rdata[i-1,2]
    x2=rdata[i,1]
    y2=rdata[i,2]
    se[i]=se[i-1]+sqrt((x2-x1)^2+(y2-y1)^2)
end

#解常微分方程
using DifferentialEquations
function shapeline(du,u,β,t)
    x,z,ϕ=u
    du[1] = dx = cos(ϕ)
    du[2] = dz = sin(ϕ)
    du[3] = dϕ = 2+β*z-sin(ϕ)/x
end


#参数和求解
function R0β(R0,β)
    data=rdata ./ R0
    sexp=se ./ R0
    u0=[data[1,1],data[1,2],sexp[1]]
    tspan=(0,sexp[end])
    prob = ODEProblem(shapeline,u0,tspan,β)
    sol = solve(prob,saveat=sexp)
    u=hcat(sol.u...)
    xz=u[1:2,:]'
    return data,xz
end

#参数扫描
function scan(R0_min,R0_max,β_min,β_max)
    ΔR0=R0_max-R0_min
    Δβ =β_max - β_min
    emin=Inf
    R0_opt=0
    β_opt=0
    for R0=R0_min:0.01ΔR0:R0_max
        for β=β_min:0.01Δβ:β_max
            d,xz=R0β(R0,β)
            e=sum((d-xz).^2)
            if e<emin
                emin=e
                R0_opt=R0
                β_opt=β
            end
        end
    end
    return R0_opt,β_opt
end

R0,β=scan(0.001,0.01,0.1,10)
d,xz=R0β(R0,β)
plot(d[:,1],d[:,2])
plot!(xz[:,1],xz[:,2])

g=9.8
ρ=1000
γ=ρ*g*R0^2/β