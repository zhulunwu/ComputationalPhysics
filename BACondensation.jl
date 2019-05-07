# 能量单位为h^2/(2mL^2)，化学势和能量单位相同。特征温度Ts=h^2/(2mkL^2),温度单位取Ts
# 则玻色爱因斯坦统计的能级粒子数表达式为a=1/(exp((ϵ-μ)/T)-1)

# 态计算函数
function states(nx::Int,ny::Int,nz::Int)
    x=(nx==0 ? 1 : 2)
    y=(ny==0 ? 1 : 2)
    z=(nz==0 ? 1 : 2)
    return x*y*z
end

# 非零能级粒子数分布
function npe(T::Float64,μ::Float64)
    nmax=Int(ceil(sqrt((T*log(1+6/eps())+μ)/3)))
    ϵmax=3*nmax^2
    ϵn=zeros(ϵmax+1)  
     
    for nx=0:nmax
        for ny=0:nmax
            for nz=0:nmax
                ϵ=nx^2+ny^2+nz^2; #能级
                s=states(nx,ny,nz)       
                a=s/(exp((ϵ-μ)/T)-1) #能级对应的粒子数
                ϵn[ϵ+1]=ϵn[ϵ+1]+a
            end
        end
    end
    return ϵn
end

# 粒子数能量密度公式
nvse(T,μ,ϵ)=2π*sqrt(ϵ)/(exp((ϵ-μ)/T)-1)
# 临界温度
Tc(n)=n^(2/3)/((2.612)^(2/3)*π)

# 非零能级积分结果
nnz(T,n)=n*(T/Tc(n))^(3/2)

# 求和与连续函数的比较
#=
using GR
T=50.0;μ=-0.1;nmax=4;
ϵ=0:3*nmax^2;
ϵn=npe(nmax,T,μ);
en=nvse.(T,μ,ϵ)
stem(ϵn)
hold(true)
plot(en)
=#