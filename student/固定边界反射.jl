using GR
include("参数以及定义.jl")

#行波传播，右边边界固定，因此为波节，波将完全反射.
E=Float64[]
function reflect()      
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
        
        v=y1 .- y0
        Ek=sum(v.^2)    
        Ep=c*sum(diff(y).^2) #势能计算  

        push!(E,Ek+Ep) #波的无量纲能量
    end
end
reflect()
ylim([0,maximum(E)])
plot(E)
入射波能量=E[length(E)÷2]
反射波能量=E[end]
反射率=反射波能量/入射波能量
println("反射率=",反射率)
