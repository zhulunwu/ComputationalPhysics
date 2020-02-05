using Statistics
using GR

# 计算翻转前后的能量差异
function ΔE(σ, i, j, L, H, J)
    top = i>1 ? i-1 : L
    bottom = i<L ? i+1 : 1
    left = j>1 ? j-1 : L
    right = j<L ? j+1 : 1
    E_former = -1.0*σ[i,j]*(J*(σ[top,j]+σ[bottom,j]+σ[i,left]+σ[i,right])+H)
    E_next = σ[i,j]*(J*(σ[top,j]+σ[bottom,j]+σ[i,left]+σ[i,right])+H)
    delta_E = E_next - E_former
    return delta_E
end

function simulation(L,K,T,J,H,nstep,relax)
    σ = rand([1,-1],L,L) 
    mag=Float64[] #某个温度对应的序列
    for step=1:(nstep+relax)
        for i=1:L
            for j=1:L
                delta_E = ΔE(σ, i, j, L, H, J)
                if delta_E <= 0
                    σ[i,j] = - σ[i,j]
                elseif exp((- delta_E)/(K * T)) > rand() #如果能量增加，按一定概率取
                    σ[i,j] = - σ[i,j]
                end
            end
        end
    
        M=abs(sum(σ))/(L*L) # 状态平均，正负等效
        push!(mag,M) #记录
    end # sweeps 结束
    result = mean(mag[relax+1:end]) # 统计平均值，去掉relax部分。
    return result
end

sweeps = 1000 
relax = 1000
K = 1
J = 1
H = 0 # 无外加磁场
L=40

results = Float64[]
for T=0.1:0.05:4.0 # 温度扫描
    result=simulation(L,K,T,J,H,sweeps,relax)
    push!(results,result)
end
plot(collect(0.1:0.05:4.0),results)
