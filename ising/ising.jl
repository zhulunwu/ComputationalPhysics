using Statistics
using GR

# 计算翻转前后的能量差异
function deltaE(σ, i, j, J, H)
    L = size(σ)[1]
    t = i>1 ? i-1 : L
    b = i<L ? i+1 : 1
    l = j>1 ? j-1 : L
    r = j<L ? j+1 : 1
    ΔE = 2σ[i,j]*(J*(σ[t,j]+σ[b,j]+σ[i,l]+σ[i,r])+H)
    return ΔE
end

function flip(σ,K,T,J,H)
    L = size(σ)[1]
    for i=1:L
        for j=1:L
            ΔE = deltaE(σ,i,j,J,H)
            if ΔE <= 0
                σ[i,j] = - σ[i,j]
            elseif rand()<exp(-ΔE/(K*T)) #如果能量增加，按一定概率取
                σ[i,j] = - σ[i,j]
            end
        end
    end
end

function sweep(σ,K,T,J,H,relax,nstep)
    L = size(σ)[1]
    for r=1:relax
        flip(σ,K,T,J,H)
    end
    mag=zeros(nstep)
    for step=1:nstep
        flip(σ,K,T,J,H)
        mag[step]=abs(mean(σ)) # 状态平均，正负等效
    end
    return mean(mag)
end

function update()
    L = 40
    K = 1
    J = 1
    H = 0 # 无外加磁场
    T = 3.5
    σ = rand([1,-1],L,L) 
    for i=1:100
        flip(σ,K,T,J,H)
        heatmap(σ)
    end
    println(mean(σ))
end

function simulation()
    L = 10
    K = 1
    J = 1
    H = 0 # 无外加磁场

    relax = 100
    nstep = 1000

    σ = rand([1,-1],L,L) 
    results = Float64[]
    for T=0.1:0.05:4.0 # 温度扫描
        result=sweep(σ,K,T,J,H,relax,nstep)
        push!(results,result)
    end
    plot(collect(0.1:0.05:4.0),results)  
end