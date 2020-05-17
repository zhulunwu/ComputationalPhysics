# 波的参数
const λ=0.2     #波长
const f=10      #频率
const u=λ*f     #波速
const ω=2π*f    #角频率
const k=2π/λ    #波数
const T=1/f     #周期
const w=2T      #脉冲宽度

const h=λ/25
const α=0.8
const c=α^2
const τ=α*h/u #如此设置是为了方便后面的计算，可以看出，一个步长需要两个时间步来完成。
const L=16λ
const X=0:h:L
const N=length(X)

# 脉冲定义
gsp(t)=sin(ω*t)exp(-(t-6T)^2/w^2)