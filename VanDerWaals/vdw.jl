# 范德瓦耳斯方程
# 等温曲线
using GR

const R=8.314472
const a_H2O=0.5535
const b_H2O=0.0305
const P0=1.01325e5

T=300
V=0.02:0.001:0.024
p=(R*T)./V
pp=p./P0
plot(pp,v)


p=(R*T)./(V.-b_H2O)-a_H2O./(V.^2)
pp=p./P0
plot(p,v)