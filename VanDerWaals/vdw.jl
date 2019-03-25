# 范德瓦耳斯方程
# 等温曲线
using GR

const R=8.314472
const a_H2O=0.5535
const b_H2O=0.03049e-3
const P0=1.01325e5
const Tc=8a_H2O/(27R*b_H2O)
const Pc=a_H2O/(27*b_H2O^2)
const Vc=3b_H2O

for T=0.5Tc:0.1Tc:1.2Tc
    V=0.388Vc:0.01Vc:2Vc
    p=(R*T)./(V.-b_H2O)-a_H2O./(V.^2)
    plot(V,p./P0)
    hold(true)
end