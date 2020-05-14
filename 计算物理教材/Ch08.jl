# 布朗运动模拟
using Distributions

function sim01(sigma = 2,N = 101)
    W1=0
    data=[]
    for n=1:N
        Z = rand(Normal(0,sigma))
        W1 = W1 + Z
        push!(data,W1)
    end
    return data
end


irt = 0:1:99;
plot(irt,W1Value)
hold on;
W2 = 0;
n = 1;
sigma = 2;
N = 101;
while (n < N)
Z = normrnd(0,sigma);
W2 = W2 + Z;
W2Value(n) = W2;
n = n + 1;
end
irt = 0:1:99;
plot(irt,W2Value)
hold on;
W3 = 0;n = 1;
sigma = 2;
N = 101;
while (n < N)
Z = normrnd(0,sigma);
W3 = W3 + Z;
W3Value(n) = W3;
n = n + 1;
end
irt = 0:1:99;
plot(irt,W3Value)
axis([0 150 -40 40])
title('布朗运动过程仿真结果')
xlabel('时间 t')
ylabel('布朗运动过程 w(t)')
legend('样本函数 1','样本函数 2','样本函数 3')