
function f=moments_linear_model_justidentified(para,num,y,x,z,W)
[N,q]=size(y);
beta=para(1);pi=para(2);

a=(y-beta*x).*z;
b=(x-pi*z).*z;

m_1N=a;   
m_2N=b;   

m_1=mean(m_1N);
m_2=mean(m_2N);
m=[m_1;m_2];
mN=[m_1N;m_2N];
obj=m'*W*m;
if num==1
    f=obj;
    elseif num==2
    f=mN;
    elseif num==3
    f=m;
end
end 