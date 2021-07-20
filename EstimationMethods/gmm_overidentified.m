function [paraest]=gmmestimation(moment,para0,data,number,K)
%This program is for  GMM replication of Lewbel
%input:
%para0:initial value for estimated parameters
%y,w:data used to estimate parameters
%number: maximum convergence number when choosing optimal weighting matrix
%K:number of moment conditions
%output:   
%paraest:parameters estimated

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Created by Monika Avila Marquez using as reference code gmmestimation.m of Cao Zhiguang
%E-mail:monika.avila@unige.ch
%Date:2021/07/08

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%main code for GMM


    
y=data(:,1);
x=data(:,2);
z=data(:,3:4);
W(:,:,1)=eye(K);
options = optimset('MaxFunEvals',2500);
[para(:,1),fv(:,1)]=fminsearch(moment,para0,options,1,y,x,z,eye(K));

for i=2:number
mom=feval(moment,para(:,i-1),2,y,x,z,eye(K));
W(:,:,i)=weighting_matrix(mom);
[para(:,i),fv(:,i)]=fminsearch(moment,para(:,i-1),options,1,y,x,z,W(:,:,i));

paraest=para(:,i);
if abs(fv(:,i)-fv(:,i-1))/abs(fv(:,i-1))<1e-4|fv(i)<=1e-15
  break  
elseif i==number     
    error('number of iteration exceeds defined maximium number')
end
 
paraest=para(:,i);
f0=feval(moment,paraest,3,y,x,z,W(:,:,i));

end 

end 


function [W]=weighting_matrix(mom)
s=size(mom);
W=inv(mom'*mom/s(1,1));
end 