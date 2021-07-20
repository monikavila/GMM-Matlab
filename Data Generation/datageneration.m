%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Monika Avila Marquez 
% Code: Data Generation
% Goal: GMM simulation and estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data]=data_generation(DGP,N,parameters)


    
if DGP==1 % Linear triangular model with just identification
    
    
    beta=parameters(1);
    pi=parameters(2);
    mu=[0;0];
    sigma=[1,0.5;0.5,1];
    U=mvnrnd(mu,sigma,N);
    eps=U(:,1);
    u=U(:,2);
    z=normrnd(0,1,N,1);
    x=pi*z+u;
    y=x*beta+eps;
    
    
    data.dependent=y;
    data.ivs=z;
    data.endogenousregressor=x;
    data.parameters=[beta,pi];
    data.errors=[eps,u];
    
    save('Data/data_linear_triangular_model_justidentified.mat')
        
elseif DGP==2
    
    beta=parameters(1);
    pi=parameters(2:3);
    mu=[0;0];
    sigma=[1,0.5;0.5,1];
    U=mvnrnd(mu,sigma,N);
    eps=U(:,1);
    u=U(:,2);
    z=mvnrnd([0,0],[1,0;0,2],N);
    x=z*pi+u;
    y=x*beta+eps;
        
    data.dependent=y;
    data.ivs=z;
    data.endogenousregressor=x;
    data.parameters=[beta;pi];
    data.errors=[eps,u];
    
    save('Data/data_linear_triangular_model_overidentified.mat')
    
    
end

end 
