%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Monika Avila Marquez 
% Code: GMM implementation
% Goal: GMM simulation and estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/Users/moka/Dropbox/1. PHD/P3/7. Replications/7.1 Lewbel/1. Matlab/Data Generation')
addpath('/Users/moka/Dropbox/1. PHD/P3/7. Replications/7.1 Lewbel/1. Matlab/Moments')
addpath('/Users/moka/Dropbox/1. PHD/P3/7. Replications/7.1 Lewbel/1. Matlab/EstimationMethods')

clear all 
close all
% Set the random seed 
rng('default')
s = rng

% Number iterations
maxiter=100;

% Set up 
% Number of observations 
N=400;
%DGP Options: 1:'Lewbel_nocovariates';2:'Lewbel_covariates';3:'TriangularModelJustIdentified';4:'TriangularModelOverIdentified'
DGP=1;
overidentification=1;



parameters_estimates=[]
loop=1
for loops=1:maxiter
    
    % Data Generation
    if DGP==1
    
        % DGP 3 Linear justidentified 
        beta=121; % 
        pi=5;
        parameters=[beta,pi];
        data=datageneration(DGP,N,parameters);
        y=data.dependent;
        x=data.endogenousregressor;
        z=data.ivs;
    
    elseif DGP==2
    
        % DGP 4 Linear overidentified 
        beta=121; % 
        pi=[5;8];
        parameters=[beta;pi];
        data=datageneration(DGP,N,parameters);
        y=data.dependent;
        x=data.endogenousregressor;
        z=data.ivs;
        
    
    end 


    
    
    number=100;
    

    if DGP==1 & overidentification==0
        para0=[1.5;1.6];
        [paraest]=gmm_justidentified(DGP,'moments_linear_model_justidentified',para0,[y,x,z],number,2);
    elseif DGP==2 & overidentification==1
        para0=[1.5;1.6;4];
        [paraest]=gmm_overidentified(DGP,'moments_linear_model_overidentified',para0,[y,x,z],number,4);            
    end 

    parameters_estimates=[parameters_estimates,paraest];
    loop=loop+1
end



if DGP==1 & overidentification==0        
    
    parameters_eval=[1,1,1.72,1.64,1.64]';
    sp=size(parameters_estimates);
    mean_par=mean(parameters_estimates,2);
    sd=std(parameters_estimates,0,2);
    lq=quantile(parameters_estimates,0.25,2);
    med=quantile(parameters_estimates,0.5,2);
    uq=quantile(parameters_estimates,0.75,2);
    bias=mean(parameters_estimates-kron(ones(1,sp(1,2)),parameters_eval),2);
    rmse=sqrt(mean(((parameters_estimates-kron(ones(1,sp(1,2)),parameters_eval)).^2),2));
    mae=mean(abs((parameters_estimates-kron(ones(1,sp(1,2)),parameters_eval))),2);
    mdae=median(abs((parameters_estimates-kron(ones(1,sp(1,2)),parameters_eval))),2);
    
    results_parameters_levels.parameters=parameters_estimates;
    results_parameters_levels.mean_par=mean_par;
    results_parameters_levels.sd=sd;
    results_parameters_levels.lq=lq;
    results_parameters_levels.med=med;
    results_parameters_levels.uq=uq;
    results_parameters_levels.bias=bias;
    results_parameters_levels.rmse=rmse;
    results_parameters_levels.mae=mae;
    results_parameters_levels.mdae=mdae;
    
    
    

elseif DGP==2 & overidentification==1
    
    parameters_eval=[1,1,1.72,1.64,1.64,10.17]';
    sp=size(parameters_estimates);
    mean_par=mean(parameters_estimates,2)
    sd=std(parameters_estimates,0,2)
    lq=quantile(parameters_estimates,0.25,2)
    med=quantile(parameters_estimates,0.5,2)
    uq=quantile(parameters_estimates,0.75,2)
    bias=mean(parameters_estimates-kron(ones(1,sp(1,2)),parameters_eval),2)
    rmse=sqrt(mean(((parameters_estimates-kron(ones(1,sp(1,2)),parameters_eval)).^2),2))
    mae=mean(abs((parameters_estimates-kron(ones(1,sp(1,2)),parameters_eval))),2)
    mdae=median(abs((parameters_estimates-kron(ones(1,sp(1,2)),parameters_eval))),2)
    
    results_parameters_levels.parameters=parameters_estimates;
    results_parameters_levels.mean_par=mean_par;
    results_parameters_levels.sd=sd;
    results_parameters_levels.lq=lq;
    results_parameters_levels.med=med;
    results_parameters_levels.uq=uq;
    results_parameters_levels.bias=bias;
    results_parameters_levels.rmse=rmse;
    results_parameters_levels.mae=mae;
    results_parameters_levels.mdae=mdae;
    
    
   
end 
    
save(sprintf('Results/Results_parameters_lev_DGP%d_N%d_overid%d.mat',DGP,N,overidentification),'results_parameters_levels')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
