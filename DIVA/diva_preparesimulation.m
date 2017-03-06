function timeseries=diva_preparesimulation(option,timeseries)

params=diva_vocaltract;
if nargin<2 || isempty(timeseries), timeseries=diva_targets('timeseries',option,'header'); end
N_productions=1;
N_samplesperproduction=size(timeseries(1).Art,1);

% defines simulation target variables
Target_production=zeros(N_productions,1);
Target_production(1)=1;
Target_production=[0,Target_production(:)'];
assignin('base','Target_production',Target_production);

% defines simulation weight matrices
W=kron(eye(N_productions),ones(1,N_samplesperproduction));
save diva_weights_SSM.mat W
W=timeseries.Art;
save diva_weights_SSM2FF.mat W
W=timeseries.Aud_min*diag(1./params.Output(1).Scale);
W0=min(-1,params.Output(1).Range(:,1)')./params.Output(1).Scale';
save diva_weights_SSM2amin.mat W W0
W=timeseries.Aud_max*diag(1./params.Output(1).Scale);
W0=max(1,params.Output(1).Range(:,2)')./params.Output(1).Scale';
save diva_weights_SSM2amax.mat W W0
W=timeseries.Som_min*diag(1./params.Output(2).Scale);
W0=min(-1,params.Output(2).Range(:,1)')./params.Output(2).Scale';
save diva_weights_SSM2smin.mat W W0
W=timeseries.Som_max*diag(1./params.Output(2).Scale);
W0=max(1,params.Output(2).Range(:,2)')./params.Output(2).Scale';
save diva_weights_SSM2smax.mat W W0

clear diva_weightsadaptive
clear diva_weights
