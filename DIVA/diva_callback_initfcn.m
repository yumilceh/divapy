function diva_callback_initfcn(forceinit)
global DIVA_x;
if ~isfield(DIVA_x,'gui'), DIVA_x.gui=0; end
%if ~isfield(DIVA_x,'init'),DIVA_x.init=1;end

if DIVA_x.gui, 
    diva_gui softinit; 
end

try
    evalin('base','Target_production;');
catch
    load('diva_weights_SSM.mat','W');
    N_productions=size(W,1);
    Target_production=[0,1,zeros(1,N_productions-1)]; % if target not defined, try running the first production available
    assignin('base','Target_production',Target_production);
end
% if (nargin>0&&forceinit) || DIVA_x.init
%     diva_preparesimulation('random');
% end

