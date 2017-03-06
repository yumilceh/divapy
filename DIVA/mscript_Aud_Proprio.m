%Created on: 15th February 2016
%by: Juan Manuel Aceveo Valle
%
%artStates and outputScale exit

samples=size(artStates,1);
auditoryStates=zeros(4,samples);
minaf=zeros(1,samples);
for k=1:samples
    [auditoryStates(:,k), som, ~, af]=diva_synth(artStates(k,1:13)');
    auditoryStates(:,k)=auditoryStates(:,k)./outputScale;
    minaf(k)=min(af);
    clear af;
end
auditoryStates=auditoryStates';
minaf=minaf';
%save test1.mat samples auditoryStates minaf artStates outputScale 

