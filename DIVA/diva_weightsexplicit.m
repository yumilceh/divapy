function diva_weightsexplicit(block)
% Level-2 M file S-Function.
  setup(block);  
end

%% Initialization   
function setup(block)
  
  % Register number of dialog parameters
  block.NumDialogPrms = 3;
  block.DialogPrmsTunable = {'Nontunable','Nontunable','Nontunable'};

  % Register number of input and output ports
  block.NumInputPorts  = 2;
  block.NumOutputPorts = 1;

  % Setup functional port properties to dynamically inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
 
  block.InputPort(1).Dimensions        = -1;
  block.InputPort(1).DirectFeedthrough = true;  
  block.InputPort(2).Dimensions        = -1;
  block.InputPort(2).DirectFeedthrough = true;  
  block.OutputPort(1).Dimensions       = -1;
  
  % Set block sample time to discrete
  block.SampleTimes = [-1 0];
  
  % Register methods
  block.RegBlockMethod('SetInputPortDimensions',  @SetInputDims);
  block.RegBlockMethod('Outputs',                 @Output);  
  
end

 
function SetInputDims(block, port, dm)
    block.InputPort(port).Dimensions = dm;
    if port==2, block.OutputPort(1).Dimensions = dm; end
end

%% Output & Update equations   
function Output(block)

    dy=block.InputPort(1).Data;
    x=block.InputPort(2).Data;
    nout=block.DialogPrm(1).Data;
%     switch(nout)
%         case 1, nout='auditory';
%         case 2, nout='somatosensory';
%     end
%     switch(lower(nout)),
%         case 'auditory', nout=1;
%         case 'somatosensory', nout=2;
%     end
    EPS=block.DialogPrm(2).Data;
    LAMBDA=block.DialogPrm(3).Data;
    N=numel(x);
    M=numel(dy);

    Ix=eye(N);
    Iy=eye(M);
    Q=diva_vocaltract('base');%Ix;%orth(randn(N));
    DX=EPS*Q; % direction of articulatory change
    DY=zeros([M,N]); % direction of auditory/somatosensory change
    y=diva_vocaltract(nout,x);
    for ndim=1:N, % computes jacobian
        xt=x+DX(:,ndim);
        DY(:,ndim)=diva_vocaltract(nout,xt)-y;
    end
    JJ=DY*DY';
    iJ=EPS*Q*DY'*pinv(JJ+LAMBDA*EPS^2*Iy); % computes pseudoinverse
    dx=-iJ*dy;
    block.OutputPort(1).Data = dx;
%     K=.05;
%     dx0=-K*(Ix-iJ*DY*Q'/EPS)*x;
%     block.OutputPort(1).Data = dx+dx0;

end

