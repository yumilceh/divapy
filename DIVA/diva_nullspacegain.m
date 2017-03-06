function diva_nullspacegain(block)

% Level-2 M file S-Function.
  setup(block);  
end

%% Initialization   
function setup(block)

%   % Register number of dialog parameters   
  block.NumDialogPrms = 4;
  block.DialogPrmsTunable = {'Nontunable','Nontunable','Nontunable','Nontunable'};

  % Register number of input and output ports
  block.NumInputPorts  = 1;
  block.NumOutputPorts = 1;

  % Setup functional port properties to dynamically inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
 
  block.InputPort(1).Dimensions        = -1;
  block.InputPort(1).DirectFeedthrough = true;
  block.OutputPort(1).Dimensions       = -1;
  
  % Set block sample time to discrete
  block.SampleTimes = [-1 0];
  
  % Register methods
  block.RegBlockMethod('SetInputPortDimensions',  @SetInputDims);
  block.RegBlockMethod('Outputs',                 @Output);  
  
end


function SetInputDims(block, port, dm)
    block.InputPort(port).Dimensions = dm;
    if port==1, block.OutputPort(1).Dimensions = dm; end
end



%% Output & Update equations   
function Output(block)

  % system output
  K=block.DialogPrm(1).Data;
  nout=block.DialogPrm(2).Data;
  EPS=block.DialogPrm(3).Data;
  LAMBDA=block.DialogPrm(4).Data;
%   switch(nout)
%       case 1, nout='auditory';
%       case 2, nout='somatosensory';
%       case 3, nout='auditory&somatosensory';
%   end
%   switch(lower(nout))
%       case 'auditory', nout=1;
%       case 'somtosensory', nout=2;
%       case 'auditory&somatosensory', nout=[1,2]; 
%   end
  
  x=block.InputPort(1).Data;
  N=numel(x);

  Ix=eye(N);
  Q=diva_vocaltract('base');%Ix;%orth(randn(N));
  %Q=orth(randn(N));
  %I=eye(N);
  DX=EPS*Q;
  y=diva_vocaltract(nout,x);
  M=size(y,1);
  Iy=eye(M);
  DY=zeros([M,N]);
  for ndim=1:N,
      xt=x+DX(:,ndim);
      DY(:,ndim)=diva_vocaltract(nout,xt);-y;
  end
  JJ=DY*DY';
  N=Ix-EPS*Q*DY'*pinv(JJ+LAMBDA*EPS^2*Iy)*DY*Q'; % computes null space projector

  block.OutputPort(1).Data = K*N*x;
end

