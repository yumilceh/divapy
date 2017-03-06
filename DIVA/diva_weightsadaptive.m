function W=diva_weightsadaptive(block,opt)
% Level-2 M file S-Function.
  if ischar(block)
      switch(lower(block))
          case 'weights',
              W=getsetWeights([], 4, opt);
      end
  else
      setup(block);
  end
end

%% Initialization   
function setup(block)
  
  ndims=getsetWeights(block, 0);
  
  % Register number of dialog parameters
  block.NumDialogPrms = 2;
  block.DialogPrmsTunable = {'Nontunable','Nontunable'};

  % Register number of input and output ports
  block.NumInputPorts  = 3;
  block.NumOutputPorts = 1;

  % Setup functional port properties to dynamically inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
 
  block.InputPort(1).Dimensions        = ndims(1);
  block.InputPort(1).DirectFeedthrough = true;  
  block.InputPort(2).Dimensions        = -1;
  block.InputPort(2).DirectFeedthrough = true;  
  block.InputPort(3).Dimensions        = -1;
  block.InputPort(3).DirectFeedthrough = true;  
  block.OutputPort(1).Dimensions       = ndims(2);
  
  % Set block sample time to discrete
  block.SampleTimes = [-1 0];
  
  % Register methods
  block.RegBlockMethod('SetInputPortDimensions',  @SetInputDims);
  block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
  block.RegBlockMethod('InitializeConditions',    @InitConditions);  
  block.RegBlockMethod('Outputs',                 @Output);  
  
end

 
function DoPostPropSetup(block)
  % Setup Dwork
  block.NumDworks = 1;
  block.Dwork(1).Name = 'WeightsMatrix'; 
  block.Dwork(1).Dimensions      = 1;
  block.Dwork(1).DatatypeID      = 0;
  block.Dwork(1).Complexity      = 'Real';
  block.Dwork(1).UsedAsDiscState = false;
end

function SetInputDims(block, port, dm)
    block.InputPort(port).Dimensions = dm;
end

function InitConditions(block)
    idx=getsetWeights(block, 1);
    block.Dwork(1).Data=idx;
end


%% Output & Update equations   
function Output(block)

    W=getsetWeights(block, 2);
    t=block.InputPort(1).Data'*W;
    block.OutputPort(1).Data = t;

    %disp(['learning ',num2str([max(block.InputPort(2).Data),max(abs(block.InputPort(3).Data))])]);
    if any(block.InputPort(2).Data)||any(block.InputPort(3).Data)
        EPS=block.DialogPrm(2).Data;
        u=reshape(block.InputPort(2).Data,size(W,1),[]);
        v=reshape(block.InputPort(3).Data,size(W,2),[]);
        dW=EPS*u*v';
        W=W+dW;
        getsetWeights(block, 3, W);
    end
end

function [out]=getsetWeights(block, opt, Wnew)
  persistent Weights Changes Filenames; %InpDims OutDims;
  if opt==0, % load filename & return weight matrix size
      if isempty(Weights),
          idx=1;
          Filenames={block.DialogPrm(1).Data};
          load(block.DialogPrm(1).Data,'W');%,'inpdims','outdims');
          Weights{1}=W;
          Changes{1}=zeros(size(W));
      else
          idx=strmatch(block.DialogPrm(1).Data,Filenames,'exact');
          if isempty(idx), 
              idx=numel(Filenames)+1; 
              Filenames{idx}=block.DialogPrm(1).Data;
              load(block.DialogPrm(1).Data,'W');%,'inpdims','outdims');
              Weights{idx}=W;
              Changes{idx}=zeros(size(W));
          end
      end
      out=size(Weights{idx});
  elseif opt==1, % return index to weight matrix
      idx=strmatch(block.DialogPrm(1).Data,Filenames,'exact');
      out=idx;
  elseif opt==2 % return weight matrix
      idx=block.Dwork(1).Data;
      out=Weights{idx};
  elseif opt==3 % set weight matrix
      idx=block.Dwork(1).Data;
      Weights{idx}=Wnew;
      %dW=Wnew-Weights{idx};
      %Weights{idx}=Weights{idx}+(Changes{idx}>=0).*(max(0,dW)+min(0,dW+Changes{idx}))+(Changes{idx}<=0).*(min(0,dW)+max(0,dW+Changes{idx}));
      %Changes{idx}=.9*Changes{idx}+dW;
      %figure(11);imagesc(Changes{1}(1:120,:));colorbar;drawnow;
      %figure(10);imagesc(Wnew);drawnow
  elseif opt==4 % return weight matrix from filename
      idx=strmatch(Wnew,Filenames,'exact');
      if isempty(idx), out=[];
      else out=Weights{idx}; end
  end
end

