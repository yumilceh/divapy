function diva_weights(block)
% Level-2 M file S-Function.
  setup(block);  
end

%% Initialization   
function setup(block)
  
  ndims=getWeights(block, 0);
  
  % Register number of dialog parameters
  block.NumDialogPrms = 1;
  block.DialogPrmsTunable = {'Nontunable'};

  % Register number of input and output ports
  block.NumInputPorts  = 1;
  block.NumOutputPorts = 1;

  % Setup functional port properties to dynamically inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
 
  block.InputPort(1).Dimensions        = ndims(1);
  block.InputPort(1).DirectFeedthrough = true;  
  block.OutputPort(1).Dimensions       = ndims(2);
  
  % Set block sample time to discrete
  block.SampleTimes = [-1 0];
  
  % Register methods
%   block.RegBlockMethod('SetInputPortDimensions',  @SetInputDims);
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

% function SetInputDims(block, port, dm)
%     block.InputPort(port).Dimensions = dm;
% end

function InitConditions(block)
    idx=getWeights(block, 1);
    block.Dwork(1).Data=idx;
end


%% Output & Update equations   
function Output(block)

    [W,W0]=getWeights(block, 2);
    x=block.InputPort(1).Data;
    if max(abs(x))<1e-10, % bias term
        t=W0;
    else
        t=x'*W;
    end
    block.OutputPort(1).Data = t;

end

%function [out,inpdims,outdims]=getWeights(block, opt)
function [out,bias]=getWeights(block, opt)
  persistent Weights Bias Filenames; %InpDims OutDims;
  if opt==0, % load filename & return weight matrix size
      if isempty(Weights),
          idx=1;
          Filenames={block.DialogPrm(1).Data};
          data=load(block.DialogPrm(1).Data);%,'inpdims','outdims');
          Weights{1}=data.W;
          if isfield(data,'W0'), Bias{1}=data.W0;
          else                   Bias{1}=zeros(1,size(data.W,2));
          end
          %InpDims{1}=inpdims;
          %OutDims{1}=outdims;
      else
          idx=strmatch(block.DialogPrm(1).Data,Filenames,'exact');
          if isempty(idx), 
              idx=numel(Filenames)+1; 
              Filenames{idx}=block.DialogPrm(1).Data;
              data=load(block.DialogPrm(1).Data);%,'inpdims','outdims');
              Weights{idx}=data.W;
              if isfield(data,'W0'), Bias{idx}=data.W0;
              else                   Bias{idx}=zeros(1,size(data.W,2));
              end
          end
      end
      out=size(Weights{idx});
  elseif opt==1, % return index to weight matrix
      idx=strmatch(block.DialogPrm(1).Data,Filenames,'exact');
      out=idx;
  else % return weight matrix
      idx=block.Dwork(1).Data;
      out=Weights{idx};
      bias=Bias{idx};
  end
end

