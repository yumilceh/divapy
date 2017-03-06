function diva_sequencer(block)

% Level-2 M file S-Function.
  setup(block);  
end

%% Initialization   
function setup(block)

%   % Register number of dialog parameters   
%   block.NumDialogPrms = 4;
%   block.DialogPrmsTunable = {'Nontunable','Nontunable','Nontunable','Nontunable'};

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
  block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
  block.RegBlockMethod('InitializeConditions',    @InitConditions);  
  block.RegBlockMethod('SetInputPortDimensions',  @SetInputDims);
  block.RegBlockMethod('Outputs',                 @Output);  
%   block.RegBlockMethod('Update',                  @Update);  
  
end


function SetInputDims(block, port, dm)
    block.InputPort(port).Dimensions = dm;
    if port==1, block.OutputPort(1).Dimensions = dm; end
end

function DoPostPropSetup(block)
  % Setup Dwork
  block.NumDworks = 2;
  block.Dwork(1).Name = 'Input'; 
  block.Dwork(1).Dimensions      = block.InputPort(1).Dimensions;
  block.Dwork(1).DatatypeID      = 0;
  block.Dwork(1).Complexity      = 'Real';
  block.Dwork(1).UsedAsDiscState = false;
  block.Dwork(2).Name = 'Sample'; 
  block.Dwork(2).Dimensions      = 1;
  block.Dwork(2).DatatypeID      = 0;
  block.Dwork(2).Complexity      = 'Real';
  block.Dwork(2).UsedAsDiscState = false;
end

function InitConditions(block)
    block.Dwork(1).Data=zeros(block.Dwork(1).Dimensions,1);
    block.Dwork(2).Data=0;
end



%% Output & Update equations   
function Output(block)

  % system output
    EPS=.01;
    N=block.Dwork(1).Dimensions;
    out=zeros(N,1);
    if any(block.InputPort(1).Data>EPS) % loads input 
        idx=find(block.InputPort(1).Data>EPS);
        block.Dwork(1).Data=cat(1,idx,zeros(N-numel(idx),1));
        block.Dwork(2).Data=1;
    end
    if block.InputPort(2).Data>EPS&&block.Dwork(2).Data>0, % outputs input sequentially
        idx=block.Dwork(1).Data;
        n=block.Dwork(2).Data;
        if n==1,n=n+block.InputPort(2).Data; end % add sub-sample jitter
        n0=floor(n);
        n1=n-n0;
        if n0+1>N||~idx(n0+1), % end of sequence
            block.Dwork(2).Data=0;
        else               % new sample point
            if n0<=1, out(idx(1))=1;
            else
                out(idx(n0))=1-n1;
                out(idx(n0+1))=n1;
            end
            block.Dwork(2).Data=n+block.InputPort(2).Data;
        end
    else
        %out(1)=1;
        %disp('no input to diva_sequencer');
    end
  block.OutputPort(1).Data = out;

end


