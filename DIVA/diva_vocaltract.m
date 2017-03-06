function varargout=diva_vocaltract(block,varargin)
global DIVA_x
DIVA_x.debug=0;

  if ~nargin,
      if ~DIVA_x.debug
          Nart=13;
          varargout{1}=struct(...
              'dosound',1,...
              'Input',struct(...
              'Name','Articulatory',...
              'Dimensions',Nart,...
              'Range',repmat([-1,1],[Nart,1]),...
              'Scale',ones(Nart,1),...
              'Default',[repmat([-1,1],[11,1]);-1,0;-1,1],...
              'BlockDiagonal',bsxfun(@eq,[ones(1,10),2,3,4],[ones(1,10),2,3,4]'),...
              'Plots_dim',{{1:10,11,12,13}},...
              'Plots_label',{{'VocalTract','tension','pressure','voicing'}} ),...
              'Output',[...
              struct(...
              'Name','Auditory',...
              'Dimensions',4,...
              'Range',[0,200;0,1000;500,3000;2000,4000],...
              'Scale',[100,500,1500,3000]',...
              'Default',[0,200;0,1000;500,3000;2000,4000],...
              'Plots_dim',{{2:4,1,2,3,4}},...
              'Plots_label',{{'Formants','F0','F1','F2','F3'}} ),...
              struct(...
              'Name','Somatosensory',...
              'Dimensions',8,...
              'Range',repmat([-1,1],[8,1]),...
              'Scale',ones(8,1),...
              'Default',[repmat([-1,-.25],[6,1]);.75,1;.75,1],...
              'Plots_dim',{{1:6,7,8,1,2,3,4,5,6}},...
              'Plots_label',{{'PlaceofArt','pressure','voicing','PA_pharyngeal','PA_uvular','PA_velar','PA_palatal','PA_alveolardental','PA_labial'}} )]);
      else
          Nart=4;
          varargout{1}=struct(...
              'dosound',0,...
              'Input',struct(...
              'Name','Motor',... %'Articulatory',...
              'Dimensions',Nart,...
              'Range',repmat([-2,2],[Nart,1]),...
              'Scale',ones(Nart,1),...
              'Plots_dim',{cat(2,{1:Nart},mat2cell(1:Nart,1,ones(1,Nart)))},...
              'Plots_label',{cat(1,{'Motor'},cellstr([repmat('motor_',[Nart,1]),num2str((1:Nart)')]))} ),...
              'Output',[...
              struct(...
              'Name','Spatial',...
              'Dimensions',2,...
              'Range',repmat([-5,5],[2,1]),...
              'Scale',ones(2,1),...
              'Default',repmat([-5,5],[2,1]),...
              'Plots_dim',{{1:2,1,2}},...
              'Plots_label',{{'Spatial','spatial_x','spatial_y'}} ),...
              struct(...
              'Name','Somatosensory',...
              'Dimensions',Nart,...
              'Range',repmat([-5,5],[Nart,1]),...
              'Scale',ones(Nart,1),...
              'Default',repmat([-5,5],[Nart,1]),...
              'Plots_dim',{cat(2,{1:Nart},mat2cell(1:Nart,1,ones(1,Nart)))},...
              'Plots_label',{cat(1,{'Somatosensory'},cellstr([repmat('somatosensory_',[Nart,1]),num2str((1:Nart)')]))} )]);
      end
  elseif ischar(block)
      switch(lower(block)),
          case 'auditory'
              varargout{1}=diva_vocaltractcompute(varargin{:});
          case 'somatosensory'
              [out{1},out{2}]=diva_vocaltractcompute(varargin{:});
              varargout{1}=out{2};
          case 'auditory&somatosensory'
              [out{1},out{2}]=diva_vocaltractcompute(varargin{:});
              varargout{1}=cat(1,out{:});
          case 'output'
              [varargout{1:nargout}]=diva_vocaltractcompute(varargin{:});
          case 'base'
              varargout{1}=eye(DIVA_x.params.Input.Dimensions);
%               Q=randn(DIVA_x.params.Input.Dimensions);
%               if isfield(DIVA_x.params.Input,'BlockDiagonal')
%                   Q=Q.*DIVA_x.params.Input.BlockDiagonal;
%               end
%               varargout{1}=orth(Q);
          otherwise
              error(['unrecognized option ',lower(block)]);
              
      end
  else
      % Level-2 M file S-Function.
      setup(block);
  end
end

%% Initialization   
function setup(block)

  params=diva_vocaltract;
%   % Register number of dialog parameters   
%   block.NumDialogPrms = 4;
%   block.DialogPrmsTunable = {'Nontunable','Nontunable','Nontunable','Nontunable'};

  % Register number of input and output ports
  block.NumInputPorts  = 1;
  block.NumOutputPorts = 2;

  % Setup functional port properties to dynamically inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
 
  block.InputPort(1).Dimensions        = params.Input.Dimensions;
  block.InputPort(1).DirectFeedthrough = true;
  block.OutputPort(1).Dimensions       = params.Output(1).Dimensions;
  block.OutputPort(2).Dimensions       = params.Output(2).Dimensions;
  
  % Set block sample time to discrete
  block.SampleTimes = [-1 0];
  
  % Register methods
  block.RegBlockMethod('SetInputPortSamplingMode',@SetInputSampling);
  %block.RegBlockMethod('SetInputPortDimensions',  @SetInputDims);
  block.RegBlockMethod('Outputs',                 @Output);  
  
end


function SetInputSampling(block, port, dm)
    block.InputPort(port).SamplingMode= dm;
    block.OutputPort(1).SamplingMode= dm;
    block.OutputPort(2).SamplingMode= dm;
end
% function SetInputDims(block, port, dm)
%     block.InputPort(port).Dimensions = dm;
%     %if port==1, block.OutputPort(1).Dimensions = dm; end
% end



%% Output & Update equations   
function Output(block)
persistent t
if isempty(t),t=0;end
t=t+1;

  % system output
  [y,z]=diva_vocaltractcompute(block.InputPort(1).Data,~rem(t,2));

  block.OutputPort(1).Data = y;
  block.OutputPort(2).Data = z;

end


function [y,z]=diva_vocaltractcompute(x,dodisp)
global DIVA_x;
if nargin<2, dodisp=0; end

  if ~DIVA_x.debug
      x=x.*DIVA_x.params.Input.Scale;
      if nargout==1&&~dodisp
          [y,z,Outline,af]=diva_synth(x);
          y=y./DIVA_x.params.Output(1).Scale;
          DIVA_x.af=af;
          %display(af);
          if(sum(af<=0)<0)
              return;
          end
      else
          [y,z,Outline,af]=diva_synth(x);
          y=y./DIVA_x.params.Output(1).Scale;
          z=z./DIVA_x.params.Output(2).Scale;
          DIVA_x.af=af;
          %display(af);
          if(sum(af<=0)<0)
              return;
          end
      end
      
      if DIVA_x.gui&&nargin>1&&dodisp % display vt
          if dodisp==-1, % initialize display
              DIVA_x.figure.handles.h1=plot(nan,nan,'k','color',.75*DIVA_x.color(2,:),'linewidth',3);
              %hold on; DIVA_x.figure.handles.h1=[plot(nan,nan,'k','color',.75*DIVA_x.color(2,:),'linewidth',4);plot(nan,nan,'w-','linewidth',2)]; hold off
              axis equal; 
              set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[],'xdir','reverse');
          end
          set(DIVA_x.figure.handles.h1,'xdata',real(Outline),'ydata',imag(Outline));
          drawnow;
      end
  else
      L=ones(1,numel(x))*4/max(1,sum(sin(linspace(0,pi,numel(x)))));
      %L=(numel(x):-1:1)/numel(x);L=L*4/sum(L.*sin(linspace(0,pi,numel(x))));
      %L=[1,1];for n1=3:numel(x),L(n1)=L(n1-2)+L(n1-1);end;L=fliplr(L);L=L*sqrt(2*4^2)/sum(L);
      %L=ones(1,numel(x))*sqrt(2*5^2)/numel(x);
      ang=cumsum(x)'+linspace(0,pi,numel(x));
      px=[0,cumsum(cos(ang).*L)]+0;
      py=[0,cumsum(sin(ang).*L)]-0*4;
      yz_x=px(end);
      yz_y=py(end);
      
      y = [yz_x;yz_y];
      z = x(:);
      
      if DIVA_x.gui&&nargin>1&&dodisp % display vt
          if dodisp==-1, % initialize display
              DIVA_x.figure.data.x0=diag([1,.3])*[-.2,1.2,1.2,-.2;-.5,-.5,.5,.5];%[.5+.6*cos(linspace(0,2*pi,64));sin(linspace(0,2*pi,64))];
              for n1=1:numel(ang),DIVA_x.figure.handles.h1(n1)=patch(nan,nan,'k','edgecolor','k','facecolor','w'); hold on; end
              DIVA_x.figure.handles.h2=plot(nan,nan,'ko','markersize',4,'markeredgecolor','k','markerfacecolor','k');
              DIVA_x.figure.handles.h3=plot(nan,nan,'k-','color',1*[0,0,1],'linewidth',2);
              hold off;
              set(gca,'xlim',[-5,5],'ylim',[-5,5]);
              set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[]);
          end
          for n1=1:numel(ang),set(DIVA_x.figure.handles.h1(n1),'xdata',px(n1)+L(n1)*[cos(ang(n1)),-sin(ang(n1))]*DIVA_x.figure.data.x0,'ydata',py(n1)+L(n1)*[sin(ang(n1)),cos(ang(n1))]*DIVA_x.figure.data.x0);end
          set(DIVA_x.figure.handles.h2,'xdata',px,'ydata',py);
          set(DIVA_x.figure.handles.h3,'xdata',[get(DIVA_x.figure.handles.h3,'xdata'),px(end)],'ydata',[get(DIVA_x.figure.handles.h3,'ydata'),py(end)]);
          drawnow;
      end
  end
end
