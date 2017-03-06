function diva_gui(option,varargin)
global DIVA_x;
if ~nargin, option='init'; end
%if ~isfield(DIVA_x,'guiinit'),DIVA_x.guiinit=1;end

switch(lower(option)),
    case {'init','simulation'}
        DIVA_x.gui=1;
        DIVA_x.model=gcs;
        idx=find(DIVA_x.model=='/');
        if ~isempty(idx),DIVA_x.model=DIVA_x.model(1:idx(1)-1);end
        if ~strncmpi(DIVA_x.model,'diva',4)
            open diva;
            DIVA_x.model='diva';
        end
        DIVA_x.params=diva_vocaltract;
        DIVA_x.color=[.37,.74,1;1,1,1];
        DIVA_x.production_list=diva_targets('list');
        if ~isfield(DIVA_x,'production'), DIVA_x.production='new@random'; end
        if ~isfield(DIVA_x,'cycles'),DIVA_x.cycles=10; end
        if ~isfield(DIVA_x,'dosound'), DIVA_x.dosound=1; end
        DIVA_x.changed=0;

        if ~isfield(DIVA_x,'figure')||isempty(DIVA_x.figure.handles.figure)||~ishandle(DIVA_x.figure.handles.figure),
            DIVA_x.figure.handles.figure=figure('units','norm','position',[.1,.55,.8,.4],'menubar','none','name',['diva_gui (model: ',DIVA_x.model,')'],'numbertitle','off','color',DIVA_x.color(2,:),'tag','diva_gui','closerequestfcn','diva_gui(''close'')');
        else
            figure(DIVA_x.figure.handles.figure);
            clf;
        end
        
        % gui initialization
        uicontrol('units','norm','position',[0,.9,1,.1],'style','frame','backgroundcolor',DIVA_x.color(1,:));
        uicontrol('units','norm','position',[0,0,1,.1],'style','frame','backgroundcolor',DIVA_x.color(1,:));
        DIVA_x.figure.handles.button1=uicontrol('units','norm','position',[.1,.9,.15,.05],'style','togglebutton','string','simulation','value',1,'backgroundcolor','w','foregroundcolor','k','fontweight','bold','callback','diva_gui(''simulation'')');
        DIVA_x.figure.handles.button2=uicontrol('units','norm','position',[.25,.9,.15,.05],'style','togglebutton','string','targets','value',0,'backgroundcolor',DIVA_x.color(1,:),'foregroundcolor','k','fontweight','normal','callback','diva_gui(''targets'')');
        DIVA_x.figure.handles.list1=uicontrol('units','norm','position',[.75,.8,.175,.07],'style','popupmenu','string',{'Target trajectories','Target errors'},'foregroundcolor',.75*DIVA_x.color(2,:),'backgroundcolor',DIVA_x.color(2,:),'callback','diva_gui(''update_inputoutputplots'')');
        DIVA_x.figure.handles.list2=uicontrol('units','norm','position',[.1,.8,.175,.07],'style','popupmenu','string',{'Motor Command','FeedForward Command','FeedBack Command'},'foregroundcolor',.75*DIVA_x.color(2,:),'backgroundcolor',DIVA_x.color(2,:),'fontweight','normal','callback','diva_gui(''update_inputoutputplots'')');
        DIVA_x.figure.handles.buttonlist1=uicontrol('units','norm','position',[.75+.175,.8,.025,.07],'style','pushbutton','string','+','foregroundcolor',.75*DIVA_x.color(2,:),'backgroundcolor',DIVA_x.color(2,:),'fontweight','normal','callback','diva_gui(''extend_outputplots'')');
        DIVA_x.figure.handles.buttonlist2=uicontrol('units','norm','position',[.1+.175,.8,.025,.07],'style','pushbutton','string','+','foregroundcolor',.75*DIVA_x.color(2,:),'backgroundcolor',DIVA_x.color(2,:),'fontweight','normal','callback','diva_gui(''extend_inputplots'')');
        
        DIVA_x.figure.handles.button4=uicontrol('units','norm','position',[.4,.01,.07,.08],'style','pushbutton','string','Start','foregroundcolor','k','backgroundcolor',1-.5*(1-DIVA_x.color(1,:)),'fontweight','bold','callback','diva_gui(''start'')');
        DIVA_x.figure.handles.field1=uicontrol('units','norm','position',[.49,.05,.07,.04],'style','edit','string',num2str(DIVA_x.cycles),'backgroundcolor',DIVA_x.color(1,:));
        DIVA_x.figure.handles.text1=uicontrol('units','norm','position',[.56,.05,.1,.04],'style','text','string','learning cycles','backgroundcolor',DIVA_x.color(1,:),'foregroundcolor','k','horizontalalignment','left');
        DIVA_x.figure.handles.field2=uicontrol('units','norm','position',[.49,.01,.07,.04],'style','edit','string','1000','backgroundcolor',DIVA_x.color(1,:));
        DIVA_x.figure.handles.text2=uicontrol('units','norm','position',[.56,.01,.1,.04],'style','text','string','ms','backgroundcolor',DIVA_x.color(1,:),'foregroundcolor','k','horizontalalignment','left');
        list=cat(1,{'new@default';'new@random'},DIVA_x.production_list);
        idxvalue=strmatch(DIVA_x.production,list,'exact'); if isempty(idxvalue), idxvalue=1; DIVA_x.production=list{1}; end
        DIVA_x.figure.handles.list3=uicontrol('units','norm','position',[.05,.02,.15,.06],'style','popupmenu','string',list,'value',idxvalue,'backgroundcolor',DIVA_x.color(1,:),'callback','diva_gui(''load'')');
        DIVA_x.figure.handles.button5=uicontrol('units','norm','position',[.21,.02,.09,.06],'style','pushbutton','string','Save changes','foregroundcolor','k','backgroundcolor',1-.5*(1-DIVA_x.color(1,:)),'fontweight','normal','callback','diva_gui(''save'')','enable','off');
        DIVA_x.figure.handles.button6=uicontrol('units','norm','position',[.67,.02,.1,.06],'style','checkbox','string','Play sound','foregroundcolor','k','backgroundcolor',DIVA_x.color(1,:),'fontweight','normal','callback','diva_gui(''soundonoff'')','value',DIVA_x.dosound);
        
        % center plot
        DIVA_x.figure.handles.ax1=axes('units','norm','position',[.35,.2,.15,.6],'color',DIVA_x.color(2,:));
        diva_vocaltract('output',zeros(DIVA_x.params.Input.Dimensions,1),-1); 
        DIVA_x.figure.handles.ax0=axes('units','norm','position',[.48,.5,.12,.08],'color',DIVA_x.color(2,:));
        DIVA_x.figure.handles.pl0=plot([1,2],[0,0],'-','color',.75*[0,0,1]);axis off;
        set(DIVA_x.figure.handles.ax0,'xcolor',DIVA_x.color(2,:),'ycolor',DIVA_x.color(2,:));

        % Input/Output plots
        DIVA_x.params.Plots.Output={1,[]};%{1:numel(DIVA_x.params.Output(1).Plots_dim),[]};
        DIVA_x.params.Plots.Input={1};%{1:numel(DIVA_x.params.Input.Plots_dim)};
        diva_gui('init_inputoutputplots');
        if isfield(DIVA_x,'production'),
            diva_gui('load',DIVA_x.production);
        end

    case 'init_inputoutputplots',
        Vars={'Output','Input'};
        Coords={.75, .1};
        if nargin<2, vars=1:2; else vars=varargin{1}; end
        if nargin<3, ylim=[]; else ylim=varargin{2}; end
        for nvar=vars,
            var=Vars{nvar};
            coords=Coords{nvar};
            
            cs=[];N0=[];N1=[];
            for n0=1:numel(DIVA_x.params.Plots.(var)),
                for n1=1:numel(DIVA_x.params.Plots.(var){n0}),
                    cs=[cs,numel(DIVA_x.params.(var)(n0).Plots_dim{DIVA_x.params.Plots.(var){n0}(n1)})];
                    N0=[N0,n0];
                    N1=[N1,DIVA_x.params.Plots.(var){n0}(n1)];
                end
            end
            DIVA_x.params.Plots_.(var).setindex=N0;
            DIVA_x.params.Plots_.(var).plotindex=N1;
            cs=[0,cumsum(cs)/sum(cs)];
            if isfield(DIVA_x.figure.handles,'ax2')&&numel(DIVA_x.figure.handles.ax2)>=nvar&&any(ishandle(DIVA_x.figure.handles.ax2{nvar})),
                delete(DIVA_x.figure.handles.ax2{nvar}(ishandle(DIVA_x.figure.handles.ax2{nvar})));
            end
            if numel(N0)>0
                for n2=1:numel(N0),
                    n0=N0(n2);
                    n1=N1(n2);
                    idxdims=DIVA_x.params.(var)(n0).Plots_dim{n1};
                    DIVA_x.figure.handles.ax2{nvar}(n2)=axes('units','norm','position',[coords,.225+.575*cs(n2),.2,.5*(cs(n2+1)-cs(n2))],'color',DIVA_x.color(2,:));
                    set(gca,'xcolor',.75*DIVA_x.color(2,:),'ycolor',.75*DIVA_x.color(2,:));
                    if ~isempty(ylim)&&isfield(DIVA_x.params.(var)(n0),'Range'), set(gca,'ylim',[min(DIVA_x.params.(var)(n0).Range(idxdims,1)),max(DIVA_x.params.(var)(n0).Range(idxdims,2))]); end
                    ylabel(DIVA_x.params.(var)(n0).Plots_label{n1},'rotation',0,'horizontalalignment','right','interpreter','none');
                    if n2>1,set(gca,'xticklabel',[]);
                    else xlabel('time (ms)'); end
                end
            end
        end
            
    case 'update_inputoutputplots',
        for n2=1:numel(DIVA_x.params.Plots_.Output.plotindex),
            axes(DIVA_x.figure.handles.ax2{1}(n2));cla;
            hold on;
            n0=DIVA_x.params.Plots_.Output.setindex(n2);
            n1=DIVA_x.params.Plots_.Output.plotindex(n2);
            idxdims=DIVA_x.params.Output(n0).Plots_dim{n1};
            if n0>1, idxdims=idxdims+sum(cat(2,DIVA_x.params.Output(1:n0-1).Dimensions)); end
            plottype=get(DIVA_x.figure.handles.list1,'value');
            scale=cat(1,DIVA_x.params.Output(:).Scale);
            switch(plottype)
                case 1,
                    for n3=1:numel(idxdims),
                        patch([DIVA_x.logs.time;flipud(DIVA_x.logs.time)]*1000,[DIVA_x.logs.AuditorySomatosensoryTargetMax(:,idxdims(n3))*scale(idxdims(n3));flipud(DIVA_x.logs.AuditorySomatosensoryTargetMin(:,idxdims(n3))*scale(idxdims(n3)))],'k','facecolor',.9*[1,1,1],'edgecolor',.8*[1,1,1]);
                    end
                    plot(DIVA_x.logs.time*1000,DIVA_x.logs.AuditorySomatosensoryState(:,idxdims)*diag(scale(idxdims)),'-','linewidth',2);
                case 2,
                    plot(DIVA_x.logs.time*1000,DIVA_x.logs.AuditorySomatosensoryError(:,idxdims)*diag(scale(idxdims)),'-','linewidth',2);
            end
            if ~isempty(DIVA_x.logs.time)&&DIVA_x.logs.time(end)>0, set(gca,'xlim',[0,DIVA_x.logs.time(end)*1000]); end
            hold off;
            axis tight;
        end
        
        for n2=1:numel(DIVA_x.params.Plots_.Input.plotindex),
            axes(DIVA_x.figure.handles.ax2{2}(n2));cla;
            hold on;
            n0=DIVA_x.params.Plots_.Input.setindex(n2);
            n1=DIVA_x.params.Plots_.Input.plotindex(n2);
            idxdims=DIVA_x.params.Input(n0).Plots_dim{n1};
            if n0>1, idxdims=idxdims+sum(cat(2,DIVA_x.params.Input(1:n0-1).Dimensions)); end
            plottype=get(DIVA_x.figure.handles.list2,'value');
            scale=cat(1,DIVA_x.params.Input(:).Scale);
            switch(plottype)
                case 1,
                    plot(DIVA_x.logs.time*1000,DIVA_x.logs.ArticulatoryPosition(:,idxdims)*diag(scale(idxdims)),'-','linewidth',1);
                case 2,
                    plot(DIVA_x.logs.time*1000,DIVA_x.logs.FeedForwardMotorCommand(:,idxdims)*diag(scale(idxdims)),'-','linewidth',1);
                case 3,
                    plot(DIVA_x.logs.time*1000,DIVA_x.logs.FeedBackMotorCommand(:,idxdims)*diag(scale(idxdims)),'-','linewidth',1);
            end
            set(gca,'xlim',[0,DIVA_x.logs.time(end)*1000]);
            hold off;
            axis tight;
        end
        drawnow

    case 'extend_outputplots'
        txt=cat(2,reshape(DIVA_x.params.Output(1).Plots_label,1,[]),reshape(DIVA_x.params.Output(2).Plots_label,1,[]));
        idx=cat(2,1+zeros(1,numel(DIVA_x.params.Output(1).Plots_label)),2+zeros(1,numel(DIVA_x.params.Output(2).Plots_label)));
        value=zeros(size(idx));
        value(DIVA_x.params.Plots.Output{1})=1;
        value(numel(DIVA_x.params.Output(1).Plots_label)+DIVA_x.params.Plots.Output{2})=1;
        [s,ok]=listdlg('PromptString','Select signals to display:','SelectionMode','multiple','ListString',txt,'listsize',[150,200],'initialvalue',find(value));
        if ok,%&&~isempty(s),
            value(:)=0;
            value(s)=1;
            DIVA_x.params.Plots.Output={find(value(idx==1)),find(value(idx==2))};
        end
        diva_gui('init_inputoutputplots');
        diva_gui('update_inputoutputplots');
        
    case 'extend_inputplots'
        txt=reshape(DIVA_x.params.Input(1).Plots_label,1,[]);
        idx=1+zeros(1,numel(DIVA_x.params.Input(1).Plots_label));
        value=zeros(size(idx));
        value(DIVA_x.params.Plots.Input{1})=1;
        [s,ok]=listdlg('PromptString','Select signals to display:','SelectionMode','multiple','ListString',txt,'listsize',[150,200],'initialvalue',find(value));
        if ok,%&&~isempty(s),
            value(:)=0;
            value(s)=1;
            DIVA_x.params.Plots.Input={find(value(idx==1))};
        end
        diva_gui('init_inputoutputplots');
        diva_gui('update_inputoutputplots');

    case 'soundonoff'
        DIVA_x.dosound=get(DIVA_x.figure.handles.button6,'value');

    case 'softinit'
        DIVA_x.gui=1;
        if ~isfield(DIVA_x,'figure')||isempty(DIVA_x.figure.handles.figure)||~ishandle(DIVA_x.figure.handles.figure),
            diva_gui('init');
        else
            figure(DIVA_x.figure.handles.figure);
            if DIVA_x.debug
                set(DIVA_x.figure.handles.h3,'xdata',[],'ydata',[]);
            end
        end
        
    case 'close'
        if DIVA_x.changed&&~strcmp(DIVA_x.production,'new@default')&&~strcmp(DIVA_x.production,'new@random')
            answ=questdlg('Disregard changes (FeedForward learning)?','','Yes (Continue)','No (Save)','No (Save)');
            if strcmp(answ,'No (Save)')
                diva_gui save
            end
        end
        DIVA_x.gui=0;
        delete(gcbf);
        
    case 'start'
        nruns=str2double(get(DIVA_x.figure.handles.field1,'string'));
        DIVA_x.cycles=nruns;
        DIVA_x.changed=1;
        time=str2double(get(DIVA_x.figure.handles.field2,'string'));
        set_param(DIVA_x.model,'StopTime',num2str(time/1000))
        set(DIVA_x.figure.handles.button4,'string','Stop','callback','diva_gui(''stop'')');
        set([DIVA_x.figure.handles.button2,DIVA_x.figure.handles.field1,DIVA_x.figure.handles.field2,DIVA_x.figure.handles.list3,DIVA_x.figure.handles.button5],'enable','off');
        evalin('base','global DIVA_x');
        for nrun=1:nruns
            DIVA_x.simopt=simset('OutputVariables','t');
            simout=evalin('base','sim(DIVA_x.model,[],DIVA_x.simopt)');
            if simout(end)<time/1000 || ~strcmp(get(DIVA_x.figure.handles.button4,'string'),'Stop'), break; end
        end
        set(DIVA_x.figure.handles.button4,'string','Start','callback','diva_gui(''start'')');
        set([DIVA_x.figure.handles.button2,DIVA_x.figure.handles.field1,DIVA_x.figure.handles.field2,DIVA_x.figure.handles.list3,DIVA_x.figure.handles.button5],'enable','on');
        
    case 'stop'
        set(DIVA_x.figure.handles.button4,'string','Start','callback','diva_gui(''start'')');
        set([DIVA_x.figure.handles.button2,DIVA_x.figure.handles.field1,DIVA_x.figure.handles.field2,DIVA_x.figure.handles.list3,DIVA_x.figure.handles.button5],'enable','on');
        set_param(DIVA_x.model,'simulationcommand','stop') ;
        
    case 'load'
        if nargin>1,
            DIVA_x.production=varargin{1};
        else
            str=get(DIVA_x.figure.handles.list3,'string');
            value=get(DIVA_x.figure.handles.list3,'value');
            DIVA_x.production=str{value};
        end
        if strcmp(DIVA_x.production,'new@default')
            DIVA_x.production_info=diva_targets('new','txt');
            DIVA_x.production_art=[];
        elseif strcmp(DIVA_x.production,'new@random')
            DIVA_x.production_info=diva_targets('random','txt');
            DIVA_x.production_art=[];
        else
            DIVA_x.production_info=diva_targets('load','txt',DIVA_x.production);
            DIVA_x.production_art=diva_targets('load','mat',DIVA_x.production);
        end
        timeseries=diva_preparesimulation(DIVA_x.production_info,DIVA_x.production_art);
        maxt=max(timeseries.time)/1000;
        %DT=min(diff(timeseries.time))/1000;
        set_param(DIVA_x.model,'StopTime',num2str(.10+maxt));%.15+DT*size(Wart,1)));
        set(DIVA_x.figure.handles.field2,'string',num2str(1000*(.10+maxt)));
        DIVA_x.logs.time=timeseries.time/1000;
        DIVA_x.logs.AuditorySomatosensoryTargetMax=cat(2,timeseries.Aud_max*diag(1./DIVA_x.params.Output(1).Scale),timeseries.Som_max*diag(1./DIVA_x.params.Output(2).Scale));
        DIVA_x.logs.AuditorySomatosensoryTargetMin=cat(2,timeseries.Aud_min*diag(1./DIVA_x.params.Output(1).Scale),timeseries.Som_min*diag(1./DIVA_x.params.Output(2).Scale));
        DIVA_x.logs.AuditorySomatosensoryState=nan(size(DIVA_x.logs.AuditorySomatosensoryTargetMax));
        DIVA_x.logs.AuditorySomatosensoryError=nan(size(DIVA_x.logs.AuditorySomatosensoryTargetMax));
        DIVA_x.logs.FeedForwardMotorCommand=timeseries.Art*diag(1./DIVA_x.params.Input(1).Scale);
        DIVA_x.logs.ArticulatoryPosition=nan(size(DIVA_x.logs.FeedForwardMotorCommand));
        DIVA_x.logs.FeedBackMotorCommand=nan(size(DIVA_x.logs.FeedForwardMotorCommand));
        DIVA_x.changed=0;
        %set(DIVA_x.figure.handles.list1,'value',1);
        %set(DIVA_x.figure.handles.list2,'value',2);
        set(DIVA_x.figure.handles.button5,'enable','off');
        diva_gui update_inputoutputplots;
        
    case 'save'
        Art=diva_weightsadaptive('weights','diva_weights_SSM2FF.mat');
        if ~isempty(Art)
            if strcmpi(DIVA_x.production,'new@default')||strcmpi(DIVA_x.production,'new@random'),
                answ=inputdlg('Enter new target name','',1,{''});
                if isempty(answ), return;  end
                DIVA_x.production=answ{1};
                diva_targets('save','txt',DIVA_x.production,DIVA_x.production_info);
                diva_targets('save','mat',DIVA_x.production,DIVA_x.production_info,1);
                clear diva_weightsadaptive;
                reload=1;
            else reload=0; end
            diva_targets('save','art',DIVA_x.production,Art);
            DIVA_x.production_list=diva_targets('list');
            list=cat(1,{'new@default';'new@random'},DIVA_x.production_list);
            idxvalue=strmatch(DIVA_x.production,list,'exact'); if isempty(idxvalue), idxvalue=1; DIVA_x.production=list{1}; end
            set(DIVA_x.figure.handles.list3,'string',list,'value',idxvalue);
            if reload diva_gui load; end
        else 
            disp(['warning: weight matrix not loaded yet. Run simulation first']);
        end
        DIVA_x.changed=0;
        
    case 'targets'
        if DIVA_x.changed&&~strcmp(DIVA_x.production,'new@default')&&~strcmp(DIVA_x.production,'new@random')
            answ=questdlg('Disregard changes (FeedForward learning)?','','Yes (Continue)','No (Save)','No (Save)');
            if strcmp(answ,'No (Save)')
                diva_gui save
            end
        end
        % gui initialization
        idxremove=[DIVA_x.figure.handles.list1,DIVA_x.figure.handles.list2,DIVA_x.figure.handles.buttonlist1,DIVA_x.figure.handles.buttonlist2,DIVA_x.figure.handles.button4,DIVA_x.figure.handles.field1,DIVA_x.figure.handles.field2,DIVA_x.figure.handles.ax1,DIVA_x.figure.handles.ax0,DIVA_x.figure.handles.pl0,DIVA_x.figure.handles.ax2{2},DIVA_x.figure.handles.text1,DIVA_x.figure.handles.text2,DIVA_x.figure.handles.button6];
        delete(idxremove(ishandle(idxremove)));
        if isfield(DIVA_x.figure,'thandles')
            idxremove=[DIVA_x.figure.thandles.button6,DIVA_x.figure.thandles.button7,DIVA_x.figure.thandles.text1,DIVA_x.figure.thandles.field1,DIVA_x.figure.thandles.text2,DIVA_x.figure.thandles.field2,DIVA_x.figure.thandles.text3,DIVA_x.figure.thandles.field3,DIVA_x.figure.thandles.text4,DIVA_x.figure.thandles.field4,DIVA_x.figure.thandles.text5,DIVA_x.figure.thandles.field5,DIVA_x.figure.thandles.text6,DIVA_x.figure.thandles.field6,DIVA_x.figure.thandles.text7,DIVA_x.figure.thandles.field7];
            delete(idxremove(ishandle(idxremove)));
        end
        if ~isfield(DIVA_x,'resetff'), DIVA_x.resetff=1; end
        set(DIVA_x.figure.handles.button1,'value',0,'backgroundcolor',DIVA_x.color(1,:),'foregroundcolor','k','fontweight','normal');
        set(DIVA_x.figure.handles.button2,'value',1,'backgroundcolor','w','foregroundcolor','k','fontweight','bold');
        DIVA_x.production_list=diva_targets('list');
        value=DIVA_x.production;
        list=cat(1,{'new@default';'new@random'},DIVA_x.production_list);
        idxvalue=strmatch(value,list,'exact'); if isempty(idxvalue), idxvalue=1; else idxvalue=idxvalue(1); end
        set(DIVA_x.figure.handles.list3,'callback','diva_gui(''targetsgui_load'')','string',list,'value',idxvalue);
        set(DIVA_x.figure.handles.button5,'callback','diva_gui(''targetsgui_save'')');
        DIVA_x.figure.thandles.button6=uicontrol('units','norm','position',[.56,.02,.09,.06],'style','pushbutton','string','Delete target','foregroundcolor','k','backgroundcolor',1-.5*(1-DIVA_x.color(1,:)),'fontweight','normal','callback','diva_gui(''targetsgui_delete'')');
        DIVA_x.figure.thandles.button7=uicontrol('units','norm','position',[.31,.02,.24,.06],'style','checkbox','string','Reset learned feedforward sequence','foregroundcolor','k','backgroundcolor',DIVA_x.color(1,:),'fontweight','normal','callback','diva_gui(''targetsgui_reset'')','value',DIVA_x.resetff);
        
        DIVA_x.figure.thandles.text1=uicontrol('units','norm','position',[.05,.8,.1,.05],'style','text','string','Target name','backgroundcolor',DIVA_x.color(2,:),'horizontalalignment','left');
        DIVA_x.figure.thandles.field1=uicontrol('units','norm','position',[.05,.75,.1,.05],'style','edit','string','','backgroundcolor',DIVA_x.color(2,:),'horizontalalignment','left','callback','diva_gui(''update_targetdown'');');
        DIVA_x.figure.thandles.text2=uicontrol('units','norm','position',[.05,.65,.1,.05],'style','text','string','Target length (ms)','backgroundcolor',DIVA_x.color(2,:),'horizontalalignment','left');
        DIVA_x.figure.thandles.field2=uicontrol('units','norm','position',[.05,.6,.1,.05],'style','edit','string','','backgroundcolor',DIVA_x.color(2,:),'horizontalalignment','left','callback','diva_gui(''update_targetdown'');');
        DIVA_x.figure.thandles.text3=uicontrol('units','norm','position',[.45,.8,.09,.05],'style','text','string','Interpolation','backgroundcolor',DIVA_x.color(2,:),'horizontalalignment','right');
        DIVA_x.figure.thandles.field3=uicontrol('units','norm','position',[.55,.81,.1,.05],'style','popupmenu','string',{'nearest','linear','spline'},'value',3,'backgroundcolor',DIVA_x.color(2,:),'horizontalalignment','left','callback','diva_gui(''update_targetdown'');');

        str={};fieldidx=[];value=[];
        for n0=1:numel(DIVA_x.params.Output),
            for n1=1:numel(DIVA_x.params.Output(n0).Plots_label),
                if numel(DIVA_x.params.Output(n0).Plots_dim{n1})==1
                    str{end+1}=DIVA_x.params.Output(n0).Plots_label{n1};
                    value(end+1)=(n0==1 & any(DIVA_x.params.Output(n0).Plots_dim{1}==DIVA_x.params.Output(n0).Plots_dim{n1})); %default to first plot in auditory representation
                    fieldidx=cat(2,fieldidx,[n0;DIVA_x.params.Output(n0).Plots_dim{n1}]);
                end
            end
        end
        DIVA_x.figure.thandles.text4=uicontrol('units','norm','position',[.2,.8,.1,.05],'style','text','string','Target field(s)','backgroundcolor',DIVA_x.color(2,:),'horizontalalignment','left');
        DIVA_x.figure.thandles.field4=uicontrol('units','norm','position',[.2,.15,.1,.65],'style','listbox','string',str,'value',find(value),'backgroundcolor',DIVA_x.color(2,:),'horizontalalignment','left','max',2,'callback','diva_gui(''update_targetup'');');
        DIVA_x.figure.thandles.text5=uicontrol('units','norm','position',[.35,.8,.1,.05],'style','text','string','Control points (ms)','backgroundcolor',DIVA_x.color(2,:),'horizontalalignment','left');
        DIVA_x.figure.thandles.field5=uicontrol('units','norm','position',[.35,.75,.3,.05],'style','edit','string','','backgroundcolor',DIVA_x.color(2,:),'horizontalalignment','left','callback','diva_gui(''update_targetdown'');');
        DIVA_x.figure.thandles.text6=uicontrol('units','norm','position',[.35,.65,.3,.05],'style','text','string','Minimum value within target (for each control point)','backgroundcolor',DIVA_x.color(2,:),'horizontalalignment','left');
        DIVA_x.figure.thandles.field6=uicontrol('units','norm','position',[.35,.45,.3,.2],'style','edit','string','','backgroundcolor',DIVA_x.color(2,:),'horizontalalignment','left','max',2,'callback','diva_gui(''update_targetdown'');');
        DIVA_x.figure.thandles.text7=uicontrol('units','norm','position',[.35,.35,.3,.05],'style','text','string','Maximum value within target (for each control point)','backgroundcolor',DIVA_x.color(2,:),'horizontalalignment','left');
        DIVA_x.figure.thandles.field7=uicontrol('units','norm','position',[.35,.15,.3,.2],'style','edit','string','','backgroundcolor',DIVA_x.color(2,:),'horizontalalignment','left','max',2,'callback','diva_gui(''update_targetdown'');');
        
        if ~isfield(DIVA_x,'production'), DIVA_x.production='new@default'; end
        diva_gui('targetsgui_load',DIVA_x.production);
        
    case 'targetsgui_load'
        if nargin>1,
            DIVA_x.production=varargin{1};
        else
            str=get(DIVA_x.figure.handles.list3,'string');
            value=get(DIVA_x.figure.handles.list3,'value');
            DIVA_x.production=str{value};
        end
        if strcmp(DIVA_x.production,'new@default')
            DIVA_x.production_info=diva_targets('new','txt');
            names=get(DIVA_x.figure.thandles.field3,'string');
            DIVA_x.production_info.interpolation=names{get(DIVA_x.figure.thandles.field3,'value')};
            set(DIVA_x.figure.thandles.button6,'enable','off');
        elseif strcmp(DIVA_x.production,'new@random')
            DIVA_x.production_info=diva_targets('random','txt');
            names=get(DIVA_x.figure.thandles.field3,'string');
            DIVA_x.production_info.interpolation=names{get(DIVA_x.figure.thandles.field3,'value')};
            set(DIVA_x.figure.thandles.button6,'enable','off');
        else
            DIVA_x.production_info=diva_targets('load','txt',DIVA_x.production);
            set(DIVA_x.figure.thandles.button6,'enable','on');
        end
        diva_gui update_targetup;
        
    case 'targetsgui_save'
        production_name=get(DIVA_x.figure.thandles.field1,'string');
        if isempty(production_name), errordlg('Enter production name field',''); return; end
        DIVA_x.resetff=get(DIVA_x.figure.thandles.button7,'value');
        diva_targets('save','txt',production_name,DIVA_x.production_info);
        diva_targets('save','mat',production_name,DIVA_x.production_info,DIVA_x.resetff);
        DIVA_x.production_list=diva_targets('list');
        DIVA_x.production=production_name;
        value=DIVA_x.production;
        idxvalue=strmatch(value,DIVA_x.production_list,'exact'); if isempty(idxvalue), idxvalue=1; else idxvalue=idxvalue+2; end
        list=cat(1,{'new@default';'new@random'},DIVA_x.production_list);
        set(DIVA_x.figure.handles.list3,'string',list,'value',idxvalue);
        set(DIVA_x.figure.thandles.button6,'enable','on');
        
    case 'targetsgui_reset'
        DIVA_x.resetff=get(DIVA_x.figure.thandles.button7,'value');
        
    case 'targetsgui_delete'
        diva_targets('delete',DIVA_x.production);
        DIVA_x.production_list=diva_targets('list');
        list=cat(1,{'new@default';'new@random'},DIVA_x.production_list);
        set(DIVA_x.figure.handles.list3,'string',list,'value',min(numel(list),get(DIVA_x.figure.handles.list3,'value')));
        diva_gui targetsgui_load
        
    case 'update_targetup'
        set(DIVA_x.figure.thandles.field1,'string',DIVA_x.production_info.name);
        set(DIVA_x.figure.thandles.field2,'string',num2str(DIVA_x.production_info.length));
        set(DIVA_x.figure.thandles.field3,'value',strmatch(DIVA_x.production_info.interpolation,get(DIVA_x.figure.thandles.field3,'string'),'exact'));
        fieldnames=get(DIVA_x.figure.thandles.field4,'string');
        fieldname={fieldnames{get(DIVA_x.figure.thandles.field4,'value')}};
        a=[]; b={}; c={}; ok=1;
        set(DIVA_x.figure.thandles.field5,'foregroundcolor','k');
        set(DIVA_x.figure.thandles.field6,'foregroundcolor','k');
        set(DIVA_x.figure.thandles.field7,'foregroundcolor','k');
        for n1=1:numel(fieldname)
            if n1==1, a=DIVA_x.production_info.([fieldname{n1},'_control']); elseif n1>1&&~isequal(a,DIVA_x.production_info.([fieldname{n1},'_control'])), ok=0; end
            b{end+1}=num2str(DIVA_x.production_info.([fieldname{n1},'_min']));
            c{end+1}=num2str(DIVA_x.production_info.([fieldname{n1},'_max']));
            if size(DIVA_x.production_info.([fieldname{n1},'_min']),2)~=size(a,2), set(DIVA_x.figure.thandles.field6,'foregroundcolor','r'); end
            if size(DIVA_x.production_info.([fieldname{n1},'_max']),2)~=size(a,2), set(DIVA_x.figure.thandles.field7,'foregroundcolor','r'); end
        end
        if ~ok, disp(['warning: unequal control points across fields, displaying control points for first field only']); end
        set(DIVA_x.figure.thandles.field5,'string',num2str(a));
        set(DIVA_x.figure.thandles.field6,'string',b);
        set(DIVA_x.figure.thandles.field7,'string',c);
        if ok, 
            DIVA_x.params.Plots.Output={[],[]};
            DIVA_x.params.Plots.Input={[]};
            dims=cell(1,numel(DIVA_x.params.Output)); ndims=dims;
            for n2=1:numel(DIVA_x.params.Output),
                ndims{n2}=cellfun('length',DIVA_x.params.Output(n2).Plots_dim);
                dims{n2}=zeros([numel(fieldname),numel(DIVA_x.params.Output(n2).Plots_dim)]);
                for n1=1:numel(fieldname)
                    idx=strmatch(fieldname{n1},DIVA_x.params.Output(n2).Plots_label,'exact');
                    if ~isempty(idx), 
                        idx=idx(1);
                        dimsidx=DIVA_x.params.Output(n2).Plots_dim{idx};
                        for n3=1:numel(DIVA_x.params.Output(n2).Plots_dim),
                            if all(ismember(dimsidx,DIVA_x.params.Output(n2).Plots_dim{n3})),
                                dims{n2}(n1,n3)=1;
                            end
                        end
                    end
                end
                mdims=(ndims{n2}==sum(dims{n2},1));
                [sndims,idx]=sort(ndims{n2}.*mdims,'descend');
                for n1=1:numel(idx),
                    if mdims(idx(n1)),
                        DIVA_x.params.Plots.Output{n2}=cat(2,DIVA_x.params.Plots.Output{n2},idx(n1)); 
                        dims{n2}(dims{n2}(:,idx(n1))>0,:)=0;
                        mdims=(ndims{n2}==sum(dims{n2},1));
                    end
                end
            end
            diva_gui('init_inputoutputplots',1,'ylim');
            diva_gui update_targetplot; 
        end
        if isempty(DIVA_x.production_info.name), set([DIVA_x.figure.handles.button5,DIVA_x.figure.thandles.button7],'enable','off'); else set([DIVA_x.figure.handles.button5,DIVA_x.figure.thandles.button7],'enable','on'); end

    case 'update_targetdown'
        DIVA_x.production_info.name=get(DIVA_x.figure.thandles.field1,'string');
        value=str2num(get(DIVA_x.figure.thandles.field2,'string')); if isempty(value), value=DIVA_x.production_info.length; set(DIVA_x.figure.thandles.field2,'string',num2str(value)); end
        DIVA_x.production_info.length=value;
        names=get(DIVA_x.figure.thandles.field3,'string');
        DIVA_x.production_info.interpolation=names{get(DIVA_x.figure.thandles.field3,'value')};
        fieldnames=get(DIVA_x.figure.thandles.field4,'string');
        fieldname={fieldnames{get(DIVA_x.figure.thandles.field4,'value')}};
        valuea=str2num(char(get(DIVA_x.figure.thandles.field5,'string'))); 
        valueb=str2num(char(get(DIVA_x.figure.thandles.field6,'string'))); 
        valuec=str2num(char(get(DIVA_x.figure.thandles.field7,'string'))); 
        ok=1;
        if isempty(valuea), 
            ok=0;set(DIVA_x.figure.thandles.field5,'foregroundcolor','r');
        else
            set(DIVA_x.figure.thandles.field5,'foregroundcolor','k');
            for n1=1:numel(fieldname),
                DIVA_x.production_info.([fieldname{n1},'_control'])=valuea;
            end
        end
        if isempty(valueb)||size(valueb,2)~=size(valuea,2)||size(valueb,1)~=numel(fieldname), 
            ok=0;set(DIVA_x.figure.thandles.field6,'foregroundcolor','r');
        else
            set(DIVA_x.figure.thandles.field6,'foregroundcolor','k');
            for n1=1:numel(fieldname),
                DIVA_x.production_info.([fieldname{n1},'_min'])=valueb(n1,:);
            end
        end
        if isempty(valuec)||size(valuec,2)~=size(valuea,2)||size(valuec,1)~=numel(fieldname), 
            ok=0;set(DIVA_x.figure.thandles.field7,'foregroundcolor','r');
        else
            set(DIVA_x.figure.thandles.field7,'foregroundcolor','k');
            for n1=1:numel(fieldname),
                DIVA_x.production_info.([fieldname{n1},'_max'])=valuec(n1,:);
            end
        end
        if ok, diva_gui update_targetplot; end
        if isempty(DIVA_x.production_info.name), set([DIVA_x.figure.handles.button5,DIVA_x.figure.thandles.button7],'enable','off'); else set([DIVA_x.figure.handles.button5,DIVA_x.figure.thandles.button7],'enable','on'); end

    case 'update_targetplot'
        timeseries=diva_targets('timeseries',DIVA_x.production_info);
        timeseries.AudSom_min=cat(2,timeseries.Aud_min,timeseries.Som_min);
        timeseries.AudSom_max=cat(2,timeseries.Aud_max,timeseries.Som_max);
        for n2=1:numel(DIVA_x.params.Plots_.Output.plotindex),
            axes(DIVA_x.figure.handles.ax2{1}(n2));cla;
            hold on;
            n0=DIVA_x.params.Plots_.Output.setindex(n2);
            n1=DIVA_x.params.Plots_.Output.plotindex(n2);
            idxdims=DIVA_x.params.Output(n0).Plots_dim{n1};
            if n0>1, idxdims=idxdims+sum(cat(2,DIVA_x.params.Output(1:n0-1).Dimensions)); end
            for n3=1:numel(idxdims),
                patch([timeseries.time;flipud(timeseries.time)],[timeseries.AudSom_max(:,idxdims(n3));flipud(timeseries.AudSom_min(:,idxdims(n3)))],'k','facecolor',.9*[1,1,1],'edgecolor',.8*[1,1,1]);
            end
            if ~isempty(timeseries.time)&&timeseries.time(end)>0, set(gca,'xlim',[0,timeseries.time(end)]); end
            hold off;
        end
        
end
        
        
end

