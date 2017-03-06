function diva_callback_stopfcn(block)
global DIVA_x;

% import all log variables
logsout=evalin('base','logsout');
names1=logsout.who;
for n1=1:numel(names1)
    names2=eval(['logsout.',names1{n1},'.who;']);
    for n2=1:numel(names2),
        eval(['DIVA_x.logs.',names2{n2},'=logsout.',names1{n1},'.',names2{n2},'.Data;']);
        eval(['DIVA_x.logs.time=logsout.',names1{n1},'.',names2{n2},'.Time;']);
    end
end


% display
if DIVA_x.gui
    % collapse AuditorySomatosensory logs
    names=fieldnames(DIVA_x.logs);
    txt={'Auditory','Somatosensory'};
    idxtxt={};remtxt={};
    for n0=1:numel(txt),
        for n1=1:numel(names)
            idxtxt{n0,n1}=strfind(names{n1},txt{n0});
            if ~isempty(idxtxt{n0,n1}),remtxt{n0,n1}=names{n1}([1:idxtxt{n0,n1}(1)-1,idxtxt{n0,n1}(1)+numel(txt{n0}):numel(names{n1})]);
            else remtxt{n0,n1}=''; end
        end
    end
    for n1=1:numel(names),
        if ~isempty(remtxt{1,n1}),
            idx=strmatch(remtxt{1,n1},remtxt(2,:),'exact');
            if ~isempty(idx)
                newname=[names{n1}(1:idxtxt{1,n1}(1)-1),cat(2,txt{:}),names{n1}(idxtxt{1,n1}(1)+numel(txt{1}):numel(names{n1}))];
                DIVA_x.logs.(newname)=cat(2,DIVA_x.logs.(names{n1}),DIVA_x.logs.(names{idx(1)}));
            end
        end
    end    
    if DIVA_x.params.dosound&&DIVA_x.dosound
        s=diva_synth(diag(DIVA_x.params.Input.Scale)*DIVA_x.logs.ArticulatoryPosition','sound');
        sound(s,11025);
        set(DIVA_x.figure.handles.pl0,'xdata',1:numel(s),'ydata',s);
        set(DIVA_x.figure.handles.ax0,'xlim',[1,numel(s)],'ylim',[-1,1]);
    end
    diva_gui update_inputoutputplots;
end

