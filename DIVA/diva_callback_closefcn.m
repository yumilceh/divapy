function diva_callback_closefcn
global DIVA_x;
if isfield(DIVA_x,'figure')&&~isempty(DIVA_x.figure.handles.figure)&&ishandle(DIVA_x.figure.handles.figure),close(DIVA_x.figure.handles.figure);end
DIVA_x=[];
