function varargout = plot(f, varargin)
%PLOT   Plot a TREEFUN3.
%   PLOT(F) gives a 3D plot of the TREEFUN3 F and shows the tree on
%   which F is defined.
%
%   See also SURF, MESH.

func = varargin{1};

holdState = ishold();
nplotpts = 51;

h = instantiateSlice3GUI();
handles = guihandles(h);

% for k = 1:length(ids)
%     id = ids(k);
%     [x, y] = meshgrid(linspace(f.domain(1,id), f.domain(2,id), nplotpts), ...
%                       linspace(f.domain(3,id), f.domain(4,id), nplotpts));
%     u = coeffs2plotvals(f.coeffs{id});
%     hk = surface(x, y, 0*u, u, 'EdgeAlpha', 1, varargin{:});
%     if ( nargout > 0 )
%         h(k) = hk; %#ok<AGROW>
%     end
% end

[xx, yy, zz] = meshgrid(linspace(f.domain(1), f.domain(2), nplotpts), ...
                        linspace(f.domain(3), f.domain(4), nplotpts), ...
                        linspace(f.domain(5), f.domain(6), nplotpts));
% v = func( xx, yy, zz);
nd = numel(func);
v = zeros(size(xx));
for k = 1:nd
  v = v + func{k}(xx,yy,zz);
end
if isreal(v)
    [row,col,tube] = ind2sub(size(v), find(v(:) == max(v(:)), 1, 'last'));
else
    [row, col, tube] = ind2sub(size(v), find(abs(v(:)) == max(abs(v(:))), 1, 'last'));    
end
xslice = xx(row,col,tube); 
yslice = yy(row,col,tube); 
zslice = zz(row,col,tube); 

set(handles.xSlider, 'Min', f.domain(1));
set(handles.xSlider, 'Max', f.domain(2));
set(handles.xSlider, 'Value', xslice);

set(handles.ySlider, 'Min', f.domain(3));
set(handles.ySlider, 'Max', f.domain(4));
set(handles.ySlider, 'Value', yslice);

set(handles.zSlider, 'Min', f.domain(5));
set(handles.zSlider, 'Max', f.domain(6));
set(handles.zSlider, 'Value', zslice);

nSteps = 15; % number of slices allowed
set(handles.xSlider, 'SliderStep', [1/nSteps , 1 ]);
set(handles.ySlider, 'SliderStep', [1/nSteps , 1 ]);
set(handles.zSlider, 'SliderStep', [1/nSteps , 1 ]);

% Choose default command line output for the slice command:
handles.xx = xx;
handles.yy = yy;
handles.zz = zz;
handles.xslice = xslice;
handles.yslice = yslice;
handles.zslice = zslice;
handles.v = v;

if isreal(v)
    slice(xx, yy, zz, v, xslice, yslice, zslice)
    shading interp
    colorbar
else
    hh = slice(xx, yy, zz, angle(-v), xslice, yslice, zslice); 
    set(hh, 'EdgeColor','none')
    caxis([-pi pi]),
    colormap('hsv')
    axis('equal')     
end

hold on

% Plot the boxes
ids = leaves(f);
xdata = [f.domain([1 2 2 1 1 1 2 2 1 1], ids) ; nan(1, length(ids)); ... % bottom & top
         f.domain([2 2 2 2 1 1], ids) ; nan(1, length(ids));]; % the rest to complete the box
ydata = [f.domain([3 3 4 4 3 3 3 4 4 3], ids) ; nan(1, length(ids)); ...
         f.domain([3 3 4 4 4 4], ids) ; nan(1, length(ids));];
zdata = [f.domain([5 5 5 5 5 6 6 6 6 6], ids) ; nan(1, length(ids)); ...
         f.domain([5 6 6 5 5 6], ids) ; nan(1, length(ids));];
line('XData', xdata(:), 'YData', ydata(:), 'ZData', zdata(:), 'LineWidth', 1)
hold off

axis equal
xlim(f.domain(1:2, 1))
ylim(f.domain(3:4, 1))
zlim(f.domain(5:6, 1))

%
handles.xdata = xdata(:);
handles.ydata = ydata(:);
handles.zdata = zdata(:);

% Update handles structure
guidata(h, handles);
handles.output = handles.xSlider;

if ( ~holdState )
    hold off
end

if ( nargout > 0 )
    varargout = {h};
end

% Force the figure to clear when another plot is drawn on it so that GUI
% widgets don't linger.  (NB:  This property needs to be reset to 'add' every
% time we change the plot using a slider; otherwise, the slider movement will
% itself clear the figure, which is not what we want.)
set(h, 'NextPlot', 'replacechildren');

end

function h = instantiateSlice3GUI()

% Load up the GUI from the *.fig file.
installDir = treefunroot();
h = openFigInCurrentFigure([installDir '/@treefun3/slice.fig']);

% Do any required initialization of the handle graphics objects.
G = get(h, 'Children');
for (i = 1:1:length(G))
    if ( isa(G(i), 'matlab.ui.control.UIControl') )
        % Adjust the background colors of the sliders.
        if ( strcmp(G(i).Style, 'slider') )
            if ( isequal(get(G(i), 'BackgroundColor'), ...
                    get(0, 'defaultUicontrolBackgroundColor')) )
                set(G(i), 'BackgroundColor', [.9 .9 .9]);
            end
        end
        % Register callbacks.
        switch ( G(i).Tag )
            case 'xSlider'
                G(i).Callback = @(hObj, data) ...
                    xSlider_Callback(hObj, data, guidata(hObj));
            case 'ySlider'
                G(i).Callback = @(hObj, data) ...
                    ySlider_Callback(hObj, data, guidata(hObj));
            case 'zSlider'
                G(i).Callback = @(hObj, data) ...
                    zSlider_Callback(hObj, data, guidata(hObj));
        end
    end
end

% Store handles to GUI objects so that the callbacks can access them. 
guidata(h, guihandles(h));

end

% --- Executes on xSlider movement.
function xSlider_Callback(hObject, eventdata, handles)
% hObject    handle to xSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nextPlot = get(hObject.Parent, 'NextPlot');
set(hObject.Parent, 'NextPlot', 'add');

xslice = get(hObject, 'Value');         %returns position of slider
yslice = get(handles.ySlider, 'Value'); %returns position of slider
zslice = get(handles.zSlider, 'Value'); %returns position of slider

% The next slice command clears the title, if there was any. So, get that
% and put it again afterwards.
tit = get(gca(), 'Title');
titText = tit.String;

if ( isreal(handles.v) )
    handles.slice = slice(handles.xx, handles.yy, handles.zz, handles.v, ...
        xslice, yslice, zslice);
    shading interp
    colorbar, 
else
    handles.slice = slice(handles.xx, handles.yy, handles.zz, angle(-handles.v), ...
        xslice, yslice, zslice); 
    set(handles.slice, 'EdgeColor','none')
    caxis([-pi pi]),
    colormap('hsv')
    axis('equal')
end
handles.line = line('XData', handles.xdata(:), 'YData', handles.ydata(:), 'ZData', handles.zdata(:), 'LineWidth', 1);
title(titText)
handles.output = hObject;

set(hObject.Parent, 'NextPlot', nextPlot);

end

function ySlider_Callback(hObject, eventdata, handles)
% --- Executes on ySlider movement.
% hObject    handle to ySlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nextPlot = get(hObject.Parent, 'NextPlot');
set(hObject.Parent, 'NextPlot', 'add');

yslice = get(hObject, 'Value');         %returns position of slider
xslice = get(handles.xSlider, 'Value'); %returns position of slider
zslice = get(handles.zSlider, 'Value'); %returns position of slider

% The next slice command clears the title, if there was any. So, get that
% and put it again afterwards.
tit = get(gca(), 'Title');
titText = tit.String;

if ( isreal(handles.v) )
    handles.slice = slice(handles.xx, handles.yy, handles.zz, handles.v, ...
        xslice, yslice, zslice);
    shading interp
    colorbar, 
else
    handles.slice = slice(handles.xx, handles.yy, handles.zz, angle(-handles.v), ...
        xslice, yslice, zslice); 
    set(handles.slice, 'EdgeColor','none')
    caxis([-pi pi]),
    colormap('hsv')
    axis('equal')    
end
handles.line = line('XData', handles.xdata(:), 'YData', handles.ydata(:), 'ZData', handles.zdata(:), 'LineWidth', 1);
title(titText)
handles.output = hObject;

set(hObject.Parent, 'NextPlot', nextPlot);

end

function zSlider_Callback(hObject, eventdata, handles)
% --- Executes on zSlider movement.
% hObject    handle to zSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nextPlot = get(hObject.Parent, 'NextPlot');
set(hObject.Parent, 'NextPlot', 'add');

zslice = get(hObject, 'Value');         %returns position of slider
xslice = get(handles.xSlider, 'Value'); %returns position of slider
yslice = get(handles.ySlider, 'Value'); %returns position of slider

% The next slice command clears the title, if there was any. So, get that
% and put it again afterwards.
tit = get(gca(), 'Title');
titText = tit.String;

if ( isreal(handles.v) )
    handles.slice = slice(handles.xx, handles.yy, handles.zz, handles.v, ...
        xslice, yslice, zslice);
    shading interp
    colorbar
else
    handles.slice = slice(handles.xx, handles.yy, handles.zz, angle(-handles.v), ...
        xslice, yslice, zslice); 
    set(handles.slice, 'EdgeColor','none')
    caxis([-pi pi]),
    colormap('hsv')
    axis('equal')    
end
handles.line = line('XData', handles.xdata, 'YData', handles.ydata, 'ZData', handles.zdata, 'LineWidth', 1);
title(titText)
handles.output = hObject;

set(hObject.Parent, 'NextPlot', nextPlot);

end

