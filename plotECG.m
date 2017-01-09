function [ varargout ] = plotECG( varargin )
% -------------------------------------------------------
%
%    plotECG( X,Y, varargin ) - Plot Very Long, Multichannel Signals - Zoom & Scroll via Slider
%
%    Enables you to plot and zoom & scroll through signals with millions of samples. 
%    There is one slider for scrolling and one for zooming. 
%    To move forward or backward by exactly one screen width, 
%    use the scroll wheel or click on the slider trough. 
%    Clicking on the arrow buttons moves the signal by 1/100 
%    the screen width. Scrolling by mouse wheel works only if  
%    no interactive mode (Datacursor, Zoom, Pan, Edit Plot, etc.) is active. 
%
%    If you launch plotECG() without any arguments, 
%    one of two demos is shown. 
%
%    Ver. 1.1.1
%
%    Created:         Daniel Frisch        (15.02.2015)
%    Last modified:   Daniel Frisch        (09.01.2017)
%
%
% -------------------------------------------------------
%
%
%   Input Arguments
%     X                     Vector of timesteps, or scalar sample rate. scalar or [N x 1] double  
%     Y                     Signal vector, or signal matrix where each column represents one signal. [N x m] double   
%     LineSpec              LineSpec string, max. 4 characters, see Matlab plot() documentation. string (optional, default: '.') 
%                           If you define a marker here, set 'AutoMarkers' to 'none' so it won't be overridden. 
%
%     Key-Value Parameters  (optional):
%       Filter              Change parameters like cutoff frequencies via slider or text edit and see the results on the fly.
%                           Filter function handle for Y. function_handle or string (default: 'none') 
%                           Use one of the builtin filter functions as string, 'filter_FFT', or 'filter_bandpass_notch',  
%                           or define an own filter function and pass its name as string or function_handle (see below). 
%       mmPerSec            Initial zoom in screen millimeters per second. scalar double (default: 50) 
%                           If you change the figure size in this mode, the axis XLim will change such that mmPerSec stays the same. 
%       secPerScreenWidth   Initial zoom setting in seconds that are displayed on screen at a time. scalar double 
%                           If you change the figure size in this mode, the signal scale will change too, such that XLim stays the same.  
%       ShowAxisTicks       Shows the axis ticks and labels and a grid. string, 'on' or 'off' (default: 'on') 
%                           You can also change axis properties on the returned handles: hs.ax.XLabel.String = 'Seconds';
%       ShowInformationList Shows information about the location of the last clicks in the signal. function_handle or string 
%                           Specify one of 'none' (default), 'std_InformationList', 'ecg_InformationList', as string, 
%                           or define an own function like those with two inputs and cell output and pass its name or function_handle. 
%                           Mouse clicks are captured only if no interactive mode (pan, zoom etc.) is active in the figure. 
%       AutoStackSignals    Cell array of strings with signal names. Length must be equal to number of columns of Y, or 0. (default: {}) 
%                           Stacks the signals vertically, so for example a multipolar ECG can be shown.  
%                           Example: {'I','II','III', 'aVR','aVL','aVF', 'V1','V2','V3','V4','V5','V6'} 
%       SecondXAxisFunction Function handle that maps from the provided X values to a different x axis scale 
%                           that will be displayed above the plotted signal. function_handle or string (default: 'none') 
%                           Example: @(x)x/60^2, shows the time also in hours.  
%       SecondXAxisLabel    Label of second x axis. Example: 'Time in h'. string (default: '') 
%       YLimMode            'dynamic': dynamic y axis limits according to minimum and maximum of currently visible signal (default) 
%                           'fixed'  : fixed y axis limits according to minimum and maximum of entire signal. 
%                           In 'fixed' mode, you shouldn't apply filters that change the signal's mean much because YLim won't be updated. 
%                           You can change the fixed interval afterwards using the returned axes handle: hs.ax.YLim = [-10,10]; 
%       AutoMarkers         Automatically shows the specified marker at the sample locations if you zoom in such that there are 
%                           more than 5 pixes per sample. Use 'none' to disable this behaviour. string (default: '.')
%       ColorOrder          Custom ColorOrder for axes properties. [N x 3] double 
%       ScrollCallback      Function handle that will be called when the signal is scrolled or zoomed. function_handle or string (default: 'none') 
%                           An information structure will be passed to the function, with the fields 
%                           - XLim (containing the left and right x-axis value) 
%                           - ContinuousValueChange (whether it was called during or at the end of sliding) 
%                           To keep the sliding smooth, you should not execute expensive code if ContinuousValueChange is true.  
%       Parent              Parent figure. (default: new figure is created)
%                           Use delete(findall(hs.panel)) or close the figure to delete the old plot and its signals. 
%                           The 'HandleVisibility' of the figure created by default is set to 'Callback' to prevent you from  
%                           plotting into it accidentially from the command line. 
%                           To close the created figures from the command line, you have to use "close all hidden" instead of "close all".  
%       Units               Units for positioning plotECG in the parent figure. (default: normalized)
%       Position            Position of plotECG in the parent figure. (default: [0,0,1,1]) 
%
%
%   Output Arguments
%     h                     Returns the chart line objects as a vector. Use h to modify a chart line after it is created. 
%     hs                    Returns a struct with handles to some more GUI objects for later modifiaction  
%
%     You can set Name-Value parameters on the returned chart line handles: 
%     set(h, 'LineWidth',3) etc.
%
%     Furthermore, you might want to change the X/YLabel: 
%     hs.ax.XLabel.String = 'Seconds'; or: xlabel(hs.ax,'Seconds')
%       
%
%
%   Example
%     X = 0:0.001:100-0.001;
%     Y = sin(2*pi*0.1*X) + sin(2*pi*X) + .1*sin(2*pi*50*X);
%     [h,hs] = plotECG(X,Y);
%     [h,hs] = plotECG(1000,Y, 'Filter','filter_FFT');
%
%
%   Filter Function
%     A filter function (@filter) can be specified that modifies Y. 
%     It must have at least two inputs - X and Y. X can be the scalar 
%     sampling rate or a time vector; Y can be a vector containing one signal
%     or a matrix with multiple signals along its columns. 
%     In addition, @filter can have any number of scalar inputs that specify 
%     alterable parameters like cutoff frequencies of filters. 
%     For each of these parameters, a slider and a text edit appears in the 
%     plot figure, where the value can be adjusted with direct feedback. 
%     The filter function file can be modified while the plot window is open. 
%     A refiltering can be triggered by hitting enter in a text edit. 
%     @filter must have three output arguments: the modified X and Y and a
%     description struct containing information about the filter function
%     and its parameters. An example filter function with two parameters 
%     and its description struct is given below: 
%
%         function [ time_out, sig_out, description ] = filter_bandpass( time_in, sig_in, highpass, lowpass )
%             % --- Description struct ---
%             % Find samplerate
%             if isscalar(time_in), samplerate=time_in; else samplerate=(length(time_in)-1)/(time_in(end)-time_in(1)); end
%             % Filter title
%             description.string = 'Butterworth bandpass';
%             % Slider for lower cutoff frequency
%             description.slider(1).Label  = 'highpass';
%             description.slider(1).Min    = 0.001;
%             description.slider(1).Value  = 0.2;
%             description.slider(1).Max    = 4;
%             % Slider for higher cutoff frequency
%             description.slider(2).Label  = 'lowpass';
%             description.slider(2).Min    = min(  5,samplerate*0.01);
%             description.slider(2).Value  = min( 50,samplerate/2.01);
%             description.slider(2).Max    = min(200,samplerate/2.01);
%             % If an empty or too short signal is given, just return the description struct
%             if size(sig_in,1)<=12; time_out=time_in; sig_out=NaN*sig_in; return; end
% 
%             % --- Signal processing ---
%             bpFilt = designfilt('bandpassiir', ...
%                 'FilterOrder',2, ...
%                 'HalfPowerFrequency1',highpass, ...
%                 'HalfPowerFrequency2',lowpass, ...
%                 'SampleRate',samplerate);
%             sig_out = filtfilt(bpFilt,sig_in);
%             time_out = time_in;
%         end
%
%
%
%
%   Changes
%     Ver. 1.1.1
%        The Matlab Data Cursor works now too, if you switch it on. Of course you can still use the built-in InformationList. 
%        (However the scroll wheel does not work as long as Data Cursor Mode is on.) 
%        Replaced the "<HTML>" tags with "<html>" because MATLAB Online uicontrol listbox doesn't support uppercase HTML tags. 
%        Replaced the dash with a hyphen in the file's documentation because MATLAB editor on Mac doesn't support unicode. 
%        
%
%
%
%   TODO
%     Export tight images (PDF, PNG) via menu bar
%     Multiple filters for signal
%     Add scaling information bars into figure
%     Currently no line if panned by hand instead of slider
%     Enable signal data update during runtime (return function handle) 
%     Test thoroughly with MATLAB Online and Linux 
%
%
%




    

    
%% Dependencies
% [flist,plist] = matlab.codetools.requiredFilesAndProducts('plotECG.m'); [flist'; {plist.Name}']

% plotECG(): no dependencies
% built-in local filter function 'filter_bandpass_notch': Signal Processing Toolbox




%% Parse Inputs

args = varargin;

iSig = 1; 
SIG = struct;

defaultXLabel = 'Time in s';
mm_default = 50;
if length(args)>=2 && isnumeric(args{1}) && isnumeric(args{2})
    % plotECG(X,Y)
    SIG.X0 = args{1};
    SIG.Y0 = args{2};
    args = args(3:end);
elseif length(args)>=1 && isnumeric(args{1})
    % plotECG(Y)
    defaultXLabel = '';
    SIG.Y0 = args{1};
    SIG.X0 = 1;
    args = args(2:end);
    mm_default = 0.001;
else
    % plotECG()
    show_demo();
    return;
end

% Check X and Y and change to column-oriented data
if numel(SIG.X0)==1 % X = sample rate
    validateattributes(SIG.X0,{'double'},{'nonempty','real','finite','scalar','positive'     }, 'plotECG','scalar X',1)
else % X = timestamps
    validateattributes(SIG.X0,{'double'},{'nonempty','real','finite','vector','nondecreasing'}, 'plotECG','vector X',1)
end
    validateattributes(SIG.Y0,{'double'},{'nonempty','real','2d'}, 'plotECG','Y')

if size(SIG.X0,1)==1, SIG.X0=SIG.X0'; end % change to column vector 
if isrow(SIG.Y0) || ~isscalar(SIG.X0) && size(SIG.Y0,1)~=size(SIG.X0,1)
    SIG.Y0=SIG.Y0';  % each column must be one signal
end

assert(isscalar(SIG.X0) || size(SIG.Y0,1)==size(SIG.X0,1),'Y must have dimensions such that one of its dimensions equals length(X) if X is not a scalar')
assert(size(SIG.Y0,1)>1, 'Signal must have at least two samples')
assert(nnz(~isfinite(SIG.Y0))<numel(SIG.Y0),'Y completely consists of infinite values')
if nnz(~isfinite(SIG.Y0))>.5*numel(SIG.Y0), warning('Y contains %.2f %% infinite values\n',nnz(~isfinite(SIG.Y0))/numel(SIG.Y0)*100); end
assert(nargout<=2,'Maximum 2 (instead of %u) output arguments are possible',nargout)

% % Works already with piecewise equally-spaced signals
% % TODO test with non-equally spaced signals
% if ~isscalar(SIG.X0)
%     x0_diff = diff(SIG.X0);
%     variation = (max(x0_diff)-min(x0_diff)) / max(abs(SIG.X0));
%     assert(variation<eps*1000, 'X must be equally spaced.')
% end

if isscalar(SIG.X0)
    period = 1/SIG.X0;
else
    period = median(diff(SIG.X0),'omitnan');
    %period = (SIG.X0(end)-SIG.X0(1))/(length(SIG.X0)-1);
end
assert(period>0,'The sampling period must be > 0')
if isscalar(SIG.X0)
    timeBoundary = [0, period*(size(SIG.Y0,1)-1)];
else
    xFinite = isfinite(SIG.X0);
    timeBoundary = [SIG.X0(find(xFinite,1,'first')),SIG.X0(find(xFinite,1,'last'))];
end
duration = diff(timeBoundary);
assert(duration>0,'The duration must be > 0')

lineSpec = '-';
if size(SIG.Y0,2)==1
    lineSpec = '-k'; 
end    
if length(args)>=1 && isLineSpec(args{1})
    lineSpec = args{1};
    args = args(2:end);
end

parser = inputParser;
parser.FunctionName = 'plotECG';
%parser.KeepUnmatched = true; % additional Matlab plot arguments 

% Classes: 
%    numeric; double; single; int8; int16; int32; int64; uint8; uint16; uint32; uint64; 
%    logical; char; struct; cell; function_handle;  
% Attributes:  
%    binary; integer; even; odd; real; finite; nonnan; 
%    positive; nonnegative; nonzero; >,N; >=,N; <,N; <=,N;  
%    decreasing; increasing; nondecreasing; nonincreasing; 
%    nonempty, numel,N; ncols,N; nrows,N; ndims,N; square, nonsparse; 
%    scalar; vector; column; row; 2d; 3d; size,[d1,...,dN]; 
%
parser.addParameter('Filter'             , 'none'      , @(x)validateattributes(x,{'char','function_handle'},{'vector','nonempty'}))
parser.addParameter('mmPerSec'           , mm_default  , @(x)validateattributes(x,{'double'},{'real','finite','positive','scalar'}))
parser.addParameter('secPerScreenWidth'  , 1           , @(x)validateattributes(x,{'double'},{'real','finite','positive','scalar'}))
parser.addParameter('ShowAxisTicks'      , 'on'        , @(x)any(validatestring(x,{'on','off'})))
parser.addParameter('ShowInformationList', 'none'      , @(x)validateattributes(x,{'char','function_handle'},{'vector','nonempty'}))
parser.addParameter('AutoStackSignals'   , {}          , @(x)iscellstr(x))
parser.addParameter('SecondXAxisFunction', 'none'      , @(x)validateattributes(x,{'char','function_handle'},{'vector','nonempty'}))
parser.addParameter('SecondXAxisLabel'   , ''          , @(x)validateattributes(x,{'char'},{})) 
parser.addParameter('YLimMode'           , 'dynamic'   , @(x)any(validatestring(x,{'dynamic','fixed'})))
parser.addParameter('AutoMarkers'        , '.'         , @(x)any(validatestring(x,{'+','o','*','.','x','square','diamond','v','^','>','<','pentagram','hexagram','none'})))
parser.addParameter('ColorOrder'         , []          , @(x)validateattributes(x,{'double'},{'real','finite','nonnegative', '<=',1, 'size',[NaN,3]}))
parser.addParameter('ScrollCallback'     , 'none'      , @(x)validateattributes(x,{'char','function_handle'},{'vector','nonempty'}))
parser.addParameter('Parent'             , 0           , @(x)isscalar(x) && isgraphics(x) && x~=0)
parser.addParameter('Units'              , 'normalized', @(x)any(validatestring(x,{'pixels','normalized','inches','centimeters','points','characters'})))
parser.addParameter('Position'           , [0,0,1,1]   , @(x)validateattributes(x,{'double'},{'real','finite','nonnegative', 'size',[1 4]}))

parser.parse(args{:})

if ~strcmp(parser.Results.Filter,'none')
    % Check if function works with empty signal input
    feval(parser.Results.Filter,1,[]);
    % Check if function has three output arguments
    assert(nargout(parser.Results.Filter)==3,'Filter %s must have 3 output arguments (instead of %u)',func2str2(parser.Results.Filter),nargout(parser.Results.Filter))
end

SIG.AutoStackSignals = parser.Results.AutoStackSignals;
if ~isempty(SIG.AutoStackSignals)
    nStr = numel(SIG.AutoStackSignals);
    nSig = size(SIG.Y0,2);
    assert(nStr==nSig, 'You specified %u Strings in ''AutoStackSignals'', but the number of signals in Y is %u.',nStr,nSig);
end
        

%% Add the GUI components

mmPerSec = [];
secPerScreenWidth = [];
if ismember('secPerScreenWidth',parser.UsingDefaults)
    mmPerSec = parser.Results.mmPerSec;
else
    secPerScreenWidth = parser.Results.secPerScreenWidth;
end
mmPerSec_slider = mmPerSec;
% N = size(SIG.Y0,1); % number of data points
N = round(duration / period)+1; 
axWidth_cm = 100;
axWidth_px = 100;

% Layout constants
fontSize = 11;
units = 'centimeters';
space = 0.05;
sliderHeight = .5;
checkboxHeight = 1;
editWidth = 2;

% Add components, save handles in a struct
hs = struct;
if parser.Results.Parent==0
    % Create new figure
    % TODO make initial figure height dependent on number of filters; voltage scaling  
    hs.parent = figure('Units','normalized', 'OuterPosition',[0.3,0.53,0.65,0.45], 'HandleVisibility','Callback');
else
    hs.parent = parser.Results.Parent;
end

% Find parent figure
hs.fig = hs.parent;
while ~isempty(hs.fig) && ~strcmp('figure', get(hs.fig,'type'))
    hs.fig = get(hs.fig,'parent');
end

% Disable all interactive modes. 
% Only then the WindowScrollWheelFcn can be set. 
rotate3d(hs.fig,'off')
zoom(hs.fig,'off')
pan(hs.fig,'off')
% Now set custom WindowScrollWheelFcn
hs.fig.WindowScrollWheelFcn = @figScroll;

hs.panel = uipanel(... % This uipanel can be put into another GUI
    'Parent',hs.parent,...
    'Units',parser.Results.Units,...
    'Position',parser.Results.Position,...
    'BorderWidth',0,...
    'SizeChangedFcn',@resizegui,... 
    'Visible','off');

hs.ax2 = axes(... 
    'Parent',hs.panel,...
    'ActivePositionProperty','Position',...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'YTickLabel',{''},...
    'Color','none');
if strcmp(parser.Results.SecondXAxisFunction,'none')
    set(hs.ax2,'Visible','off')
end

hs.ax = axes(...
    'Parent',hs.panel,...
    'TickLabelInterpreter','none',...
    'ActivePositionProperty','Position');
if ~isempty(parser.Results.ColorOrder)
    hs.ax.ColorOrder = parser.Results.ColorOrder;
    hs.ax.NextPlot = 'replacechildren';
end

hs.scroll = uicontrol(...
    'Parent',hs.panel,...
    'Style','slider',...
    'Min',0,...
    'Value',0,...
    'Max',1,...
    'SliderStep',[1e-4,.07],... 
    'TooltipString','Click the slider trough or ouse mouse wheel to page forward one screen width.',...
    'Interruptible','on',...
    'Callback',{@redraw,true});
    hListener = addlistener(hs.scroll,'ContinuousValueChange',@redraw); 
    setappdata(hs.scroll,'sliderListener',hListener);
    
hs.zoom = uicontrol(...
    'Parent',hs.panel,...
    'Style','slider',...
    'Min',0,...
    'Value',.5,...
    'Max',1,...
    'SliderStep',[1e-4,.07],...
    'TooltipString','Zoom',...
    'Interruptible','on',...
    'Callback',{@redraw,true});
    hListener = addlistener(hs.zoom,'ContinuousValueChange',@zoom_callback); 
    setappdata(hs.zoom,'sliderListener',hListener);
    
    function zoom_callback(varargin)
        redraw(varargin{:});
        mmPerSec_slider = (axWidth_cm*10)/(numPoints*period);
    end
d = -log(N/7);

    
hs.list = uicontrol(...
    'Parent',hs.panel,...
    'Style','listbox',...
    'Max',2,...
    'FontSize',fontSize,...
    'Value',[]);
if strcmpi(parser.Results.ShowInformationList,'none')
    hs.list.Visible = 'off';
end

    


%% Add Filter Controls

if ~strcmp(parser.Results.Filter,'none')  % a filter function was specified 
    
    SIG.filter.hFilter = parser.Results.Filter;
    
    % Get the filter description by passing it an empty signal
    dummy = cell(1,nargin(SIG.filter.hFilter));
    dummy{1} = SIG.X0; % for samplerate detection
    [~,~,description] = feval(SIG.filter.hFilter, dummy{:}); 
    
    % Find the function's file
    filterFile = function_file(SIG.filter.hFilter);
    SIG.filter.panel = uipanel(... 
        'Parent',hs.panel,...
        'FontSize',fontSize,...
        'BorderType','none',...
        'Title',filterFile);
    nControls = nargin(SIG.filter.hFilter)-2;
    SIG.filter.mainCheck = uicontrol(...
        'Parent',SIG.filter.panel,...
        'Style','checkbox',...
        'Value',0,...
        'FontSize',fontSize,...
        'Callback',@updateEnabled,...
        'String',description.string);
    
    % Format of description
    % description.slider(1).Label  = 'f_low';
    % description.slider(1).Min    = 0.01;
    % description.slider(1).Value  = 0.5;
    % description.slider(1).Max    = 3;
    
    % Create the sliders
    % TODO check 'description' better
    assert(length(description.slider)==nControls,'Error in filter function ''%s'': 3rd output argument description.slider has %u entries where the function has %u user control inputs.',func2str2(SIG.filter.hFilter),length(description),nControls)
    for iControl = 1:nControls
        assert(description.slider(iControl).Min>0  , 'Due to logarithmic slider scaling, the slider ranges must be geater than zero.')
        assert(description.slider(iControl).Value>0, 'Due to logarithmic slider scaling, the slider ranges must be geater than zero.')
        assert(description.slider(iControl).Max>0  , 'Due to logarithmic slider scaling, the slider ranges must be geater than zero.')
        SIG.filter.slider(iControl).text = uicontrol(...
            'Parent',SIG.filter.panel,...
            'Style','text',...
            'FontSize',fontSize,...
            'String',description.slider(iControl).Label);
        SIG.filter.slider(iControl).edit = uicontrol(...
            'Parent',SIG.filter.panel,...
            'FontSize',fontSize,...
            'String',num2str(description.slider(iControl).Value,4),...
            'Style','edit'); 
        SIG.filter.slider(iControl).slider = uicontrol(...
            'Parent',SIG.filter.panel,...
            'Style','slider',...
            'FontSize',fontSize,...
            'Interruptible','on',...
            'Min',description.slider(iControl).Min,...
            'Max',description.slider(iControl).Max,...
            'SliderStep',[1e-4,.1]);
        editToSlider(SIG.filter.slider(iControl).edit,SIG.filter.slider(iControl).slider);
        
        % Assemble callbacks
        % Edit to slider; updateFilter
        set(SIG.filter.slider(iControl).edit, 'Callback', {@editToSliderUpdateFilter, SIG.filter.slider(iControl).edit, SIG.filter.slider(iControl).slider, iSig}); 
        % Slider to edit as continuous update; no updateFilter
        sliderToEdit2 = @(varargin) sliderToEdit(SIG.filter.slider(iControl).edit,SIG.filter.slider(iControl).slider);
        hListener = addlistener(SIG.filter.slider(iControl).slider,'ContinuousValueChange',sliderToEdit2); 
        setappdata(SIG.filter.slider(iControl).slider,'sliderListener',hListener);
        sliderToEdit2();
        % Slider to edit on mouse release; updateFilter
        set(SIG.filter.slider(iControl).slider, 'Callback', {@sliderToEditUpdateFilter, SIG.filter.slider(iControl).edit, SIG.filter.slider(iControl).slider, iSig});        
    end
    
    % Use default value for first calculation
    sliderToEditUpdateFilter([],[], SIG.filter.slider(iControl).edit, SIG.filter.slider(iControl).slider, iSig);
end
    
    

    
%% Create chart line handle

click_info.x = [];
click_info.y = [];

%hs.line = plot(hs.ax,1,1:size(SIG.Y0,2),lineSpec,parser.Unmatched); % Unmatched name-value pairs as plot parameters 
hs.line = plot(hs.ax,1,1:size(SIG.Y0,2),lineSpec); % Safer: Don't allow unmatched name-value pairs. Plot can still be modified by handles. 

% Remove 10^x axes factor
hs.ax2.XAxis.Exponent = 0;
hs.ax2.XAxis.ExponentMode = 'manual';
hs.ax2.YAxis.Exponent = 0;
hs.ax2.YAxis.ExponentMode = 'manual';
hs.ax.XAxis.Exponent = 0;
hs.ax.XAxis.ExponentMode = 'manual';
hs.ax.YAxis.Exponent = 0;
hs.ax.YAxis.ExponentMode = 'manual';

if strcmpi(parser.Results.ShowAxisTicks,'on')
    hs.ax.XLabel.String = defaultXLabel;
    hs.ax.TickLength = [0.001,0.001];
    hs.ax2.TickLength = [0.001,0.001];
    hs.ax2.XLabel.String = parser.Results.SecondXAxisLabel;
    hs.ax.XMinorGrid = 'on';
    if isempty(parser.Results.AutoStackSignals)
        hs.ax.YLabel.String = 'Voltage in mV';
    else
        hs.ax.YLabel.String = 'Channel';
    end        
else
    set(hs.ax ,'XTick',[], 'YTick',[])
    set(hs.ax2,'XTick',[], 'YTick',[])
end

if ~isempty(SIG.AutoStackSignals) && strcmp(parser.Results.YLimMode,'fixed')
    % Stack signals horizontally. 
    [sigPosVec,sigAddVec] = auto_stack_nooverlap(SIG.Y0);
    hs.ax.YTick = flip(sigPosVec);
    hs.ax.YTickLabel = flip(SIG.AutoStackSignals(:));
    hs.ax.TickLabelInterpreter = 'none';
else
    sigPosVec = zeros(1,size(SIG.Y0,2));
    sigAddVec = zeros(1,size(SIG.Y0,2));
end

Y0pos = bsxfun(@plus,SIG.Y0,sigAddVec);
range = [min(Y0pos(:)),max(Y0pos(:))];
dlt = diff(range)/50;
range(1) = range(1)-dlt;
range(2) = range(2)+dlt;
if strcmp(parser.Results.YLimMode,'fixed') && nnz(isnan(range))==0 && range(2)>range(1)
    hs.ax .YLimMode = 'manual';
    hs.ax2.YLimMode = 'manual';
    hs.ax. YLim = range; 
    hs.ax2.YLim = range; 
end

if ~strcmp(parser.Results.Filter,'none')
    updateEnabled()
end

% Make figure visible after adding components
btnDown()
redraw(true);
try getframe(hs.fig); catch, end % update system queue
resizegui

if ~isempty(mmPerSec)
    numPoints = axWidth_cm*10/(mmPerSec*period);
else
    numPoints = secPerScreenWidth/period;
end
zoomValue = log(numPoints/N)/d;
zoomValue = max(zoomValue,0);
zoomValue = min(zoomValue,1);
set(hs.zoom,'Value',zoomValue)
zoomValue = hs.zoom.Value;

resizegui
redraw(true)
hs.panel.Visible = 'on';




%% Updatefunctions

    function value = editToSliderUpdateFilter(ho,cbd, edit,slider,iSig)
        value = editToSlider(edit,slider);
        updateFilter(ho,cbd, iSig);
    end

    function value = sliderToEditUpdateFilter(ho,cbd, edit,slider,iSig)
        value = sliderToEdit(edit,slider);
        updateFilter(ho,cbd, iSig);
    end

    function updateFilter(ho,~, iSig)
        nCtrl = length(SIG.filter.slider);
        parameters = cell(nCtrl,1);
        for iCtrl = 1:nCtrl
            parameters{iCtrl} = str2double(SIG.filter.slider(iCtrl).edit.String);
        end
        function arrowset(fig)
            fig.Pointer = 'arrow';
        end
        if ~strcmp(hs.fig.Pointer,'watch')
            cleaner = onCleanup(@() arrowset(hs.fig));
            hs.fig.Pointer = 'watch'; 
            drawnow; 
        end
%         Filter debugging with "dbstop if error" works better without a try/rethrow in the stack. 
%         But for end users you can uncomment this try/catch block
%         so if they type improper values in the filter text boxes, the filter panel will turn red. 
%         try
            [SIG.X, SIG.Y, ~] = feval(SIG.filter.hFilter, SIG.X0, SIG.Y0, parameters{:});
            if nnz(isfinite(SIG.X))<nnz(SIG.X); warning('Filtered X contains %u nonfinite values\n',nnz(SIG.X)-nnz(isfinite(SIG.X))); end
            if nnz(isfinite(SIG.Y))<nnz(SIG.Y); warning('Filtered Y contains %u nonfinite values\n',nnz(SIG.Y)-nnz(isfinite(SIG.Y))); end
            assert(isreal(SIG.Y), 'Signal procuded by filter must be real')
            SIG.filter.panel.BackgroundColor = [.94 .94 .94];
            SIG.filter.mainCheck.BackgroundColor = [.94 .94 .94];
%         catch ex
%             fprintf('The filter produced an error:\n') 
%             SIG.filter.panel.BackgroundColor = 'red';
%             SIG.filter.mainCheck.BackgroundColor = 'red';
%             hs.fig.Pointer = 'arrow';
%             rethrow(ex)
%         end
        
        hs.fig.Pointer = 'arrow';
        if ~isempty(ho)
            redraw(true);
        end
    end

    function updateEnabled(varargin)
        if SIG.filter.mainCheck.Value == SIG.filter.mainCheck.Max
            enable = 'on';
        else
            enable = 'off';
        end
        nCtrl = length(SIG.filter.slider);
        for iCtrl = 1:nCtrl
            SIG.filter.slider(iCtrl).text.Enable = enable;
            SIG.filter.slider(iCtrl).edit.Enable = enable;
            SIG.filter.slider(iCtrl).slider.Enable = enable;
        end
        redraw(true);
    end

    function redraw(varargin)
        if ~ishandle(hs.line(1)) && length(varargin)>1
            % figure overplotted by normal plot()
            return
        end
        
        scrollValue = get(hs.scroll,'Value');
        %fprintf('scrollValue: %f\n',scrollValue)
        zoomValue   = get(hs.zoom,'Value');
        
        % zoomValue==0: numPoints=N
        % zoomValue==1: numPoints=7
        
        % N * exp(d*x)
        numPoints = N*exp(d*zoomValue);
        numPoints = round(numPoints);
        numPoints = max(numPoints,2);
        
        % scrollValue==0: startIndex=1;
        % scrollValue==1: startIndex=N-numPoints+1;
        startIndex = (N-numPoints)*scrollValue+1; % m*x+b
        endIndex = startIndex+numPoints;
        
        if ~isscalar(SIG.X0)
            startTime = timeBoundary(1)+period*(startIndex-1);
            [~,startIndex] = min(abs(startTime-SIG.X0));
            endTime = timeBoundary(1)+period*(endIndex-1);
            [~,endIndex] = min(abs(endTime-SIG.X0));
        end
        
        startIndex = round(startIndex);
        endIndex = round(endIndex);
        startIndex = max(startIndex,1);
        endIndex = min(endIndex,N);
        
        % Maximum factor_max values per pixel, 
        % so very long signals don't hang Matlab
        factor = round(numPoints/max(1,axWidth_px));
        factor_max = 1000; % increase this if you want to find the "needle in the haystack" (single outlier sample) 
        if factor>factor_max
            spc = floor(factor/factor_max);
        else
            spc = 1;
        end
        ind = startIndex:spc:endIndex;
        %fprintf('spc: %f\n',spc)
        
        % Use original or filtered signal according to checkbox state
        if isfield(SIG,'filter') && SIG.filter.mainCheck.Value == SIG.filter.mainCheck.Max
            if isscalar(SIG.X), XData=1/SIG.X*(ind-1); else, XData=SIG.X(ind); end
            YData = SIG.Y(ind,:);
        else
            if isscalar(SIG.X0), XData=period*(ind-1); else, XData = SIG.X0(ind); end
            YData = SIG.Y0(ind,:);
        end
        
        if isscalar(SIG.X0)
            startTime = XData(1);
            endTime = XData(end);
        end
        
        % Don't show much more samples than pixels. 
        % Make sure that minimum and maximum data is shown anyway
        % (except if factor is > factor_max, for responsiveness)
        maxSamplesPerPixel = 2;
        if size(YData,1)/(axWidth_px*maxSamplesPerPixel) > 2
            factor = ceil(size(YData,1)/(max(1,axWidth_px)*maxSamplesPerPixel));
            remove = mod(size(YData,1),factor);
            XData = XData(1:(end-remove));
            YData = YData(1:(end-remove),:);
            XData = reshape(XData,factor,[]);
            XData = [min(XData,[],1);max(XData,[],1)];
            XData = XData(:)';
            YData = permute(YData,[3,1,2]);
            YData = reshape(YData,factor,[],size(YData,3));
            YData = [min(YData,[],1,'includenan');max(YData,[],1,'includenan')]; % preserves 'NaN' separations
            YData = [min(YData,[],1');max(YData,[],1)];
            YData = reshape(YData,[],size(YData,3));
        end
        
        % On the other hand, if there are much less samples than pixels, 
        % show additional dots
        factor = size(YData,1)/max(1,axWidth_px);
        if ~strcmp(parser.Results.AutoMarkers,'none')
            if factor<0.2
                set(hs.line, 'Marker',parser.Results.AutoMarkers);
            else
                set(hs.line, 'Marker','none');
            end
        end
        
        if ~isempty(SIG.AutoStackSignals) && ~strcmp(parser.Results.YLimMode,'fixed')
            % Stack signals horizontally dynamically
            [sigPosVec,sigAddVec] = auto_stack(YData);
            hs.ax.YTick = flip(sigPosVec);
            hs.ax.YTickLabel = flip(SIG.AutoStackSignals(:));
            hs.ax.TickLabelInterpreter = 'none';
        end
        YData = bsxfun(@plus,YData,sigAddVec);

        
        % hs.ax2 limits
        if ~strcmp(parser.Results.SecondXAxisFunction,'none')
            ax2Limits = [feval(parser.Results.SecondXAxisFunction,startTime), feval(parser.Results.SecondXAxisFunction,endTime)];
            if ax2Limits(1)<ax2Limits(2)
                set(hs.ax2,'XDir','normal')
            else
                ax2Limits = flip(ax2Limits);
                set(hs.ax2,'XDir','reverse')
            end
            set(hs.ax2,'XLim',ax2Limits)
        end

        set(hs.line,'XData',XData);
        for iLine = 1:size(YData,2)
            set(hs.line(iLine),'YData',YData(:,iLine));
        end
        set(hs.ax,'XLim',[startTime,endTime])
        minY = min(YData(:));
        maxY = max(YData(:));
        delta = (maxY-minY)/50;
        minY = minY-delta;
        maxY = maxY+delta;
        if nnz(sigPosVec)>1
            minY = min([minY;sigPosVec(:)]);
            maxY = max([maxY;sigPosVec(:)]);
        end
        
        if strcmp(parser.Results.YLimMode,'dynamic') && nnz(isnan([minY,maxY]))==0 && maxY>minY
            set(hs.ax2,'YLim',[minY,maxY])
            set(hs.ax,'YLim',[minY,maxY])
        end
        
        % Big scrollbar if there is nothing to scroll
        if N<=numPoints
            majorStep = Inf; 
            minorStep = .1;
        else % N > numPoints
            majorStep = max(1e-6,numPoints/(N-numPoints)); 
            % 100 steps per screen width
            minorStep = max(1e-6,(endTime-startTime)/(100*duration));
        end
        set(hs.scroll,'SliderStep',[minorStep,majorStep]);      
        
        if ~strcmp(parser.Results.ScrollCallback,'none')
            arg.XLim = [startTime;endTime];
            arg.ContinuousValueChange = ~isempty(varargin) && ~(islogical(varargin{end}) && varargin{end});
            arg.hs = hs;
            feval(parser.Results.ScrollCallback,arg);
        end
    end



    function resizegui(varargin)
        
        panelUnits = hs.panel.Units;
        
        % Centimeter layout
        set(hs.panel ,'Units',units);
        set(hs.ax    ,'Units',units);
        set(hs.ax2   ,'Units',units);
        set(hs.list  ,'Units',units);
        set(hs.scroll,'Units',units);
        set(hs.zoom  ,'Units',units);
        
        width = hs.panel.Position(3);
        height = hs.panel.Position(4);
        
        yPos = space;
        
        % Filter layout
        if isfield(SIG,'filter')
            filterHeight = layoutFilter(iSig,1,width-2*space);
            pos = [space,yPos,width-2*space,filterHeight];
            pos = [pos(1), pos(2), max(0,pos(3)), max(0,pos(4))];
            set(SIG.filter.panel, 'Units',units, 'Position',pos);
            yPos = pos(2)+pos(4)+4*space;
        end
        
        % Zoom slider
        pos = [space,yPos,width-2*space,sliderHeight];
        pos = [pos(1), pos(2), max(0,pos(3)), max(0,pos(4))];
        set(hs.zoom, 'Units',units, 'Position',pos)
        yPos = pos(2)+pos(4)+space;
        
        % Scroll slider
        pos = [space,yPos,width-2*space,sliderHeight];
        pos = [pos(1), pos(2), max(0,pos(3)), max(0,pos(4))];
        set(hs.scroll, 'Units',units, 'Position',pos)
        yPos = pos(2)+pos(4)+space;
        
        % List (I)
        if strcmpi(parser.Results.ShowInformationList,'none')
            listWidth = 0;
        else
            listWidth = 3;
            listWidth = min(listWidth,width/5);
        end
        
        % Axis
        if strcmpi(parser.Results.ShowAxisTicks,'on')
            insets = get(hs.ax,'TightInset');
            if isequal(hs.ax2.Visible,'on')
                insets = insets + get(hs.ax2,'TightInset');
            end
        else
            insets = [0,0,0,0];
        end
        pos = [space,yPos,max(1,width-3*space-listWidth),max(1,height-yPos)];
        pos = [pos(1)+insets(1), pos(2)+insets(2), pos(3)-insets(1)-insets(3), pos(4)-insets(2)-insets(4)];
        pos = [pos(1), pos(2), max(0,pos(3)), max(0,pos(4))];
        set(hs.ax, 'Units',units, 'Position',pos)
        set(hs.ax2, 'Units',units, 'Position',pos)
        set(hs.panel,'Units',panelUnits);
        
        % List (II)
        if ~strcmpi(parser.Results.ShowInformationList,'none')
            axPos = hs.ax.Position;
            pos = [width-space-listWidth,axPos(2),listWidth,axPos(4)];
            set(hs.list, 'Units',units, 'Position',pos)
        end
        
        % Update axWidth_cm and axWidth_px
        set(hs.ax,'Units','centimeters');
        axWidth_cm=get(hs.ax,'Position'); axWidth_cm=axWidth_cm(3); axWidth_cm=max(axWidth_cm,0);
        set(hs.ax,'Units','pixels');
        axWidth_px=get(hs.ax,'Position'); axWidth_px=round(axWidth_px(3)); axWidth_px=max(axWidth_px,0);
        
        % Change zooming such that mmPerSec stays the same
        if ~isempty(varargin) && ~isempty(mmPerSec) % do this only for calls by GUI, not during the initialization call
            numPoints = axWidth_cm*10/(mmPerSec_slider*period);
            zoomValue = log(numPoints/N)/d;
            zoomValue = max(zoomValue,0);
            zoomValue = min(zoomValue,1);
            set(hs.zoom,'Value',zoomValue)
            redraw(true)
        end
    end

    function height = layoutFilter(~,~,width)
        yPos = space;
        % Layout the sliders and their text fields from bottom to top
        nCtrl = length(SIG.filter.slider);
        % Find maximum label width
        labelWidth = 1;
        for iCtrl = 1:nCtrl
            ext = SIG.filter.slider(iCtrl).text.Extent(3);
            if ext > labelWidth
                labelWidth = ext;
            end
        end
        for iCtrl = flip(1:nCtrl)
            xPos = space;
            % Text label
            pos = [xPos,yPos,labelWidth,sliderHeight];
            pos = [pos(1), pos(2), max(0,pos(3)), max(0,pos(4))];
            set(SIG.filter.slider(iCtrl).text, 'Units',units, 'Position',pos);
            xPos=pos(1)+pos(3)+space;
            % Edit box
            pos = [xPos,yPos,editWidth,sliderHeight];
            pos = [pos(1), pos(2), max(0,pos(3)), max(0,pos(4))];
            set(SIG.filter.slider(iCtrl).edit, 'Units',units, 'Position',pos);
            xPos=pos(1)+pos(3)+space;
            % Slider
            pos = [xPos,yPos,width-xPos-space,sliderHeight];
            pos = [pos(1), pos(2), max(0,pos(3)), max(0,pos(4))];
            set(SIG.filter.slider(iCtrl).slider, 'Units',units, 'Position',pos);
            % New line
            yPos=pos(2)+pos(4)+space;
        end
        % Layout the checkbox
        yPos = yPos-3*space;
        pos = [space,yPos,width-2*space,checkboxHeight];
        pos = [pos(1), pos(2), max(0,pos(3)), max(0,pos(4))];
        set(SIG.filter.mainCheck, 'Units',units, 'Position',pos);
        yPos = pos(2)+pos(4)+space;
        height = yPos+3*space;
    end


    function figScroll(~,callbackdata)
        scrollCount = callbackdata.VerticalScrollCount;
        val = hs.scroll.Value + scrollCount*hs.scroll.SliderStep(2);
        val = max(hs.scroll.Min,val);
        val = min(hs.scroll.Max,val);
        hs.scroll.Value = val;
        redraw(true);
    end

    function btnDown(varargin)
        % x and y location where button was pressed is stored
        % ShowInformationList function is called with these locations
        % and the returned strings  are displayed. 
        % Numbers are converted to strings. 
        x = hs.ax.CurrentPoint(1,1);
        y = hs.ax.CurrentPoint(1,2);
        click_info.x = [x;click_info.x];
        click_info.y = [y;click_info.y];
        
        %std_InformationList(hs,click_info)
        
        if strcmpi(parser.Results.ShowInformationList,'none')
            str = {};
        else
            str = feval(parser.Results.ShowInformationList, hs,click_info);
        end
        assert(iscell(str),'ShowInformationList function must return a cell array.')
        assert(isempty(str) || isvector(str), 'ShowInformationList function must return a cell vector.')
        
        % Convert numbers to strings
        for k = 1:length(str)
            if ~ischar(str{k})
                str{k} = num2str(str{k});
            end
        end
        
        % Update uicontrol
        hs.list.String = str; 
    end

if ~strcmpi(parser.Results.ShowInformationList,'none')
    hs.ax.ButtonDownFcn = @btnDown;
    %set(hs.line,'PickableParts','none')  % DataCursor doesn't work anymore if you use this. 
    set(hs.line,'ButtonDownFcn',@btnDown) % DataCursor mode puts its own callback here as long as it is enabled. 
end

    function delete_plotECG(varargin)
        %fprintf('delete_plotECG\n')
        SIG = [];
        click_info = [];
        if isvalid(hs.fig)
            hs.fig.WindowScrollWheelFcn = '';
        end
    end
    
hs.ax.DeleteFcn    = @delete_plotECG;
hs.panel.DeleteFcn = @delete_plotECG;
hs.zoom.DeleteFcn  = @delete_plotECG;




%% Return output arguments

if nargout>=1
    varargout{1} = hs.line;
end
if nargout>=2
    varargout{2} = hs;
end


end











%% Helper Functions

function y = editToSlider(edit,slider)
y = str2double(get(edit,'String'));
if y<=0
    % negative value provided: use minimum possible value
    y = slider.Min;
    edit.String = num2str(y,4);
end
% Logarithmic Conversion
a = slider.Min;
b = slider.Max;
r = log(a/b)/(a-b);
p = a*exp(-log(a/b)*a/(a-b));
if isnan(y)
    % no valid number string typed: restore old value 
    y = p*exp(r*slider.Value);
    edit.String = num2str(y,4);
else
    % convert to logarithmic scale
    x = log(y/p)/r;
    if x<slider.Min
        slider.Value = slider.Min;
    elseif x>slider.Max
        slider.Value = slider.Max;
    else
        slider.Value = x;
    end
end
drawnow
end


function y = sliderToEdit(edit,slider)
x = slider.Value;
% convert to exponential scale
% y = p*exp(r*x), x:[a,b], y:[a,b]
a = slider.Min; 
b = slider.Max;
r = log(a/b)/(a-b);
p = a*exp(-log(a/b)*a/(a-b));
y = p*exp(r*x);
edit.String = num2str(y,4); 
%edit.String = sprintf('%0.3g',y); 
end


function ls = isLineSpec(str)
ls = ischar(str) && length(str)<=4;
allowed = '-:.+o*xsd^v><phrgbcmykw';
for pos = 1:length(str)
    ls = ls && any(str(pos)==allowed);
end
end


function str = func2str2(func)
if ischar(func)
    str = func;
else
    str = func2str(func);
end
end


function str = function_file(func)
if ischar(func)
    funH = str2func(func);
    funS = func;
else
    funH = func;
    funS = func2str(func);
end
S = functions(funH);
str = S.file;
if isempty(str)
    str = funS;
else 
    str = sprintf('%s()   %s',funS,str);
end
end


function [sigPosVec,sigAddVec] = auto_stack(YData)
% Stacks Signals Horizontally with little overlap
% Used after each scroll/zoom action for 'AutoStackSignals' 
% with 'YLimMode' set to 'dynamic'.
% You might want to adjust this for your specific needs. 

%YData = bsxfun(@minus,YData,YData(1,:));
signalMed = median(YData,1,'omitnan');
YData = bsxfun(@minus,YData,signalMed);
overlap = YData;
overlap = diff(overlap,1,2); % positive values are overlap
overlap(isnan(overlap)) = 0;
overlapS = sort(overlap,1,'descend');
index = max(1,round(size(overlapS,1)*.007));
signalSpacing = overlapS(index,:)*1.1;
stdd = std(YData,1,1);
stdd = min(stdd(1:end-1),stdd(2:end));
signalSpacing = max(signalSpacing, median(signalSpacing,'omitnan')*.5 );
signalSpacing = max(signalSpacing, stdd*4);
% Increas very small spacings
signalSpacing = max(signalSpacing, max(signalSpacing)./1000.*ones(size(signalSpacing)));
signalSpacing = max(eps,signalSpacing);
sigPosVec = -cumsum([0 signalSpacing]);
sigAddVec = sigPosVec-signalMed;
end


function [sigPosVec,sigAddVec] = auto_stack_nooverlap(YData)
% Stacks Signals Horizontally with strictly no overlap. 
% Used for 'AutoStackSignals' with 'YLimMode' set to 'fixed'.
% You might want to adjust this for your specific needs. 

signalMed = median(YData,1,'omitnan');
YData = bsxfun(@minus,YData,signalMed);
overlap = YData;
overlap = min(overlap(:,1:end-1),[],1) - max(overlap(:,2:end),[],1);
overlap(isnan(overlap)) = 0;
signalSpacing = -overlap*1.01;
% Increas very small spacings
signalSpacing = max(signalSpacing, max(signalSpacing)./1000.*ones(size(signalSpacing)));
signalSpacing = max(eps,signalSpacing);
sigPosVec = -cumsum([0 signalSpacing]);
sigAddVec = sigPosVec-signalMed;
end





%% Built-in Functions for 'ShowInformationList' 


function str = ecg_InformationList(hs,click)
% click.x: [n x 1] array, first value is last click
% click.y: [n x 1] array, first value is last click

if size(click.x,1)<2
    click.x = [click.x;0;0];
    click.y = [click.y;0;0];
end

deltaX = abs(click.x(1)-click.x(2));
freq = 1/deltaX;

str = {
    
    'BPM'
    freq*60
    ''
    
    'RR'
    deltaX
    ''
    
    };
end


function str = std_InformationList(hs,click)
% click.x: [n x 1] array, first value is last click
% click.y: [n x 1] array, first value is last click

if size(click.x,1)<2
    click.x = [click.x;0;0];
    click.y = [click.y;0;0];
end

deltaX = abs(click.x(1)-click.x(2));
freq = 1/deltaX;
deltaY = abs(click.y(1)-click.y(2));
%fprintf('dx=%f; freq=%f, dy=%f\n',deltaX,freq,deltaY)

str = {
    '<html>( &#916x <tt>&#8629</tt> 1/&#916x )</html>'
    deltaX
    freq
    ''
    
    '<html>( &#916y )</html>'
    deltaY
    ''
    
    '<html>( x<sub>1</sub> <tt>&#8629</tt> x<sub>2</sub> )</html>'
    click.x(2)
    click.x(1)
    ''
    
    '<html>( y<sub>1</sub> <tt>&#8629</tt> y<sub>2</sub> )</html>'
    click.y(2)
    click.y(1)
    ''
    
    };
end





%% Built-in Filter Functions

function [ time_out, sig_out, description ] = filter_FFT( time_in, sig_in, f1, f2 )
    % --- Description struct ---
    % Find samplerate
    if isscalar(time_in), samplerate=time_in; else samplerate=(length(time_in)-1)/(time_in(end)-time_in(1)); end
    % Filter title
    description.string = 'FFT Bandpass';
    k = 0;
    % Slider for lower cutoff frequency
    k = k+1;
    description.slider(k).Label  = 'f1';
    description.slider(k).Min    = 0.001;
    description.slider(k).Value  = 0.2;
    description.slider(k).Max    = 4;
    % Slider for higher cutoff frequency
    k = k+1;
    description.slider(k).Label  = 'f2';
    description.slider(k).Min    = min(  5,samplerate*0.01);
    description.slider(k).Value  = min( 49,samplerate/2.01);
    description.slider(k).Max    = min(200,samplerate/2.01);
    % If an empty signal [] is given, just return the description struct
    if isempty(sig_in); time_out = []; sig_out = []; return; end

    % --- Signal processing ---
    sig = sig_in;
    L = size(sig,1);
    % N = 2^nextpow2(L); sometimes slower
    N = L;
    SIG = fft(sig,N);
    freq = 0:samplerate/N:samplerate*(N-1)/N;
    freq(freq>samplerate/2) = freq(freq>samplerate/2) - samplerate;
    SIG(abs(freq)<f1 | abs(freq)>f2,:) = 0;
    sig = ifft(SIG,'symmetric');
    sig = sig(1:L,:);
    sig_out = sig;
    time_out = time_in;
end






function [ time_out, sig_out, description ] = filter_bandpass_notch( time_in, sig_in, f_high, f_low, notch_width, f_notch )
    % FILTER_BANDPASS_NOTCH - ECG Filter
    %tic

    % --- Define Filter Title and Labels and Ranges for the Parameters --- 

    % Find samplerate
    if isscalar(time_in), samplerate=time_in; else samplerate=1/(time_in(2)-time_in(1)); end

    % Filter title
    description.string     = 'Zero phase IIR bandpass & notch';

    % Slider for highpass cutoff frequency
    description.slider(1).Label  = 'highpass';
    description.slider(1).Min    = 0.00001;
    description.slider(1).Value  = 0.3; 
    description.slider(1).Max    = 2; 

    % Slider for lowpass cutoff frequency
    description.slider(2).Label  = 'lowpass';
    description.slider(2).Min    = min(  3,samplerate*0.01);
    description.slider(2).Value  = min( 60,samplerate/2.01); 
    description.slider(2).Max    = min(200,samplerate/2.01); 

    % Slider for notch width
    description.slider(3).Label  = 'notch width';
    description.slider(3).Min    =  0.000001;
    description.slider(3).Value  =  .3;
    description.slider(3).Max    =  99; 

    % Slider for notch frequency
    description.slider(4).Label  = 'notch';
    description.slider(4).Min    =  min(40,samplerate/4);
    description.slider(4).Value  =  min(50,samplerate/3);
    description.slider(4).Max    =  min(70,samplerate/2.01);

    % If an empty signal [] is given, just return the description struct
    if size(sig_in,1)<=12; time_out=time_in; sig_out=NaN*sig_in; return; end


    % --- Signal Processing --- 
    sig = sig_in;

    % Removing deviation from zero not necessary as filtfilt already minimizes 
    % start-up and ending transients by matching initial conditions
    %sig = detrend(sig);
    %sig = bsxfun(@minus,sig,mean(sig));
    %sig = bsxfun(@minus,sig,sig(1,:));

    % Bandpass filter
    bpFilt = designfilt('bandpassiir','FilterOrder',4,'HalfPowerFrequency1',f_high,'HalfPowerFrequency2',f_low,'SampleRate',samplerate);
    %fvtool(bpFilt)
    sig = filtfilt(bpFilt,sig);
    %sos = bpFilt.Coefficients; % Works only for IIR
    %[z,p,k]=bpFilt.zpk; [sos,g]=zp2sos(z,p,k);

    % Bandstop filter
    if notch_width>0
        bsFilt = designfilt('bandstopiir','FilterOrder',2,'HalfPowerFrequency1',f_notch-notch_width/2,'HalfPowerFrequency2',f_notch+notch_width/2,'SampleRate',samplerate);
        %fvtool(bsFilt)
        sig = filtfilt(bsFilt,sig);
        %sos = [sos;bsFilt.Coefficients];
        %[z,p,k]=bsFilt.zpk; [sos2,g2]=zp2sos(z,p,k); sos=[sos;sos2]; g=g*g2;
    end

    %sig = filtfilt(sos,g,sig); % no speed improvement with concatenating filters first

    sig_out = sig;
    time_out = time_in;

    %toc
end








%% Show Demo


function show_demo()

    persistent demotype
    if isempty(demotype)
        demotype = 0;
    else
        demotype = mod(demotype+1,5);
    end
    
    switch demotype
        
        case 0
            % minimalistic view
            fs = 200;
            X = (0:1/fs:10000)';
            f = .5;
            Y = [sin(2*pi*f*X)+rand(size(X))*.1, 0.5*cos(2*pi*f*X)+rand(size(X))*.3];
            plotECG(X,Y, 'AutoStackSignals',{'Sinus','Cosinus'}, 'ShowAxisTicks','off', 'YLimMode','fixed');
        
        otherwise
            % real ECG with random electrical heart axis
            dt = 0.01;
            Mx = [-0.0408 -0.0395 -0.0384 -0.0364 -0.0306 -0.0213 -0.0126 -0.00819 -0.0108 -0.0179 -0.0264 -0.0313 -0.0382 -0.0478 -0.0485 -0.0465 -0.0464 -0.0171 0.0289 0.206 0.311 -0.167 -0.249 -0.0988 -0.0255 -0.00753 -0.0147 -0.00495 0.00099 -0.000707 0.00394 0.0127 0.0228 0.0288 0.0321 0.0389 0.0497 0.0661 0.0876 0.118 0.159 0.194 0.221 0.251 0.251 0.21 0.153 0.0951 0.049 0.0125 -0.0135 -0.0276 -0.0353 -0.0401 -0.0408 -0.0406 -0.0422 -0.0438 -0.0444 -0.0429 -0.0405 -0.0381 -0.0364 -0.0341 -0.0299 -0.035 -0.0428 -0.0338 -0.0281 -0.0323 -0.0346 -0.0358 -0.0367 -0.0375 -0.0384 -0.0395 -0.0408 -0.0419 -0.0426 -0.043 -0.0431 -0.0431 -0.0429 -0.0426 -0.0422 -0.0419 -0.0415 -0.041 -0.0404 -0.0399 -0.0389 -0.0368 -0.0345];
            My = [0.071 0.0696 0.0671 0.0574 0.0567 0.0651 0.0427 0.0222 0.036 0.0528 0.0681 0.0851 0.0994 0.109 0.102 0.0963 0.0909 0.131 0.0352 -0.14 -0.63 -0.524 0.282 0.398 0.185 0.0532 0.0342 0.0199 0.0163 0.00534 0.00251 -0.00119 -0.0131 -0.0237 -0.0333 -0.0469 -0.0609 -0.0757 -0.096 -0.126 -0.169 -0.218 -0.265 -0.298 -0.327 -0.341 -0.299 -0.219 -0.125 -0.0572 -0.0185 0.0113 0.0299 0.0398 0.0418 0.041 0.0418 0.0419 0.041 0.0387 0.0362 0.0347 0.0344 0.0359 0.0394 0.0376 0.0344 0.0415 0.0476 0.0487 0.0513 0.0546 0.0578 0.0609 0.0639 0.0667 0.0691 0.0711 0.0725 0.0734 0.0739 0.0741 0.0741 0.0742 0.0751 0.0769 0.0788 0.0796 0.0786 0.0769 0.0756 0.0761 0.0765];
            P  = [-0.089 -0.107 -0.12 -0.129 -0.136 -0.14 -0.141 -0.141 -0.14 -0.138 -0.134 -0.129 -0.123 -0.116 -0.11 -0.104 -0.0977 -0.0916 -0.0857 -0.0799 -0.0741 -0.0684 -0.0632 -0.0586 -0.0549 -0.0522 -0.0508 -0.0506 -0.0517 -0.0538 -0.0571 -0.0611 -0.0654 -0.0693 -0.0728 -0.0764 -0.079 -0.0809 -0.0836 -0.09 -0.102 -0.113 -0.119 -0.115 -0.102 -0.081 -0.0542 -0.0186 0.0287 0.0943 0.173 0.252 0.322 0.378 0.418 0.438 0.434 0.413 0.377 0.331 0.277 0.22 0.166 0.114 0.0634 0.0131 -0.0354 -0.0783 -0.114 -0.141 -0.159 -0.165 -0.159 -0.143 -0.12 -0.0927 -0.062 -0.0291 0.00392 0.0356 0.0645 0.0909 0.113 0.126 0.128 0.117 0.0976 0.0728 0.045 0.0154 -0.0132 -0.0377 -0.0591];
            MXP = [Mx',My',P']; % [93 x 3]
            L = size(MXP,1);
            
            % Interpolation 
            dt_interp = 0.001;
            xq = 0:dt_interp:L*dt;
            % add some values for periodic continuation
            MXP_ex = [MXP;MXP(1:10,:)]; % [103 x 3]
            L_ex = size(MXP_ex,1);
            t = 0:dt:(L_ex-1)*dt;
            ind_interp = true(size(t));
            ind_interp(L-2:L) = false; % leave out for better periodic continuation
            MP_interp = interp1(t(ind_interp), MXP_ex(ind_interp,:), xq, 'spline');
            L_interp = size(MP_interp,1);
            
            % Lead Field Matrix for {I,II,III,aVR,aVlL,aVF} 
            A = [1,0; 1/2,-sqrt(3)/2; -1/2,-sqrt(3)/2; -sqrt(3)/2,1/2; sqrt(3)/2,1/2; 0,-1];
            % Rotation Matrix for changing electrical heart axis
            R = @(alpha) [cos(alpha),-sin(alpha); sin(alpha),cos(alpha)];
            alpha = 60*randn()+45;
            % Calculate ECG Leads: Lead Field Matrix & Rotation
            ECG = A * R(alpha/180*pi) * MP_interp(:,1:2)';
            ECG = ECG';
            ECGP = [ECG,MP_interp(:,3)];
            
            % Repeat single heartbeat multiple times
            N = ceil(3*60/(L_interp*dt_interp));
            ECGP_rep = repmat(ECGP,N,1);
            t_rep = (0:dt_interp:(L_interp*N-1)*dt_interp)';
            
            % Add 50Hz Interfering Signal
            noise = 0.05*sin(2*pi*50*t_rep);
            ECGP_rep = bsxfun(@plus,ECGP_rep,noise);
            
            % Use plotECG
            fig = figure(...
                'Units','normalized', ...
                'OuterPosition',[0.05,0.1,.9,.88],...
                'HandleVisibility','Callback');
            
            [h,hs] = plotECG(t_rep,ECGP_rep, ...
                'Parent',fig, ...
                'AutoStackSignals',{'I','II','III','aVR','aVL','aVF','Pleth'}, ...
                'YLimMode','fixed',...
                'SecondXAxisFunction',@(x) x/60, ...
                'SecondXAxisLabel', 'Time in min', ...
                'Filter', 'filter_FFT');
            
            set(h,'LineWidth',1.5)
            set(hs.fig, 'Name','Demo plotECG', 'NumberTitle','off')
            hs.ax.FontWeight='bold';
            hs.ax.YLabel.String = 'ECG Channel';
            
    end
    
end






