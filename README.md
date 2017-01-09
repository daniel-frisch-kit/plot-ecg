# plot-ecg
plotECG( X,Y ) — Plot Very Long, Multichannel Signals — Zoom &amp; Scroll/Pan via Slider

-------------------------------------------------------------------------------------------- 
   
Enables you to plot and zoom & scroll/pan through signals with millions of samples. 
There is one slider for panning and one for zooming. 
To move forward or backward by exactly one screen width, 
use the scroll wheel or click on the slider trough. 
  
If you launch plotECG() without any arguments, 
one of two demos is shown. 
  
-------------------------------------------------------------------------------------------- 
  
Input Arguments 
  
- X: 
    Vector of timesteps, or scalar sample rate. scalar or [N x 1] double 
      
- Y: 
    Signal vector, or signal matrix where each column represents one signal. [N x m] double 
       
- LineSpec: 
    LineSpec string, max. 4 characters, see Matlab plot() documentation. 
    string (optional, default: '.') 
    If you define a marker here, set 'AutoMarkers' to 'none' so it won't be overridden. 
  
  
Important Key-Value Parameters: 
(For more parameters, look into plotECG.) 
  
- 'Filter' 
    Change parameters like cutoff frequencies via slider, or text edit, and see the results on the fly. 
    Use the builtin 'filter_FFT', or 'filter_bandpass_notch', or define an own 
    filter function and pass its name as string or function_handle (see plotECG help). 
     
- 'mmPerSec' 
    Initial zoom in screen millimeters per second. scalar double (default: 50) 
    If you change the figure size in this mode, the axis XLim will change 
    such that mmPerSec stays the same. 
     
- 'secPerScreenWidth' 
    Initial zoom setting in seconds that are displayed on screen at a time. scalar double 
    If you change the figure size in this mode, the signal scale will change too, 
    such that XLim stays the same. 
       
- 'ShowAxisTicks' 
    Shows the axis ticks and labels and a grid. string, 'on' or 'off' (default: 'on') 
    You can also change axis properties on the returned handles: hs.ax.XLabel.String = 'Seconds'; 
     
- 'ShowInformationList' 
    Shows information about the location of the last clicks in the signal. function_handle or string 
    Specify one of 'none' (default), 'std_InformationList', 'ecg_InformationList', as string, 
    or define an own function like those with two inputs and cell output 
    and pass its name or function_handle. 
    Mouse clicks are captured only if no interactive mode (pan, zoom etc.) is active in the figure. 
     
- 'AutoStackSignals' 
    Cell array of strings with signal names. Length must be equal to number of columns of Y, or 0. 
    (default: {}) Stacks the signals vertically, so for example a multipolar ECG can be shown. 
    Example: {'I','II','III', 'aVR','aVL','aVF', 'V1','V2','V3','V4','V5','V6'} 
     
- 'SecondXAxisFunction' 
    Function handle that maps from the provided X values to a different x axis scale 
    that will be displayed above the plotted signal. function_handle or string (default: 'none') 
    Example: @(x)x/60^2, shows the time also in hours. 
     
- 'SecondXAxisLabel' 
    Label of second x axis. Example: 'Time in h'. string (default: '') 
     
- 'YLimMode' 
    'dynamic': dynamic y axis limits according to minimum and maximum of currently visible signal 
    'fixed' : fixed y axis limits according to minimum and maximum of entire signal. 
    (default: 'dynamic') In 'fixed' mode, you shouldn't apply filters that change 
    the signal's mean much because YLim won't be updated. 
    You can change the fixed interval afterwards using the returned axes handle: 
    hs.ax.YLim = [-10,10]; 
     
- 'Parent' 
    Parent figure. (default: new figure is created) 
    Use delete(findall(hs.panel)) or close the figure to delete the old plot and its signals. 
    The 'HandleVisibility' of the figure created by default is set to 'Callback' to prevent you from 
    plotting into it accidentially from the command line. 
    To close the created figures from the command line, you have to use "close all hidden" 
    instead of "close all". 
    
-------------------------------------------------------------------------------------------- 
    
Output Arguments 
- h: 
    Returns the chart line objects as a vector. Use h to modify a chart line after it is created. 
     
- hs: 
    Returns a struct with handles to some more GUI objects for later modifiaction 
  
  
You can set Name-Value parameters on the returned chart line handles: 
set(h, 'LineWidth',3) etc. 
  
Furthermore, you might want to change the X/YLabel: 
hs.ax.XLabel.String = 'Seconds'; or: xlabel(hs.ax,'Seconds') 
   
-------------------------------------------------------------------------------------------- 
  
Example

X = 0:0.001:100-0.001;

Y = sin(2*pi*0.1*X) + sin(2*pi*X) + .1*sin(2*pi*50*X);

[h,hs] = plotECG(X,Y);

[h,hs] = plotECG(1000,Y, 'Filter','filter_FFT'); 
   
--------------------------------------------------------------------------------------------
