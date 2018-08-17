function varargout = sineSweepGui_sub(varargin)
% SINESWEEPGUI_SUB MATLAB code for sineSweepGui_sub.fig
%      SINESWEEPGUI_SUB, by itself, creates a new SINESWEEPGUI_SUB or raises the existing
%      singleton*.
%
%      H = SINESWEEPGUI_SUB returns the handle to a new SINESWEEPGUI_SUB or the handle to
%      the existing singleton*.
%
%      SINESWEEPGUI_SUB('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SINESWEEPGUI_SUB.M with the given input arguments.
%
%      SINESWEEPGUI_SUB('Property','Value',...) creates a new SINESWEEPGUI_SUB or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sineSweepGui_sub_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sineSweepGui_sub_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sineSweepGui_sub

% Last Modified by GUIDE v2.5 24-Dec-2017 14:53:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sineSweepGui_sub_OpeningFcn, ...
                   'gui_OutputFcn',  @sineSweepGui_sub_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before sineSweepGui_sub is made visible.
function sineSweepGui_sub_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sineSweepGui_sub (see VARARGIN)

% Choose default command line output for sineSweepGui_sub
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes sineSweepGui_sub wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = sineSweepGui_sub_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
