
function varargout = Interface_Tandoori_TRG(varargin)

% INTERFACE_TANDOORI_TRG MATLAB code for Interface_Tandoori_TRG.fig
%      INTERFACE_TANDOORI_TRG, by itself, creates a new INTERFACE_TANDOORI_TRG or raises the existing
%      singleton*.
%
%      H = INTERFACE_TANDOORI_TRG returns the handle to a new INTERFACE_TANDOORI_TRG or the handle to
%      the existing singleton*.
%
%      INTERFACE_TANDOORI_TRG('CALLBACK',hObject,eventData,handles,...) calls the local

%      function named CALLBACK in INTERFACE_TANDOORI_TRG.M with the given input arguments.
%
%      INTERFACE_TANDOORI_TRG('Property','Value',...) creates a new INTERFACE_TANDOORI_TRG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Interface_Tandoori_TRG_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Interface_Tandoori_TRG_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help Interface_Tandoori_TRG
% Last Modified by GUIDE v2.5 29-Sep-2017 14:25:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Interface_Tandoori_TRG_OpeningFcn, ...
                   'gui_OutputFcn',  @Interface_Tandoori_TRG_OutputFcn, ...
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

% --- Executes just before Interface_Tandoori_TRG is made visible.
function Interface_Tandoori_TRG_OpeningFcn(hObject, eventdata, handles, varargin)
clc
%% ############################################################################
%% ####################### INITIALISATION DES VARIABLES #######################
%% ############################################################################
Fit696=0; Fit706=0; Fit727=0; Fit738=0; Fit750=0; Fit751=0; Fit763=0; Fit794=0; Fit800=0;
Fit801=0; Fit811=0; Fit826=0; Fit840=0; Fit842=0; Fit852=0; Fit866=0; Fit810=0; Fit810Kr=0;
Fit758=0; Fit760=0; Fit768=0; Fit769=0; Fit785=0; Fit819=0; Fit811Kr=0; Fit829=0; Fit850=0; Fit877=0; 
Fit788=0; Fit823=0; Fit828=0; Fit881=0; Fit904=0; Fit820=0; Fit834=0; Fit640=0; 
Fit585=0; Fit667=0; Fit714=0; Fit892=0; Fit895=0;   


filetype=0; GraphExp=0; TeGraph1=0; TeGraph2=0; TeGrapheFinaux=0; graphePopDepop=0;
P=0.01; Tg=300; 
TeMin=0.1; TeMax=2; TeStep=0.005;
NeMin=10; NeMax=16; NePts=200;

% Choose default command line output for Interface_Tandoori_TRG
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Interface_Tandoori_TRG wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = Interface_Tandoori_TRG_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% ################################################################################
%% ########### DÉBUT DES INPUTS VIA L'INTERFACE_TANDOORI_TRG ###########
%% ################################################################################

%% =============== Selection du répertoire où se trouvent les spectres ===============
function Repertoire_Callback(hObject, eventdata, handles)
    Repertoire=uigetdir('/Users/hackintosh/Documents/Physique/Maitrise_2016/Experimental Data');
    set(handles.RepertoireSelect, 'String', Repertoire);


%% =============== Sélection de l'extension des fichiers ===============
function ExtensionFichier_Callback(Extension, eventdata, handles)
function ExtensionFichier_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end    

function ExtensionFichier_ButtonDownFcn(hObject, eventdata, handles)
function ExtensionFichier_KeyPressFcn(hObject, eventdata, handles)

%% =============== Paramètres Température électronique ===============
%% Te Min
function TeMin_Callback(InputTeMin, eventdata, handles)
function TeMin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Te Max
function TeMax_Callback(InputTeMax, eventdata, handles)
function TeMax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Step pour Te
function TeStep_Callback(InputTeStep, eventdata, handles)
function TeStep_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% ============================== Paramètres Densité électronique ==============================
%% Ne MIN
function NmMin_Callback(InputNmMin, eventdata, handles)
function NmMin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Ne MAX
function NmMax_Callback(InputNmMax, eventdata, handles)
function NmMax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Nombre de valeurs de Ne
function NmPts_Callback(InputNmPts, eventdata, handles)
function NmPts_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% ============================== Température du gaz ==============================
function Temperature_Callback(InputTg, eventdata, handles)
function Temperature_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% ============================== Pression d'opération ==============================
function Pression_Callback(InputPression, eventdata, handles)
function Pression_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ChoixHautePression_Callback(hObject, eventdata, handles)
function choixAutoabs_Callback(hObject, eventdata, handles)
function exposant_Callback(hObject, eventdata, handles)
function exposant_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% ============================== Longueur d'absorption ==============================
function Absorption_Callback(InputAbsorption, eventdata, handles)
function Absorption_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% ============================== Flow de gas ==============================
function MainGas_Callback(hObject, eventdata, handles)
function MainGas_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function TRGflow_Callback(hObject, eventdata, handles)
function TRGflow_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Mainflow_Callback(hObject, eventdata, handles)
function Mainflow_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% ============================== Sélection des graphiques à afficher ==============================
function GraphExp_Callback(hObject, eventdata, handles)
function Graphe1_Callback(hObject, eventdata, handles)
function GrapheFinaux_Callback(hObject, eventdata, handles)
function Graphe2_Callback(hObject, eventdata, handles)       
     
%% ======================= Selection de l'overwrite des raies =======================
% Ar
function Raie696_Callback(hObject, eventdata, handles)
function Raie727_Callback(hObject, eventdata, handles)
function Raie826_Callback(hObject, eventdata, handles)
function Raie738_Callback(hObject, eventdata, handles)
function Raie840_Callback(hObject, eventdata, handles)
function Raie794_Callback(hObject, eventdata, handles)
function Raie852_Callback(hObject, eventdata, handles)
function Raie751_Callback(hObject, eventdata, handles)
function Raie763_Callback(hObject, eventdata, handles)
function Raie800_Callback(hObject, eventdata, handles)
function Raie810_Callback(hObject, eventdata, handles)
function Raie866_Callback(hObject, eventdata, handles)
function Raie811_Callback(hObject, eventdata, handles)
function Raie801_Callback(hObject, eventdata, handles)
function Raie842_Callback(hObject, eventdata, handles)
function Raie706_Callback(hObject, eventdata, handles)
function Raie750_Callback(hObject, eventdata, handles)
function Raie667_Callback(hObject, eventdata, handles)
function Raie714_Callback(hObject, eventdata, handles)
 %Kr
function Raie758_Callback(hObject, eventdata, handles)
function Raie760_Callback(hObject, eventdata, handles)
function Raie768_Callback(hObject, eventdata, handles)
function Raie769_Callback(hObject, eventdata, handles)
function Raie785_Callback(hObject, eventdata, handles)
function Raie810Kr_Callback(hObject, eventdata, handles)
function Raie811Kr_Callback(hObject, eventdata, handles)
function Raie819_Callback(hObject, eventdata, handles)
function Raie829_Callback(hObject, eventdata, handles)
function Raie850_Callback(hObject, eventdata, handles)
function Raie877_Callback(hObject, eventdata, handles) 
function Raie892_Callback(hObject, eventdata, handles)  
 %Xe   
function Raie788_Callback(hObject, eventdata, handles)
function Raie823_Callback(hObject, eventdata, handles)
function Raie828_Callback(hObject, eventdata, handles)
function Raie881_Callback(hObject, eventdata, handles)
function Raie904_Callback(hObject, eventdata, handles)
function Raie820_Callback(hObject, eventdata, handles)
function Raie834_Callback(hObject, eventdata, handles)
function Raie895_Callback(hObject, eventdata, handles)  
 %Ne  
function Raie585_Callback(hObject, eventdata, handles)
function Raie640_Callback(hObject, eventdata, handles)
    

    




    
%% ======================= Choix de correction par une fonction de réponse =======================    
function CheckFctRep_Callback(hObject, eventdata, handles)
ChoixFctRep=get(handles.CheckFctRep,'Value');
if ChoixFctRep==1
    FctRep=uigetfile('Complements/Fonction_Reponse');
    set(handles.FileFctRep, 'String', FctRep);
elseif ChoixFctRep==0
    set(handles.FileFctRep, 'String', 'Aucun fichier sélectionné');
end


%% ================================ Démarrer la comparaison ================================
    function pushbutton1_Callback(hObject, eventdata, handles)
         %% Raies
         %Ne  
         Fit640=get(handles.Raie640,'Value');
         Fit585=get(handles.Raie585,'Value');
         %Ar
         Fit750=get(handles.Raie750,'Value');
         Fit696=get(handles.Raie696,'Value');
         Fit727=get(handles.Raie727,'Value');
         Fit826=get(handles.Raie826,'Value');
         Fit706=get(handles.Raie706,'Value');
         Fit738=get(handles.Raie738,'Value');
         Fit840=get(handles.Raie840,'Value');
         Fit794=get(handles.Raie794,'Value');
         Fit852=get(handles.Raie852,'Value');
         Fit751=get(handles.Raie751,'Value');
         Fit763=get(handles.Raie763,'Value');
         Fit800=get(handles.Raie800,'Value');
         Fit810=get(handles.Raie810,'Value');
         Fit866=get(handles.Raie866,'Value');
         Fit801=get(handles.Raie801,'Value');
         Fit842=get(handles.Raie842,'Value');
         Fit811=get(handles.Raie811,'Value');
         Fit667=get(handles.Raie667,'Value');
         Fit714=get(handles.Raie714,'Value');
         
         %Kr
         Fit758 = get(handles.Raie758,'Value');
         Fit760 = get(handles.Raie760,'Value');
         Fit768 = get(handles.Raie768,'Value');
         Fit769 = get(handles.Raie769,'Value');
         Fit785 = get(handles.Raie785,'Value');
         Fit810Kr=get(handles.Raie810Kr,'Value');
         Fit811Kr=get(handles.Raie811Kr,'Value');
         Fit819 = get(handles.Raie819,'Value');
         Fit829 = get(handles.Raie829,'Value');
         Fit850 = get(handles.Raie850,'Value');
         Fit877 = get(handles.Raie877,'Value');
         Fit892 = get(handles.Raie892,'Value');
         %Xe
         Fit788 = get(handles.Raie788,'Value');
         Fit823 = get(handles.Raie823,'Value');
         Fit828 = get(handles.Raie828,'Value');
         Fit881 = get(handles.Raie881,'Value');
         Fit904 = get(handles.Raie904,'Value');
         Fit820 = get(handles.Raie820,'Value');
         Fit834 = get(handles.Raie834,'Value');
         Fit895 = get(handles.Raie895,'Value');
         
         Overwrite=[Fit585 Fit640 ... %Ne
                    Fit667 Fit696 Fit706 Fit714 Fit727 Fit738 Fit750 Fit751 Fit763 Fit794 Fit800 Fit801 Fit810 Fit811 Fit826 Fit840 Fit842 Fit852 Fit866 ...%Ar
                    Fit758 Fit760 Fit768 Fit769 Fit785 Fit810Kr Fit811Kr Fit819 Fit829 Fit850 Fit877 Fit892 ... %Kr
                    Fit788 Fit820  Fit823 Fit828 Fit834 Fit881 Fit895 Fit904]; %Xe
         Overwrite=abs(1-Overwrite);
         
         %% rejette le loop pour un gaz qu'on ne check pas!!
         global gaz_i gaz_f
         gaz_f = 5 ;
         gaz_i = 2 ;
         if  sum(Overwrite(1:2)) == 0 %Ne
            gaz_i = 3 ;        
         end
         if  sum(Overwrite(34:41)) == 0 % pas de Xe
             if sum(Overwrite(22:33)) == 0 %pas de Kr
                 gaz_f = 3 ; 
             else 
                 gaz_f = 4 ;
             end
           
         end
         % si j'ai pas de Xe, on va jusquau Kr
         % si je n'ai pas de Xe et de Kr, on va jusqu'au Ar
         % si j'ai pas de Kr et jai du Xe, on fait toute.

         %% Calcul d'erreur
         
         ChoixErreur = get(handles.choixerreur,'value'); % 1= 1-sig 2= 1=1/sig 3= 1
         
         %% Graphes
         GraphExp=get(handles.GraphExp,'Value');    
         TeGraph1=get(handles.Graphe1,'Value');    
         TeGraph2=get(handles.Graphe2,'Value');
         TeGrapheFinaux=get(handles.GrapheFinaux,'Value');
         graphePopDepop = get(handles.graphePopDepop,'Value');
         
         InfoGraphes=[GraphExp TeGraph1 TeGraph2 TeGrapheFinaux graphePopDepop];
         
         %% Te
         %pour la STD
         TeMin=str2double(get(handles.TeMin,'String'));
         TeMax=str2double(get(handles.TeMax,'String'));
         TeStep=str2double(get(handles.TeStep,'String'));
         
         InfoTe=[TeMin TeMax TeStep];
                
         %% Ne
         %pour la STD normal sur un range décidé
         NeMin=str2double(get(handles.NmMin,'String'));
         NeMax=str2double(get(handles.NmMax,'String'));
         NePts=str2double(get(handles.NmPts,'String')); 
         
         InfoNe=[NeMin NeMax NePts];
         
         %% Exposant 1=maxwell, 2=Druyvesteyn
         exposant=str2double(get(handles.exposant,'String')); 
         
         %% Extension fichier
         contentsExtension = cellstr(get(handles.ExtensionFichier,'String'));
         if strcmp(contentsExtension{get(handles.ExtensionFichier,'Value')},'.csv')==1;
            filetype=0;
         elseif strcmp(contentsExtension{get(handles.ExtensionFichier,'Value')},'.csv (Morgane)')==1;
             filetype=0.1;
         elseif strcmp(contentsExtension{get(handles.ExtensionFichier,'Value')},'.trt')==1;
            filetype=1;
         elseif strcmp(contentsExtension{get(handles.ExtensionFichier,'Value')},'.trtx')==1;
            filetype=2;
         elseif strcmp(contentsExtension{get(handles.ExtensionFichier,'Value')},'.TXT')==1;
            filetype=4;   
         elseif strcmp(contentsExtension{get(handles.ExtensionFichier,'Value')},'.txt')==1;
            filetype=3;   
         end
                  
        %% Paramètres d'operation
         Tg=str2double(get(handles.Temperature,'String'));
         l=str2double(get(handles.Absorption,'String')); 
         sig_l=str2double(get(handles.sig_Absorption,'String'));       
        %% Choix High Pressure Cross Section
            ChoixHautePression =  get(handles.ChoixHautePression,'Value'); 
        %% Choix résonant et Autoabsorption
            ChoixAutoabs =  get(handles.choixAutoabs,'Value');    
        %% Commentaire à mettre en output
        Commentaire = get(handles.Commentaire,'String')
        %% Choix des dimension du réacteur
        if get(handles.ChoixDimension,'value')==1 %
            ChoixDimension='Dimension_Garofano.txt';
        elseif get(handles.ChoixDimension,'value')==2 %
            ChoixDimension='Dimension_C400.txt';
        elseif get(handles.ChoixDimension,'value')==3 %
            ChoixDimension='Simon 2.45GhZ.txt';
        elseif get(handles.ChoixDimension,'value')==4 %
            ChoixDimension='Dimension_HiPIMS.txt';
        end
         
        %% pressions partielles
         flow = zeros(1,25); % les 25 gaz de donnelly
         P_tot=str2double(get(handles.Pression,'String'));
         %Flow de TRG
         flow(1) = 0; % 0 % Helium
         flow(2) = 0.4*str2double(get(handles.TRGflow,'String'));% 40% Ne
         flow(3) = 0.2*str2double(get(handles.TRGflow,'String'));% 20% Ar
         flow(4) = 0.2*str2double(get(handles.TRGflow,'String'));% 20% Kr
         flow(5) = 0.2*str2double(get(handles.TRGflow,'String'));% 20% Xe
           
         %Ajout du Main flow
         if get(handles.MainGas,'value')==1 %Argon
         flow(3) = flow(3) + str2double(get(handles.Mainflow,'String'));
         end
         if get(handles.MainGas,'value')==2 %Helium
         flow(1) = flow(1) + str2double(get(handles.Mainflow,'String'));
         end
         if get(handles.MainGas,'value')==3 %O2
         flow(13) =  str2double(get(handles.Mainflow,'String'));
         end
         if get(handles.MainGas,'value')==4 %N2
         flow(14) = str2double(get(handles.Mainflow,'String'));
         end
         
         %% Flow de la trace de HMDSO 
         global nu_hmdso;
         nu_hmdso =str2double(get(handles.FluxHMDSO,'String'));
            %% Correction Pour la Needle Valve!!
%             needle_correction = sum(flow)/0.1833; %plus petit flow que ce qu'on envoit.
%             flow = flow/needle_correction;
%            clear needle_correction
%             
            
         %% Emplacements des dossiers
         current_folder = pwd; %dossier du programme
         Repertoire=get(handles.RepertoireSelect,'String');
         old=cd(Repertoire);
         addpath(strcat(current_folder,'/Complements')) %ajoute le dossier complements
         addpath(strcat(current_folder,'/Complements/Taux de reaction')) %ajout pour TeTauxDeReaction_TRG.m
         addpath(strcat(current_folder,'/Complements/Fonction_Reponse'))
         %% Fonction de réponse
         ChoixFctRep=get(handles.CheckFctRep,'Value');
         if ChoixFctRep==1
             FctRep=get(handles.FileFctRep,'String');
             FctRep=load(FctRep);
             FctRepName=fieldnames(FctRep);
             FctRep=FctRep.(FctRepName{1});
         elseif ChoixFctRep==0
             FctRep=0;
         end
        
          %% constante
          global c 
          c = 299792458; %vitesse de la lumière
          
         
         %% Execution du code
         tic %timer
         Calcul_Te_TRG_ViaInterface(Overwrite,InfoGraphes,InfoTe,InfoNe,filetype,Tg,P_tot,l,FctRep,ChoixErreur,ChoixHautePression,ChoixAutoabs,Commentaire,flow,exposant,ChoixDimension,sig_l); 
         cd(old);
         toc %fin du timer
         beep % fait un son pour la fin du calcul
%          load handel % fait un son pour la fin du calcul
%          sound(y,Fs) % fait un son pour la fin du calcul

% ================================================================
%% Boutons qui coches les bonne raies pour Te_High, Low et Tail
% ================================================================
% --- Executes on button press in Te_High.
function Te_High_Callback(hObject, eventdata, handles) %PressButton
    
%Ne  
set(handles.Raie640,'Value',1);
set(handles.Raie585,'Value',1);
%Ar
set(handles.Raie667,'Value',1);
set(handles.Raie696,'Value',1);
set(handles.Raie706,'Value',1);
set(handles.Raie714,'Value',1);
set(handles.Raie727,'Value',1);
set(handles.Raie738,'Value',1);
set(handles.Raie750,'Value',0);
set(handles.Raie751,'Value',0);
set(handles.Raie763,'Value',1);
set(handles.Raie794,'Value',1);
set(handles.Raie800,'Value',1);
set(handles.Raie801,'Value',1);
set(handles.Raie810,'Value',1);
set(handles.Raie811,'Value',1);
set(handles.Raie826,'Value',1);
set(handles.Raie840,'Value',1);
set(handles.Raie842,'Value',1);
set(handles.Raie852,'Value',1);
set(handles.Raie866,'Value',1);
%Kr
set(handles.Raie758,'Value',0);
set(handles.Raie760,'Value',1);
set(handles.Raie768,'Value',0);
set(handles.Raie769,'Value',0);
set(handles.Raie785,'Value',1);
set(handles.Raie810Kr,'Value',1);
set(handles.Raie811Kr,'Value',1);
set(handles.Raie819,'Value',1);
set(handles.Raie829,'Value',0);
set(handles.Raie850,'Value',1);
set(handles.Raie877,'Value',1);
set(handles.Raie892,'Value',1);
%Xe
set(handles.Raie788,'Value',0);
set(handles.Raie820,'Value',1);
set(handles.Raie823,'Value',1);
set(handles.Raie828,'Value',0);
set(handles.Raie834,'Value',0);
set(handles.Raie881,'Value',1);
set(handles.Raie895,'Value',1);
set(handles.Raie904,'Value',1);

% --- Executes on button press in Te_Tail.
function Te_Tail_Callback(hObject, eventdata, handles)
%coche les raies appropriées
%Ne  
set(handles.Raie640,'Value',1);
set(handles.Raie585,'Value',0);
%Ar
set(handles.Raie667,'Value',1);
set(handles.Raie696,'Value',1);
set(handles.Raie706,'Value',1);
set(handles.Raie714,'Value',1);
set(handles.Raie727,'Value',1);
set(handles.Raie738,'Value',1);
set(handles.Raie750,'Value',0);
set(handles.Raie751,'Value',0);
set(handles.Raie763,'Value',1);
set(handles.Raie794,'Value',1);
set(handles.Raie800,'Value',1);
set(handles.Raie801,'Value',1);
set(handles.Raie810,'Value',1);
set(handles.Raie811,'Value',1);
set(handles.Raie826,'Value',1);
set(handles.Raie840,'Value',1);
set(handles.Raie842,'Value',1);
set(handles.Raie852,'Value',1);
set(handles.Raie866,'Value',1);
%Kr
set(handles.Raie758,'Value',1);
set(handles.Raie760,'Value',1);
set(handles.Raie768,'Value',1);
set(handles.Raie769,'Value',1);
set(handles.Raie785,'Value',1);
set(handles.Raie810Kr,'Value',1);
set(handles.Raie811Kr,'Value',1);
set(handles.Raie819,'Value',1);
set(handles.Raie829,'Value',1);
set(handles.Raie850,'Value',1);
set(handles.Raie877,'Value',1);
set(handles.Raie892,'Value',1);
%Xe
set(handles.Raie788,'Value',1);
set(handles.Raie820,'Value',1);
set(handles.Raie823,'Value',1);
set(handles.Raie828,'Value',1);
set(handles.Raie834,'Value',1);
set(handles.Raie881,'Value',1);
set(handles.Raie895,'Value',1);
set(handles.Raie904,'Value',1);

% --- Executes on button press in Te_Low.
function Te_Low_Callback(hObject, eventdata, handles)
%coche les raies appropriées
%Ne  
set(handles.Raie640,'Value',1);
set(handles.Raie585,'Value',1);
%Ar
set(handles.Raie667,'Value',1);
set(handles.Raie696,'Value',1);
set(handles.Raie706,'Value',1);
set(handles.Raie714,'Value',1);
set(handles.Raie727,'Value',1);
set(handles.Raie738,'Value',1);
set(handles.Raie750,'Value',1);
set(handles.Raie751,'Value',1);
set(handles.Raie763,'Value',1);
set(handles.Raie794,'Value',1);
set(handles.Raie800,'Value',1);
set(handles.Raie801,'Value',1);
set(handles.Raie810,'Value',1);
set(handles.Raie811,'Value',0);
set(handles.Raie826,'Value',1);
set(handles.Raie840,'Value',1);
set(handles.Raie842,'Value',1);
set(handles.Raie852,'Value',1);
set(handles.Raie866,'Value',1);
%Kr
set(handles.Raie758,'Value',1);
set(handles.Raie760,'Value',0);
set(handles.Raie768,'Value',1);
set(handles.Raie769,'Value',1);
set(handles.Raie785,'Value',1);
set(handles.Raie810Kr,'Value',1);
set(handles.Raie811Kr,'Value',1);
set(handles.Raie819,'Value',0);
set(handles.Raie829,'Value',1);
set(handles.Raie850,'Value',1);
set(handles.Raie877,'Value',1);
set(handles.Raie892,'Value',1);
%Xe
set(handles.Raie788,'Value',1);
set(handles.Raie820,'Value',1);
set(handles.Raie823,'Value',0);
set(handles.Raie828,'Value',1);
set(handles.Raie834,'Value',1);
set(handles.Raie881,'Value',0);
set(handles.Raie895,'Value',1);
set(handles.Raie904,'Value',1);

% --- Executes on button press in Reset.
function Reset_Callback(hObject, eventdata, handles)
%Ne  
set(handles.Raie640,'Value',1);
set(handles.Raie585,'Value',1);
%Ar
set(handles.Raie667,'Value',1);
set(handles.Raie696,'Value',1);
set(handles.Raie706,'Value',1);
set(handles.Raie714,'Value',1);
set(handles.Raie727,'Value',1);
set(handles.Raie738,'Value',1);
set(handles.Raie750,'Value',0);
set(handles.Raie751,'Value',0);
set(handles.Raie763,'Value',0);
set(handles.Raie794,'Value',0);
set(handles.Raie800,'Value',0);
set(handles.Raie801,'Value',0);
set(handles.Raie810,'Value',0);
set(handles.Raie811,'Value',0);
set(handles.Raie826,'Value',1);
set(handles.Raie840,'Value',0);
set(handles.Raie842,'Value',0);
set(handles.Raie852,'Value',0);
set(handles.Raie866,'Value',0);
%Kr
set(handles.Raie758,'Value',0);
set(handles.Raie760,'Value',0);
set(handles.Raie768,'Value',1);
set(handles.Raie769,'Value',1);
set(handles.Raie785,'Value',0);
set(handles.Raie810Kr,'Value',1);
set(handles.Raie811Kr,'Value',1);
set(handles.Raie819,'Value',0);
set(handles.Raie829,'Value',0);
set(handles.Raie850,'Value',1);
set(handles.Raie877,'Value',1);
set(handles.Raie892,'Value',1);
%Xe
set(handles.Raie788,'Value',0);
set(handles.Raie820,'Value',0);
set(handles.Raie823,'Value',0);
set(handles.Raie828,'Value',0);
set(handles.Raie834,'Value',0);
set(handles.Raie881,'Value',1);
set(handles.Raie895,'Value',1);
set(handles.Raie904,'Value',1);



function Pression_KeyPressFcn(hObject, eventdata, handles)

function Temperature_KeyPressFcn(hObject, eventdata, handles)

function graphePopDepop_Callback(hObject, eventdata, handles)

function Commentaire_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function Commentaire_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ChoixDimension.
function ChoixDimension_Callback(hObject, eventdata, handles)
% hObject    handle to ChoixDimension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ChoixDimension contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ChoixDimension


% --- Executes during object creation, after setting all properties.
function ChoixDimension_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ChoixDimension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CHECK_w_1.
function CHECK_w_1_Callback(hObject, eventdata, handles)
% hObject    handle to CHECK_w_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CHECK_w_1



function sig_Absorption_Callback(hObject, eventdata, handles)
% hObject    handle to sig_Absorption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sig_Absorption as text
%        str2double(get(hObject,'String')) returns contents of sig_Absorption as a double


% --- Executes during object creation, after setting all properties.
function sig_Absorption_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sig_Absorption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in choixerreur.
function choixerreur_Callback(hObject, eventdata, handles)
% hObject    handle to choixerreur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns choixerreur contents as cell array
%        contents{get(hObject,'Value')} returns selected item from choixerreur


% --- Executes during object creation, after setting all properties.
function choixerreur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to choixerreur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PointHMDSO_Callback(hObject, eventdata, handles)
% hObject    handle to PointHMDSO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PointHMDSO as text
%        str2double(get(hObject,'String')) returns contents of PointHMDSO as a double


% --- Executes during object creation, after setting all properties.
function PointHMDSO_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PointHMDSO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FacteurHMDSO_Callback(hObject, eventdata, handles)
% hObject    handle to FacteurHMDSO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FacteurHMDSO as text
%        str2double(get(hObject,'String')) returns contents of FacteurHMDSO as a double


% --- Executes during object creation, after setting all properties.
function FacteurHMDSO_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FacteurHMDSO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FluxHMDSO_Callback(hObject, eventdata, handles)
% hObject    handle to FluxHMDSO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FluxHMDSO as text
%        str2double(get(hObject,'String')) returns contents of FluxHMDSO as a double


% --- Executes during object creation, after setting all properties.
function FluxHMDSO_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FluxHMDSO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
