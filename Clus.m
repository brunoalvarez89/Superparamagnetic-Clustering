function varargout = GUI(varargin)
% GUI M-file for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 19-Apr-2013 20:13:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_OutputFcn, ...
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


% --- Executes just before GUI is made visible.
function GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI (see VARARGIN)

% Choose default command line output for GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%p = uiextras.TabPanel( 'Parent', handles.uipanel1, 'Padding', 5 );
%c = uicontrol('Parent', p);
%axes('Parent', c);
%axes( 'Parent', p );
% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


a = zeros(1,3);
a(:) = 1:3;
set(handles.popup_distribuciones,'String', a);

function distribucion_1()

global tresde;
tresde = 0;

cluster_1 = 360;
cluster_2 = 720;
cluster_3 = 1080;

r1 = 1;
r2 = 2;
r3 = 3;

desvio_r1 = 0.15;
desvio_r2 = 0.20;
desvio_r3 = 0.25;

global cant_puntos;
cant_puntos = cluster_1 + cluster_2 + cluster_3;

global Distribucion;
Distribucion = zeros(cant_puntos, 2);

tita = 2 * pi / 360;

for i=1:cant_puntos
    
    if i < cluster_1
        
        R1 = normrnd(r1, desvio_r1);
    
        Distribucion(i,1) = R1 * cos(tita * i); 
        Distribucion(i,2) = R1 * sin(tita * i);
    
    end
    
    
    if i >= cluster_1 && i < (cluster_2 + cluster_1)
    
        R2 = normrnd(r2, desvio_r2);
    
        Distribucion(i,1) = R2 * cos(tita * i);
        Distribucion(i,2) = R2 * sin(tita * i);
        
    end
    
    
    if i >= (cluster_2 + cluster_1) && i <= cant_puntos
    
        R3 = normrnd(r3, desvio_r3);
    
        Distribucion(i,1) = R3 * cos(tita * i);
        Distribucion(i,2) = R3 * sin(tita * i);
        
    end
    
end

handles = guidata(gcf); 
axes(handles.Grafico_Puntos);
scatter(handles.Grafico_Puntos, Distribucion(:,1), Distribucion(:,2));
axis normal

calcular();

function distribucion_2()

global tresde;
tresde = 0;

cluster_1 = 360;
cluster_2 = 360;
cluster_3 = 360;
cluster_4 = 360;
cluster_5 = 200;
cluster_6 = 1080;

r1 = 3;
r2 = 3;
r3 = 3;
r4 = 3;
r5 = 3;
r6 = 40;

desvio_r1 = 4;
desvio_r2 = 4;
desvio_r3 = 4;
desvio_r4 = 4;
desvio_r5 = 4;
desvio_r6 = 15;

global cant_puntos;
cant_puntos = cluster_1 + cluster_2 + cluster_3 + cluster_4 + cluster_5 + cluster_6;

global Distribucion;
Distribucion = zeros(cant_puntos, 2);

tita = 2 * pi / 360;

for i=1:cluster_1
   
    R1 = normrnd(r1, desvio_r1);
    
    Distribucion(i, 1) = R1 * cos(tita * i) + 30;
    Distribucion(i, 2) = R1 * sin(tita * i) + 30;
    
end

for i=1:cluster_2
   
    R2 = normrnd(r2, desvio_r2);
    
    Distribucion(i + cluster_1, 1) = R2 * cos(tita * i) + 30;
    Distribucion(i + cluster_1, 2) = R2 * sin(tita * i) - 30;
    
end

for i=1:cluster_3
   
    R3 = normrnd(r3, desvio_r3);
    
    Distribucion(i + cluster_1 + cluster_2, 1) = R3 * cos(tita * i) - 30;
    Distribucion(i + cluster_1 + cluster_2, 2) = R3 * sin(tita * i) + 30;
    
end

for i=1:cluster_4
   
    R4 = normrnd(r4, desvio_r4);
    
    Distribucion(i + cluster_1 + cluster_2 + cluster_3, 1) = R4 * cos(tita * i) - 30;
    Distribucion(i + cluster_1 + cluster_2 + cluster_3, 2) = R4 * sin(tita * i) - 30;
    
end

for i=1:cluster_5
   
    R5 = normrnd(r5, desvio_r5);
    
    Distribucion(i + cluster_1 + cluster_2 + cluster_3 + cluster_4, 1) = R5 * cos(tita * i);
    Distribucion(i + cluster_1 + cluster_2 + cluster_3 + cluster_4, 2) = R5 * sin(tita * i);
    
end

for i=1:cluster_6
   
    R6 = normrnd(r6, desvio_r6);
    
    Distribucion(i + cluster_1 + cluster_2 + cluster_3 + cluster_4 + cluster_5, 1) = R6 * cos(tita * i);
    Distribucion(i + cluster_1 + cluster_2 + cluster_3 + cluster_4 + cluster_5, 2) = R6 * sin(tita * i);
    
end

handles = guidata(gcf); 
axes(handles.Grafico_Puntos);
scatter(handles.Grafico_Puntos, Distribucion(:,1), Distribucion(:,2));
axis normal

calcular();

function distribucion_3()

global tresde;
tresde = 1;

cluster_1 = 300;
cluster_2 = 300;
cluster_3 = 300;

r1 = 3;
r2 = 3;
r3 = 3;

desvio_r1 = 10;
desvio_r2 = 10;
desvio_r3 = 10;

global cant_puntos;
cant_puntos = cluster_1 + cluster_2 + cluster_3;

global Distribucion;
Distribucion = zeros(cant_puntos, 3);

tita = 2 * pi / 180;
fi = 2 * pi / 180;

for i=1:cluster_1
   
    R1 = normrnd(r1, desvio_r1);
        
        Distribucion(i, 1) = R1 * cos(tita * i) * sin(fi * i) + 15;
        Distribucion(i, 2) = R1 * sin(tita * i) * sin(fi * i) + 15;
        Distribucion(i, 3) = R1 * cos(fi * i) + 15;
    
end

for i=1:cluster_2
   
    R2 = normrnd(r2, desvio_r2);
    
    Distribucion(i + cluster_1, 1) = R2 * cos(tita * i) * sin(fi * i) - 15;
    Distribucion(i + cluster_1, 2) = R2 * sin(tita * i) * sin(fi * i) + 15;
    Distribucion(i + cluster_1, 3) = R2 * cos(fi * i) - 15;

end

for i=1:cluster_3
   
    R3 = normrnd(r3, desvio_r3);

    Distribucion(i + cluster_1 + cluster_2, 1) = R3 * cos(tita * i) * sin(fi * i);
    Distribucion(i + cluster_1 + cluster_2, 2) = R3 * sin(tita * i) * sin(fi * i);
    Distribucion(i + cluster_1 + cluster_2, 3) = R3 * cos(fi * i);
    
end

handles = guidata(gcf); 
axes(handles.Grafico_Puntos);
scatter3(handles.Grafico_Puntos, Distribucion(:,1), Distribucion(:,2), Distribucion(:,3));
axis normal
h=rotate3d;
set(h,'Enable','on'); 

calcular();

function calcular()

global Distribucion;

global cant_puntos;

global cant_vecinos_cercanos;
cant_vecinos_cercanos = 10;

global Vecinos_Cercanos;
Vecinos_Cercanos = zeros(cant_puntos, cant_vecinos_cercanos);

global Distancias;
Distancias = zeros(cant_puntos, cant_puntos);

global J;
J = zeros(cant_puntos, cant_vecinos_cercanos);

% Funcion moduladora de a
f = 0.45;

% Calculo distancias
for i=1:cant_puntos
    
    for j=1:cant_puntos
               
        Distancias(i,j) = (Distribucion(i,1) - Distribucion(j,1))^2 + (Distribucion(i,2) - Distribucion(j,2))^2;
        
        Distancias(i,j) = sqrt(Distancias(i,j));
                            
    end
    
end

clear Distribucion;

% Calculo vecinos cercanos.
for i=1:cant_puntos
    
    % Almaceno el valor C(j) y la coordenada I(j) de los m�nimos de
    % distancia
  
    [C,I] = sort(Distancias(i,:));      
    
    for j=1:cant_vecinos_cercanos
        
        % Salteo la distancia de Spike(i) consigo misma poniendo j+1
        Vecinos_Cercanos(i,j) = I(j+1);
               
    end   

end

% Contador que suma si y sólo si Spike(i) es vecino cercano de Spike(j) y 
% Spike(j) es vecino cercano de Spike(i).
aux1 = 0;

% Contador para calcular a _Promedio(almacena distancias entre vecinos cercanos).
aux2 = 0;
    
for i=1:cant_puntos
        
    for j=1:cant_vecinos_cercanos
        
        flag = 0;
        
        for l=1:cant_vecinos_cercanos
            
            if Vecinos_Cercanos(Vecinos_Cercanos(i,j),l) == i
                
                aux1 = aux1 + 1;
                                
                aux2 = aux2 + Distancias(i, Vecinos_Cercanos(i,j));
                   
                flag = 1;
                
            end
            
            if flag == 1
                
                break;
            
            end
            
        end
        
        if flag == 0
        
            Vecinos_Cercanos(i,j) = -1;
        
        end
        
    end
    
end

% Vecinos cercanos PROMEDIO 
K = aux1 / cant_puntos;

% Distancia promedio entre vecinos cercanos
a = aux2 / aux1;

% Calculo J.
for i=1:cant_puntos
       
    for j=1:cant_vecinos_cercanos
        
        if Vecinos_Cercanos(i,j) == -1
            
            break;
            
        end
        
        J(i,j) = 1/K * exp( - 1/2 * 1/(f * a)^2 * Distancias(i,Vecinos_Cercanos(i,j))^2);                         
        
    end
        
end

clear Distancias;

% --- Executes on button press in boton_monte_carlo.
function boton_monte_carlo_Callback(hObject, eventdata, handles)
% hObject    handle to boton_monte_carlo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
global Vecinos_Cercanos;
global J;
global cant_vecinos_cercanos;
global cant_puntos;

% Estados posibles de spin.
q = str2double(get(handles.edit_q,'String'));

% Incremento de temperatura.
global delta_T;
delta_T = 0.01;

% Temperatura inicial (T = 0 es semilla, no itera).
T_inicial = delta_T;

% Temperatura final
T_final = 0.2;

% Temperatura de transición ferromagnética-superparamagnética.
T_spm = 0;

% Temperatura de transición superparamagnética-paramagnética.
T_pm = 0;

% Temperatura a evaluar.
global T_deseada;
T_deseada = 0;

% Vector de temperaturas.
global vector_T;
vector_T = 0:delta_T:T_final;

% Número de iteración de temperatura.
step_T = 1;

% Número de iteraciones de temperatura que hará el algoritmo.
global step_T_max;
step_T_max = floor(T_final / delta_T);

% Cantidad máxima de iteraciones de Monte Carlo para cada temperatura.
m_max = str2double(get(handles.edit_m_max,'String'));

% Hamiltoniano
H = zeros(1, m_max);
global H_prom;
H_prom = zeros(1, step_T_max);

% Vector de spines.
% s(i) = spin de Spike(i).
s = zeros(cant_puntos, 1);

% Vector de Clústeres.
% Cluster(i) = A que cluster pertenece Spike(i).
global Cluster;
Cluster = zeros(cant_puntos, step_T_max);

% Vector de cantidad de spines.
% N(i) = Cantidad de espines en el estado q = i.
N = zeros(q, 1);

% Vector de magnetización.
% M(i,j) = magnetizaci�n del sistema en la iteraci�n i para la temperatura j * delta_T.
M = zeros(1, m_max);

% Vector de Magnetización promedio.
global M_prom;
M_prom = zeros(1, step_T_max);

% Vector de densidad de Susceptibilidad.
global S;
S = zeros(1, step_T_max);

% Tamaño en funcion de T de los 5 clusteres mas grandes
global tam_clusters
tam_clusters = zeros(5, step_T_max);

% Matriz de correlación espín-espin
G = zeros(cant_puntos, cant_vecinos_cercanos, step_T_max + 1);

% Umbral de correlación
tita = (1 - 1/q) / 2;

% Vector de punto ya visitado
Punto_Visitado = zeros(cant_puntos, 1);

% Variables en T = 0
M_prom(1) = 1;
S(1) = 0;
s(:) = 1;
S_max = 0;
flag = 0;
tam_clusters(1,1) = cant_puntos;
Cluster(:,1) = -1;
G(:, :, 1) = 1;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop principal de Monte Carlo
%

for T=T_inicial:delta_T:T_final
    
    step_T = step_T + 1;
    
    for i=cant_puntos
        
        s(i) = floor(1 + q*rand(1));
    
    end

    for m=1:m_max
        
        % Reseteo clusters.
        Cluster(:, step_T) = -1;

        % Reseteo Hamiltoniano
        H(:) = 0;
        
        % Reseteo cantidad de clusters.
        cant_clusters = 0;
            
        % Recorro todas las spikes.
        for i=1:cant_puntos
           
            for j=1:cant_vecinos_cercanos
                 
                % ¿ Es Spike(i) vecino cercano de Spike(j) ?
                if Vecinos_Cercanos(i,j) ~= -1
                    
                    % Calculo Hamiltoniano
                    if s(i) ~= s(Vecinos_Cercanos(i,j))
                        
                        H(m) = H(m) + J(i,j);
                        
                    end
                    
                    % ¿ Tiene Spike(i) el mismo valor de spin que Spike(j) ?
                    if s(i) == s(Vecinos_Cercanos(i,j)) 
        
                        % ¿ Se genera enlace congelado ?
                        if rand < (1 - exp(-1/T * J(i,j)))
                            
                            % Calculo conectividad
                            G(i, j, step_T) = G(i, j, step_T) + 1;
                                       
                            % ¿ Spike(j) no tiene clúster asignado ...
                            if Cluster(Vecinos_Cercanos(i,j), step_T) == -1
                
                                % ... y Spike(i) no tiene clúster asignado ?
                                if Cluster(i, step_T) == -1
                    
                                    % Se crea cluster nuevo.
                                    cant_clusters = cant_clusters + 1;
                    
                                    % Se lo asigno a ambas Spikes.
                                    Cluster(i, step_T) = cant_clusters;
                                    Cluster(Vecinos_Cercanos(i,j), step_T) = cant_clusters;
                    
                                else
                
                                % ... y Spike(i) tiene clúster asignado ?
                                % Le asigno a Spike(j) el mismo cluster que Spike(i)
                                Cluster(Vecinos_Cercanos(i,j), step_T) = Cluster(i, step_T);
                    
                                end
                
                            else
           
                            % ¿ Spike(j) tiene clúster asignado ...
                                % ... y Spike(i) no tiene clúster asignado ?
                                if Cluster(i, step_T) == -1
                    
                                % Le asigno a Spike(i) el mismo cluster que Spike(j).
                                Cluster(i,step_T) = Cluster(Vecinos_Cercanos(i,j), step_T);
                    
                                else
                
                                % ... y Spike(i) tiene clúster asignado ?                    
                                    % Chequeo cual es el cluster más chico (por lo tanto el que se creó primero)
                                    % Le asigno Cluster(j) a todas las Spikes cuyo cluster sea Cluster(i) y reduzco la cantidad de clústeres
                                    if Cluster(i, step_T) > Cluster(Vecinos_Cercanos(i,j), step_T)
                        
                                        for k=1:cant_puntos
                            
                                            if Cluster(k, step_T) == Cluster(i, step_T)
                                
                                                Cluster(k, step_T) = Cluster(Vecinos_Cercanos(i,j), step_T);
                            
                                            end
                            
                                            if Cluster(k, step_T) > Cluster (i, step_T)
                                
                                                Cluster(k, step_T) = Cluster(k, step_T) - 1;
                                
                                            end
                           
                                        end
                        
                                        cant_clusters = cant_clusters - 1;
                        
                                    end
                                    
                                    % Le asigno Cluster(i) a todas las Spikes cuyo cluster sea Cluster(j) y reduzco la cantidad de cl�steres;
                                    if Cluster(Vecinos_Cercanos(i,j), step_T) > Cluster(i, step_T) 
                        
                                        for k=1:cant_puntos
                            
                                            if Cluster(k, step_T) == Cluster(Vecinos_Cercanos(i,j), step_T)
                                
                                                Cluster(k, step_T) = Cluster(i, step_T);
                            
                                            end
                            
                                            if Cluster(k, step_T) > Cluster (Vecinos_Cercanos(i,j), step_T)
                                
                                                Cluster(k, step_T) = Cluster(k, step_T) - 1;
                                
                                            end
                            
                                        end
                        
                                        cant_clusters = cant_clusters - 1;
                        
                                    end
                    
                                end
                
                            end
          
                        end
                    
                    end
                
                end
                
            end
    
        end
              
        % Borro N.
        N(:) = 0;

        % Asigno los nuevos valores de N.
        for i=1:cant_puntos
    
            N(s(i)) = N(s(i)) + 1;

        end
              
        % Flipeo spines.
        for i=1:cant_clusters
    
            aux = 1 + q*rand(1);
    
            aux = floor(aux);
        
            for j=1:cant_puntos
        
                if Cluster(j, step_T) == i
        
                    s(j) = aux;
            
                end
    
            end
    
        end
        
        % Calculo N máximo.
        N_max = max(N);
        
        % Calculo magnetización.
        M(m) = ((q * (N_max / cant_puntos)) - 1) / (q - 1);
               
    end
       
    % Calculo densidad de Susceptibilidad (X * T / N).
    S(step_T) = var(M);
    
    % Calculo Magnetización promedio.
    M_prom(step_T) = mean(M);
       
    % Calculo Hamiltoniano promedio.
    H_prom(step_T) = mean(H);
    
    % Normalizo Hamiltoniano
    H_prom(:) = H_prom(:)/max(H_prom);
    
    % Calculo temperatura de transición superparamagnética
    if S(step_T) > S_max
        
        S_max = S(step_T);
    
        T_spm = (step_T - 1) * delta_T;
        
    end
    
    % Calculo temperatura de transicion paramagnetica.
    if step_T > 2
        
        if S(step_T) < S_max * 0.1 && S(step_T) < S(step_T - 2) && flag == 0
            
            T_pm = (step_T - 1) * delta_T;
        
            flag = 1;
        
            % break;
            
        end
        
    end
    
    % Calculo matriz de correlación
    for i=1:cant_puntos
        
        for j=1:cant_vecinos_cercanos
   
            G(i, j, step_T) =((G(i, j, step_T) / m_max * (q - 1)) + 1) / q;
        
        end
        
    end
    
    % Genero clusters
    for i=1:cant_puntos
    
        for j=1:cant_vecinos_cercanos
            
            if Vecinos_Cercanos(i,j) ~= -1
            
                if G(i, j, step_T) > tita && Cluster(i) == Cluster(Vecinos_Cercanos(i,j))
            
                    Punto_Visitado(i) = 1; 
                    
                    Punto_Visitado(Vecinos_Cercanos(i,j)) = 1;
                
                end
            
                if G(i, j, step_T) < tita && Punto_Visitado(i) ~= 1 && Punto_Visitado(Vecinos_Cercanos(i,j)) ~= 1
                        
                    Cluster(i, step_T) = -1;
                        
                    Cluster(Vecinos_Cercanos(i,j), step_T) = -1;
                    
                end
                             
            end
            
        end
    
    end
    
    % Número de clústeres en T = Tclus
    cluster_max = max(Cluster(:, step_T));

    % cant_elementos(i) = cuantos puntos tiene el cluster i 
    cant_elementos = zeros( cluster_max, 1);

    % En la última posición de cant_elementos se almacena la cantidad de puntos
    % no clusterizados
    %cant_elementos(cluster_max + 1) = sum(Cluster(:,step_T) == -1);

    for i=1:cluster_max
            
        % Lleno con la cantidad de puntos de cada cluster
        cant_elementos(i) = sum(Cluster(:,step_T) == i);
        
    end
    
    % Almaceno los tamaños de los 5 clusteres mas grandes
    for v=1:5
        
         
         [C, I] = max(cant_elementos);
         
         if I ~= (cluster_max + 1)
             
            tam_clusters(v, step_T) = C;
            
         end
         
         cant_elementos(I) = -1;
    
    end
        
    % Muestro el paso de temperatura
    set(handles.Texto_T,'String', T)
    drawnow
    
end


% Calculo temperatura de clustering
T_deseada = (T_spm + T_pm) / 2;
set(handles.Texto_Tclus,'String', T_deseada);

Graficar_Curvas();
Graficar_Tamano_Clusters();
Graficar_Clusters();

function Graficar_Curvas()

global vector_T;
global M_prom;
global S;
global H_prom;
global T_deseada;

handles = guidata(gcf); 

% Ploteo Magnetización + Hamiltoniano
axes(handles.Grafico_Magnetizacion);
plot(handles.Grafico_Magnetizacion, vector_T, M_prom, vector_T, H_prom)
xlabel('T');
title('Magnetizacion & Energia')
line([T_deseada T_deseada], [0 1], 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
axis tight;

% Ploteo Susceptibilidad.
axes(handles.Grafico_Susceptibilidad);
plot(handles.Grafico_Susceptibilidad, vector_T, S)
xlabel('T')
title('Susceptibilidad')
line([T_deseada T_deseada], [0 max(S)], 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
axis tight;


function Graficar_Tamano_Clusters()

global vector_T;
global tam_clusters;
global T_deseada;
global cant_puntos;


handles = guidata(gcf); 

% Ploteo tamaño de clusters
axes(handles.Grafico_Hamiltoniano);
plot(handles.Grafico_Hamiltoniano, vector_T, tam_clusters(1,:), vector_T, tam_clusters(2,:), vector_T, tam_clusters(3,:), vector_T, tam_clusters(4,:), vector_T, tam_clusters(5,:));
xlabel('T')
title('Tamano de clusters')
line([T_deseada T_deseada], [0 cant_puntos], 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
axis tight;

tam_clusters(:,floor(T_deseada/0.01)+1)
function Graficar_Clusters()

global T_deseada;
global Cluster;
global delta_T;
global Distribucion;
global cant_puntos;
global tresde;
global step_T_max;

Tclus = T_deseada;

% Array con colores
colores = {'m', 'g', 'c', 'r', 'b', 'y', 'k'};
iterador_colores = 1;

% Array con simbolos
simbolos = {'o', '.', '*', '+', 's', 'd', '^', 'v', '>', '<', 'p', 'h', 'x'};
iterador_simbolos = 1;

% Step de temperatura correspondiente a T = Tclus
step_T_Tclus = floor((Tclus / delta_T)) + 1;

% Número de clústeres en T = Tclus
cluster_max = max(Cluster(:, step_T_Tclus));

% cant_elementos(i) = cuantos puntos tiene el cluster i 
cant_elementos = zeros( cluster_max + 1, 1);

% En la última posición de cant_elementos se almacena la cantidad de puntos
% no clusterizados
cant_elementos(cluster_max + 1) = sum(Cluster(:,step_T_Tclus) == -1);

for i=1:cluster_max
            
    % Lleno con la cantidad de puntos de cada cluster
    cant_elementos(i) = sum(Cluster(:,step_T_Tclus) == i);
        
end

% Obtengo cuantos puntos tiene el cluster mas grande (esto sirve para tener
% un criterio para los puntos sin clasificar, que pertenecen a clusters muy
% pequeños)
[valor_max, ~] = max(cant_elementos);

% Almaceno en la última posición de cant_elementos todos los puntos que
% respondan al criterio de "no clasificados". En este caso el criterio es
% que el cluster que conforman tiene que tener menos puntos que el 10% del
% cluster mas grande. Se coloca un cero en la entrada que responda ante dicho
% criterio para después eliminar dichos clusteres y condensarlos todos en
% un unico cluster sin clasificar.
for i=1:cluster_max
    
    if cant_elementos(i) < valor_max * 0.1
        
        cant_elementos(cluster_max + 1) = cant_elementos(cluster_max + 1) + cant_elementos(i);
        
        cant_elementos(i) = 0;
        
    end
        
end

% Obtengo la cantidad de elementos NO NULOS de cant_elementos (o sea,
% correspondientes a clusteres clasificados)
array_clusters = zeros(nnz(cant_elementos), 1);

% Obtengo los indices de los elementos NO NULOS
ind = find(cant_elementos);

% Genero un array del tamaño de la cantidad de elementos no nulos de
% cant_elementos y en cada entrada coloco el número de cluster clasificado.
for i=1:size(array_clusters)
    
    array_clusters(i) = ind(i);
    
end

handles = guidata(gcf); 

if tresde == 0
% Genero dos vectores con los puntos de x e y de la distribucion para poder
% ir encontrando  y ploteando los clusteres.
x = Distribucion(:, 1);
y = Distribucion(:, 2);

% Vector auxiliar para facilitar el ploteo de los puntos sin clasificar
auxi = zeros(cant_puntos, 1);
auxi(:) = -1;

axes(handles.Grafico_Clusters);

for i=1:(size(array_clusters))
   
    % Vector logico que contiene los indices de los puntos correspondientes
    % al numero de cluster de array_clusters(i).
    ind = ( Cluster(:, step_T_Tclus) == array_clusters(i));
                
    % Almaceno en dos vectores nuevos (correspondientes a x e y) los puntos
    % que corresponden al cluster
    xx = x(ind);
    yy = y(ind);
        
    % Coloco un 1 en el auxiliar indicando que el punto ESTA CLASIFICADO
    auxi(ind) = 1;
        
    % Ploteo
    
    scatter(handles.Grafico_Clusters, xx, yy, simbolos{iterador_simbolos}, 'MarkerEdgeColor', colores{iterador_colores}); 
    hold on;
    axis normal

    
    iterador_colores = iterador_colores + 1;
    
    if iterador_colores == 7
        
        iterador_simbolos = iterador_simbolos + 1;
        
        iterador_colores = 1;
        
    end
   
end

% Vector logico que contiene los indices de los puntos correspondientes
% al numero de cluster -1, o sea, SIN CLASIFICAR.
ind = ( auxi == -1) ;     

% Almaceno en dos vectores nuevos (correspondientes a x e y) los puntos
% que corresponden al cluster -1
xx = x(ind);
yy = y(ind);
        
% Ploteo cluster -1, siempre con color negro y cruces.
scatter(handles.Grafico_Clusters, xx, yy, 'x','MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); 
           
hold off;

zoom on;
end

if tresde == 1

x = Distribucion(:, 1);
y = Distribucion(:, 2);
z = Distribucion(:, 3);

% Vector auxiliar para facilitar el ploteo de los puntos sin clasificar
auxi = zeros(cant_puntos, 1);
auxi(:) = -1;

axes(handles.Grafico_Clusters);

for i=1:(size(array_clusters))
   
    % Vector logico que contiene los indices de los puntos correspondientes
    % al numero de cluster de array_clusters(i).
    ind = ( Cluster(:, step_T_Tclus) == array_clusters(i));
                
    % Almaceno en dos vectores nuevos (correspondientes a x e y) los puntos
    % que corresponden al cluster
    xx = x(ind);
    yy = y(ind);
    zz = z(ind);
        
    % Coloco un 1 en el auxiliar indicando que el punto ESTA CLASIFICADO
    auxi(ind) = 1;
        
    % Ploteo
    
    scatter3(handles.Grafico_Clusters, xx, yy, zz, simbolos{iterador_simbolos}, 'MarkerEdgeColor', colores{iterador_colores}); 
    hold on;
    axis normal

    
    iterador_colores = iterador_colores + 1;
    
    if iterador_colores == 7
        
        iterador_simbolos = iterador_simbolos + 1;
        
        iterador_colores = 1;
        
    end
   
end

% Vector logico que contiene los indices de los puntos correspondientes
% al numero de cluster -1, o sea, SIN CLASIFICAR.
ind = ( auxi == -1) ;     

% Almaceno en dos vectores nuevos (correspondientes a x e y) los puntos
% que corresponden al cluster -1
xx = x(ind);
yy = y(ind);
zz = z(ind);

% Ploteo cluster -1, siempre con color negro y cruces.
scatter3(handles.Grafico_Clusters, xx, yy, zz, 'x','MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); 
           
hold off;
end

% Seteo slider de temperaturas
slider_Min = 1;
slider_Max = step_T_max;
slider_Step = [1, 1] / (slider_Max - slider_Min);
slider_T_Tclus = step_T_Tclus;
set(handles.slider_T, 'Min', slider_Min);
set(handles.slider_T, 'Max', slider_Max);
set(handles.slider_T, 'SliderStep', slider_Step);
set(handles.slider_T, 'Value', slider_T_Tclus);


% --- Executes on selection change in popup_distribuciones.
function popup_distribuciones_Callback(hObject, eventdata, handles)
% hObject    handle to popup_distribuciones (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_distribuciones contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_distribuciones

valor = get(hObject,'Value');

if valor == 1
    
    distribucion_1();
    
end

if valor == 2
    
    distribucion_2();
    
end

if valor == 3
    
    distribucion_3();
    
end


% --- Executes during object creation, after setting all properties.
function popup_distribuciones_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_distribuciones (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_m_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_m_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_m_max as text
%        str2double(get(hObject,'String')) returns contents of edit_m_max as a double

% --- Executes during object creation, after setting all properties.
function edit_m_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_m_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_q_Callback(hObject, eventdata, handles)
% hObject    handle to edit_q (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_q as text
%        str2double(get(hObject,'String')) returns contents of edit_q as a double


% --- Executes during object creation, after setting all properties.
function edit_q_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_q (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_delta_T_Callback(hObject, eventdata, handles)
% hObject    handle to edit_delta_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_delta_T as text
%        str2double(get(hObject,'String')) returns contents of edit_delta_T as a double


% --- Executes during object creation, after setting all properties.
function edit_delta_T_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_delta_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function slider_T_Callback(hObject, eventdata, handles)
% hObject    handle to slider_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global delta_T;
global T_deseada;

step_T = get(hObject,'Value');

T_deseada = (step_T - 1) * delta_T;

set(handles.Texto_Tclus,'String', T_deseada);

if T_deseada ~= 0
    
    Graficar_Curvas();
    
    Graficar_Tamano_Clusters();
    
    Graficar_Clusters();

end

% --- Executes during object creation, after setting all properties.
function slider_T_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
