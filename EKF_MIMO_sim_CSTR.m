clear all
close all
clc

import casadi.*
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x3 = SX.sym('x3');

u1 = SX.sym('u1');
u2 = SX.sym('u2');
u3 = SX.sym('u3');

qi = SX.sym('qi');
Ti = SX.sym('Ti');

d = 5;
delay_real     = [d d d];
delay_modelado = delay_real+10;
% delay_modelado = [3 5 3];
d_max          = max(delay_real);
dmodelado_max  = max(delay_modelado);
delay_total    = d_max+dmodelado_max;
%% CSTR MIMO instável
global Ts Ac dH p cp ER k0 V
%---- Constantes ----
Ac = 0.05;        %Area do reator
dH = -50e6;       %Calor da reação
p  = 1e6;         %Densidade dos líquidos
cp = 1;           %Calor específico
ER = 1000;        %Termo de energia de ativação
k0 = 4.872997584; %Constante da taxa de reação
V  = 0.05;        %Volume do reator

Ts = 3;           %Período de amostragem

%---- Ponto de Operação ----

% Entradas
qo0  = 5e-3;
Caf0 = 5;
Qh0  = 0.75;

u0 = [qo0; Caf0; Qh0];

% Perturbações
qi0 = 5e-3;
Ti0 = 350;

q0  = [qi0; Ti0];

% Saídas
h0  = 1;
Ca0 = 1;
T0 = 400;

x0  = [h0; Ca0; T0];

%Estados / Variáveis Controladas
% x1 = h - nível dentro do tanque; 
% x2 = Ca - Concentração de saída do produto A;
% x3 = T  - Temperatura dentro do reator; 

%Entradas / Variáveis Manipuladas
% u1 = q0     - Vazão de saída;
% u2 = Caf    - Concentração do produto A na alimentação do tanque; 
% u3 = Qh/pCp - taxa de remoção de calor normalizada;

%Perturbações
% qi = qi - Vazão de entrada; 
% Ti = Ti - Temperatura externa;


%---- Modelo ----

Rt = (k0*exp(-ER/x3));

dx1_discreto = x1 + Ts*((qi-u1)/Ac);
dx2_discreto = x2 + Ts*(qi*((u2-x2)/V) - Rt*x2);
dx3_discreto = x3 + Ts*(qi*((Ti-x3)/V) - (dH/p/cp)*Rt*x2 - u3/V);

fun_ax_ext = [dx1_discreto;dx2_discreto;dx3_discreto;qi;Ti];     %Sistema aumentado discreto
fun_x      = Function('fun_x',{x1,x2,x3,qi,Ti,u1,u2,u3},{fun_ax_ext});

fun_yx_ext = [x1;x2;x3];
fun_y      = Function('fun_y',{x1,x2,x3,qi,Ti,u1,u2,u3},{fun_yx_ext});

na = size(fun_ax_ext,1);                      %Dimensão vetor de estados aumentado
m  = size(fun_yx_ext,1);                      %Dimensão vetor de saídas
%% --------------- Inicialização das variáveis --------------------
iteracoes   = 900;

%---- Saídas ----
saidas       = zeros(3,iteracoes);            %inicia vetor de saídas

saidas(1,:)  = x0(1); %1   h0
saidas(2,:)  = x0(2); %1   Ca0
saidas(3,:)  = x0(3); %400 T0

%---- Entradas ----
entradas      = zeros(3,iteracoes);              %inicia vetor de entradas

entradas(1,:) = u0(1); %0.005  q0
entradas(2,:) = u0(2); %5      Caf0
entradas(3,:) = u0(3); %0.75   Qh/pcp 0

%---- Perturbações ----

perturbacoes      = zeros(2,iteracoes);              %inicia vetor de perturbações

perturbacoes(1,:) = q0(1);%0.005  qi0                    
perturbacoes(2,:) = q0(2);%350    Ti0 

%---- Inicilização das referências de controle ----
ref_h            = zeros(1,iteracoes);
ref_h(1:end)     = saidas(1,1);

ref_ca           = zeros(1,iteracoes);
ref_ca(1:end)    = saidas(2,1);

ref_T            = zeros(1,iteracoes);
ref_T(1:end)     = saidas(3,1);

%---- Ruido - ainda não sendo utilizado ----
ruido = zeros(3,iteracoes);

%---- Variáveis De Estimação -----
x_estimado  = [saidas(1,1);saidas(2,1);saidas(3,1)];     %Cria vetor de estados estimados

d_estimado  = [perturbacoes(1,1);perturbacoes(2,1)];     %Cria vetor de perturbações estimadas

y_estimado  = [saidas(1,1);saidas(2,1);saidas(3,1)];     %Cria vetor de saídas estimadas

x_a_estim   = [x_estimado;d_estimado];                   %Compila as estimações iniciais de estados e perturbações

y_estimado_vect      = y_estimado;
x_estimado_vect      = zeros(na,iteracoes);
x_estimado_vect(1,:) = x_estimado(1);
x_estimado_vect(2,:) = x_estimado(2);
x_estimado_vect(3,:) = x_estimado(3);
x_estimado_vect(4,:) = d_estimado(1);
x_estimado_vect(5,:) = d_estimado(2);

x_a_pred              = zeros(na,iteracoes);
x_pred_vect           = zeros(na,iteracoes);
x_pred_vect(1,:)      = x_estimado(1);
x_pred_vect(2,:)      = x_estimado(2);
x_pred_vect(3,:)      = x_estimado(3);
x_pred_vect(4,:)      = d_estimado(1);
x_pred_vect(5,:)      = d_estimado(2);

%---- Parâmetros do Estimador ----   
%h - nível dentro do tanque (0~2); Ca - Concentração de saída do produto A (1~3); T  - Temperatura dentro do reator (390~410); 
%qi - Vazão de entrada(0.001~0.007); %Ti - Temperatura externa(350~370);
h_range=[0 2];       Ca_range=[1 3];         T_range=[390 410];       qi_range=[0.001 0.007];      Ti_range=[350 370];
h_med=mean(h_range); Ca_med=mean(Ca_range);  T_med=mean(T_range);     qi_med=mean(qi_range);       Ti_med=mean(Ti_range);
h_dp=abs(h_med-2)/2; Ca_dp=abs(Ca_med-3)/2;  T_dp=abs(T_med-410)/2;   qi_dp=abs(qi_med-0.007)/2;   Ti_dp=abs(Ti_med-370)/2;

% Q = diag([1*ones(1,na-2),1,1]);               %Variável da ponderação dos estados
% R = diag(ones(1,m));                          %Variável da ponderação da saída       
P_estimado_at = diag([1*ones(1,na-2),1,1]); 

%Normalização = xi-media/desvio
% x1_norm = (x_pred_vect(1,1) - h_med)/h_dp;
% x2_norm = (x_pred_vect(2,1) - Ca_med)/Ca_dp;
% x3_norm = (x_pred_vect(3,1) - T_med)/T_dp;
% d1_norm = (x_pred_vect(4,1) - qi_med)/qi_dp;
% d2_norm = (x_pred_vect(5,1) - Ti_med)/Ti_dp;

% Q = diag([h_dp^2, Ca_dp^2, T_dp^2, qi_dp^2, Ti_dp^2]);
% R = diag([((h_dp*0.9)^2), ((Ca_dp*1.5)^2), ((T_dp*0.9)^2)]);

% Testes
Q = diag([(h_dp*0.7)^2, (Ca_dp*0.25)^2, (T_dp*0.3)^2, (qi_dp*300)^2, (Ti_dp*0.45)^2]);
R = diag([((h_dp*8)^2), ((Ca_dp*2)^2), ((T_dp*4)^2)]);

%---- Funções CasADi ----
A_jacobian = jacobian(fun_ax_ext,[x1 x2 x3 qi Ti]);  %Cálculo do jacobiano para matriz A
B_jacobian = jacobian(fun_ax_ext,[u1 u2 u3]);        %Cálculo do jacobiano para matriz B
C_jacobian = jacobian(fun_yx_ext,[x1 x2 x3 qi Ti]);  %Cálculo do jacobiano para matriz C
fun_a = Function('fun_a',{x1,x2,x3,qi,Ti,u1,u2,u3},{A_jacobian});
fun_b = Function('fun_b',{x1,x2,x3,qi,Ti,u1,u2,u3},{B_jacobian});
fun_c = Function('fun_a',{x1,x2,x3,qi,Ti,u1,u2,u3},{C_jacobian});

%% --------------- Projeto de Controle ---------------
for lqr=1
%---- Linearização no ponto de operação ----
Aa = full(fun_a(saidas(1,1),saidas(2,1),saidas(3,1),perturbacoes(1,1),perturbacoes(2,1),entradas(1,1),entradas(2,1),entradas(3,1)));
Ba = full(fun_b(saidas(1,1),saidas(2,1),saidas(3,1),perturbacoes(1,1),perturbacoes(2,1),entradas(1,1),entradas(2,1),entradas(3,1)));
Ca = full(fun_c(saidas(1,1),saidas(2,1),saidas(3,1),perturbacoes(1,1),perturbacoes(2,1),entradas(1,1),entradas(2,1),entradas(3,1)));
Ba_x = Ba(1:na-2,1:3);
Aa_x = Aa(1:na-2,1:3);
Ca_x = Ca(1:na-2,1:3);

%---- Insere dinâmica integradora no controlador ----
A_exp = [Aa_x zeros(size(Aa_x,1),size(Ca_x,1));
         -Ca_x         [1 0 0; 0 1 0;0 0 1]];
B_exp = [Ba_x;zeros(3,3)];
C_exp = [Ca_x zeros(3,3)];
D_exp = 0;
sistema = ss(A_exp,B_exp,C_exp,D_exp,Ts);                            %Sistema expandido

Co = ctrb(sistema.A,sistema.B);                                      %Matriz de Controlabilidade
rank(Co);
Q_lqr = diag([1,1,1/400^2,1,1,1/400^2]);                             %Matrizes de ponderação
R_lqr = 30*diag([1/0.005^2,1/5^2,1/0.75^2]);
[K,P,poles] = dlqr(sistema.A,sistema.B,Q_lqr,R_lqr);                 %Ganho do controlador via dlqr

%---- Inicializa integrador no ponto de operação ---- 
integrador_anterior = -pinv(K(:,4:6))*(entradas(:,1)+K(:,1:3)*x_estimado);

integrador_anterior_1 = integrador_anterior(1);
integrador_atual_1 = 0;

% integrador_anterior_2 = (entradas(:,1)+K(1:3,1:3)*x_estimado)/(-K(1:3,5));
integrador_anterior_2 = integrador_anterior(2);
integrador_atual_2 = 0;
% end

integrador_anterior_3 = integrador_anterior(3);
integrador_atual_3 = 0;
end

%Variável para fechar a malha de controle; 0 = malha aberta; 1 = malha fechada
controle = 0;

%% --------------- Definição das referências ---------------
if controle == 1 % --- MALHA FECHADA ---
    
%     ref_h(200:end)     = saidas(1,1)*1.02;          
    ref_ca(100:end)    = saidas(2,1)*1.02;            
%     ref_T(600:end)     = saidas(3,1)*1.02;
    perturbacoes(1,200:end)      = q0(1)*0.98;    
    perturbacoes(2,400:end)      = q0(2)*0.98;

else % --- MALHA ABERTA ---
%     entradas(2,600:end) = u0(2)*1.02;               %Degrau na entrada de 2%
%     entradas(1,200:end) = u0(1)*0.98;             
%     entradas(3,500:end) = u0(3)*0.98;    
%     perturbacoes(1,100:end)      = q0(1)*1.02;    
%     perturbacoes(2,700:end)      = q0(2)*1.02;
end

%% --------------- Simulação --------------------
tic
for k = 2+delay_total:iteracoes

    %Compila vetor de saídas e perturbações
    xd = [saidas(:,k-1)' perturbacoes(:,k-1)']; 
    
    %Compila as entradas com atraso real
    entrada_atrasada_real = [entradas(1,k-1-delay_real(1)) entradas(2,k-1-delay_real(2)) entradas(3,k-1-delay_real(3))]';
    
    %Compila as entradas com atraso modelado
    entrada_atrasada_mod  = [entradas(1,k-1-delay_modelado(1)) entradas(2,k-1-delay_modelado(2)) entradas(3,k-1-delay_modelado(3))]';

%     entrada_atrasada_real_ruido = entrada_atrasada_real + ruido(k);
    
    %---- Simulação do processo ----
    saidas(:,k) = modeloCSTR_naoIsotermico(xd,entrada_atrasada_real);

    %---- Linearização a cada iteração para EKF ----
    Aa = full(fun_a(x_a_estim(1),x_a_estim(2),x_a_estim(3),x_a_estim(4),x_a_estim(5),...
                    entrada_atrasada_mod(1),entrada_atrasada_mod(2),entrada_atrasada_mod(3)));

    Ca = full(fun_c(x_a_estim(1),x_a_estim(2),x_a_estim(3),x_a_estim(4),x_a_estim(5),...
                    entrada_atrasada_mod(1),entrada_atrasada_mod(2),entrada_atrasada_mod(3)));

    %---- Estimação ---- 
    [x_a_estim,P_estimado] = ekfMIMO(fun_x,fun_y,saidas(:,k),...
                                     entrada_atrasada_mod,x_a_estim,P_estimado_at,Aa,Ca,Q,R);

    %---- Predição ----
    %Primeira iteração
    x_a_pred = x_a_estim;   
    if(dmodelado_max>0)
       %Segunda iteração em diante:
        for kk=1:dmodelado_max
            %Tratamento das entradas
            if(k-1-delay_modelado(1)+kk>k) 
                entradas(1,k-1-delay_modelado(1)+kk) = entradas(1,k-1-delay_modelado(1)+kk-1);
            end
            if(k-1-delay_modelado(2)+kk>k)
                entradas(2,k-1-delay_modelado(2)+kk) = entradas(2,k-1-delay_modelado(2)+kk-1);
            end
            if(k-1-delay_modelado(3)+kk>k)
                entradas(3,k-1-delay_modelado(3)+kk) = entradas(3,k-1-delay_modelado(3)+kk-1);
            end
            x_a_pred = modeloCSTR_naoIsotermico([x_a_pred(1:3);x_a_estim(4);x_a_estim(5)]',...
                                                [entradas(1,k-1-delay_modelado(1)+kk) entradas(2,k-1-delay_modelado(2)+kk) entradas(3,k-1-delay_modelado(3)+kk)]')';
        end
    end
    
    x_a_pred = [x_a_pred(1:3);x_a_estim(4);x_a_estim(5)]'; %Compila estados e perturbações preditas em um vetor
    
    %---- Controle ----
    if controle == 1
        %Integrador para altura do tanque - h
        erro_1                = ref_h(k) - x_a_pred(1);
        integrador_atual_1    = integrador_anterior_1 + erro_1;
        integrador_anterior_1 = integrador_atual_1;
                
        %Integrador para concentração do produto A - Ca
        erro_2                = ref_ca(k) - x_a_pred(2);
        integrador_atual_2    = integrador_anterior_2 + erro_2;
        integrador_anterior_2 = integrador_atual_2;
        
        %Integrador para temperatura interna - T
        erro_3                = ref_T(k) - x_a_pred(3);
        integrador_atual_3    = integrador_anterior_3 + erro_3;
        integrador_anterior_3 = integrador_atual_3;
      
        nova_entrada = -K*[x_a_pred(1);x_a_pred(2);x_a_pred(3); integrador_atual_1; integrador_atual_2; integrador_atual_3];
        entradas(:,k) = nova_entrada;
    end
    
    %---- Atualização das variáveis ----
    x_estimado_vect(:,k) = x_a_estim(1:end);
    x_pred_vect(:,k)     = x_a_pred';
    P_estimado_at        = P_estimado;
    
end

Elapsed_time = toc;
disp('Tempo gasto na simulação:');
disp(Elapsed_time);

%---- Plot de gráficos ----
plot_cstrNaoIsotermico
%% Função do modelo

function dxdt = odefcn(t, estados, entradas)
    global Ac dH p cp ER k0 V;
    
    Rt = (k0*exp(-ER/estados(3)));

    dxdt = zeros(5,1);
    dxdt(1) = ((estados(4)-entradas(1))/Ac);
    dxdt(2) = (estados(4)*((entradas(2)-estados(2))/V) - Rt*estados(2));
    dxdt(3) = (estados(4)*((estados(5)-estados(3))/V) - (dH/p/cp)*Rt*estados(2) - entradas(3)/V);
    dxdt(4) = 0;
    dxdt(5) = 0;
end

function x = modeloCSTR_naoIsotermico(estados,entradas)
    global Ts Ac dH p cp ER k0 V;
    %estados(1) = h; estados(2) = Ca; estados(3) = T ;
    %estados(4) = qi; estados(5) = Ti;
    %entradas(1) = qo; entradas(2) = Caf; entradas(3) = Qh/pCp;
   
%     tspan = [0 Ts];
%     [t, xfun] = ode45(@(t,x) odefcn(t, x, entradas), tspan, estados);
%     
%     x = xfun(end,:)';
   
    Rt = (k0*exp(-ER/estados(3)));

    x(1) = estados(1) + Ts*((estados(4)-entradas(1))/Ac);
    x(2) = estados(2) + Ts*(estados(4)*((entradas(2)-estados(2))/V) - Rt*estados(2));
    x(3) = estados(3) + Ts*(estados(4)*((estados(5)-estados(3))/V) - (dH/p/cp)*Rt*estados(2) - entradas(3)/V);
   
end

%% Função ekf MIMO

function[x_estimado,P_prox] = ekfMIMO(x_sistema,y_sistema,y_t,u_t,x_a_estimado,P_estimado,Aa,Ca,Q,R)

    %Estrutura do filtro de Kalman, perturbação apenas nos estados x1, x2 e x3
    
    G = zeros(size(Aa,1),size(Aa,1));
    G(end-1,end-1) = 1;
    G(end,end)     = 1;
    
    %Substitui os valores estimados (estados, perturbação e entrada) no sistema real
    x_hat_a_atual = full(x_sistema(x_a_estimado(1),x_a_estimado(2),x_a_estimado(3),x_a_estimado(4),x_a_estimado(5),...
                                    u_t(1),u_t(2),u_t(3)));
                                
    %Calcula nova matriz P
    P_atual = Aa*P_estimado*Aa' + G*Q*G';
   
    %Saída estimada
    z_hat_atual = full(y_sistema(x_hat_a_atual(1),x_hat_a_atual(2),x_hat_a_atual(3),x_hat_a_atual(4),x_hat_a_atual(5),...
                                 u_t(1),u_t(2),u_t(3)));
                             
    %Calcula matriz de ganhos L do filtro
    L = P_atual*Ca'/((Ca*P_atual*Ca' + R));
    
    %Atualização dos estados = EstadosEstimados + GanhoFiltro*(SaidaReal-SaidaEstimada)
    x_estimado = x_hat_a_atual + L*(y_t-z_hat_atual);
    
    %Atualização do erro de covariância
    P_prox = (eye(length(L))-L*Ca)*P_atual;
    
    %Matriz de inovação
%     Mx = Aa*L;
end