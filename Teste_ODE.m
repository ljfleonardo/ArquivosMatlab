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

delay_real     = [0 0 0];
delay_modelado = delay_real;
% delay_modelado = [0 1 0];
d_max          = max(delay_real);
dmodelado_max  = max(delay_modelado);
delay_total    = d_max+dmodelado_max;
%% CSTR MIMO instável
global Ts Ac dH p cp ER k0 V
%---- Constantes ----
Ac = 0.05;        %Area do reator
dH = -50e6;       %Calor da reação
p  = 1e6;          %Densidade dos líquidos
cp = 1;           %Calor específico
ER = 1000;        %Termo de energia de ativação
k0 = 4.872997584; %Constante da taxa de reação
V  = 0.05;         %Volume do reator

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
% Ca0 = 1;
Ca0 = 1.000000000069293;
% T0 = 400;
T0  = 3.999999999965353e+02;

x0  = [h0; Ca0; T0];

%Estados / Variáveis Controladas
% x1 = h - nível dentro do tanque; 
% x2 = Ca - Concentração de saída do produto A;
% x3 = T  - Temperatura dentro do reator; 
% x4 = R  - Velocidade de reação

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

%     tspan = [0 Ts];
%     [t, xfun] = ode45(@(t,x) odefcn(t, x, entradas), tspan, estados);
%     
%     x = xfun(end,:)';

% tspan = [0 Ts];
% x_q_0 = [x0;q0];
% 
% [t,dx]=ode45(@(t,dx) odefcn(t,dx,u0),tspan,x_q_0);

% dx1_discreto = dx(1);

fun_ax_ext = [dx1_discreto;dx2_discreto;dx3_discreto;qi;Ti];     %Sistema aumentado discreto
fun_x      = Function('fun_x',{x1,x2,x3,qi,Ti,u1,u2,u3},{fun_ax_ext});

fun_yx_ext = [x1;x2;x3];
fun_y      = Function('fun_y',{x1,x2,x3,qi,Ti,u1,u2,u3},{fun_yx_ext});

na = size(fun_ax_ext,1);                      %Dimensão vetor de estados aumentado
m  = size(fun_yx_ext,1);                      %Dimensão vetor de saídas
%% --------------- Inicialização das variáveis --------------------
iteracoes   = 900;

%---- Saídas ----
saidas       = zeros(5,iteracoes+delay_total);            %inicia vetor de saídas

saidas(1,:)  = x0(1); %1   h0
saidas(2,:)  = x0(2); %1   Ca0
saidas(3,:)  = x0(3); %400 T0
saidas(4,:)  = q0(1);
saidas(5,:)  = q0(2);

%---- Entradas ----
entradas      = zeros(3,iteracoes+delay_total);              %inicia vetor de entradas

entradas(1,:) = u0(1); %0.005  q0
entradas(2,:) = u0(2); %5      Caf0
entradas(3,:) = u0(3); %0.75   Qh/pcp 0

% entradas(2,400:end) = u0(2)*1.02;             %Degrau de 2%
% entradas(1,200:end) = u0(1)*1.02;             %Degrau de 2%
% entradas(3,700:end) = u0(3)*1.02;             %Degrau de 2%


%---- Perturbações ----

perturbacoes      = zeros(2,iteracoes+delay_total);              %inicia vetor de perturbações

perturbacoes(1,:) = q0(1);%0.005  qi0                    
perturbacoes(2,:) = q0(2);%350    Ti0 

% perturbacoes(1,200:end)      = q0(1)*1.02;           %degrau na perturbação de 2% a partir de 500 it
% perturbacoes(2,800:end)      = q0(2)*1.02;

%---- Referências de controle ----

ref_h            = zeros(1,iteracoes+delay_total);
ref_h(1:end)     = saidas(1,1);

ref_ca           = zeros(1,iteracoes+delay_total);
ref_ca(1:end)    = saidas(2,1);
ref_ca(400:end)  = saidas(2,1)*1.02;           %Degrau de 2%
% ref_ca(700:end)   = saidas(2,1)*1.1;

ref_T            = zeros(1,iteracoes+delay_total);
ref_T(1:end)     = saidas(3,1);
% ref_T(700:end)   = saidas(3,1)*1.02;           %Degrau de 2%

%---- Ruido - ainda não sendo utilizado ----
ruido = zeros(3,iteracoes+delay_total);

%---- Variáveis De Estimação -----
x_estimado  = [saidas(1,1);saidas(2,1);saidas(3,1)];     %Cria vetor de estados estimados

d_estimado  = [perturbacoes(1,1);perturbacoes(2,1)];     %Cria vetor de perturbações estimadas

y_estimado  = [saidas(1,1);saidas(2,1);saidas(3,1)];     %Cria vetor de saídas estimadas

x_a_estim   = [x_estimado;d_estimado];                     %Compila as estimações de estados e perturbações

y_estimado_vect      = y_estimado;
x_estimado_vect      = zeros(na,iteracoes+delay_total);
x_estimado_vect(1,:) = x_estimado(1);
x_estimado_vect(2,:) = x_estimado(2);
x_estimado_vect(3,:) = x_estimado(3);
x_estimado_vect(4,:) = d_estimado(1);
x_estimado_vect(5,:) = d_estimado(2);

x_a_pred              = zeros(na,iteracoes+delay_total);
x_pred_vect           = zeros(na,iteracoes+delay_total);
x_pred_vect(1,:)      = x_estimado(1);
x_pred_vect(2,:)      = x_estimado(2);
x_pred_vect(3,:)      = x_estimado(3);
x_pred_vect(4,:)      = d_estimado(1);
x_pred_vect(5,:)      = d_estimado(2);


%---- Parâmetros do Estimador ----   

Q = diag([1*ones(1,na-2),1,1]);               %Variável da ponderação dos estados
% Q = 0.1*Q;
% Q = diag([1 1]);
R = diag(ones(1,m));                          %Variável da ponderação da saída                    
P_estimado_at = diag([1*ones(1,na-2),1,1]);   

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
Q_lqr = diag([1,1,1,1,1,1]);                                          %Matrizes de ponderação
R_lqr = 10*eye(size(sistema.B,2));
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
%Variável para fechar a malha; 0 = malha aberta; 1 = malha fechada
controle = 1;

%% --------------- Simulação --------------------
tic
for k = 2+delay_total:iteracoes

    xd = [saidas(1:3,k-1)' perturbacoes(:,k-1)'];

    entradas_comp = [entradas(1,k-1-delay_real(1)) entradas(2,k-1-delay_real(2)) entradas(3,k-1-delay_real(3))];
    
%     entradas_ruidosas = entradas_comp + ruido(k);
    
    %---- Simulação do processo ----
%     var_aux = modeloCSTR_naoIsotermico(xd,entradas_comp);
    saidas(:,k) = modeloCSTR_naoIsotermico(xd,entradas_comp);
    
%     saidas(:,k) = var_aux;

    
    %---- Linearização a cada iteração para EKF ----
    
    Aa = full(fun_a(x_a_estim(1),x_a_estim(2),x_a_estim(3),x_a_estim(4),x_a_estim(5),...
                    entradas(1,k-1-delay_modelado(1)),entradas(2,k-1-delay_modelado(2)),entradas(3,k-1-delay_modelado(3))));                %substitui os valores reais em A_jacobian
 
    Ca = full(fun_c(x_a_estim(1),x_a_estim(2),x_a_estim(3),x_a_estim(4),x_a_estim(5),...
                    entradas(1,k-1-delay_modelado(1)),entradas(2,k-1-delay_modelado(2)),entradas(3,k-1-delay_modelado(3)))); 

    %---- Estimação ---- 
    entrada_atrasada  = [entradas(1,k-1-delay_modelado(1)) entradas(2,k-1-delay_modelado(2)) entradas(3,k-1-delay_modelado(3))]';
    
    [x_a_estim,P_estimado] = ekfMIMO(saidas(1:3,k),entrada_atrasada,x_a_estim,P_estimado_at,Aa,Ca,Q,R);

    %---- Predição ----
    %Primeira iteração
    x_a_pred = x_a_estim;   
    if dmodelado_max>0
        %Segunda iteração em diante:
        for kk=1:dmodelado_max
    %com atrasos diferentes, o indice k-delay+kk-1 ultrapassa o número de iterações
    %vetores aumentados como iteracoes+delay_total
            x_a_pred = modeloCSTR_naoIsotermico(x_a_pred,[entradas(1,k-delay_modelado(1)+kk-1) entradas(2,k-delay_modelado(2)+kk-1) entradas(3,k-delay_modelado(3)+kk-1)]');

        end
    end
    
    %---- Controle ----
    if controle == 1
        %Integrador para altura do tanque - h
        erro_1                = ref_h(k) - x_a_pred(1);
%         erro_1                = ref_h(k) - x_a_estim(1);
%         erro_1                = ref_h(k) - saidas(1,k);
        integrador_atual_1    = integrador_anterior_1 + erro_1;
        integrador_anterior_1 = integrador_atual_1;
                
        %Integrador para concentração do produto A - Ca
        erro_2                = ref_ca(k) - x_a_pred(2);
%         erro_2                = ref_ca(k) - x_a_estim(2);
%         erro_2                = ref_ca(k) - saidas(2,k);
        integrador_atual_2    = integrador_anterior_2 + erro_2;
        integrador_anterior_2 = integrador_atual_2;
        
        %Integrador para temperatura interna - T
        erro_3                = ref_T(k) - x_a_pred(3);
%         erro_3                = ref_T(k) - x_a_estim(3);
%         erro_3                = ref_T(k) - saidas(3,k);
        integrador_atual_3    = integrador_anterior_3 + erro_3;
        integrador_anterior_3 = integrador_atual_3;
      
      nova_entrada = -K*[x_a_pred(1);x_a_pred(2);x_a_pred(3); integrador_atual_1; integrador_atual_2; integrador_atual_3];
%       
%       nova_entrada = -K*[x_a_estim(1);x_a_estim(2);x_a_estim(3); integrador_atual_1; integrador_atual_2; integrador_atual_3];
%       nova_entrada = -K*[saidas(1,k);saidas(2,k);saidas(3,k); integrador_atual_1; integrador_atual_2; integrador_atual_3];

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
    global Ts % Ac dH p cp ER k0 V;
    %estados(1) = h; estados(2) = Ca; estados(3) = T ;
    %estados(4) = qi; estados(5) = Ti;
    %entradas(1) = qo; entradas(2) = Caf; entradas(3) = Qh/pCp;
   
    tspan = [0 Ts];
    [t, xfun] = ode45(@(t,x) odefcn(t, x, entradas), tspan, estados);
    
    x = xfun(end,:)';
   
%     Rt = (k0*exp(-ER/estados(3)));
% 
%     x(1) = estados(1) + Ts*((estados(4)-entradas(1))/Ac);
%     x(2) = estados(2) + Ts*(estados(4)*((entradas(2)-estados(2))/V) - Rt*estados(2));
%     x(3) = estados(3) + Ts*(estados(4)*((estados(5)-estados(3))/V) - (dH/p/cp)*Rt*estados(2) - entradas(3)/V);
   
end

%% Função ekf MIMO

function[x_estimado,P_prox] = ekfMIMO(y_t,u_t,x_a_estimado,P_estimado,Aa,Ca,Q,R)

    %Estrutura do filtro de Kalman, perturbação apenas nos estados x1, x2 e x3
    
    G = zeros(size(Aa,1),size(Aa,1));
    G(end-1,end-1) = 1;
    G(end,end)     = 1;
    
    %Substitui os valores estimados (estados, perturbação e entrada) no sistema real
%     x_hat_a_atual = full(x_sistema(x_a_estimado(1),x_a_estimado(2),x_a_estimado(3),x_a_estimado(4),x_a_estimado(5),...
%                                     u_t(1),u_t(2),u_t(3)));
    x_hat_a_atual = modeloCSTR_naoIsotermico([x_a_estimado(1),x_a_estimado(2),x_a_estimado(3),x_a_estimado(4),x_a_estimado(5)],...
                                             [u_t(1),u_t(2),u_t(3)]);
                                
    %Calcula nova matriz P
    P_atual = Aa*P_estimado*Aa' + G*Q*G';
   
    %Saída estimada
%     z_hat_atual = full(y_sistema(x_hat_a_atual(1),x_hat_a_atual(2),x_hat_a_atual(3),x_hat_a_atual(4),x_hat_a_atual(5),...
%                                  u_t(1),u_t(2),u_t(3)));
    z_hat_atual_aux = modeloCSTR_naoIsotermico([x_hat_a_atual(1),x_hat_a_atual(2),x_hat_a_atual(3),x_hat_a_atual(4),x_hat_a_atual(5)],...
                                               [u_t(1),u_t(2),u_t(3)]);
    z_hat_atual = z_hat_atual_aux(1:3);
    %Calcula matriz de ganhos L do filtro
    L = P_atual*Ca'/((Ca*P_atual*Ca' + R));
    
    %Atualização dos estados = EstadosEstimados + GanhoFiltro*(SaidaReal-SaidaEstimada)
    x_estimado = x_hat_a_atual + L*(y_t-z_hat_atual);
    
    %Atualização do erro de covariância
    P_prox = (eye(length(L))-L*Ca)*P_atual;
    
    %Matriz de inovação
%     Mx = Aa*L;
end