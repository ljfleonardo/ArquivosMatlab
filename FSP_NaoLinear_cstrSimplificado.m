% clear; 
% close all; 
clc
s = tf('s');

import casadi.*
x1 = SX.sym('x1');
x2 = SX.sym('x2');

u1 = SX.sym('u1');
u2 = SX.sym('u2');

qi = SX.sym('qi');
Ti = SX.sym('Ti');

d = 10;                               %Atraso real
d_nom = d;                          %Atraso nominal das entradas

delay_real     = [d d];               %Vetor de atraso real
delay_modelado = delay_real;          %Vetor de atraso modelado
d_max          = max(delay_real);     %Máximo valor do atraso (caso sejam diferentes)
dmodelado_max  = max(delay_modelado); %Máximo valor do atraso (caso sejam diferentes)
delay_total    = d_max+d_nom; %Soma dos maiores atrasos para predição

ekf            = 0;                   %Variável usada somente no plot
%% --------------- CSTR MIMO instável ---------------
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
%Entradas / Variáveis Manipuladas
% u1 = q0     - Vazão de saída;
% u2 = Caf    - Concentração do produto A na alimentação do tanque; 
% u3 = Qh/pCp - taxa de remoção de calor normalizada;
% qo0  = 5e-3;
Caf0 = 5;
Qh0  = 0.75;

u0 = [Caf0; Qh0];

%Perturbações
% qi = qi - Vazão de entrada; 
% Ti = Ti - Temperatura externa;
qi0 = 5e-3;
Ti0 = 350;

q0  = [qi0; Ti0];

% Saídas
%Estados / Variáveis Controladas
% x1 = h - nível dentro do tanque; 
% x2 = Ca - Concentração de saída do produto A;
% x3 = T  - Temperatura dentro do reator; 
% h0  = 1;
Ca0 = 1;
T0 = 400;

x0  = [Ca0; T0];

%---- Modelo ----
Rt = (k0*exp(-ER/x2));

dx1_discreto = x1 + Ts*(qi*((u1-x1)/V) - Rt*x1);
dx2_discreto = x2 + Ts*(qi*((Ti-x2)/V) - (dH/p/cp)*Rt*x1 - u2/V);

fun_ax_ext = [dx1_discreto;dx2_discreto;qi;Ti];     %Sistema aumentado discreto
fun_x      = Function('fun_x',{x1,x2,qi,Ti,u1,u2},{fun_ax_ext});

fun_yx_ext = [x1;x2];
fun_y      = Function('fun_y',{x1,x2,qi,Ti,u1,u2},{fun_yx_ext});

na = size(fun_ax_ext,1);                      %Dimensão vetor de estados aumentado
m  = size(fun_yx_ext,1);                      %Dimensão vetor de saídas
%% --------------- Inicialização das variáveis --------------------
iteracoes   = 400;

%% ---- Saídas ----
saidas       = zeros(2,iteracoes);          %inicia vetor de saídas
saidas(1,:)  = x0(1);
saidas(2,:)  = x0(2);
% x(3,:)  = x0(3);

%---- Entradas ----
entradas      = zeros(2,iteracoes);    %inicia vetor de entradas
entradas(1,:) = u0(1);
entradas(2,:) = u0(2);
% entradas(3,:) = u0(3); 

%---- Perturbações ----
perturbacoes      = zeros(2,iteracoes); %inicia vetor de perturbações
perturbacoes(1,:) = q0(1);                   
perturbacoes(2,:) = q0(2);

%---- Inicilização das referências de controle ----
% ref_h            = zeros(1,iteracoes);
% ref_h(1:end)     = x(1,1);

ref_ca           = zeros(1,iteracoes);
ref_ca(1:end)    = saidas(1,1);

ref_T            = zeros(1,iteracoes);
ref_T(1:end)     = saidas(2,1);

%---- Variáveis De Estimação -----
e           = zeros(na-2,iteracoes); %Erro de estimação
ef          = zeros(na-2,iteracoes); %Erro de estimação filtrado
x_n         = zeros(na-2,iteracoes); %Saídas nominais
x_n(1,:) = x0(1);
x_n(2,:) = x0(2);
% x_n(3,:) = x0(3);

x_n(:,1)    = saidas(:,1);                %Saidas nominais anteriores
x_pred_vect = saidas;                     %Vetor de saídas preditas usado para plot

%---- Funções CasADi ----
A_jacobian = jacobian(fun_ax_ext,[x1 x2 qi Ti]);  %Cálculo do jacobiano para matriz A
B_jacobian = jacobian(fun_ax_ext,[u1 u2]);        %Cálculo do jacobiano para matriz B
C_jacobian = jacobian(fun_yx_ext,[x1 x2 qi Ti]);  %Cálculo do jacobiano para matriz C
fun_a = Function('fun_a',{x1,x2,qi,Ti,u1,u2},{A_jacobian});
fun_b = Function('fun_b',{x1,x2,qi,Ti,u1,u2},{B_jacobian});
fun_c = Function('fun_a',{x1,x2,qi,Ti,u1,u2},{C_jacobian});

%% --------------- Projeto de Controle ---------------
for lqr=1
%---- Linearização no ponto de operação ----
Aa = full(fun_a(saidas(1,1),saidas(2,1),perturbacoes(1,1),perturbacoes(2,1),entradas(1,1),entradas(2,1)));
Ba = full(fun_b(saidas(1,1),saidas(2,1),perturbacoes(1,1),perturbacoes(2,1),entradas(1,1),entradas(2,1)));
Ca = full(fun_c(saidas(1,1),saidas(2,1),perturbacoes(1,1),perturbacoes(2,1),entradas(1,1),entradas(2,1)));
Aa_x = Aa(1:na-2,1:2);
Ba_x = Ba(1:na-2,1:2);
Ca_x = Ca(1:na-2,1:2);
Da_x = zeros(size(Ca_x,1),size(Ba_x,2));

%---- Insere dinâmica integradora no controlador ----
A_exp = [Aa_x zeros(size(Aa_x,1),size(Ca_x,1));
         -Ca_x         [1 0; 0 1]];
B_exp = [Ba_x;zeros(2,2)];
C_exp = [Ca_x zeros(2,2)];
D_exp = 0;
sistema = ss(A_exp,B_exp,C_exp,D_exp,Ts);                            %Sistema expandido

Co = ctrb(sistema.A,sistema.B);                                      %Matriz de Controlabilidade
rank(Co);
Q_lqr = diag([1,1/400^2,1,1/400^2]);                             %Matrizes de ponderação
R_lqr = 100*diag([1/5^2,1/0.75^2]);
[K,P,poles] = dlqr(sistema.A,sistema.B,Q_lqr,R_lqr);                 %Ganho do controlador via dlqr

%---- Inicializa integrador no ponto de operação ---- 
integrador_anterior = -pinv(K(:,3:4))*(entradas(:,1)+K(:,1:2)*saidas);

integrador_anterior_1 = integrador_anterior(1);
integrador_atual_1    = 0;

% integrador_anterior_2 = (entradas(:,1)+K(1:3,1:3)*x_estimado)/(-K(1:3,5));
integrador_anterior_2 = integrador_anterior(2);
integrador_atual_2    = 0;
% end
end

%Variável para fechar a malha de controle; 0 = malha aberta; 1 = malha fechada
controle = 1;

%% --------------- Definição das referências ---------------
if controle == 1 % --- MALHA FECHADA ---
    
    ref_ca((round(50/Ts)):end)               = saidas(1,1)*1.6;
%     ref_T((round(100/Ts)):end)               = saidas(1,1)*1.1;
%     ref_ca((round(100/Ts)):end)                = saidas(1,1)*1.02;
%     perturbacoes(1,(round(400/Ts)):end)       = q0(1)*1.02;       %200s
%     perturbacoes(2,(round(700/Ts)):end)       = q0(2)*1.02;       %400s
    perturbacoes(1,(round(450/Ts)):end)       = q0(1)*1.2;       %200s
    perturbacoes(2,(round(750/Ts)):end)       = q0(2)*1.1;       %400s
 

else % --- MALHA ABERTA ---
    entradas(1,(round(100/Ts)):end)           = u0(1)*1.1;             
%     entradas(2,(round(250/Ts)):end)           = u0(2)*1.02;       %100s
    perturbacoes(1,(round(400/Ts)):end)       = q0(1)*1.02;       %200s
    perturbacoes(2,(round(700/Ts)):end)      = q0(2)*1.02;       %400s

end
%% --------------- Projeto do Filtro ---------------
%Passa baixas - maior robustez
%Passa altas  - melhor rejeição de perturbação
%Variável para escolha do filtro; 0 = filtro unitário; 1 = primeira ordem; 2 = segunda ordem;
f = 1;
for filtro=1
    z = tf('z',Ts);
    if f == 2
        af          = 0.06895;               %zero a - calculado para cancelar polo lento do sistema no ponto de operação
        bf          = 0.6873;                %polo b - deve ser posicionado (projetado para zero lento de G22)
        Kf          = 1/((1-af)/((1-bf)^2)); %K para garantir ganho unit
        Fd          = (Kf*((z-af)*z))/(z-bf)^2;
        [num1,den1] = ss2tf(Aa_x,Ba_x,Ca_x,Da_x,1);
        [num2,den2] = ss2tf(Aa_x,Ba_x,Ca_x,Da_x,2);
        G11         = tf(num1(1,:),den1);
        G12         = tf(num1(2,:),den1);
        G21         = tf(num2(1,:),den2);
        G22         = tf(num2(2,:),den2);
        G           = [G11 G12;G21 G22];
        Gd          = zpk(c2d(G,Ts));        %polo mais lento (mais próximo do eixo imag) deve ser polo do filtro 'af'
    else if f == 1
            lambda      = 0.95;
            Fd          = (z*(1-lambda))/((z-lambda));
        else if f == 0                 %Sem filtro / Filtro unitário
                af = 0;
                bf = 0;
                Kf = 1/((1-af)/((1-bf)^2));
                Fd = (Kf*((z-af)*z))/(z-bf)^2;
            end
        end
    end
    den_f       = Fd.den{1,1};
    num_f       = Fd.num{1,1};
    syms z
    
    den_s = poly2sym(den_f,z);
    num_s = poly2sym(num_f,z);
    filtro_s = num_s/den_s;
    filtro_sym = vpa(filtro_s);
    tf_title = latex(filtro_sym);
    n           = length(den_f(2:end));      %Ordem do filtro
end
%% --------------- Simulação --------------------
tic
for k = 2+delay_total:iteracoes
    %% ---- Simulação do processo ----
    %Compila as entradas com atraso real
    entrada_atrasada_real   = [entradas(1,k-1-delay_real(1)) entradas(2,k-1-delay_real(2))];
    
%     saidas(:,k) = modeloCSTR_naoIsotermico(saidas(:,k-1)', entrada_atrasada_real, perturbacoes(:,k-1)');
        tspan = [0 Ts];
    [t, xfun] = ode45(@(t,x) odefcn(t, x, entrada_atrasada_real), tspan, [saidas(:,k-1)' perturbacoes(:,k-1)']);
%     
    saidas(:,k) = xfun(end,1:2)';
    
    %% ---- Predição ----
    %Compila as entradas com atraso pequeno
    entrada_sem_atraso = [entradas(1,k-1) entradas(2,k-1)];
    
    %---- Simulação do sistema sem atraso (atraso min) ----
    %estados nominais, entrada sem atraso e perturbação no ponto de operação
    x_n(:,k) = modeloCSTR_naoIsotermico(x_n(:,k-1), entrada_sem_atraso, q0);
    
    %---- Filtro ----
    e(:,k)   = saidas(:,k) - x_n(:,k-d_nom);                                      %Erro de predição = saida predita - saída real
    ef(:,k)  = -den_f(2:end)*ef(:,k-1:-1:k-n)' + num_f*e(:,k:-1:k-n)'; %Filtragem do sinal do erro
    %---- Saída do preditor ----
    x_pred   = ef(:,k) + x_n(:,k);    
    
    %% ---- Controle ----
    if controle == 1
        %Integrador para concentração do produto A - Ca
        erro_1                = ref_ca(k) - x_pred(1);
        integrador_atual_1    = integrador_anterior_1 + erro_1;
        integrador_anterior_1 = integrador_atual_1;
            
        %Integrador para temperatura interna - T
        erro_2                = ref_T(k) - x_pred(2);
        integrador_atual_2    = integrador_anterior_2 + erro_2;
        integrador_anterior_2 = integrador_atual_2;
              
        nova_entrada  = -K*[x_pred(1);x_pred(2);integrador_atual_1;integrador_atual_2];
        entradas(:,k) = nova_entrada;
    end
    
    %---- Atualização das variáveis ----
    x_pred_vect(:,k)     = x_pred';              %Vetor de plot;
    
end

Elapsed_time = toc;
text1 = ['Tempo de simulação: ',num2str(Elapsed_time),' s'];
disp(text1)
%% --------------- Cálculo de erro --------------------
err_MSE_1_fsp = MSE(ref_ca,saidas(1,:),iteracoes);
err_MSE_2_fsp = MSE(ref_T,saidas(2,:),iteracoes);

err_RMSE_1_fsp = RMSE(ref_ca,saidas(1,:),iteracoes);
err_RMSE_2_fsp = RMSE(ref_T,saidas(2,:),iteracoes);

text2 = ['MSE FSNPNL C_a: ',num2str(err_MSE_1_fsp(end)),newline,'MSE FSNPNL T: ',num2str(err_MSE_2_fsp(end))];
disp(text2)

% err_MAPE_1_fsp = MAPE(ref_ca,saidas(1,:),iteracoes);
% err_MAPE_2_fsp = MAPE(ref_T,saidas(2,:),iteracoes);

%% --------------- Plot de Gráficos --------------------
plot_cstr
% plot_erros
%% --------------- Função do modelo ---------------
function dxdt = odefcn(t, estados, entradas)
    global dH p cp ER k0 V;
    
    Rt = (k0*exp(-ER/estados(2)));

    dxdt = zeros(4,1);
    dxdt(1) = (estados(3)*((entradas(1)-estados(1))/V) - Rt*estados(1));
    dxdt(2) = (estados(3)*((estados(4)-estados(2))/V) - (dH/p/cp)*Rt*estados(1) - entradas(2)/V);
%     dxdt(3) = (estados(4)*((estados(5)-estados(3))/V) - (dH/p/cp)*Rt*estados(2) - entradas(3)/V);
    dxdt(3) = 0;
    dxdt(4) = 0;
end
function x = modeloCSTR_naoIsotermico(estados,entradas,pert) 
    global Ts dH p cp ER k0 V;
    
    Rt = (k0*exp(-ER/estados(2)));
    
    x(1) = estados(1) + Ts*(pert(1)*((entradas(1)-estados(1))/V) - Rt*estados(1));
    x(2) = estados(2) + Ts*(pert(1)*((pert(2)-estados(2))/V) - (dH/p/cp)*Rt*estados(1) - entradas(2)/V);
    
end
%% --------------- Função erro MSE ---------------
function erro = MSE(referencia,saida,N) 
    erro = zeros(1,N);
    erro(1,1)=((referencia(1) - saida(1))^2);
    for i=2:N
        aux     = (referencia(i) - saida(i))^2;
        erro(i) = ((erro(i-1) + aux)); %integral do erro = anterior + atual
%         erro(i) = (predito(i)-real(i))^2;
    end    
%     erro = erro/N;
end
%% --------------- Função erro RMSE ---------------
function erro = RMSE(real,predito,N) %erro na mesma unidade da variável
    erro = zeros(1,N);
    erro(1)=((real(1) - predito(1))^2);
    for i=2:N
        aux     = (real(i) - predito(i))^2;
        erro(i) = ((erro(i-1) + aux)); %integral do erro = anterior + atual
    end    
    erro = sqrt(erro/N);
end
%% --------------- Função erro MAPE ---------------
function erro = MAPE(real,predito,N) %erro em porcentagem
    erro = zeros(1,N);
    erro(1)=((real(1) - predito(1)))/real(1);
    for i=2:N
        aux     = ((real(i) - predito(i)))/real(i);
        erro(i) = ((erro(i-1) + aux)); %integral do erro = anterior + atual
    end    
    erro = (erro/N);
end