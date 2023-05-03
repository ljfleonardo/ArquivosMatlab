clear; %close all; clc
clc
import casadi.*
x1 = SX.sym('x1');
x2 = SX.sym('x2');

u1 = SX.sym('u1');
u2 = SX.sym('u2');

q_1 = SX.sym('q_1');
q_2 = SX.sym('q_2');

d = 10;                               %Atraso real
delay_real     = [d d];               %Vetor de atraso real
delay_modelado = delay_real;          %Vetor de atraso modelado
d_max          = max(delay_real);     %M�ximo valor do atraso (caso sejam diferentes)
dmodelado_max  = max(delay_modelado); %M�ximo valor do atraso (caso sejam diferentes)
delay_total    = d_max+dmodelado_max; %Soma dos maiores atrasos para predi��o

ekf            = 1;                   %Vari�vel usada somente no plot
%% --------------- CSTR MIMO inst�vel ---------------
global Ts Ac dH p cp ER k0 V q0
%---- Constantes ----
Ac = 0.05;        %Area do reator
dH = -50e6;       %Calor da rea��o
p  = 1e6;         %Densidade dos l�quidos
cp = 1;           %Calor espec�fico
ER = 1000;        %Termo de energia de ativa��o
k0 = 4.872997584; %Constante da taxa de rea��o
V  = 0.05;        %Volume do reator
Ts = 3;           %Per�odo de amostragem

%---- Ponto de Opera��o ----
%Entradas / Vari�veis Manipuladas
% u1 = Caf    - Concentra��o do produto A na alimenta��o do tanque; 
% u2 = Qh/pCp - taxa de remo��o de calor normalizada;
Caf0 = 5;
Qh0  = 0.75;

u0 = [Caf0; Qh0];

%Perturba��es
% qi = qi - Vaz�o de entrada; 
% Ti = Ti - Temperatura externa;
qi0 = 5e-3;
Ti0 = 350;

q0  = [qi0; Ti0];

% Sa�das
%Estados / Vari�veis Controladas
% x1 = Ca - Concentra��o de sa�da do produto A;
% x2 = T  - Temperatura dentro do reator; 
Ca0 = 1;
T0  = 400;

x0  = [Ca0; T0];

%---- Modelo ----
Rt = (k0*exp(-ER/x2));

u_q1 = u1+q_1;
u_q2 = u2+q_2;
dx1_discreto = x1 + Ts*(q0(1)*((u_q1-x1)/V) - Rt*x1);
dx2_discreto = x2 + Ts*(q0(1)*((q0(2)-x2)/V) - (dH/p/cp)*Rt*x1 - u_q2/V);

fun_ax_ext = [dx1_discreto;dx2_discreto;q_1;q_2];     %Sistema aumentado discreto
fun_x      = Function('fun_x',{x1,x2,q_1,q_2,u1,u2},{fun_ax_ext});

fun_yx_ext = [x1;x2];
fun_y      = Function('fun_y',{x1,x2,q_1,q_2,u1,u2},{fun_yx_ext});

na = size(fun_ax_ext,1);                      %Dimens�o vetor de estados aumentado
m  = size(fun_yx_ext,1);                      %Dimens�o vetor de sa�das
%% --------------- Inicializa��o das vari�veis --------------------
iteracoes   = 400;
tspan       = [0 Ts];

%---- Sa�das ----
saidas       = zeros(2,iteracoes);            %inicia vetor de sa�das

saidas(1,:)  = x0(1); %1   h0
saidas(2,:)  = x0(2); %1   Ca0

%---- Entradas ----
entradas      = zeros(2,iteracoes);           %inicia vetor de entradas

entradas(1,:) = u0(1); %0.005  q0
entradas(2,:) = u0(2); %5      Caf0

%---- Perturba��es ----
perturbacoes      = zeros(2,iteracoes);       %inicia vetor de perturba��es

perturbacoes(1,:) = q0(1);%0.005  qi0                    
perturbacoes(2,:) = q0(2);%350    Ti0 

%---- Iniciliza��o Das Refer�ncias De Controle ----
ref_ca           = zeros(1,iteracoes);
ref_ca(1:end)    = saidas(1,1);

ref_T            = zeros(1,iteracoes);
ref_T(1:end)     = saidas(2,1);

%---- Vari�veis De Estima��o -----
x_estimado  = [saidas(1,1);saidas(2,1)];     %Cria vetor de estados estimados

q_estimado  = [0;0];                         %Cria vetor de perturba��es estimadas

x_a_estim   = [x_estimado;q_estimado];       %Compila as estima��es iniciais de estados e perturba��es

x_estimado_vect      = zeros(na,iteracoes);  %Vetor para plot
x_estimado_vect(1,:) = x_estimado(1);
x_estimado_vect(2,:) = x_estimado(2);
x_estimado_vect(3,:) = q_estimado(1);
x_estimado_vect(4,:) = q_estimado(2);

x_a_pred              = zeros(na,iteracoes); 
x_pred_vect           = zeros(na,iteracoes); %Vetor para plot
x_pred_vect(1,:)      = x_estimado(1);
x_pred_vect(2,:)      = x_estimado(2);
x_pred_vect(3,:)      = q_estimado(1);
x_pred_vect(4,:)      = q_estimado(2);

%---- Par�metros do Estimador ----   
%qi - Vaz�o de entrada(0.001~0.007); %Ti - Temperatura externa(350~370);
qi_range=[0.001 0.007];      Ti_range=[350 370];
qi_med=mean(qi_range);       Ti_med=mean(Ti_range);
qi_dp=abs(qi_med-0.007)/2;   Ti_dp=abs(Ti_med-370)/2;
   
P_estimado_at = diag([1*ones(1,na-2),1,1]); 

%Vari�vel para ajuste da sintonia; 0 = filtro unit�rio; 1 = melhor rejei��o de perturba��o; 2 = maior robustez
filtro = 2;
for f = 1
    if filtro==0
        text1 = 'Filtro unit�rio';
        disp(text1)
        Q = eye(2,2);
        R = eye(2,2);
    else if filtro == 1
            %Sintonia para melhorar rejei��o de perturba��o
            text1 = 'Filtro com melhor rejei��o de perturba��o';
            disp(text1)
            Q = diag([(qi_dp*1e1)^2, (Ti_dp*1e1)^2]);
            R = diag([((Ca_dp*1)^2), ((T_dp*1)^2)]);
        else if filtro == 2
                text1 = 'Filtro com maior robustez';
                disp(text1)
                %Sintonia boa para erro de modelagem 10+5 (real+erro)
%                 Q = diag([(qi_dp*1e-3)^2, (Ti_dp*1e-5)^2]);
%                 R = diag([((Ca_dp*1e2)^2), ((T_dp*1e2)^2)]);
                Q = diag([(qi_dp*1e-2)^2, (Ti_dp*1e-2)^2]);
                R = diag([((Ca_dp*1e0)^2), ((T_dp*1e0)^2)]);
            end
        end
    end
end

%---- Fun��es CasADi ----
A_jacobian = jacobian(fun_ax_ext,[x1 x2 q_1 q_2]);  %C�lculo do jacobiano para matriz A
B_jacobian = jacobian(fun_ax_ext,[u1 u2]);          %C�lculo do jacobiano para matriz B
C_jacobian = jacobian(fun_yx_ext,[x1 x2 q_1 q_2]);  %C�lculo do jacobiano para matriz C
fun_a = Function('fun_a',{x1,x2,q_1,q_2,u1,u2},{A_jacobian});
fun_b = Function('fun_b',{x1,x2,q_1,q_2,u1,u2},{B_jacobian});
fun_c = Function('fun_a',{x1,x2,q_1,q_2,u1,u2},{C_jacobian});

%% --------------- Projeto de Controle ---------------
for lqr=1
%---- Lineariza��o no ponto de opera��o ----
Aa = full(fun_a(saidas(1,1),saidas(2,1),0,0,entradas(1,1),entradas(2,1)));
Ba = full(fun_b(saidas(1,1),saidas(2,1),0,0,entradas(1,1),entradas(2,1)));
Ca = full(fun_c(saidas(1,1),saidas(2,1),0,0,entradas(1,1),entradas(2,1)));

Ob = obsv(Aa,Ca);
rank(Ob);
if rank(Ob)==size(Ob,2)
    text1 = ['� observ�vel, posto ',num2str(rank(Ob))];
    disp(text1)
else
    text1 = ['N�O � observ�vel, posto ',num2str(rank(Ob)),' e Ob tem tamanho ',num2str(size(Ob,1))];
    disp(text1)
end

Aa_x = Aa(1:na-2,1:2);
Ba_x = Ba(1:na-2,1:2);
Ca_x = Ca(1:na-2,1:2);

%---- Insere din�mica integradora no controlador ----
A_exp = [Aa_x zeros(size(Aa_x,1),size(Ca_x,1));
         -Ca_x         [1 0; 0 1]];
B_exp = [Ba_x;zeros(2,2)];
C_exp = [Ca_x zeros(2,2)];
D_exp = 0;
sistema = ss(A_exp,B_exp,C_exp,D_exp,Ts);                            %Sistema expandido

Co = ctrb(sistema.A,sistema.B);                                      %Matriz de Controlabilidade
rank(Co);
Q_lqr = diag([1,1/400^2,1,1/400^2]);                             %Matrizes de pondera��o
R_lqr = 100*diag([1/5^2,1/0.75^2]);
[K,P,poles] = dlqr(sistema.A,sistema.B,Q_lqr,R_lqr);                 %Ganho do controlador via dlqr

%---- Inicializa integrador no ponto de opera��o ---- 
integrador_anterior = -pinv(K(:,3:4))*(entradas(:,1)+K(:,1:2)*saidas);

integrador_anterior_1 = integrador_anterior(1);
integrador_atual_1    = 0;

integrador_anterior_2 = integrador_anterior(2);
integrador_atual_2    = 0;
% end
end
%Vari�vel para fechar a malha de controle; 0 = malha aberta; 1 = malha fechada
controle = 1;

%% --------------- Defini��o das refer�ncias ---------------
if controle == 1 % --- MALHA FECHADA ---
    text1 = 'Sistema em malha fechada';
    disp(text1)
    ref_ca((round(50/Ts)):end)                = saidas(1,1)*1.6;
%     ref_T((round(100/Ts)):end)                = saidas(1,1)*1.1;
    perturbacoes(1,(round(400/Ts)):end)       = q0(1)*1.02;       %200s
    perturbacoes(2,(round(700/Ts)):end)       = q0(2)*1.02;       %400s

else % --- MALHA ABERTA ---
    text1 = 'Sistema em malha aberta';
    disp(text1)
    entradas(1,(round(100/Ts)):end)           = u0(1)*1.1;             
%     entradas(2,(round(250/Ts)):end)           = u0(2)*1.02;       %100s
    perturbacoes(1,(round(400/Ts)):end)       = q0(1)*1.02;       %200s
    perturbacoes(2,(round(700/Ts)):end)       = q0(2)*1.02;       %400s

end
%% --------------- Simula��o --------------------
tic
for k = 2+delay_total:iteracoes

    %Compila as entradas com atraso real
    entrada_atrasada_real = [entradas(1,k-1-delay_real(1)) entradas(2,k-1-delay_real(2))]';
    
    %Compila as entradas com atraso modelado
    entrada_atrasada_mod  = [entradas(1,k-1-delay_modelado(1)) entradas(2,k-1-delay_modelado(2))]';
    
    %---- Simula��o do processo ----
    [t, xfun]   = ode45(@(t,x) odeCSTR_nominal(t, x, entrada_atrasada_real, perturbacoes(:,k-1)),...
                                               tspan, saidas(:,k-1));
    saidas(:,k) = xfun(end,1:2)';

    %---- Lineariza��o a cada itera��o para EKF ----
    Aa = full(fun_a(x_a_estim(1),x_a_estim(2),x_a_estim(3),x_a_estim(4),...
                    entrada_atrasada_mod(1),entrada_atrasada_mod(2)));

    Ca = full(fun_c(x_a_estim(1),x_a_estim(2),x_a_estim(3),x_a_estim(4),...
                    entrada_atrasada_mod(1),entrada_atrasada_mod(2)));

    %---- Estima��o ---- 
    [x_a_estim,P_estimado] = ekfMIMO(fun_x,fun_y,saidas(:,k),...
                                     entrada_atrasada_mod,x_a_estim,P_estimado_at,Aa,Ca,Q,R);
    %---- Predi��o ----
    %Primeira itera��o
    x_a_pred = x_a_estim; 
    if(dmodelado_max>0)
       %Segunda itera��o em diante:
        for kk=1:dmodelado_max
            %Tratamento das entradas
            if(k-1-delay_modelado(1)+kk>k) 
                entradas(1,k-1-delay_modelado(1)+kk) = entradas(1,k-1-delay_modelado(1)+kk-1);
            end
            if(k-1-delay_modelado(2)+kk>k)
                entradas(2,k-1-delay_modelado(2)+kk) = entradas(2,k-1-delay_modelado(2)+kk-1);
            end
            x_a_pred = modeloCSTR_naoIsotermico(x_a_pred,...
                                      [entradas(1,k-1-delay_modelado(1)+kk) entradas(2,k-1-delay_modelado(2)+kk)]);
        end
    end
        
    %---- Controle ----
    if controle == 1
        %Integrador para concentra��o do produto A - Ca
        erro_1                = ref_ca(k) - x_a_pred(1);
        integrador_atual_1    = integrador_anterior_1 + erro_1;
        integrador_anterior_1 = integrador_atual_1;
                
        %Integrador para temperatura interna - T
        erro_2                = ref_T(k) - x_a_pred(2);
        integrador_atual_2    = integrador_anterior_2 + erro_2;
        integrador_anterior_2 = integrador_atual_2;
        
        %Lei de Controle
        nova_entrada  = -K*[x_a_pred(1);x_a_pred(2);integrador_atual_1;integrador_atual_2];
        entradas(:,k) = nova_entrada;
    end
    
    %---- Atualiza��o das vari�veis ----
    x_estimado_vect(:,k) = x_a_estim(1:end);
    x_pred_vect(:,k)     = x_a_pred';
    P_estimado_at        = P_estimado;
    
end

Elapsed_time = toc;
text1 = ['Tempo de simula��o: ',num2str(Elapsed_time),' s'];
disp(text1)
%% --------------- C�lculo de erro --------------------

err_MSE_1_ekf = MSE(ref_ca,saidas(1,:),iteracoes);
err_MSE_2_ekf = MSE(ref_T,saidas(2,:),iteracoes);

text2 = ['MSE EKF C_a: ',num2str(err_MSE_1_ekf(end)),newline,'MSE EKF T: ',num2str(err_MSE_2_ekf(end)) ];
disp(text2)

err_RMSE_1_ekf = RMSE(ref_ca,saidas(1,:),iteracoes);
err_RMSE_2_ekf = RMSE(ref_T,saidas(2,:),iteracoes);

%% --------------- Plot de Gr�ficos --------------------
plot_cstr
% plot_erros
%% --------------- Fun��o do modelo ---------------
function x = modeloCSTR_naoIsotermico(estados,entradas) 
    global Ts dH p cp ER k0 V q0;
    
    Rt         = (k0*exp(-ER/estados(2)));
    
    entrada_q1 = entradas(1)+estados(3);
    entrada_q2 = entradas(2)+estados(4);
    x(1)       = estados(1) + Ts*(q0(1)*((entrada_q1-estados(1))/V) - Rt*estados(1));
    x(2)       = estados(2) + Ts*(q0(1)*((q0(2)-estados(2))/V) - (dH/p/cp)*Rt*estados(1) - entrada_q2/V);
    x(3)       = estados(3);
    x(4)       = estados(4);
        
end

%% --------------- Fun��o processo nominal ---------------
function dxdt = odeCSTR_nominal(t, estados, entradas,perturbacoes)
    global dH p cp ER k0 V;
    
    Rt = (k0*exp(-ER/estados(2)));

    dxdt    = zeros(2,1);
    dxdt(1) = (perturbacoes(1)*((entradas(1)-estados(1))/V) - Rt*estados(1));
    dxdt(2) = (perturbacoes(1)*((perturbacoes(2)-estados(2))/V) - (dH/p/cp)*Rt*estados(1) - entradas(2)/V);

end

%% --------------- Fun��o EKF ---------------
function[x_estimado,P_prox] = ekfMIMO(x_sistema,y_sistema,y_t,u_t,x_a_estimado,P_estimado,Aa,Ca,Q,R)

    %Estrutura do filtro de Kalman, perturba��o apenas nos estados x1, x2 e x3
    
    G = zeros(size(Aa,1),2);
    G(end-1,end-1) = 1;
    G(end,end)     = 1;
    
    %Substitui os valores estimados (estados, perturba��o e entrada) no sistema real
    x_hat_a_atual = full(x_sistema(x_a_estimado(1),x_a_estimado(2),x_a_estimado(3),x_a_estimado(4),...
                                    u_t(1),u_t(2)));
                                
    %Calcula nova matriz P
    P_atual = Aa*P_estimado*Aa' + G*Q*G';
   
    %Sa�da estimada
    z_hat_atual = full(y_sistema(x_hat_a_atual(1),x_hat_a_atual(2),x_hat_a_atual(3),x_hat_a_atual(4),...
                                 u_t(1),u_t(2)));
                             
    %Calcula matriz de ganhos L do filtro
    L = P_atual*Ca'/((Ca*P_atual*Ca' + R));
    
    %Atualiza��o dos estados = EstadosEstimados + GanhoFiltro*(SaidaReal-SaidaEstimada)
    x_estimado = x_hat_a_atual + L*(y_t-z_hat_atual);
    
    %Atualiza��o do erro de covari�ncia
    P_prox = (eye(length(L))-L*Ca)*P_atual;
    
    %Matriz de inova��o
%     Mx = Aa*L;
end
%% --------------- Fun��o erro MSE ---------------
function erro = MSE(real,predito,N)   %Integral do erro
    erro = zeros(1,N);
    erro(1)=((real(1) - predito(1))^2);
    for i=2:N
        aux     = (real(i) - predito(i))^2;
        erro(i) = ((erro(i-1) + aux)); %integral do erro = anterior + atual
%         erro(i) = (predito(i)-real(i))^2;
    end    
%     erro = erro/N;
end
%% --------------- Fun��o erro RMSE ---------------
function erro = RMSE(real,predito,N)  %Erro na mesma unidade da vari�vel
    erro = zeros(1,N);
    erro(1)=((real(1) - predito(1))^2);
    for i=2:N
        aux     = (real(i) - predito(i))^2;
        erro(i) = ((erro(i-1) + aux)); %integral do erro = anterior + atual
    end    
    erro = sqrt(erro/N);
end
%% --------------- Fun��o erro MAPE ---------------
function erro = MAPE(real,predito,N)  %Erro em porcentagem
    erro = zeros(1,N);
    erro(1)=((real(1) - predito(1)))/real(1);
    for i=2:N
        aux     = ((real(i) - predito(i)))/real(i);
        erro(i) = ((erro(i-1) + aux)); %integral do erro = anterior + atual
    end    
    erro = (erro/N);
end