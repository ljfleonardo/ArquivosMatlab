%Estados / Variáveis Controladas
% x1 = h - nível dentro do tanque; 
% x2 = Ca - Concentração de saída do produto A;
% x3 = T - Temperatura dentro do reator; 
% x4 = R - Velocidade de reação?
figure
tamLetra = 10;
tamTitulo = 12;
subplot(3,2,1);
plot(saidas(1,1:iteracoes),'r','linewidth',2);
hold on
plot(ref_h(1:iteracoes),'b--','linewidth',2);
plot(x_pred_vect(1,1:iteracoes),'k--','linewidth',2);
legend({'Real','Referência','Predição'},'FontSize',tamLetra);
xlabel('Iterações [-]','FontSize',tamLetra);
title('$h$ - Altura do tanque','interpreter','latex','FontSize',tamTitulo)
grid


subplot(3,2,3);
plot(saidas(2,1:iteracoes),'r','linewidth',2);
hold on
plot(ref_ca(1:iteracoes),'b--','linewidth',2);
plot(x_pred_vect(2,1:iteracoes),'k--','linewidth',2)%,'marker','x');
legend({'Real','Referência','Predição'},'Location','best','FontSize',tamLetra);
xlabel('Iterações [-]','FontSize',tamLetra);
title('$C_a$ Concentra\c{c}\~{a}o do produto A','interpreter','latex');
grid


subplot(3,2,5);
plot(saidas(3,1:iteracoes),'r','linewidth',2);
hold on
plot(ref_T(1:iteracoes),'b--','linewidth',2);
plot(x_pred_vect(3,1:iteracoes),'k--','linewidth',2);
legend({'Real','Referência','Predição'},'Location','best','FontSize',tamLetra);
xlabel('Iterações [-]','FontSize',tamLetra);
title('$T$ Temperatura dentro do tanque','interpreter','latex');
grid


% figure
%----- Entradas / Variáveis Manipuladas -----
% u1 = qo - Vazão de saída;
% u2 = Caf - Concentração do produto A na alimentação do tanque; 
% u3 = Qh/pCp - taxa de remoção de calor normalizada;
subplot(3,2,2);
% plot(entradas_atrasadas_vect(1,1:iteracoes),'k--','linewidth',2);
% hold on
% plot(entradas_comp_vect(1,1:iteracoes),'b','linewidth',2);
plot(entradas(1,1:iteracoes),'r','linewidth',2);
title('Entradas','FontSize',tamTitulo)
legend({'$q_0$ Vaz\~{a}o de sa\''ida'},'interpreter','latex','Location','best','FontSize',tamLetra);
xlabel('Iterações [-]','FontSize',tamLetra);
grid

subplot(3,2,4);
% plot(entradas_atrasadas_vect(2,1:iteracoes),'k--','linewidth',2);
% hold on
% plot(entradas_comp_vect(2,1:iteracoes),'b','linewidth',2);
plot(entradas(2,1:iteracoes),'r','linewidth',2);
legend({'$C_{af}$ Concentra\c{c}\~{a}o produto A na alimenta\c{c}\~{a}o do tanque'},'interpreter','latex','Location','best','FontSize',tamLetra);
xlabel('Iterações [-]','FontSize',tamLetra);
grid

subplot(3,2,6);
% plot(entradas_atrasadas_vect(3,1:iteracoes),'k--','linewidth',2);
% hold on
% plot(entradas_comp_vect(3,1:iteracoes),'b','linewidth',2);
plot(entradas(3,1:iteracoes),'r','linewidth',2);
legend({'${Qh}/{pc_p}$  Taxa de remo\c{c}\~{a}o de calor normalizada'},'interpreter','latex','Location','best','FontSize',tamLetra);
xlabel('Iterações [-]','FontSize',tamLetra);
grid

figure
%----- Perturbações ------
% q = qi - vazao de entrada ?; 
% Ti = Ti - Temeratura externa ?;
subplot(2,1,1)
plot(perturbacoes(1,1:iteracoes),'r','linewidth',2);
hold on
plot(x_pred_vect(4,1:iteracoes),'k--','linewidth',2);
title('$q_i$ Vaz\~{a}o de entrada','interpreter','latex','FontSize',tamTitulo);
legend({'Real','Predição'},'Location','best','FontSize',tamLetra);
xlabel('Iterações [-]','FontSize',tamLetra);
grid


subplot(2,1,2);
plot(perturbacoes(2,1:iteracoes),'r','linewidth',2);
hold on
plot(x_pred_vect(5,1:iteracoes),'k--','linewidth',2);
title('$T_i$ Temperatura externa','interpreter','latex','FontSize',tamTitulo)
legend({'Real','Predição'},'Location','best','FontSize',tamLetra);
xlabel('Iterações [-]','FontSize',tamLetra);
grid


