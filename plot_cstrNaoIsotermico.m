if ekf == 0
    saidas=x;
end
tamLetra = 10;
tamTitulo = 12;
espes = 3;

h = figure();
h.WindowState = 'maximized';
p = get(0,"MonitorPositions");%pega o 2 monitor
h.Position = p(2,:); %plota no 2 monitor
%% ------------ Gráfico das saídas ------------
%Estados / Variáveis Controladas
% x1 = h - nível dentro do tanque; 
% x2 = Ca - Concentração de saída do produto A;
% x3 = T - Temperatura dentro do reator; 
% x4 = R - Velocidade de reação?

subplot(3,2,1);
plot(saidas(1,1:iteracoes),'r','linewidth',espes);
hold on
plot(x_pred_vect(1,1:iteracoes),'k--','linewidth',espes);
if controle==1 
    plot(ref_h(1:iteracoes),'b--','linewidth',espes);
    legend({'Real','Predição','Referência'},'FontSize',tamLetra);
else
   legend({'Real','Predição'},'FontSize',tamLetra); 
end
xlim([0 iteracoes])
% ylim([0.92 1.2])
xlabel('Iterações (-)','FontSize',tamLetra);
ylabel('Altura (m)','FontSize',tamLetra);
title('$h$ - Altura do tanque','interpreter','latex','FontSize',tamTitulo)
grid


subplot(3,2,3);
plot(saidas(2,1:iteracoes),'r','linewidth',espes);
hold on
plot(x_pred_vect(2,1:iteracoes),'k--','linewidth',espes)%,'marker','x');
if controle==1 
    plot(ref_ca(1:iteracoes),'b--','linewidth',espes);
    legend({'Real','Predição','Referência'},'FontSize',tamLetra);
else
   legend({'Real','Predição'},'FontSize',tamLetra); 
end
xlim([0 iteracoes])
% ylim([0.98 1.08])
xlabel('Iterações (-)','FontSize',tamLetra);
ylabel('Concentração (kmol/m^3)','FontSize',tamLetra);
title('$C_a$ - Concentra\c{c}\~{a}o do produto A','interpreter','latex','FontSize',tamTitulo);
grid


subplot(3,2,5);
plot(saidas(3,1:iteracoes),'r','linewidth',espes);
hold on
plot(x_pred_vect(3,1:iteracoes),'k--','linewidth',espes);
if controle==1 
    plot(ref_T(1:iteracoes),'b--','linewidth',espes);
    legend({'Real','Predição','Referência'},'FontSize',tamLetra);
else
   legend({'Real','Predição'},'FontSize',tamLetra); 
end
xlim([0 iteracoes])
% ylim([390 410])
xlabel('Iterações (-)','FontSize',tamLetra);
ylabel('Temperatura Interna (K)','FontSize',tamLetra);
title('$T$ - Temperatura dentro do tanque','interpreter','latex','FontSize',tamTitulo);
grid

%% ------------ Gráfico das entradas ------------
% figure
%----- Entradas / Variáveis Manipuladas -----
% u1 = qo - Vazão de saída;
% u2 = Caf - Concentração do produto A na alimentação do tanque; 
% u3 = Qh/pCp - taxa de remoção de calor normalizada;
subplot(3,2,2);
% plot(entradas_atrasadas_vect(1,1:iteracoes),'k--','linewidth',2);
% hold on
ax=gca;
xlim([0 iteracoes])
ylim([4.8e-3 5.2e-3])
ax.YAxis.Exponent = -3;
% plot(entradas_comp_vect(1,1:iteracoes),'b','linewidth',2);
plot(entradas(1,1:iteracoes),'r','linewidth',espes);
title('Entradas','FontSize',tamTitulo)
legend({'$q_0$ - Vaz\~{a}o de sa\''ida'},'interpreter','latex','Location','best','FontSize',tamLetra);
xlabel('Iterações (-)','FontSize',tamLetra);
ylabel('Vazão (m^3/s^{-1})','FontSize',tamLetra)
grid

subplot(3,2,4);
% plot(entradas_atrasadas_vect(2,1:iteracoes),'k--','linewidth',2);
% hold on
% plot(entradas_comp_vect(2,1:iteracoes),'b','linewidth',2);
plot(entradas(2,1:iteracoes),'r','linewidth',espes);
legend({'$C_{af}$ - Concentra\c{c}\~{a}o produto A na alimenta\c{c}\~{a}o do tanque'},'interpreter','latex','Location','best','FontSize',tamLetra);
xlim([0 iteracoes])
ylim([4.9 5.2])
xlabel('Iterações (-)','FontSize',tamLetra);
ylabel('Concentração (kmol/m^3)','FontSize',tamLetra)
grid

subplot(3,2,6);
% plot(entradas_atrasadas_vect(3,1:iteracoes),'k--','linewidth',2);
% hold on
% plot(entradas_comp_vect(3,1:iteracoes),'b','linewidth',2);
plot(entradas(3,1:iteracoes),'r','linewidth',espes);
legend({'${Qh}/{pc_p}$ - Taxa de remo\c{c}\~{a}o de calor normalizada'},'interpreter','latex','Location','best','FontSize',tamLetra);
xlim([0 iteracoes])
ylim([0.71 0.8])
xlabel('Iterações (-)','FontSize',tamLetra);
ylabel({'Taxa de remoção de'; 'calor (Km^3/s^{-1})'},'FontSize',tamLetra)
grid

if ekf == 0
    if controle == 0
        sgtitle('FSP Não-Linear em malha aberta');
    else 
        sgtitle('FSP Não-Linear em malha fechada');
    end
else
     if controle == 0
        sgtitle('O&P com EKF em malha aberta');
    else 
        sgtitle('O&P com EKF em malha fechada');
    end
end

%% ------------ Gráfico de erros ------------
h2 = figure();
h2.WindowState = 'maximized';
h2.Position = p(2,:);

subplot(3,1,1);
plot(err_1,'r','linewidth',espes);
xlim([0 iteracoes])
% ylim([0.92 1.2])
xlabel('Iterações (-)','FontSize',tamLetra);
ylabel('Altura (m)','FontSize',tamLetra);
title('$h$ - Altura do tanque','interpreter','latex','FontSize',tamTitulo)
grid


subplot(3,1,2);
plot(err_2,'r','linewidth',espes);
xlim([0 iteracoes])
% ylim([0.98 1.08])
xlabel('Iterações (-)','FontSize',tamLetra);
ylabel('Concentração (kmol/m^3)','FontSize',tamLetra);
title('$C_a$ - Concentra\c{c}\~{a}o do produto A','interpreter','latex','FontSize',tamTitulo);
grid

subplot(3,1,3);
plot(err_3,'r','linewidth',espes);
xlim([0 iteracoes])
% ylim([390 410])
xlabel('Iterações (-)','FontSize',tamLetra);
ylabel('Temperatura Interna (K)','FontSize',tamLetra);
title('$T$ - Temperatura dentro do tanque','interpreter','latex','FontSize',tamTitulo);
grid

if ekf == 0
    sgtitle('Erro Quadrático Médio FSP Não-Linear');
else
    sgtitle('Erro Quadrático Médio O&P com EKF');
end