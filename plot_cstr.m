tamLetra = 10;
tamTitulo = 12;
espes = 3;
global Ts
segundos = (1:iteracoes)*Ts;

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

subplot(2,2,1);
plot(segundos,saidas(1,:),'r','linewidth',espes);
hold on
plot(segundos,x_pred_vect(1,:),'k--','linewidth',espes)%,'marker','x');
if controle==1 
    plot(segundos,ref_ca(:),'b--','linewidth',espes);
    legend({'Real','Predição','Referência'},'FontSize',tamLetra);
else
   legend({'Real','Predição'},'FontSize',tamLetra); 
end
xlim([0 segundos(end)])
xticks(0:150:segundos(end))
% ylim([0.98 1.08])
xlabel('Tempo (s)','FontSize',tamLetra);
ylabel('Concentração (kmol/m^3)','FontSize',tamLetra);
title('$C_a$ - Concentra\c{c}\~{a}o do produto A','interpreter','latex','FontSize',tamTitulo);
grid


subplot(2,2,3);
plot(segundos,saidas(2,:),'r','linewidth',espes);
hold on
plot(segundos,x_pred_vect(2,:),'k--','linewidth',espes);
if controle==1 
    plot(segundos,ref_T(:),'b--','linewidth',espes);
    legend({'Real','Predição','Referência'},'FontSize',tamLetra);
else
   legend({'Real','Predição'},'FontSize',tamLetra); 
end
xlim([0 segundos(end)])
xticks(0:150:segundos(end))
% ylim([390 410])
xlabel('Tempo (s)','FontSize',tamLetra);
ylabel('Temperatura Interna (K)','FontSize',tamLetra);
title('$T$ - Temperatura dentro do tanque','interpreter','latex','FontSize',tamTitulo);
grid

%% ------------ Gráfico das entradas ------------
% figure
%----- Entradas / Variáveis Manipuladas -----
% u1 = qo - Vazão de saída;
% u2 = Caf - Concentração do produto A na alimentação do tanque; 
% u3 = Qh/pCp - taxa de remoção de calor normalizada;


subplot(2,2,2);
% plot(entradas_atrasadas_vect(2,1:iteracoes),'k--','linewidth',2);
% hold on
% plot(entradas_comp_vect(2,1:iteracoes),'b','linewidth',2);
plot(segundos,entradas(1,:),'r','linewidth',espes);
legend(['$C_{af}$ - Concentra\c{c}\~{a}o produto A' newline 'na alimenta\c{c}\~{a}o do tanque'],'interpreter','latex','Location','best','FontSize',tamLetra);
xlim([0 segundos(end)])
xticks(0:150:segundos(end))
% ylim([4.9 5.2])
xlabel('Tempo (s)','FontSize',tamLetra);
ylabel('Concentração (kmol/m^3)','FontSize',tamLetra)
grid

subplot(2,2,4);
% plot(entradas_atrasadas_vect(3,1:iteracoes),'k--','linewidth',2);
% hold on
% plot(entradas_comp_vect(3,1:iteracoes),'b','linewidth',2);
plot(segundos,entradas(2,:),'r','linewidth',espes);
legend(['${Qh}/{pc_p}$ - Taxa de remo\c{c}\~{a}o' newline 'de calor normalizada'],'interpreter','latex','Location','best','FontSize',tamLetra);
xlim([0 segundos(end)])
xticks(0:150:segundos(end))
% ylim([0.71 0.82])
xlabel('Tempo (s)','FontSize',tamLetra);
ylabel({'Taxa de remoção de'; 'calor (Km^3/s^{-1})'},'FontSize',tamLetra)
grid

if ekf == 0
    if controle == 0
%         if f==1
        sg_1 = sprintf('FSP Nao-Linear em malha aberta');
        sg_2 = sprintf('$$ F_r = %s $$ ', tf_title);
        sgtitle({sg_1,sg_2},'interpreter','latex')
    else 
        sg_1 = sprintf('FSP Nao-Linear em malha fechada');
        sg_2 = sprintf('$$ F_r = %s $$ ', tf_title);
        sgtitle({sg_1,sg_2},'interpreter','latex')
    end
else
     if controle == 0
        sgtitle('O&P com EKF em malha aberta');
    else 
        sgtitle('O&P com EKF em malha fechada');
    end
end