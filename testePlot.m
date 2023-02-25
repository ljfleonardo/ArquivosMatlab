%Estados / Variáveis Controladas
% x1 = h - nível dentro do tanque; 
% x2 = Ca - Concentração de saída do produto A;
% x3 = T - Temperatura dentro do reator; 
% x4 = R - Velocidade de reação?
figure
subplot(3,1,1);
plot(x_pred_vect(1,300:end),'b--','linewidth',2);
hold on
    plot(ref_h(300:end),'r--','linewidth',2);

plot(saidas(1,300:end),'b','linewidth',2);
grid
% legend({'Predição','Real'},'Location','best','FontSize',10);
legend({'Predição','Real','Referência'},'Location','best','FontSize',10);

title('h - Altura do tanque')

subplot(3,1,2);
plot(x_pred_vect(2,300:end),'b--','linewidth',2);
hold on
plot(saidas(2,300:end),'b','linewidth',2);
% if(ref_ca)
    plot(ref_ca(300:end),'r--','linewidth',2);
% end
legend({'Predição','Real','Referência'},'Location','best','FontSize',10);
% legend({'Predição','Real'},'Location','best','FontSize',10);
title('Ca - Concentração do produto A');
grid

subplot(3,1,3);
plot(x_pred_vect(3,300:end),'b--','linewidth',2);
hold on
    plot(ref_T(300:end),'r--','linewidth',2);
plot(saidas(3,300:end),'b','linewidth',2);
% legend({'Predição','Real'},'Location','best','FontSize',10);
legend({'Predição','Real','Referência'},'Location','best','FontSize',10);

title('T - Temperatura dentro do tanque');
grid

return
% figure
%Entradas / Variáveis Manipuladas
% u1 = qo - Vazão de saída;
% u2 = Caf - Concentração do produto A na alimentação do tanque; 
% u3 = Qh/pCp - taxa de remoção de calor normalizada;
subplot(3,2,2);
plot(entradas(1,:),'k','linewidth',2);
hold on
title('Entradas')
legend({'q0 - Vazão de saída'},'Location','best','FontSize',10);
grid

subplot(3,2,4);
plot(entradas(2,:),'k','linewidth',2);
legend({'Caf - Concentração produto A na alimentação do tanque'},'Location','best','FontSize',10);
grid

subplot(3,2,6);
plot(entradas(3,:),'k','linewidth',2);
legend({'Qh/pcp - Taxa de remoção de calor normalizada'},'Location','best','FontSize',10);
grid

figure
%Perturbações
% q = qi - vazao de entrada ?; 
% Ti = Ti - Temeratura externa ?;
subplot(2,1,1)
plot(x_pred_vect(4,:),'r--','linewidth',2);
hold on
plot(perturbacoes(1,:),'r','linewidth',2);
title('qi - Vazão de entrada')
legend({'Predição','Real'},'Location','best','FontSize',10);
grid


subplot(2,1,2);
plot(x_pred_vect(5,:),'r--','linewidth',2);
hold on
plot(perturbacoes(2,:),'r','linewidth',2);
title('Ti - Temperatura externa')
legend({'Predição','Real'},'Location','best','FontSize',10);
grid


