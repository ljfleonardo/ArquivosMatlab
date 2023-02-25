%Estados / Vari�veis Controladas
% x1 = h - n�vel dentro do tanque; 
% x2 = Ca - Concentra��o de sa�da do produto A;
% x3 = T - Temperatura dentro do reator; 
% x4 = R - Velocidade de rea��o?
figure
subplot(3,1,1);
plot(x_pred_vect(1,300:end),'b--','linewidth',2);
hold on
    plot(ref_h(300:end),'r--','linewidth',2);

plot(saidas(1,300:end),'b','linewidth',2);
grid
% legend({'Predi��o','Real'},'Location','best','FontSize',10);
legend({'Predi��o','Real','Refer�ncia'},'Location','best','FontSize',10);

title('h - Altura do tanque')

subplot(3,1,2);
plot(x_pred_vect(2,300:end),'b--','linewidth',2);
hold on
plot(saidas(2,300:end),'b','linewidth',2);
% if(ref_ca)
    plot(ref_ca(300:end),'r--','linewidth',2);
% end
legend({'Predi��o','Real','Refer�ncia'},'Location','best','FontSize',10);
% legend({'Predi��o','Real'},'Location','best','FontSize',10);
title('Ca - Concentra��o do produto A');
grid

subplot(3,1,3);
plot(x_pred_vect(3,300:end),'b--','linewidth',2);
hold on
    plot(ref_T(300:end),'r--','linewidth',2);
plot(saidas(3,300:end),'b','linewidth',2);
% legend({'Predi��o','Real'},'Location','best','FontSize',10);
legend({'Predi��o','Real','Refer�ncia'},'Location','best','FontSize',10);

title('T - Temperatura dentro do tanque');
grid

return
% figure
%Entradas / Vari�veis Manipuladas
% u1 = qo - Vaz�o de sa�da;
% u2 = Caf - Concentra��o do produto A na alimenta��o do tanque; 
% u3 = Qh/pCp - taxa de remo��o de calor normalizada;
subplot(3,2,2);
plot(entradas(1,:),'k','linewidth',2);
hold on
title('Entradas')
legend({'q0 - Vaz�o de sa�da'},'Location','best','FontSize',10);
grid

subplot(3,2,4);
plot(entradas(2,:),'k','linewidth',2);
legend({'Caf - Concentra��o produto A na alimenta��o do tanque'},'Location','best','FontSize',10);
grid

subplot(3,2,6);
plot(entradas(3,:),'k','linewidth',2);
legend({'Qh/pcp - Taxa de remo��o de calor normalizada'},'Location','best','FontSize',10);
grid

figure
%Perturba��es
% q = qi - vazao de entrada ?; 
% Ti = Ti - Temeratura externa ?;
subplot(2,1,1)
plot(x_pred_vect(4,:),'r--','linewidth',2);
hold on
plot(perturbacoes(1,:),'r','linewidth',2);
title('qi - Vaz�o de entrada')
legend({'Predi��o','Real'},'Location','best','FontSize',10);
grid


subplot(2,1,2);
plot(x_pred_vect(5,:),'r--','linewidth',2);
hold on
plot(perturbacoes(2,:),'r','linewidth',2);
title('Ti - Temperatura externa')
legend({'Predi��o','Real'},'Location','best','FontSize',10);
grid


