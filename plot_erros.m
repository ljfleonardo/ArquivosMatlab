tamLetra = 10;
tamTitulo = 12;
espes = 3;
global Ts
segundos = (1:iteracoes)*Ts;

%% ------------ Gráfico de erros Erro Quadrático Médio ------------
bool_MSE = exist ('err_MSE_1_ekf','var');
if bool_MSE==1
    h2 = figure();
    h2.WindowState = 'maximized';
    h2.Position = p(2,:);    
    subplot(2,1,1);
    if exist ('err_MSE_1_ekf','var')&& exist ('err_MSE_1_fsp','var')
        plot(segundos,err_MSE_1_ekf,'r','linewidth',espes);
        hold on
        plot(segundos,err_MSE_1_fsp,'b','linewidth',espes);
        legend({'EKF','FSP NÃO-LINEAR'},'Location','best','FontSize',tamLetra);
        
    else if exist ('err_MSE_1_ekf','var') && ~exist ('err_MSE_1_fsp','var')
            plot(segundos,err_MSE_1_ekf,'r','linewidth',espes);
            legend({'EKF'},'Location','best','FontSize',tamLetra);
            
        else if ~exist ('err_MSE_1_ekf','var') && exist ('err_MSE_1_fsp','var')
                plot(segundos,err_MSE_1_fsp,'r','linewidth',espes);
                legend({'FSP NÃO-LINEAR'},'Location','best','FontSize',tamLetra);
            end
        end
    end
    xlim([0 segundos(end)])
    xticks(0:200:segundos(end))
    % ylim([0.98 1.08])
    xlabel('Tempo (s)','FontSize',tamLetra);
    ylabel('Concentração (kmol/m^3)','FontSize',tamLetra);
    title('$C_a$ - Concentra\c{c}\~{a}o do produto A','interpreter','latex','FontSize',tamTitulo);
    grid
    
    subplot(2,1,2);
    if exist ('err_MSE_2_ekf','var')&& exist ('err_MSE_2_fsp','var')
        plot(segundos,err_MSE_2_ekf,'r','linewidth',espes);
        hold on
        plot(segundos,err_MSE_2_fsp,'b','linewidth',espes);
        legend({'EKF','FSP NÃO-LINEAR'},'Location','best','FontSize',tamLetra);
        
    else if exist ('err_MSE_2_ekf','var') && ~exist ('err_MSE_2_fsp','var')
            plot(segundos,err_MSE_2,'r','linewidth',espes);
            legend({'EKF'},'Location','best','FontSize',tamLetra);
            
        else if ~exist ('err_MSE_2_ekf','var') && exist ('err_MSE_2_fsp','var')
                plot(segundos,err_MSE_2_fsp,'r','linewidth',espes);
                legend({'FSP NÃO-LINEAR'},'Location','best','FontSize',tamLetra);
            end
        end
    end
    xlim([0 segundos(end)])
    xticks(0:200:segundos(end))
    % ylim([390 410])
    xlabel('Tempo (s)','FontSize',tamLetra);
    ylabel('Temperatura Interna (K)','FontSize',tamLetra);
    title('$T$ - Temperatura dentro do tanque','interpreter','latex','FontSize',tamTitulo);
    grid
    
    if ekf == 0
        sgtitle('MSE FSP Não-Linear');
    else
        sgtitle('MSE O&P com EKF');
    end
end
%% ------------ Gráfico de erros Raiz Quadrada do Erro Quadrático Médio  ------------
bool_RMSE = exist ('err_RMSE_1_ekf','var');
if bool_RMSE==1
    h2 = figure();
    h2.WindowState = 'maximized';
    h2.Position = p(2,:);
    
    subplot(2,1,1);
    if exist ('err_RMSE_1_ekf','var')&& exist ('err_RMSE_1_fsp','var')
        plot(segundos,err_RMSE_1_ekf,'r','linewidth',espes);
        hold on
        plot(segundos,err_RMSE_1_fsp,'b','linewidth',espes);
        legend({'EKF','FSP NÃO-LINEAR'},'Location','best','FontSize',tamLetra);
        
    else if exist ('err_RMSE_1_ekf','var') && ~exist ('err_RMSE_1_fsp','var')
            plot(segundos,err_RMSE_1_ekf,'r','linewidth',espes);
            legend({'EKF'},'Location','best','FontSize',tamLetra);
            
        else if ~exist ('err_RMSE_1_ekf','var') && exist ('err_RMSE_1_fsp','var')
                plot(segundos,err_RMSE_1_fsp,'r','linewidth',espes);
                legend({'FSP NÃO-LINEAR'},'Location','best','FontSize',tamLetra);
            end
        end
    end
    xlim([0 segundos(end)])
    xticks(0:200:segundos(end))
    % ylim([0.98 1.08])
    xlabel('Tempo (s)','FontSize',tamLetra);
    ylabel('Concentração (kmol/m^3)','FontSize',tamLetra);
    title('$C_a$ - Concentra\c{c}\~{a}o do produto A','interpreter','latex','FontSize',tamTitulo);
    grid
    
    subplot(2,1,2);
    if exist ('err_RMSE_2_ekf','var')&& exist ('err_RMSE_2_fsp','var')
        plot(segundos,err_RMSE_2_ekf,'r','linewidth',espes);
        hold on
        plot(segundos,err_RMSE_2_fsp,'b','linewidth',espes);
        legend({'EKF','FSP NÃO-LINEAR'},'Location','best','FontSize',tamLetra);
        
    else if exist ('err_RMSE_2_ekf','var') && ~exist ('err_RMSE_2_fsp','var')
            plot(segundos,err_RMSE_2_ekf,'r','linewidth',espes);
            legend({'EKF'},'Location','best','FontSize',tamLetra);
            
        else if ~exist ('err_RMSE_2_ekf','var') && exist ('err_RMSE_2_fsp','var')
                plot(segundos,err_RMSE_2_fsp,'r','linewidth',espes);
                legend({'FSP NÃO-LINEAR'},'Location','best','FontSize',tamLetra);
            end
        end
    end
    xlim([0 segundos(end)])
    xticks(0:200:segundos(end))
    % ylim([390 410])
    xlabel('Tempo (s)','FontSize',tamLetra);
    ylabel('Temperatura Interna (K)','FontSize',tamLetra);
    title('$T$ - Temperatura dentro do tanque','interpreter','latex','FontSize',tamTitulo);
    grid
    
    if ekf == 0
        sgtitle('RMSE FSP Não-Linear');
    else
        sgtitle('RMSE O&P com EKF');
    end
end
%% ------------ Gráfico de erros Porcentagem Erro Absoluto Médio ------------
bool_MAPE = exist ('err_MAPE_1_ekf','var');
if bool_MAPE==1
    h2 = figure();
    h2.WindowState = 'maximized';
    h2.Position = p(2,:);
    
    subplot(2,1,1);
    if exist ('err_MAPE_1_ekf','var')&& exist ('err_MAPE_1_fsp','var')
        plot(segundos,err_MAPE_1_ekf,'r','linewidth',espes);
        hold on
        plot(segundos,err_MAPE_1_fsp,'b','linewidth',espes);
        legend({'EKF','FSP NÃO-LINEAR'},'Location','best','FontSize',tamLetra);
        
    else if exist ('err_MAPE_1_ekf','var') && ~exist ('err_MAPE_1_fsp','var')
            plot(segundos,err_MAPE_1_ekf,'r','linewidth',espes);
            legend({'EKF'},'Location','best','FontSize',tamLetra);
            
        else if ~exist ('err_MAPE_1_ekf','var') && exist ('err_MAPE_1_fsp','var')
                plot(segundos,err_MAPE_1_fsp,'r','linewidth',espes);
                legend({'FSP NÃO-LINEAR'},'Location','best','FontSize',tamLetra);
            end
        end
    end
    xlim([0 segundos(end)])
    xticks(0:200:segundos(end))
    % ylim([0.98 1.08])
    xlabel('Tempo (s)','FontSize',tamLetra);
    ylabel('(%)','FontSize',tamLetra);
    title('$C_a$ - Concentra\c{c}\~{a}o do produto A','interpreter','latex','FontSize',tamTitulo);
    grid
    
    subplot(2,1,2);
    if exist ('err_MAPE_2_ekf','var')&& exist ('err_MAPE_2_fsp','var')
        plot(segundos,err_MAPE_2_ekf,'r','linewidth',espes);
        hold on
        plot(segundos,err_MAPE_2_fsp,'b','linewidth',espes);
        legend({'EKF','FSP NÃO-LINEAR'},'Location','best','FontSize',tamLetra);
    else if exist ('err_MAPE_2_ekf','var') && ~exist ('err_MAPE_2_fsp','var')
            plot(segundos,err_MAPE_2_ekf,'r','linewidth',espes);
            legend({'EKF'},'Location','best','FontSize',tamLetra);
        else if ~exist ('err_MAPE_2_ekf','var') && exist ('err_MAPE_2_fsp','var')
                plot(segundos,err_MAPE_2_fsp,'r','linewidth',espes);
                legend({'FSP NÃO-LINEAR'},'Location','best','FontSize',tamLetra);
            end
        end
    end
    xlim([0 segundos(end)])
    xticks(0:200:segundos(end))
    % ylim([390 410])
    xlabel('Tempo (s)','FontSize',tamLetra);
    ylabel('(%)','FontSize',tamLetra);
    title('$T$ - Temperatura dentro do tanque','interpreter','latex','FontSize',tamTitulo);
    grid
    
    if ekf == 0
        sgtitle('MAPE FSP Não-Linear');
    else
        sgtitle('MAPE O&P com EKF');
    end
end