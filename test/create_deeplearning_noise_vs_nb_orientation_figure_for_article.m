clear
close all



val_lost_15_noise05 = [0.19430694700082143, 0.17021300185521443, 0.1594601011912028, 0.15875361088116963, 0.1522194670041402, ...
                       0.15028052077293397, 0.14727029361724853, 0.1454579091310501, 0.143324094525973, 0.14253146015008292, ...
                       0.13977111157576244, 0.1377464103460312, 0.13644190277258555, 0.1363016130208969, 0.13583138999144237];   
            
plot(val_lost_15_noise05, '-', 'LineWidth', 3);
hold on

val_lost_15_noise_1 = [0.20591944970289866, 0.18217554452419282, 0.1742391024986903, 0.1713144571463267, 0.16718870964050292, ...
                       0.16415977447827657, 0.1606951448281606, 0.15982555096944173, 0.1573519329468409, 0.15661284449895224, ...
                       0.1551549522717794, 0.15142291674613953, 0.150662149922053, 0.14920144217014314, 0.1489879332224528];
                       
                       
plot(val_lost_15_noise_1, '-', 'LineWidth', 3);



val_lost_15_noise2 = [0.22001085085074107, 0.2003790401617686, 0.19329011374314625, 0.18881990302403767, 0.1848405996799469, ...
                      0.1821978758096695, 0.17952379978497823, 0.17758157347043355, 0.17558464884757996, 0.17398655341466268, ...
                      0.1723179187933604, 0.1700329360405604, 0.16886445336341857, 0.16643994994958242, 0.16681532655556996];
    
plot(val_lost_15_noise2, '-', 'LineWidth', 3);

xlabel('number of orientations')
ylabel('val loss')
title('After 10 epochs')
leg = legend('0.5%', '1%', '2%');
title(leg, 'noise')
set(gca, 'FontSize', 24, 'FontWeight','bold')

% ylim([0 0.5])
    