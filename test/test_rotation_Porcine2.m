% test rotation Porcine-2
% From folder
% P:\3015069.01\derived\Porcine-2\ses-mri01\registration\gre\highres
% from files
% Porcine-2_gre_highres_ref_2_orientation-1.fsl
clear;
close all;

rotation(:, :, 1) = [0.02961507104  0.06669085489  -0.9973341665   
0.09015299249  0.9935268762  0.06911334762  
0.9954870124  -0.09195950545  0.02341099896  ];

rotation(:, :, 2) = [0.3339029182  0.06159285721  -0.9405929357  
-0.0526306458  0.9975248315  0.0466375113  
0.9411369717  0.03393171275  0.3363180577  ];

rotation(:, :, 3) = [0.5869444611  0.07325834769  -0.8063056539  
-0.04448012003  0.9973111376  0.05823364073  
0.808404465  0.001684650553  0.5886256256 ];

rotation(:, :, 4) = [0.775384534  0.06872606886  -0.6277383055  
-0.03995800494  0.9974075737  0.05984199618  
0.6302232558  -0.02131742217  0.7761209229  ];

rotation(:, :, 5) = [0.9505737286  0.1042314627  -0.2924823991 
-0.08581687788  0.9934736334  0.07513442909  
0.298405544  -0.0463209039  0.9533141382 ];

rotation(:, :, 6) = [0.9999677585  0.007600607621  0.002543853147  
-0.007648388539  0.9997837338  0.01933185978 
-0.002396375411  -0.01935073068  0.9998097711 ];


original_direction = [0; 0; 1];
for k = 1:6
    orientations(k, :) = rotation(:, :, k) * original_direction;
end
visualize3dArrow(orientations)

return;























A(:, :, 1) = [0.02961507104  0.06669085489  -0.9973341665  143.513275  
0.09015299249  0.9935268762  0.06911334762  -12.1940855  
0.9954870124  -0.09195950545  0.02341099896  -1.946005758  
0  0  0  1 ];

A(:, :, 2) = [0.3339029182  0.06159285721  -0.9405929357  114.1786054  
-0.0526306458  0.9975248315  0.0466375113  1.1493823  
0.9411369717  0.03393171275  0.3363180577  -27.70891601  
0  0  0  1  ];

A(:, :, 3) = [0.5869444611  0.07325834769  -0.8063056539  84.34663287  
-0.04448012003  0.9973111376  0.05823364073  -0.3166022608  
0.808404465  0.001684650553  0.5886256256  -33.23315685  
0  0  0  1  ];

A(:, :, 4) = [0.775384534  0.06872606886  -0.6277383055  57.27089966  
-0.03995800494  0.9974075737  0.05984199618  -0.8324666661  
0.6302232558  -0.02131742217  0.7761209229  -31.15076582  
0  0  0  1  ];

A(:, :, 5) = [0.9505737286  0.1042314627  -0.2924823991  17.39847362  
-0.08581687788  0.9934736334  0.07513442909  1.86368529  
0.298405544  -0.0463209039  0.9533141382  -16.59416843  
0  0  0  1  ];


A(:, :, 6) = [0.9999677585  0.007600607621  0.002543853147  -0.7841131436  
-0.007648388539  0.9997837338  0.01933185978  -0.7099212638  
-0.002396375411  -0.01935073068  0.9998097711  1.378210661  
0  0  0  1  ];