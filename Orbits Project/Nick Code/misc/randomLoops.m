%% Useless Loops
% for i=1:num_spacecraft
%     for j=1:length(data.launches(i).payload)
%         for k=1:num_launches
%             satellite_list(i).name = data.launches(k).payload(1).name;
%             
%             
%           

% for i=1:num_spacecraft
%     for j=1:num_launches
%         fprintf("1")
%         for k=1:length(data.launches(j).payload)
%             satellite_list(i).name = data.launches(j).payload(k).name;
%         end
%     end
% end

% list = zeros(num_spacecraft,1)

%% Backup of loadConstellation
% fid = fopen(filename);
% raw = fread(fid,inf);
% str = char(raw'); 
% fclose(fid); 
% data = jsondecode(str);
% 
% num_launches = length(data.launches);
% num_spacecraft = 0;
% for i=1:num_launches
%     num_spacecraft = num_spacecraft + length(data.launches(i).payload);
%     
% end
% satellite_list(num_spacecraft).name = '';
% satellite_list(num_spacecraft).oe0 = NaN(6,1);
% 
% i = 1;
% j = 1;
% k = 1;
% for j=1:num_launches
%     oe0Temp1 = cell2mat(struct2cell((data.launches(j).orbit)));
%     for k=1:length(data.launches(j).payload)
%         nameTemp = data.launches(j).payload(k).name;
%         oe0Temp2 = data.launches(j).payload(k).f;
%         satellite_list(i).name = nameTemp;
%         satellite_list(i).oe0 = [oe0Temp1; oe0Temp2];
%         i = i + 1;
%     end
%     
% end

%% Broken loops to plot in 3D

% % t0 = 1;
% % t = 500;
% % tf = 24*60;
% % tstep = 1;
% % 
% % oe0 = [502512;.99;0;0;0;3.14150000000000/4];
% % x = zeros(tf,6);
% % x(1,:) = oe0;
% % 
% % for i=1:tstep:tf+1
% %     x(i+1,:) = propagateState(x(i,:),i+1,i,mu,J2,Re);
% % %     memory(i+1,:) = x';
% % end
% % hold on
% % plot3(x(:,1),x(:,2),x(:,3));

% oe0 = [502512;.99;0;0;0;3.14150000000000/4];
% t0 = 1;
% t = 500;
% tf = 24*60*60;
% tstep = 30;
% memory = zeros(30,6)
% for i=1:tstep:tf+1
%     x = propagateState(oe0,i+1,i,mu,J2,Re)
%     memory(i+1,:) = x';
% end
% hold on
% plot3(memory(:,1),memory(:,2),memory(:,3));

% for i=1:86400/tstep
%     x(i+1,:) = propagateState(oe0,i*tstep,1,mu,J2,Re)';
% %     [i*tstep i*tstep-tstep]
% end

%%
% tol = 1e-9;
% f = @(E) E-e*sin(E)-M; %Kepler's ToF
% fp = @(E) 1-e*cos(E); %Derivative of Kepler's ToF
% E0 = .0000001; %Initial E guess
% E(1) = E0;
% K = 100; %Max iterations to find E
% k = 1; 
% while (k <= K)
%     fe = f(E(k));
%     fpe = fp(E(k));
%     E(k+1) = E(k) - fe/fpe;
%     if (abs(fe) <= tol)
%         break;
%     end
%     k = k + 1;
% end
% rtemp = a*(1-e*cos(E));
% fmatrix = (acos((p./rtemp)-1)./e);
% f = fmatrix(length(fmatrix));

