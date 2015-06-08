LBScap = csvread('\\icnas2.cc.ic.ac.uk\sl5412\L_Bent_Elbow_Scapula_Plane_01-cut1-.csv');
LSCor = csvread('\\icnas2.cc.ic.ac.uk\sl5412\L Straight Arm Coronal Plane 02-cut1.c3d.csv.csv');
RSScap=csvread('H:\R Straight Arm Scapula Plane 01-cut.c3d.csv.csv');
RSCor=csvread('H:\R Straight Arm Coronal Elevation 01-cut.c3d.csv.csv');



%% Calibrating between two datasets and calculating error
numFrames=600;
RMSE=zeros(1,numFrames);

for i=1:numFrames
    % Extracting single frame
    frame1=RSScap(i,2:(size(RSScap,2)));
    frame1=reshape(frame1,[3,length(frame1)/3])';


    frame2=RSCor(i,2:(size(RSCor,2)));
    frame2=reshape(frame2,[3,length(frame2)/3])';
    
    % Filter for calibration cluster in frame 1
    
    for j=1:size(frame1,1)
        counter=0;

        x0=frame1(j,1);
        y0=frame1(j,2);
        z0=frame1(j,3);

        for k=1:size(frame1,1)


            x=frame1(k,1);
            y=frame1(k,2);
            z=frame1(k,3);


            if (sqrt((x-x0)^2 + (y-y0)^2 + (z-z0)^2)<80 && (sqrt((x-x0)^2 + (y-y0)^2 + (z-z0)^2)>30))
                counter=counter+1;
                xFrame1(counter)=x;
                yFrame1(counter)=y;
                zFrame1(counter)=z;

            end

        end
        if (counter==4)
%             hold on;
%             plotFilt=scatter3(x0,z0,y0,'r');
%             view(3);
%             grid on
            xFrame1(5)=x0;
            yFrame1(5)=y0;
            zFrame1(5)=z0;
            break;
        end
    end
    
    
    % Filter for calibration cluster in frame 2
    for l=1:size(frame2,1)
        counter=0;

        x0=frame2(l,1);
        y0=frame2(l,2);
        z0=frame2(l,3);

        for m=1:size(frame2,1)


            x=frame2(m,1);
            y=frame2(m,2);
            z=frame2(m,3);


            if (sqrt((x-x0)^2 + (y-y0)^2 + (z-z0)^2)<80 && (sqrt((x-x0)^2 + (y-y0)^2 + (z-z0)^2)>30))
                counter=counter+1;
                xFrame2(counter)=x;
                yFrame2(counter)=y;
                zFrame2(counter)=z;
            end

        end
        if (counter==4)
%             hold on;
%             plotFilt=scatter3(x0,z0,y0,'r');
%             view(3);
%             grid on
            xFrame2(5)=x0;
            yFrame2(5)=y0;
            zFrame2(5)=z0;
            break;
        end
    end   
    
    
    % Make matrices for calibration algorithm
    cluster1=[xFrame1' yFrame1' zFrame1'];
    cluster2=[xFrame2' yFrame2' zFrame2'];
    
    
    

    % Calibration
    [R,t]=rigid_transform_3D(cluster1,cluster2);
    calEstimate=(cluster1*R'+[t,t,t,t,t]');
    
    
    % Calculate RMSE
    D = abs(cluster2-calEstimate).^2;
    RMSE(i) = sqrt(sum(D(:))/numel(cluster2));
end
% Average RMSE
Average_RMSE=mean(RMSE)

%Standard deviation
Standard_deviation_RMSE=std(RMSE)

% Plot RMSE over time
figure;
plot(RMSE);
hold on;
plot(ones(length(RMSE))*Average_RMSE,'r');
grid on;
xlabel('Frame');
ylabel('mm');
legend('RMSE','Average RMSE');
