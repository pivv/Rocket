clear all
close all
clc

dimension = input('Input the dimension : ');

%% Case Dimension 2

if dimension==2
%     data1 = load('level_set_algorithm\level_set_algorithm\level1.dat');
%     axis equal, grid on, figure(1), contour(data1);figure(gcf);
%     %surf(data);figure(gcf);
% 
%     data2 = load('level_set_algorithm\level_set_algorithm\level2.dat');
%     hold on, contour(data2);figure(gcf), hold off;
% 
%     figure(2), imshow(data1>0);
%     figure(3), imshow(data2>0);
    
% 



    aviobj = avifile('level.avi','compression','None');

    fig = figure;
    axis([0 1 0 1]);
    for j=0:100
        data3 = load(['level_set_algorithm\level_set_algorithm\tree_level' num2str(j) '.dat']);
        plot(data3(1:1+4,1), data3(1:1+4,2));
        axis([0 1 0 1]); axis square;
        hold on;
        for i=6:5:size(data3,1)
           plot(data3(i:i+4,1), data3(i:i+4,2));
        end
        
        data4 = load(['level_set_algorithm\level_set_algorithm\tree_level_surface' num2str(j) '.dat']);
        plot(data4(1:1+1,1), data4(1:1+1,2), 'r');
        axis([0 1 0 1]); axis square;
        hold on;
        for i=3:2:size(data4,1)
            plot(data4(i:i+1,1), data4(i:i+1,2), 'r');
        end
        
        hold off;
        
        drawnow();
        
        f = getframe(fig);
        aviobj = addframe(aviobj, f);
    end
    close(fig);
    aviobj = close(aviobj);
    
    
    
%     aviobj = avifile('level.avi','compression','None');
% 
%     fig = figure;
%     for j=0:39
%         data3 = load(['level_set_algorithm\level_set_algorithm\tree_level' num2str(j) '.dat']);
%         plot3(data3(1:1+4,1), data3(1:1+4,2), data3(1:1+4,3));
%         axis([0 1 0 1 -1 1]); axis square;
%         hold on;
%         for i=6:5:size(data3,1)
%            plot3(data3(i:i+4,1), data3(i:i+4,2), data3(i:i+4,3));
%         end
%         
%         hold off;
%         
%         drawnow();
%         
%         f = getframe(fig);
%         aviobj = addframe(aviobj, f);
%     end
%     close(fig);
%     aviobj = close(aviobj);


    data3 = load('level_set_algorithm\level_set_algorithm\tree_reinitial.dat');
    x = linspace(0,99);
    plot(x,data3); axis([0 100 0 1]);



    

%     fig = figure;
%     axis([0 100 0 100]);
%     data4 = load(['level_set_algorithm\level_set_algorithm\tree_level_surface' num2str(30) '.dat']);
%     hold on;
%     plot(data4(1:1+1,1), data4(1:1+1,2), 'r');
%     for i=3:2:size(data4,1)
%         plot(data4(i:i+1,1), data4(i:i+1,2), 'r');
%     end
%     hold off;
%     drawnow();
    


%     data3 = load('level_set_algorithm\level_set_algorithm\tree_level0.dat');
%     figure(6), axis([0 100 0 100 -100 100]);
%     hold on;
%     for i=1:5:size(data3,1)
%         line(data3(i:i+4,1), data3(i:i+4,2), data3(i:i+4,3));
%     end
% 
%     hold off;
% 
%     drawnow();

    figure(5);
%     data31 = load('level_set_algorithm\level_set_algorithm\tree_level_initial.dat');
%     plot(data31(1:1+4,1), data31(1:1+4,2));
%     axis([0 1 0 1], 'square');
%     hold on;
%     for i=6:5:size(data31,1)
%        plot(data31(i:i+4,1), data31(i:i+4,2));
%     end
    data41 = load('level_set_algorithm\level_set_algorithm\tree_level_surface_initial.dat');
    plot(data41(1:1+1,1), data41(1:1+1,2), 'r');
    axis([0 1 0 1], 'square');
    hold on;
    for i=3:2:size(data41,1)
        plot(data41(i:i+1,1), data41(i:i+1,2), 'r');
    end
    hold off;
    drawnow();

    figure(6);
    axis([0 1 0 1], 'square');
%     data31 = load('level_set_algorithm\level_set_algorithm\tree_level_result.dat');
%     plot(data31(1:1+4,1), data31(1:1+4,2));
%     hold on;
%     for i=6:5:size(data31,1)
%         plot(data31(i:i+4,1), data31(i:i+4,2));
%     end
    hold on;
    data41 = load('level_set_algorithm\level_set_algorithm\tree_level_surface_result.dat');
    plot(data41(1:1+1,1), data41(1:1+1,2), 'r');
    for i=3:2:size(data41,1)
        plot(data41(i:i+1,1), data41(i:i+1,2), 'r');
    end
%     data41 = load(['level_set_algorithm\level_set_algorithm\tree_levelr_surface' num2str(1) '.dat']);
%     plot(data41(1:1+1,1), data41(1:1+1,2), 'g');
%     for i=3:2:size(data41,1)
%         plot(data41(i:i+1,1), data41(i:i+1,2), 'g');
%     end
    hold off;
    drawnow();
    
    
%     figure,
%     
% %     data3 = load('level_set_algorithm\level_set_algorithm\tree_level_result.dat');
% %     plot(data3(1:1+4,1), data3(1:1+4,2));
% %     hold on;
% %     for i=6:5:size(data3,1)
% %         plot(data3(i:i+4,1), data3(i:i+4,2));
% %     end
% 
%     data4 = load(['level_set_algorithm\level_set_algorithm\tree_level_surface' num2str(1) '.dat']);
%     plot(data4(1:1+1,1), data4(1:1+1,2));
%     hold on;
%     for i=3:2:size(data4,1)
%         plot(data4(i:i+1,1), data4(i:i+1,2));
%     end
    
    
%     data4 = load('level_set_algorithm\level_set_algorithm\tree_level2.dat');
%     figure(5), axis([0 100 0 100]);
%     for i=1:5:size(data4,1)
%         line(data4(i:i+3,1), data4(i:i+3,2));
%     end

%     for j=0:29
%         data3 = load(['level_set_algorithm\level_set_algorithm\tree_level' num2str(j) '.dat']);
%         figure(6), axis([0 100 0 100 -100 100]);
%         hold on;
%         for i=1:5:size(data3,1)
%             line(data3(i:i+4,1), data3(i:i+4,2), data3(i:i+4,3));
%         end
%         
%         hold off;
%         
%         drawnow();
%     end
    
    
%     figure(7), axis([0 100 0 100 -100 100]);
%     for i=1:5:size(data4,1)
%         line(data4(i:i+3,1), data4(i:i+3,2), data4(i:i+3,3));
%     end
end








%% Case Dimension 3


if dimension==3
    
    max_dim = input('maximum dimension is : ');
    size = 2^max_dim;
    
    l = linspace(0,1,size+1);
    lp = linspace(1,0,size+1);
    
%     data1 = load('level_set_algorithm\level_set_algorithm\level_3d_initial.dat');
%     data1 = reshape(data1,size+1,size+1,size+1);
%     %figure(1), axis equal, grid on, contour(data1(:,:,50));figure(gcf);
%     figure(1), isosurface(l,l,l,data1, 0), axis([0 1 0 1 0 1]);
%     
%     data3 = load('level_set_algorithm\level_set_algorithm\level_3d_final.dat');
%     data3 = reshape(data3,size+1,size+1,size+1);
%     %figure(1), axis equal, grid on, contour(data2(:,:,50));figure(gcf);
%     figure(3), isosurface(l,l,l,data3, 0), axis([0 1 0 1 0 1]);
    
    
    axis([0 1 0 1 0 1]);
    for j=0:20
%         data3 = load(['C:\Users\pivv\Documents\NCIA\hyeonuk_Rocket\level_set_algorithm\level_set_algorithm\tree_level' num2str(j) '.dat']);
%         plot(data3(1:1+4,1), data3(1:1+4,2));
%         axis([0 1 0 1]); axis square;
%         hold on;
%         for i=6:5:size(data3,1)
%            plot(data3(i:i+4,1), data3(i:i+4,2));
%         end
        
        data4 = load(['level_set_algorithm\level_set_algorithm\tree_3dlevel' num2str(j) '.dat']);
        data4 = reshape(data4,size+1,size+1,size+1);
        
        %hold on;
        %clf(figure(4));
        plot3([1;2],[1;2],[1;2]);
        figure(4), isosurface(l,lp,l,data4, 0), axis([0 1 0 1 0 1]);
        
        drawnow();
        %hold off;
        
    end

end