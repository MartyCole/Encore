
ico_grid = load('/nas/longleaf/home/mrcole/ABCD/SBCI_AVG/template_sphere_grid_ico4.mat');

lh_mesh.V = normr(ico_grid.lh_V);
lh_mesh.T = ico_grid.lh_T;
rh_mesh.V = normr(ico_grid.rh_V);
rh_mesh.T = ico_grid.rh_T;

lh_grid = SphericalGrid(lh_mesh,30);
rh_grid = SphericalGrid(rh_mesh,30);

P = length(lh_grid.V);

output_folder = "/work/users/m/r/mrcole/ScratchingTheSurface/new_reg/";
lh_displacement = zeros(P,1);
rh_displacement = zeros(P,1);

for i = [1:100,201:300]
    load(sprintf('%s/css_registered_warp%0.4i.mat', output_folder, i));
    lh_displacement = lh_displacement + sqrt((sum((100*(lh_grid.V - lh_warp.V)).^2,2)));
    rh_displacement = rh_displacement + sqrt((sum((100*(rh_grid.V - rh_warp.V)).^2,2)));
    disp(i)
end

lh_displacement = lh_displacement / 200;
rh_displacement = rh_displacement / 200;


%%

fig = figure('Units','normalized','OuterPosition',[0 .2 .6 .8]);
ax1 = subplot(2,2,1,'XLim',[-1,1],'YLim',[-1,1],'ZLim',[-1,1]);
ax2 = subplot(2,2,2,'XLim',[-1,1],'YLim',[-1,1],'ZLim',[-1,1]);
ax3 = subplot(2,2,3,'XLim',[-1,1],'YLim',[-1,1],'ZLim',[-1,1]);
ax4 = subplot(2,2,4,'XLim',[-1,1],'YLim',[-1,1],'ZLim',[-1,1]);

[lh_inflated,~] = read_vtk('/nas/longleaf/home/mrcole/ABCD/SBCI_AVG/lh_inflated_avg_ico4.vtk');

trimesh(lh_inflated.tri,lh_inflated.vtx(1,:),lh_inflated.vtx(2,:),lh_inflated.vtx(3,:),...
     'Parent',ax1,'EdgeColor','None','LineWidth',0.5,...
     'FaceLighting','phong','BackFaceLighting','unlit',...     
     'AmbientStrength',0.3,'DiffuseStrength',0.6,...
     'SpecularExponent',10,'SpecularStrength',0.9,...
     'FaceColor','Interp','FaceVertexCData',lh_displacement);

axis(ax1,'equal');
axis(ax1,'off');
clim(ax1,[0,8])
view(ax1,90,0);
colormap magma

trimesh(lh_inflated.tri,lh_inflated.vtx(1,:),lh_inflated.vtx(2,:),lh_inflated.vtx(3,:),...
     'Parent',ax3,'EdgeColor','None','LineWidth',0.5,...
     'FaceLighting','phong','BackFaceLighting','unlit',...     
     'AmbientStrength',0.3,'DiffuseStrength',0.6,...
     'SpecularExponent',10,'SpecularStrength',0.9,...
     'FaceColor','Interp','FaceVertexCData',lh_displacement);

axis(ax3,'equal');
axis(ax3,'off');
clim(ax3,[0,8])
view(ax3,-90,0);
colormap magma

rh_inflated = read_vtk('/nas/longleaf/home/mrcole/ABCD/SBCI_AVG/rh_inflated_avg_ico4.vtk');

trimesh(rh_inflated.tri,rh_inflated.vtx(1,:),rh_inflated.vtx(2,:),rh_inflated.vtx(3,:),...
     'Parent',ax2,'EdgeColor','None','LineWidth',0.5,...
     'FaceLighting','phong','BackFaceLighting','unlit',...     
     'AmbientStrength',0.3,'DiffuseStrength',0.6,...
     'SpecularExponent',10,'SpecularStrength',0.9,...
     'FaceColor','Interp','FaceVertexCData',rh_displacement);

axis(ax2,'equal');
axis(ax2,'off');
clim(ax2,[0,8])
view(ax2,-90,0);
colormap magma

trimesh(rh_inflated.tri,rh_inflated.vtx(1,:),rh_inflated.vtx(2,:),rh_inflated.vtx(3,:),...
     'Parent',ax4,'EdgeColor','None','LineWidth',0.5,...
     'FaceLighting','phong','BackFaceLighting','unlit',...     
     'AmbientStrength',0.3,'DiffuseStrength',0.6,...
     'SpecularExponent',10,'SpecularStrength',0.9,...
     'FaceColor','Interp','FaceVertexCData',rh_displacement);

axis(ax4,'equal');
axis(ax4,'off');
clim(ax4,[0,8])
view(ax4,90,0);
colormap magma

h = axes(fig,'visible','off'); 
h.Title.Visible = 'on';
h.XLabel.Visible = 'on';
h.YLabel.Visible = 'on';
ylabel(h,'','FontWeight','bold');
xlabel(h,'','FontWeight','bold');
title(h,'Mean Deformation');
c = colorbar(h,'Position',[0.93 0.168 0.022 0.7]);  % attach colorbar to h
colormap(c,'magma')
clim(h,[0,8]);             % set colorbar limits


%%

for i = 1:4
    if i == 1 || i == 3
        fig = figure('Units','normalized','OuterPosition',[0 .2 .6 .6]);
        ax1 = subplot(1,2,1,'XLim',[-1,1],'YLim',[-1,1],'ZLim',[-1,1]);
        ax2 = subplot(1,2,2,'XLim',[-1,1],'YLim',[-1,1],'ZLim',[-1,1]);
        ax3 = subplot(1,2,1,'XLim',[-1,1],'YLim',[-1,1],'ZLim',[-1,1]);
        ax4 = subplot(1,2,2,'XLim',[-1,1],'YLim',[-1,1],'ZLim',[-1,1]);
    end

    load(sprintf('%s/css_registered_warp%0.4i.mat', output_folder, i));

    trimesh(lh_grid.T,lh_grid.V(:,1)*0.99,lh_grid.V(:,2)*0.99,lh_grid.V(:,3)*0.99,...
        'Parent',eval(sprintf('ax%d',i)),'EdgeColor','None','LineWidth',0.5,...
        'FaceLighting','phong','BackFaceLighting','unlit',...
        'AmbientStrength',0.3,'DiffuseStrength',0.6,...
        'SpecularExponent',10,'SpecularStrength',0.9,...
        'FaceColor','Interp','FaceVertexCData',sqrt((sum((100*(lh_grid.V - lh_warp.V)).^2,2))));

    hold(eval(sprintf('ax%d',i)),'on')
    dir = lh_warp.V - lh_grid.V;
    quiver3(lh_grid.V(:,1),lh_grid.V(:,2),...
        lh_grid.V(:,3),dir(:,1),dir(:,2),dir(:,3),1.5, ...
        'linewidth',1.5,'Parent',eval(sprintf('ax%d',i)),'Color','white','AutoScaleFactor',10)
    hold(eval(sprintf('ax%d',i)),'off')

    axis(eval(sprintf('ax%d',i)),'equal');
    axis(eval(sprintf('ax%d',i)),'off');
    clim(eval(sprintf('ax%d',i)),[-15,15])
    view(eval(sprintf('ax%d',i)),90,0);
    colormap magma
end
