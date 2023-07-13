classdef GifFactory
    % GifFactory - Makes it easier to produce animated GIFs from figures,
    % requires export_gif: https://github.com/altmany/export_fig

    properties (Access = protected)
        delaytime {mustBeScalarOrEmpty} = 1/15
        ditheroption = 'dither'
        loopcount {mustBeScalarOrEmpty} = Inf
        frame {mustBeA(frame,'handle')} = gcf
        resolution {mustBeScalarOrEmpty} = 300    
        transparent = false
        filename
    end

    properties (Access = private)
        firstframe = true
    end

    methods
        function obj = GifFactory(filename, varargin)
            % Check for an existing .gif file by the same name: 
            obj.filename = filename;

            if exist(obj.filename, 'file') == 2
                warning(['Overwriting file: ', obj.filename]);
            end              
         
            tmp = strcmpi(varargin,'DelayTime'); 
   
            if any(tmp) 
                obj.delaytime = varargin{find(tmp)+1}; 
            end
   
            if any(strcmpi(varargin,'nodither'))
                obj.ditheroption = 'nodither'; 
            end

            if any(strcmpi(varargin,'transparent'))
                obj.transparent = true; 
            end
   
            tmp = strcmpi(varargin,'LoopCount'); 

            if any(tmp) 
                obj.loopcount = varargin{find(tmp)+1};       
            end           

            tmp = strncmpi(varargin,'resolution',3); 
 
            if any(tmp) 
                obj.resolution = varargin{find(tmp)+1};                 
            end
        end

        function obj = update(obj, frame)
            args = [];

            if obj.transparent 
                args = [args,"-transparent"];
            end

            if isgraphics(frame, 'figure')
                args = [args, "-nocrop"];    
            end

            if obj.resolution == 1
                args = [args, "-native"];
            else
                args = [args, strjoin({'-r', num2str(obj.resolution)},'')];                
            end

            warning off export_fig:exportgraphics

            fig = export_fig(args{:});

            % Convert the frame to a colormap and corresponding indices: 
            [imind,cmap] = rgb2ind(fig,256,obj.ditheroption);    
            
            % Write the file     
            if obj.firstframe
                imwrite(imind,cmap,obj.filename,'gif','WriteMode','overwrite', ...
                    'LoopCount',obj.loopcount,'DelayTime',obj.delaytime)    

                obj.firstframe = false;
            else
                imwrite(imind,cmap,obj.filename,'gif','WriteMode','append', ...
                    'DelayTime',obj.delaytime)                                 
            end
        end
    end
end