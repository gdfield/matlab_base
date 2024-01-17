function plotEiWaveforms(ei,varargin)

nbin = length(varargin);
if mod(nbin,2)==1
    err = MException('MATLAB:InvArgIn','Unexpected number of arguments');
    throw(err);
end

threshold = 0 ; %-2; 
saveFile = false; 
redElectrode = [] ;
plotRange = 1:40;
if size(ei,2)<40;
    plotRange = 1:size(ei,2);
end
% Read the optional input arguments
for j=1:(nbin/2)
    if ~ischar(varargin{j*2-1})
        err = MException('MATLAB:InvArgIn',...
            'Unexpected additional property');
        throw(err);
    end
    
    switch lower(varargin{j*2-1})
        case 'threshold'
            threshold = varargin{j*2};
        case 'savefile'
            saveFile = varargin{j*2};
        case 'redelectrode'
            redElectrode = varargin{j*2};
        case 'plotrange'
            plotRange = varargin{j*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end
[xc,yc]=getElectrodeCoords512();
positions = [xc; yc]';
scaledPos = (positions - repmat(min(positions),512,1))./repmat(max(positions)*3,512,1)+0.2;


f=figure; set(gcf,'Position',[100         100        2000         1000]); 
for e = 1:size(ei,1)
    figure(f.Number);
    
    axes('Position',[scaledPos(e,1) scaledPos(e,2) .02 .02]);
    if  any(ei(e,plotRange)< threshold)
        if isempty(redElectrode)
            if any(ei(e,plotRange)<-15)
                plot(ei(e,plotRange),'b'); axis off; % title(['ekectrode ' num2str(e)]) %axis off;
            else
                figure(f.Number);
                plot(ei(e,plotRange),'k');  axis off;
                ylim([-15 10]);
            end
        else
            if ismember(e,redElectrode)
                 plot(ei(e,plotRange),'r'); axis off;
            elseif any(ei(e,plotRange)<-15)
                plot(ei(e,plotRange),'b'); axis off; % title(['ekectrode ' num2str(e)]) %axis off;
            else
                figure(f.Number);
                plot(ei(e,plotRange),'k');  axis off;
                ylim([-15 10]);
            end
        end
    else
        figure(f.Number);
        plot(zeros(1,plotRange(end)),'k'); axis off;
    end
%     ylim([-15 10]);
end

if saveFile
    [savename,savepath] = uiputfile('*.*','save images as');
    if savename % If user does not cancel selection
        savingName = [savepath savename];
        saveas(gcf, savingName,'epsc');
    end
end