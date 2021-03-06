
allCfg = allCfg{1}
filename = allCfg.inputfile;

if strcmp(allCfg.name, 'hermes')
    muaName = dir(fullfile(filename, sprintf('*.muax')));
    load(fullfile(filename, muaName.name ), '-mat');
    label = data.label;
    label_tdt = data.label_tdt;
    rfName = dir(fullfile(filename, sprintf('*.RF')));
    load(fullfile(filename, rfName.name ), '-mat');
else
    rfName = dir(fullfile(filename, sprintf('*.RF')));
    load(fullfile(filename, rfName.name ), '-mat');
end

% Parameters
screenSize = [1680 1050];
fixPoint = screenSize/2;
orient = 0:45:360-45;
orient = orient/180*pi;

caccept = strncmp(label, 'V1', 2);
RFs = RF(caccept);


savename = allCfg.outputfile;

for ch=1:length(RFs)
    load(fullfile(savename,sprintf('rfmapping-bar_1_Ch_%2d_PSTHFit.mat', ch)), 'gaussFit');
%     sh = figure(); if allCfg.print; set(sh, 'visible', 'off'); end;
    for cnd=1:size(RFs(ch).SRF, 1) % For all conditions/orientations
        d = RFs(ch).SRF(cnd, :);
        if max(d(40:end-40))<1.6;
            RFs(ch).centerposx = screenSize(1)/2;
            RFs(ch).centerposy = screenSize(2)/2;
            RFs(ch).sigmaX = 10;
            RFs(ch).sigmaY = 10;
            RFs(ch).angle = 0;
            RFs(ch).label = label(ch);
            RFs(ch).label_tdt = label_tdt(ch);
            break;
        end;
        rfMap = imresize(RFs(ch).map, [screenSize(2) screenSize(1)]);
        x = 1/size(d, 2):1/size(d, 2):(size(d, 2))/size(d,2);
        % gaussian fit options
%         [maxVal, ind ]= max(d(40:end-40));
%         options = fitoptions('gauss1');
%         options.startPoint = [maxVal x(ind) 0.5];
%         options.Lower = [1 0.2 5/size(d,2)];
%         options.Upper = [10 0.8 0.5];
%         options.Robust = 'On';
%         % gaussian fit to zero mean
%         dSmooth = d-mean(d);
%         gaussFit{cnd} = fit(x', dSmooth', 'gauss1', options);
        fitDist{cnd} = gaussFit{cnd}.a1*exp(-((x-gaussFit{cnd}.b1)/gaussFit{cnd}.c1).^2)...
            +mean(d);
        xCutOff{cnd} = [gaussFit{cnd}.b1-gaussFit{cnd}.c1;gaussFit{cnd}.b1+gaussFit{cnd}.c1];
        xSTD{cnd} = [-gaussFit{cnd}.c1/sqrt(2)*1.645; gaussFit{cnd}.c1/sqrt(2)*1.645];
        tmp = fitDist{cnd};
        [pkMax, pkInd] = max(tmp(:));
        xPeak{cnd} = x(pkInd);
        yPeak{cnd} = pkMax;
        xSel = xPeak{cnd}+xSTD{cnd}(1) <= x & x <= xPeak{cnd}+xSTD{cnd}(2);
        gX{cnd} = (x(xSel))*size(d,2);
        gVal{cnd} = gaussFit{cnd}.a1*exp(-((x(xSel)-gaussFit{cnd}.b1)/gaussFit{cnd}.c1).^2);
        gEdgeVal{cnd} = [gVal{cnd}(1), gVal{cnd}(end)];
    end
    
    allPeak = ((vertcat(yPeak{:})));
    allCenter = (vertcat(xPeak{:})-0.5)*norm(screenSize);
    allCut = (horzcat(xSTD{:})')*norm(screenSize);
    
    allVal = cellfun(@(cIn) cIn/max([gVal{:}]), gVal, 'UniformOutput', false);
    
    %     % Find center
    meanCent = (allCenter(1:4)-allCenter(5:end))/2;
    %     [meanDegX, meanDegY] = pol2cart(sum(allPeak'.*orient),min(allCut));
    
    pairXY = [ meanCent(1:2) meanCent(3:4) ];
    newX = zeros(2,1); newY = zeros(2,1);
    for i=1:2
        z = meanCent(i).*exp(1i*orient(i)) + meanCent(i+2)*exp(1i*orient(i+2));
        [x1,y1]=pol2cart(angle(z), abs(z));
        newX(i) = x1; newY(i) = y1;
    end
    
    x1 = fixPoint(1) + mean(newX);
    y1 = fixPoint(2) - mean(newY);
    
    % smooth the image by 0.25 degree
    convIm = conv2(rfMap,ones(20,20)/400,'same');
    sc = conv2(ones(size(rfMap)),ones(20,20)/400,'same');
    convIm(sc<0.99) = NaN;
    [val,indx] = nanmax(convIm(:));
    [suby,subx] = ind2sub(size(rfMap),indx);
    rng_x = fixPoint(1) + newX > subx-40 & fixPoint(1) + newX<subx+40;
    rng_y = fixPoint(2) - newY > suby-40 & fixPoint(2) - newY<suby+40;
    if sum(rng_x & rng_y)<2
        x1 = subx; y1 = suby;
        disp('center from 2dmap')
    end
    
    z = zscore(allCut,1);
    todel = (any(abs(z)>1.98,2)) ;
    %      % Find ellipse points
    allCut(todel,:) = 0;
    cx1 = []; cy1 = [];cx2 = []; cy2 = [];
    
    cnt = 0;
    for ii = 1:4
        w = sum(~todel(ii));
        if w
            cnt = cnt + 1;
            w = sum(~todel([ii ii+4]));
            dx(cnt) = (allCut(ii,1).*cos(orient(ii)) - allCut(ii+4,2).*cos(orient(ii)))/w;
            dy(cnt) = (allCut(ii,1).*sin(orient(ii)) - allCut(ii+4,2).*sin(orient(ii)))/w;
            cx1(cnt) = x1 + (allCut(ii,1).*cos(orient(ii)) - allCut(ii+4,2).*cos(orient(ii)))/w;
            cy1(cnt) = y1 + (allCut(ii,1).*sin(orient(ii)) - allCut(ii+4,2).*sin(orient(ii)))/w;
            cx2(cnt) = x1 + (allCut(ii,2).*cos(orient(ii)) - allCut(ii+4,1).*cos(orient(ii)))/w;
            cy2(cnt) = y1 + (allCut(ii,2).*sin(orient(ii)) - allCut(ii+4,1).*sin(orient(ii)))/w;
        end
    end
    
    sh = figure(); set(sh, 'Visible', 'off');
    imagesc(rfMap); hold on;
    hold on, plot(cx1,cy1,'k.'), hold on, plot(x1,y1,'r.'), hold on, plot(cx2,cy2,'k.');
    %     [a,b, x0, y0, phi] = ellipse_fit([cx1(:) cx2(:)], [cy1(:) cy2(:)]);
    try
        [z, a, b, phi] = fitellipse([cx1 cx2;cy1 cy2], 'constraint', 'trace');
    catch
        [z(1), z(2), a, b, phi] = deal(x0, y0, 30, 30, 0);
    end
    x0 = z(1); y0 = z(2);
    hold on, ellipsedraw(a,b,x0,y0,phi,'k', [screenSize(1) screenSize(2)], 0);
    print(sh, fullfile(savename,sprintf('rfmapping-bar_1_Ch_%02d_RFEllipseFit.png', ch)), '-dpng', '-r300');
    RFs(ch).centerposx = x0;
    RFs(ch).centerposy = y0;
    RFs(ch).sigmaX = a;
    RFs(ch).sigmaY = b;
    RFs(ch).angle = phi;
    RFs(ch).label = label(ch);
    RFs(ch).label_tdt = label_tdt(ch);
    % RFs(ch).label = RFs(ch).label;
    clf(sh);close all;
    save(fullfile(savename,sprintf('GRFs_Fit_%02d.mat', ch)), 'RFs')
end