function getRFfromMUA(cfg)
% get gpsth from MUA
screenSize = [1680 1050];
fixPoint = screenSize/2;
binSize = cfg.binSize;
edges = 0:binSize:3.1-binSize;
filename = cfg.filename;
savename = cfg.savename;
orient = cfg.orient;
rfDir = dir(fullfile(filename, '*chopped.muat'));
load(fullfile(filename, rfDir.name), '-mat');
taccept = data.trialinfo(:, 2)==0;  % get correct trials
caccept = cfg.caccept;
Cond = data.trialinfo(:, 3); % get all conditions
muat.data = data; clear data;
libCond = unique(Cond);
nCond = length(unique(Cond));
assert(nCond == length(orient), 'nCond does not match orientations')
for ch=(caccept) % For each channel
    allDev = []; GPSTH = {};
    for cn =1:nCond
        cnd = libCond(cn);
        % Get all trials for a condition
        trialsChosen = find(Cond==cnd & taccept==1);
        % Get the MUA
        cfg = [];
        cfg.edges = edges;
        % z-score
        [tmpGPSTH, devGPSTH] = getGPSTHfromMUA(cfg, muat.data, ch, trialsChosen);
        GPSTH{cnd} = zscore(tmpGPSTH); allDev = [allDev; devGPSTH];
    end
    [gMap, gLatency] = getResponseLatency(GPSTH, edges(1:end-1));
    rfMap = imresize(gMap, [norm(screenSize) norm(screenSize)]);
    mapSize = size(rfMap);
    cropSize = [mapSize/2-screenSize/2+1  mapSize/2+screenSize/2];
    rfMap = flipud(fliplr(rfMap(cropSize(2):cropSize(4), cropSize(1):cropSize(3))));
    zMap = zscore(rfMap); max(abs(zMap(:)));
    % fit gaussian
    GPSTH = vertcat(GPSTH{:});
    sh = figure; set(sh, 'Visible', 'off');
    for ii = 1:size(GPSTH, 1)
        d = zscore(GPSTH(ii, :));
        %         d(abs(d)<2) = 0;
        [maxVal, ind ]= max(d(40:end-40));
        ind = ind + 40;
        x = 1/size(d, 2):1/size(d, 2):(size(d, 2))/size(d,2);
        options = fitoptions('gauss1');
        options.startPoint = [d(ind) x(ind) 0.005];
        options.Lower = [1 0.25  5/size(d,2)];
        options.Upper = [10 0.75 100/size(d,2)];
                options.Robust = 'Off';
        dSmooth = d-mean(d);
        gaussFit{ii} = fit(x', dSmooth', 'gauss1', options);
        fitDist{ii} = gaussFit{ii}.a1*exp(-((x-gaussFit{ii}.b1)/gaussFit{ii}.c1).^2)...
            +mean(d);
        xCutOff{ii} = [gaussFit{ii}.b1-gaussFit{ii}.c1;gaussFit{ii}.b1+gaussFit{ii}.c1];
        %         xSTD{ii} = [-gaussFit{ii}.c1/sqrt(2)*1.345; gaussFit{ii}.c1/sqrt(2)*1.345];
        xSTD{ii} = [-gaussFit{ii}.c1/sqrt(2)*1.645; gaussFit{ii}.c1/sqrt(2)*1.645];
        tmp = fitDist{ii};
        [pkMax, pkInd] = max(tmp(:));
        xPeak{ii} = x(pkInd);
        yPeak{ii} = pkMax;
        xSel = xPeak{ii}+xSTD{ii}(1) <= x & x <= xPeak{ii}+xSTD{ii}(2);
        %         gX{ii} = (x(xSel)-0.5)*size(d,2);
        gVal{ii} = gaussFit{ii}.a1*exp(-((x(xSel)-gaussFit{ii}.b1)/gaussFit{ii}.c1).^2);
        gEdgeVal{ii} = [gVal{ii}(1), gVal{ii}(end)];
        subplot(4,4,ii);
        plot(x, d, '.', 'MarkerSize', 5); hold on;
        plot(x, fitDist{ii}); hold on;
        plot( xCutOff{ii},gaussFit{ii}.a1*exp(-((xCutOff{ii}-gaussFit{ii}.b1)/gaussFit{ii}.c1).^2)...
            +mean(d), 'r.'); hold on;
        plot( xPeak{ii} ,yPeak{ii} , 'r.'); hold on;
        plot(x(ind), d(ind), 'g');
        [mVal, mInd] = max(GPSTH(:, 40:end-40), [], 2);
        mInd = mInd + 40; mInt = round(mInd-degDistances(1)/8+1:mInd+degDistances(1)/8);
%         mInt = true(size(allDev, 2));
        fanoWindow = allDev(:, mInt);
        fano = abs(mean(fanoWindow.^2./GPSTH(:, mInt), 2));
%         fano = abs(mean(allDev/GPSTH, 2));
%         fano = abs(diag(allDev(:, mInd)).^2./mVal);
        title(sprintf('%.2f, %.2f, %.2f, %.2f', ...
            mean(allDev(ii, :)), max(GPSTH(ii,:)), mean(allDev(ii, :))/max(GPSTH(ii,:)), fano(ii) ));
    end
    print(strcat(savename,'_RFPSTHFit.png'), '-dpng', '-r300')
    close; h = figure; set(h, 'visible', 'off');
    imagesc(rfMap);
    print(strcat(savename, '_RFEllipseFit.png'), '-dpng', '-r300');
    close;
    %% RF Criteria which is yet to be tested
    if ~( mean(allDev(:)) < 0.95 && mean(max(GPSTH, [], 2)) > 3)
        continue
    end
%     todel = abs(zscore(fano))>1.98;
    todel = mean(allDev, 2)/max(GPSTH, [], 2) > .3;
    allPeak = ((vertcat(yPeak{:})));
    allCenter = (vertcat(xPeak{:})-0.5)*norm(screenSize);
    allCut = (horzcat(xSTD{:})')*norm(screenSize);
    
    allVal = cellfun(@(cIn) cIn/max([gVal{:}]), gVal, 'UniformOutput', false);
    % Find center
    rad2Zero = 2*pi-orient;
    if size(GPSTH, 1) == 16
        meanCent = (allCenter(1:8)-allCenter(9:end))/2;
        %         [meanDegX, meanDegY] = pol2cart(sum(allPeak'.*orient),min(allCut));
        pairXY = [meanCent(1:4) meanCent(5:8)];
        newX = zeros(4,1); newY = zeros(4,1);
        for i=1:4
            z = meanCent(i).*exp(1i*orient(i)) + meanCent(i+4)*exp(1i*orient(i+4));
            [x1,y1]=pol2cart(angle(z), abs(z));
            newX(i) = x1; newY(i) = y1;
        end
    else
        meanCent = (allCenter(1:4)-allCenter(5:end))/2;
        pairXY = [meanCent(1:2) meanCent(3:4)];
        newX = zeros(2,1); newY = zeros(2,1);
        for i=1:2
            z = meanCent(i).*exp(1i*orient(i)) + meanCent(i+2)*exp(1i*orient(i+2));
            [x1,y1]=pol2cart(angle(z), abs(z));
            newX(i) = x1; newY(i) = y1;
        end
    end
    convIm = conv2(rfMap,ones(20,20)/400,'same');
    sc = conv2(ones(size(rfMap)),ones(20,20)/400,'same');
    convIm(sc<0.99) = NaN;
    [val,indx] = nanmax(convIm(:));
    [suby,subx] = ind2sub(size(rfMap),indx);
    rng_x = fixPoint(1) + newX > subx-40 & fixPoint(1) + newX<subx+40;
    rng_y = fixPoint(2) - newY > suby-40 & fixPoint(2) - newY<suby+40;
    
    %         if size(thisPSTH.GPSTH, 1) == 16
    %             for ii = 1:16
    %                 xl(ii) = allCenter(ii)*cos(orient(ii));
    %                 yl(ii) = allCenter(ii)*(-sin(orient(ii)));
    %                 d = sqrt((fixPoint(1)-subx).^2 + (fixPoint(2)-suby).^2);
    %                 xl_o(ii) = cos(orient(ii))*d;
    %                 yl_o(ii) = sin(orient(ii))*d;
    %             end
    %         else
    %         end
    
    x1 = fixPoint(1) + mean(newX);
    y1 = fixPoint(2) - mean(newY);
    
    % replace this mean estimation if not enough points are around the
    % estimated maximum in the image.
    if sum(rng_x & rng_y)<2
        x1 = subx; y1 = suby;
        disp('center from 2dmap')
    end
    
    %     z = zscore(allCut,1);
    %     todel = todel | (any(abs(z)>1.98, 2)) ;
    %      % Find ellipse points
    
    allCut(:,1) = allCut(:,1);
    allCut(:,2) = allCut(:,2);
    allCut(todel,:) = 0;
    cx1 = []; cy1 = [];cx2 = []; cy2 = [];
    if size(GPSTH, 1) == 16
        cnt = 0;
        for ii = 1:8
            cnt = cnt + 1;
            w = sum(~todel([ii ii+8]));
            dx(cnt) = (allCut(ii,1).*cos(orient(ii)) - allCut(ii+8,2).*cos(orient(ii)))/(w+1e-14);
            dy(cnt) = (allCut(ii,1).*sin(orient(ii)) - allCut(ii+8,2).*sin(orient(ii)))/(w+1e-14);
            cx1(cnt) = x1 + (allCut(ii,1).*cos(orient(ii)) - allCut(ii+8,2).*cos(orient(ii)))/(w+1e-14);
            cy1(cnt) = y1 + (allCut(ii,1).*sin(orient(ii)) - allCut(ii+8,2).*sin(orient(ii)))/(w+1e-14);
            cx2(cnt) = x1 + (allCut(ii,2).*cos(orient(ii)) - allCut(ii+8,1).*cos(orient(ii)))/(w+1e-14);
            cy2(cnt) = y1 + (allCut(ii,2).*sin(orient(ii)) - allCut(ii+8,1).*sin(orient(ii)))/(w+1e-14);
        end
    else
        cnt = 0;
        for ii = 1:4
            w = sum(~todel(ii));
            cnt = cnt + 1;
            w = sum(~todel([ii ii+4]));
            dx(cnt) = (allCut(ii,1).*cos(orient(ii)) - allCut(ii+4,2).*cos(orient(ii)))/(w+1e-14);
            dy(cnt) = (allCut(ii,1).*sin(orient(ii)) - allCut(ii+4,2).*sin(orient(ii)))/(w+1e-14);
            cx1(cnt) = x1 + (allCut(ii,1).*cos(orient(ii)) - allCut(ii+4,2).*cos(orient(ii)))/(w+1e-14);
            cy1(cnt) = y1 + (allCut(ii,1).*sin(orient(ii)) - allCut(ii+4,2).*sin(orient(ii)))/(w+1e-14);
            cx2(cnt) = x1 + (allCut(ii,2).*cos(orient(ii)) - allCut(ii+4,1).*cos(orient(ii)))/(w+1e-14);
            cy2(cnt) = y1 + (allCut(ii,2).*sin(orient(ii)) - allCut(ii+4,1).*sin(orient(ii)))/(w+1e-14);
        end
    end
    tdel = abs(cx1-x1)<1e-1 | abs(cy1-y1)<1e-1;
    tdel2 = abs(cx2-x1)<1e-1 | abs(cy2-y1)<1e-1;
    cx1(tdel) = []; cy1(tdel) = [];
    cx2(tdel2) = []; cy2(tdel2) = [];
    h = figure; set(h, 'Visible', 'off');
    imagesc(rfMap); hold on;
    hold on, plot(cx1,cy1,'k.'), hold on, plot(x1,y1,'r.'), hold on, plot(cx2,cy2,'k.');
    
    try
%         [z, a, b, phi] = fitellipse([cx1 cx2;cy1 cy2], 'constraint', 'trace');
    [z, a, b, phi] = fitellipse([cx1 cx2;cy1 cy2], 'constraint', 'bookstein');
    catch
        [a, b, x0, y0, phi] = ellipse_fit([cx1(:) cx2(:)], [cy1(:) cy2(:)]);
        z = zeros(2,1);
        z(1) = x0; z(2) = y0;
    end
    x0 = z(1); y0 = z(2);
    [a, b, x0, y0, phi]
    hold on, ellipsedraw(a,b,x0,y0,phi,'k', [128 128], 0);
    print(strcat(savename, '_RFEllipseFit.png'), '-dpng', '-r300');
    RF.centerposx = x0;
    RF.centerposy = y0;
    RF.sigmaX = a;
    RF.sigmaY = b;
    RF.angle = phi;
    save(strcat(savename, '_RFFit.mat'), 'RF')
    %     zGPSTH = zscore(GPSTH);
    %     fname = fullfile(outputDirSes, 'data', sprintf('RF_isisc420a01_Ch_%02d_gpsth.mat', ch));
    %     save(fname, 'GPSTH');
    %     fname = fullfile(outputDirSes, sprintf('RF_isisc420a01_Ch_%02d_gpsth.png', ch));
    %     print(h, fname, '-dpng', '-r300');
end
end


function [gMap, gLatency] = getResponseLatency(psth, edges)
% look for the response latency in 40-80ms with 1ms resolution

orient = 0:22.5:360-22.5;
respLatency = 0.04:0.001:0.08;
stimDuration = 3;
thisPSTH = vertcat(psth{:});
thisPSTH = [thisPSTH(9:end, :); thisPSTH(1:8, :)];
allMap = []; allMax = [];
for rind =1:length(respLatency)
    rDel = respLatency(rind);
    tsel = rDel < edges & edges < rDel+stimDuration;
    delPSTH = thisPSTH(:, tsel);
    allMap{rind} = back_project(delPSTH, orient);
    
    % smooth the image by 0.25 degree
    convIm = conv2(allMap{rind}, ones(20,20)/400, 'same');
    sc = conv2(ones(size(allMap{rind})),ones(20,20)/400,'same');
    convIm(sc<0.99) = NaN;
    [val, indx] = nanmax(convIm(:));
    allMax = [allMax val];
end
[~, ind] = max(allMax);
gMap = allMap{ind}; gLatency = respLatency(ind);
end

function [GPSTH, devGPSTH] = getGPSTHfromMUA(cfg, MUA, ch, trialsChosen)
% smooth the MUA and correct for the response latency
edges = cfg.edges;
allMUA = zeros(length(trialsChosen), length(edges));
for ii=1:numel(trialsChosen)
    allMUA(ii, :) = hist(MUA.trial{ch, trialsChosen(ii)}, edges);
end
sigma = 2.5;
cutoff = ceil(3*sigma);
GaussianK = fspecial('gaussian',[1,2*cutoff+1],sigma); % 1D filter
% GPSTH = zscore(conv(mean(allMUA(:, 1:end-1)),GaussianK,'same'));

allGPSTH = [];
for tr=1:size(allMUA, 1)
    GPSTH = zscore(conv((allMUA(tr, 1:end-1)),GaussianK,'same'));
    allGPSTH = [allGPSTH; GPSTH];
end
clear GPSTH; GPSTH = mean(allGPSTH, 1); devGPSTH = std(allGPSTH, [], 1);
%Save them
% fname = fullfile(outputDirSes, 'data', sprintf('RF_isisc420a01_Ch_%02d_gpsth.mat', ch));
% save(fname, 'GPSTH');
% fname = fullfile(outputDirSes, sprintf('RF_isisc420a01_Ch_%02d_gpsth.png', ch));
% print(h, fname, '-dpng', '-r300');
end