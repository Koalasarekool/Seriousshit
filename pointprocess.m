function Criticality = pointprocess( inputdata, varargin )

% inputdata is unthresholded data ranging from one to four dimesions  with
% the last dimension being temporal - i.e. t, (x,t), (x,y,t) or (x,y,z,t).

% e.g. pointprocess( inputdata )

%koala2

% Name-value pair input arguments:

% 'sdthresh': the threshold to be used in units of standard deviation
% along the temporal dimension of the data. sdthresh can be either positive
% or negative. Input must be scalar and numeric. Default = 2SD.

% e.g. pointprocess( inputdata, 'sdthresh', 3 )

% 'abs': if set to 'on' takes the absolute value of the input data, so that
% an input of sdthresh = 2 corresponds to thresholds of +-2SD. Input must
% be string, 'on', or 'off'. Default = 'off'.

% e.g. pointprocess( inputdata, 'sdthresh', 3, 'abs', 'on')

% Name-value pair output arguments:

% 'point': If set to 'block', calculates point process based all times at
% which the signal is above threshold. 
% If set to 'onset' calculates point process based on times at which the
% signal first crosses the threshold
% If set to 'peak' calculates point process based on times at which the
% signal reaches its peak value above threshold
% Default calculates all three versions

% ensure that last dimension is temporal for 1D case
if numel(size(inputdata)) == 2 & size(inputdata,2) == 1
    inputdata = inputdata';
end

p = inputParser;

paramName1 = 'sdthresh'; % standard deviation threshold
default1 = 1;
errorStr = 'Value must be scalar and numeric.';
validationFcn = @(x) assert(isnumeric(x) && isscalar(x) ...
    ,errorStr);
addParameter(p,paramName1,default1,validationFcn);

paramName2 = 'abs'; % take absolute value of dat
default2 = 'off';
errorStr = 'Value must be either "on" or "off"';
validationFcn = @(x) assert(max(strcmpi(x,{'on','off'})), errorStr);
addParameter(p,paramName2,default2,validationFcn);

paramName3 = 'point'; % version of point process detection
default3 = 'ALL';
errorStr = 'Value must be either "block", "onset" or "peak"';
validationFcn = @(x) assert(max(strcmpi(x,{'block','onset','peak'})), errorStr);
addParameter(p,paramName3,default3,validationFcn);

p.parse(varargin{:});

threshdata = zeros(size(inputdata));

zscoreddata = zscore(inputdata, 0, numel(size(inputdata)));

% threshold data (block)
% in your original code you do this in a loop over time,
% but that doesn't seem necessary?
if strcmp(p.Results.abs,'on')
    threshdata = abs(zscoreddata) < abs(p.Results.sdthresh);
else if p.Results.sdthresh < 0
        threshdata = zscoreddata < p.Results.sdthresh;
    else if p.Results.sdthresh > 0
            threshdata = zscoreddata > p.Results.sdthresh;
        end
    end
end

pointtype = p.Results.point;

switch pointtype
    
    case 'ALL' 

% point process (onset)
% is there a better way to generalize this to N dimensions?
if numel(size(zscoreddata)) == 2
    onset=cat(2, threshdata(:,1), and( threshdata(:,2:end), ...
        not(threshdata(:,1:(end-1)))));
else if numel(size(zscoreddata)) == 3
        onset=cat(3, threshdata(:,:,1), and( threshdata(:,:,2:end), ...
            not(threshdata(:,:,1:(end-1)))));
    else if numel(size(zscoreddata)) == 4
            onset=cat(4, threshdata(:,:,:,1), and( threshdata(:,:,:,2:end), ...
                not(threshdata(:,:,:,1:(end-1)))));
        end
    end
end

% point process (ending) - to use later in calculating peak deflection
% is there a better way to generalize this to N dimensions?
if numel(size(zscoreddata)) == 2
    ending=cat(2, and(threshdata(:,1:(end-1)), not(threshdata(:,2:end))), ...
        threshdata(:,end));
else if numel(size(zscoreddata)) == 3
        ending=cat(3, and(threshdata(:,:,1:(end-1)), not(threshdata(:,:,2:end))), ...
            threshdata(:,:,end));
    else if numel(size(zscoreddata)) == 4
            ending=cat(4, and(threshdata(:,:,:,1:(end-1)), not(threshdata(:,:,:,2:end))), ...
                threshdata(:,:,:,end));
        end
    end
end

% point process (peak) - is there a better way of doing this?
% is there a better way to generalize this to N dimensions?
peak  = zeros(size(inputdata));
if numel(size(zscoreddata)) == 2
    for j = 1:size(zscoreddata,1)
        Onset = find(onset(j,:) == 1 );
        Ending = find(ending(j,:) == 1 );
        for i = 1:numel(Onset)
            Range = zscoreddata(j,Onset(i):Ending(i));
            [~,index] = sort(Range,'descend');
            peakfinder(i) = Onset(i) + index(1) - 1;
        end
        peak(j,peakfinder) = 1;
    end
else if numel(size(zscoreddata)) == 3 
    for j = 1:size(zscoreddata,1)
        for k = 1:size(zscoreddata,2)
        Onset = find(onset(j,k,:) == 1 );
        Ending = find(ending(j,k,:) == 1 );
        for i = 1:numel(Onset)
            Range = zscoreddata(j,k,Onset(i):Ending(i));
            [~,index] = sort(Range,'descend');
            peakfinder(i) = Onset(i) + index(1) - 1;
        end
        peak(j,k,peakfinder) = 1;
        end
    end
  else if numel(size(zscoreddata)) == 4      
        for j = 1:size(zscoreddata,1)
        for k = 1:size(zscoreddata,2)
            for l = 1:size(zscoreddata,3)
        Onset = find(onset(j,k,l,:) == 1 );
        Ending = find(ending(j,k,l,:) == 1 );
        for i = 1:numel(Onset)
            Range = zscoreddata(j,k,l,Onset(i):Ending(i));
            [~,index] = sort(Range,'descend');
            peakfinder(i) = Onset(i) + index(1) - 1;
        end
        peak(j,k,l,peakfinder) = 1;
        end
        end 
        end
      end
    end
end
      
peak = logical(peak);

% populate Criticality object
Criticality.inputdata = inputdata;
Criticality.zscoreddata = zscoreddata;
Criticality.sdthresh = p.Results.sdthresh;
Criticality.abs = p.Results.abs;
Criticality.pointtype = p.Results.point;

Criticality.pointprocess.bin.block = threshdata;
Criticality.pointprocess.bin.onset = onset;
Criticality.pointprocess.bin.peak = peak;

Criticality.pointprocess.mag.zscored.block = threshdata .* zscoreddata;
Criticality.pointprocess.mag.zscored.onset = onset .* zscoreddata;
Criticality.pointprocess.mag.zscored.peak = peak .* zscoreddata;

Criticality.pointprocess.mag.input.block = threshdata .* inputdata;
Criticality.pointprocess.mag.input.onset = onset .* inputdata;
Criticality.pointprocess.mag.input.peak = peak .* inputdata;

case 'block'
    
    % populate Criticality object
Criticality.inputdata = inputdata;
Criticality.zscoreddata = zscoreddata;
Criticality.sdthresh = p.Results.sdthresh;
Criticality.abs = p.Results.abs;
Criticality.pointtype = p.Results.point;

Criticality.pointprocess.bin.block = threshdata;

Criticality.pointprocess.mag.zscored.block = threshdata .* zscoreddata;

Criticality.pointprocess.mag.input.block = threshdata .* inputdata;

case 'onset'
    
% point process (onset)
% is there a better way to generalize this to N dimensions?
if numel(size(zscoreddata)) == 2
    onset=cat(2, threshdata(:,1), and( threshdata(:,2:end), ...
        not(threshdata(:,1:(end-1)))));
else if numel(size(zscoreddata)) == 3
        onset=cat(3, threshdata(:,:,1), and( threshdata(:,:,2:end), ...
            not(threshdata(:,:,1:(end-1)))));
    else if numel(size(zscoreddata)) == 4
            onset=cat(4, threshdata(:,:,:,1), and( threshdata(:,:,:,2:end), ...
                not(threshdata(:,:,:,1:(end-1)))));
        end
    end
end

% populate Criticality object
Criticality.inputdata = inputdata;
Criticality.zscoreddata = zscoreddata;
Criticality.sdthresh = p.Results.sdthresh;
Criticality.abs = p.Results.abs;
Criticality.pointtype = p.Results.point;

Criticality.pointprocess.bin.onset = onset;

Criticality.pointprocess.mag.zscored.onset = onset .* zscoreddata;

Criticality.pointprocess.mag.input.onset = onset .* inputdata;

case 'peak'
    
   
% point process (onset)
% is there a better way to generalize this to N dimensions?
if numel(size(zscoreddata)) == 2
    onset=cat(2, threshdata(:,1), and( threshdata(:,2:end), ...
        not(threshdata(:,1:(end-1)))));
else if numel(size(zscoreddata)) == 3
        onset=cat(3, threshdata(:,:,1), and( threshdata(:,:,2:end), ...
            not(threshdata(:,:,1:(end-1)))));
    else if numel(size(zscoreddata)) == 4
            onset=cat(4, threshdata(:,:,:,1), and( threshdata(:,:,:,2:end), ...
                not(threshdata(:,:,:,1:(end-1)))));
        end
    end
end

% point process (ending) - to use later in calculating peak deflection
% is there a better way to generalize this to N dimensions?
if numel(size(zscoreddata)) == 2
    ending=cat(2, and(threshdata(:,1:(end-1)), not(threshdata(:,2:end))), ...
        threshdata(:,end));
else if numel(size(zscoreddata)) == 3
        ending=cat(3, and(threshdata(:,:,1:(end-1)), not(threshdata(:,:,2:end))), ...
            threshdata(:,:,end));
    else if numel(size(zscoreddata)) == 4
            ending=cat(4, and(threshdata(:,:,:,1:(end-1)), not(threshdata(:,:,:,2:end))), ...
                threshdata(:,:,:,end));
        end
    end
end

% point process (peak) - is there a better way of doing this?
% is there a better way to generalize this to N dimensions?
peak  = zeros(size(inputdata));
if numel(size(zscoreddata)) == 2
    for j = 1:size(zscoreddata,1)
        Onset = find(onset(j,:) == 1 );
        Ending = find(ending(j,:) == 1 );
        for i = 1:numel(Onset)
            Range = zscoreddata(j,Onset(i):Ending(i));
            [~,index] = sort(Range,'descend');
            peakfinder(i) = Onset(i) + index(1) - 1;
        end
        peak(j,peakfinder) = 1;
    end
else if numel(size(zscoreddata)) == 3 
    for j = 1:size(zscoreddata,1)
        for k = 1:size(zscoreddata,2)
        Onset = find(onset(j,k,:) == 1 );
        Ending = find(ending(j,k,:) == 1 );
        for i = 1:numel(Onset)
            Range = zscoreddata(j,k,Onset(i):Ending(i));
            [~,index] = sort(Range,'descend');
            peakfinder(i) = Onset(i) + index(1) - 1;
        end
        peak(j,k,peakfinder) = 1;
        end
    end
  else if numel(size(zscoreddata)) == 4      
        for j = 1:size(zscoreddata,1)
        for k = 1:size(zscoreddata,2)
            for l = 1:size(zscoreddata,3)
        Onset = find(onset(j,k,l,:) == 1 );
        Ending = find(ending(j,k,l,:) == 1 );
        for i = 1:numel(Onset)
            Range = zscoreddata(j,k,l,Onset(i):Ending(i));
            [~,index] = sort(Range,'descend');
            peakfinder(i) = Onset(i) + index(1) - 1;
        end
        peak(j,k,l,peakfinder) = 1;
        end
        end 
        end
      end
    end
end
      
peak = logical(peak);

% populate Criticality object
Criticality.inputdata = inputdata;
Criticality.zscoreddata = zscoreddata;
Criticality.sdthresh = p.Results.sdthresh;
Criticality.abs = p.Results.abs;
Criticality.pointtype = p.Results.point;

Criticality.pointprocess.bin.peak = peak;

Criticality.pointprocess.mag.zscored.peak = peak .* zscoreddata;

Criticality.pointprocess.mag.input.peak = peak .* inputdata; 

end
