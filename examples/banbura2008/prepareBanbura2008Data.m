function data = prepareBanbura2008Data(sourceFile, outputFile)
%PREPAREBANBURA2008DATA Prepare the official JAE BGR data archive for MAIN.
%   DATA = PREPAREBANBURA2008DATA(SOURCEFILE, OUTPUTFILE) reads the official
%   Banbura-Giannone-Reichlin replication workbook, renames the key paper
%   variables, applies the paper's level/log convention, and saves
%   banbura2008Data.mat with TT and PRIORMEAN.
%
%   The official archive is hosted by the Journal of Applied Econometrics at
%   qed.econ.queensu.ca/jae/2010-v25.1/banbura-giannone-reichlin/.
%
%   See also READTIMETABLE, BANBURA2008MODELSPECIFICATIONS.

arguments
    sourceFile (1,1) string {mustBeFile}
    outputFile (1,1) string = fullfile(fileparts(sourceFile), "banbura2008Data.mat")
end

raw = readcell(sourceFile, Sheet="Data");
codes = string(raw(1, 2:end));
transformCodes = cell2mat(raw(3, 2:end));
dates = datetime([raw{4:end, 1}]');
values = cell2mat(raw(4:end, 2:end));

aliases = localAliasMap();
names = codes;
for j = 1:numel(codes)
    if isKey(aliases, codes(j))
        names(j) = aliases(codes(j));
    end
end
names = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(names));

Y = localTransformLevels(values, transformCodes);
tt = array2timetable(Y, RowTimes=dates, VariableNames=cellstr(names));

% Stock-Watson transformation codes 5 and 6 are nonstationary log-level
% variables in the paper's setup, so their Minnesota prior mean is 1.
% Stationary quantities and rates get a white-noise prior mean of 0.
priorMean = double(transformCodes == 5 | transformCodes == 6);

save(outputFile, "tt", "priorMean", "codes", "transformCodes");
data = struct(tt=tt, priorMean=priorMean, codes=codes, ...
    transformCodes=transformCodes, outputFile=outputFile);
end

function Y = localTransformLevels(values, transformCodes)
Y = values;
logColumns = transformCodes >= 4;
Y(:, logColumns) = log(Y(:, logColumns));
end

function aliases = localAliasMap()
aliases = dictionary( ...
    ["CES002" "PUNEW" "FYFF" "PSM99Q" "FMRNBA" "FMRRA" "FM2" ...
     "a0m052" "A0M224_R" "IPS10" "A0m082" "LHUR" "HSFR" "PWFSA" ...
     "GMDC" "CES275" "FM1" "FSPCOM" "FYGT10" "EXRUS"], ...
    ["EMPL" "CPI" "FFR" "COMM_PR" "NBORR_RES" "TOT_RES" "M2" ...
     "INCOME" "CONSUM" "IP" "CAP_UTIL" "UNEMPL" "HOUS_START" "PPI" ...
     "PCE_DEFL" "HOUR_EARN" "M1" "SP" "TB_YIELD" "EXR"]);
end
