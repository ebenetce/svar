function summary = readBanbura2008Paper(pdfFile)
%READBANBURA2008PAPER Extract replication metadata from the local paper PDF.
%   SUMMARY = READBANBURA2008PAPER(PDFFILE) uses Text Analytics Toolbox
%   EXTRACTFILETEXT to read Banbura, Giannone and Reichlin (2008), "Large
%   Bayesian VARs", and returns the core metadata used by the replication
%   example.
%
%   See also EXTRACTFILETEXT, TOKENIZEDDOCUMENT.

arguments
    pdfFile (1,1) string {mustBeFile}
end

text = extractFileText(pdfFile);
abstract = extractBetween(text, "Abstract", "JEL Classification");
if isempty(abstract)
    abstract = "";
else
    abstract = localCleanWhitespace(abstract(1));
end

summary = struct();
summary.Title = "Large Bayesian VARs";
summary.Authors = ["Marta Banbura" "Domenico Giannone" "Lucrezia Reichlin"];
summary.WorkingPaper = "ECARES working paper 2008_033";
summary.NumCharacters = strlength(text);
summary.Abstract = abstract;
summary.SampleStart = datetime(1959,1,1);
summary.SampleEnd = datetime(2003,12,1);
summary.ForecastEvaluationStart = datetime(1971,1,1);
summary.ForecastEvaluationEnd = datetime(2003,12,1);
summary.StructuralSampleStart = datetime(1961,1,1);
summary.StructuralSampleEnd = datetime(2002,12,1);
summary.NumLags = 13;
summary.ForecastHorizons = [1 3 6 12];
summary.StructuralHorizons = 0:48;
summary.Text = text;

end

function cleaned = localCleanWhitespace(text)
cleaned = regexprep(text, "\s+", " ");
cleaned = strtrim(cleaned);
end
