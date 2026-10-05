function plan = buildfile
import matlab.buildtool.tasks.CodeIssuesTask
import matlab.buildtool.tasks.TestTask

% Create a plan from task functions
plan = buildplan(localfunctions);

% Add the "check" task to identify code issues
plan("check") = CodeIssuesTask('tbx', InfoThreshold = 0, WarningThreshold=0);

% Add the "test" task to run tests
plan("test") = TestTask(SourceFiles = "tbx/svar", ...
    TestResults = "public/results.html", ...
    CodeCoverageResults = ["public/coverage.html", "public/coverage.xml"]);

% Make the "archive" task the default task in the plan
plan.DefaultTasks = "test";

end

function docTask(~)

if isempty(ver('docmaker'))
    websave('MATLAB_DocMaker.mltbx','https://github.com/mathworks/docmaker/releases/latest/download/MATLAB_DocMaker.mltbx');
    cobj = onCleanup(@() delete('MATLAB_DocMaker.mltbx'));
    matlab.addons.install('MATLAB_DocMaker.mltbx', true);
end

doc = fullfile( currentProject().RootFolder, "tbx", "doc" );

docdelete(doc)

md = fullfile(doc,"**","*.md"); % Markdown documents

html = docconvert(md, Scripts = fullfile(doc, 'mathjax-config.js')); % convert to HTML

docrun(html) % run code and insert output
docindex(doc); % index

end

function archiveTask(~)

v = ver('svar').Version;
opts = matlab.addons.toolbox.ToolboxOptions('tbx', "cb828c0f-df82-4306-a94a-5b09e773b60f");
opts.ToolboxName = "SVAR";
opts.PackageName = "SVAR";
opts.AuthorCompany = 'MathWorks';
opts.AuthorEmail = 'ebenetce@mathworks.com';
opts.AuthorName = 'Eduard Benet';
opts.Description = 'SVAR extension to Econometrics Toolbox.';
opts.OutputFile = "releases/svar.mltbx";
opts.Summary = 'SVAR extension to Econometrics Toolbox.';
% opts.ToolboxGettingStartedGuide = fullfile(currentProject().RootFolder,'tbx','doc','mfiles','GettingStarted.m');
opts.ToolboxVersion = v;
matlab.addons.toolbox.packageToolbox(opts)

end