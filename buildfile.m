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