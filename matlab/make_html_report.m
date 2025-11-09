function make_html_report()
%MAKE_HTML_REPORT Run replication and publish HTML report to /docs/index.html.

% locate dirs
thisFile   = mfilename('fullpath');
matlabDir  = fileparts(thisFile);
rootDir    = fileparts(matlabDir);
scriptsDir = fullfile(matlabDir,'scripts');
funcDir    = fullfile(matlabDir,'functions');
addpath(scriptsDir, funcDir);

% 1. run the full replication to refresh results
cd(matlabDir);
main_fiscal;

% 2. publish replication_report.m into /docs/index.html
opts = struct;
opts.format        = 'html';
opts.outputDir     = fullfile(rootDir,'..','docs');
opts.showCode      = true;     % set false if you prefer code hidden
opts.evalCode      = false;    % evaluation happens inside report
opts.useNewFigure  = true;

if ~exist(opts.outputDir,'dir'); mkdir(opts.outputDir); end

htmlFile = publish('main_fiscal.m', opts);

target = fullfile(opts.outputDir,'index.html');
if ~strcmp(htmlFile, target)
    movefile(htmlFile, target, 'f');
end

fprintf('HTML report generated at %s\n', target);
end
