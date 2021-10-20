function test_install(waiting_key_press)
% This script runs all the examples of the toolbox
clc;
if nargin < 1
    waiting_key_press = 0;
end
msg_length=90;
maketitle('Performance Estimation Toolbox (PESTO) -- TESTS: run all examples of the toolbox',msg_length,2)

rest_path = pwd;
[filepath,~,~] = fileparts(which('test_install.m'));
cd(filepath);

filePattern = fullfile('Examples', '**/*.m');
theFiles = dir(filePattern);
nbFiles  = max(size(theFiles));
skip_list = {'A_Minimizer_of_a_sum.m','C_DirectAccelerationSAGA.m'}; skip_reason = {'too long','too long'};
nb_skip = length(skip_list);

waitfor;
fprintf('(PESTO TESTS): %d examples to be run\n',nbFiles)
waitfor;
for i_outer = 1:nbFiles
    tt=sprintf('%s/%s',theFiles(i_outer).folder,theFiles(i_outer).name);
    skip_case = 0;
    for j_outer = 1:nb_skip
        if strcmp(skip_list{j_outer},theFiles(i_outer).name)
            skip_case = 1; skip_c = j_outer;
            break
        end
    end
    if ~skip_case
        fprintf('(PESTO TESTS): running ''%s'' [test %d on %d]\n',theFiles(i_outer).name,i_outer,nbFiles)
        run(tt)
        fprintf('(PESTO TESTS): ''%s'' finished [test %d on %d]\n',theFiles(i_outer).name,i_outer,nbFiles)
    else
        fprintf('(PESTO TESTS): skipping ''%s'' [test %d on %d], reason: %s\n',theFiles(i_outer).name,i_outer,nbFiles,skip_reason{skip_c});
    end
    if waiting_key_press
        waitfor;
    end
end


cd(rest_path);

end

function waitfor
waitmsg=sprintf('\n\nPress any key to continue\n');
waitlength = numel(waitmsg);
disp(waitmsg)
pause
fprintf(repmat('\b', 1, waitlength));
end
function maketitle(msg,width,nblines)
if nargin<3
    nblines=1;
end
fprintf('\n');
for i=1:nblines
    fprintf([repmat('*',1,width) '\n']);
end
msg_l=numel(msg);
left_margin=floor((width-msg_l)/2);
right_margin=width-msg_l-left_margin;
for i=1:nblines-1
    fprintf(['*' repmat(' ',1,width-2) '*\n']);
end
msg_mod=['*' repmat(' ',1,left_margin-1) msg repmat(' ',1,right_margin-1) '*\n'];
fprintf(msg_mod)
for i=1:nblines-1
    fprintf(['*' repmat(' ',1,width-2) '*\n']);
end
for i=1:nblines
    fprintf([repmat('*',1,width) '\n']);
end
fprintf('\n');
end