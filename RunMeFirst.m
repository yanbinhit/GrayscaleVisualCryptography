prevDir = pwd;
[dir, dummy, dummy2] = fileparts(mfilename('fullpath'));
addpath(dir,'-begin');
cd(dir);