% Set up python

% (set pykilosort python environment)
py_version = pyenv('Version','/home/julie/anaconda3/envs/pyks2/python.py');

% (add pykilosort environment paths to windows system path)
pre_pykilosort_syspath = getenv('PATH');
py_env_paths = {
    fullfile(char(py_version.Home),'Library','bin'); ...
    fullfile(char(py_version.Home),'Scripts')};

run_pykilosort_syspath = strjoin( ...
    unique(cat(1,py_env_paths, ...
    split(pre_pykilosort_syspath,pathsep)),'stable'),pathsep);

setenv('PATH',run_pykilosort_syspath);

% Run pykilosort
pyrunfile('YOUR PYKILOSORT PYTHON SCRIPT GOES HERE');

% (optional - after running pykilosort, revert system path to pre-pykilosort)
setenv('PATH',pre_pykilosort_syspath);