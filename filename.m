%%%%%%%%%%%
%   Janez Prešern, 2014
%   
%   Function allows loading/saving .mat & .other files independent of the
%   OS
%   It checks inputDir variable for the (back)slashes concatenetes the
%   components of path into correct form, using either slashes or
%   backslashes.

%   inputDir .. starting directory as red by parent function
%   path     .. vector of strings with all the components of the file path,
%               including inputDir, Project, Results in the correct order
%               as they should be: [inputDir,Project,Project]
%   fileName .. name of the file we wish to load/save

function fn = filename (inputDir,path,fName)
    fn= [];
    if strfind(inputDir,'\');
        for i = 1:length(path)
            fn = [fn,path{i},'\'];
        end;
        fn = [fn,fName];
    else
        for i = 1:length(path)
            fn = [fn,path{i},'/'];
        end;
        fn = [fn,fName];
    end;