function util_matched_replace(dirstr, filetype, matchkey, replace, testmode) %#ok<INUSD>
%UTIL_MATCHED_REPLACE Scans through a directory, looking through all files
%	with a given extension, replacing the contents of n lines after
%	matching a string with another string.
%Usage example; util_matched_replace('./', '.m', '#key123foobar', 'newstr', 1)
%will scan through the current directory, looking through files with the
%extension '.m', looking for lines containing '#key123foobar' surrounded by
%whitespace, replacing the contents of the next line with 'newstr'.
%
%The testmode argument, if set, will not actually modify any files, only
%displaying the changes that would be made on the console.
%
%The key and replacement string must not be identical, and the key must not
%contain any whitespace.
%
%For example, given that there exists a file '1.m' with the following
%content;
%
% a = 1;
% b = 2;
% %#key123foobar match
% oldstr
% anotheroldstr
% %#key123foobarnomatch
% c = 3;
%
% After executing util_matched_replace('./', '.m', '#key123foobar', {'newstr', 'newstr2'})
% the file would then contain;

% a = 1;
% b = 2;
% %#key123foobar
% newstr
% newstr2
% %#key123foobarnomatch
% c = 3;
%

if exist('testmode', 'var')
	fprintf('Running in test mode, no files will actually be modified.\n');
end

%Convert non-cell replacement string to cell type variable
if(~iscellstr(replace))
	dummy = replace;
	replace = cell(1);
	replace{1} = dummy;
	clear dummy;
end

files = dir(dirstr);

for k = 1:numel(files)
	%Ignore all directory types
	if(files(k).isdir)
		continue;
	end
	%Check the file extension
	cand_file_ext = files(k).name(end - (numel(filetype) - 1) : end);
	if(~strcmp(cand_file_ext, filetype))
		continue;
	end
	
	%At this stage we know the k'th file is the correct type. Now read the
	%whole file, scan through it and search for the key.
	fid = fopen(files(k).name, 'r');
	if(fid < 0)
		%Couldn't open file =(
		fprintf('Warning: could not open file %s\n', files(k).name);
		continue;
	end
	
	%Read entire file to a cell array with one line per element
	lines = cell(0);
	n = 1;
	while(true)
		res = fgetl(fid);
		if(res == -1)
			break;
		end
		lines{n} = res;
		n = n + 1;
	end
	
	match = false;
	%Search through cell array for a line ending in matchkey
	for n=1:numel(lines)-1
		%Check if the line contains the string we're looking for
		linewords = strsplit(lines{n});
		if(~any(strcmp(linewords, matchkey)))
			continue;
		end
		
		fprintf('\nFound a match in ''%s'' at line %d; contains;\n''%s''\n', files(k).name, n, lines{n});
		%At this stage we know the n'th line contains the key, so replace
		%the next line(s) with the replacement line(s)
		for m=1:numel(replace)
			fprintf('Replacing ''%s'' with ''%s''\n', lines{n+m}, replace{m});
			lines{n+m} = replace{m};
		end
		fprintf('\n');
		match = true;
	end
	
	if(~match)
		fclose(fid);
		continue;
	end
	
	%At least one modification has been made, so close the file and re-open
	%it in write mode if not in test mode
	fclose(fid);
	if(~exist('testmode', 'var'))
		fopen(files(k).name, 'w');
		for n=1:numel(lines)
			fprintf(fid, '%s\n', lines{n});
		end
		fclose(fid);
	end
end

