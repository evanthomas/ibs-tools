function [name, file] = filegui
% Widget to allow the user to select the file to load
  h = axes('position', [0 0 1 1]);
  set(h, 'TickLength', [0 0]);

  set(gcf,'Name','File details', ...
	  'NumberTitle','off', ...
	  'toolbar', 'none', ...
	  'menubar', 'none');

  text(0.2, 0.6, 'Specify file and name', ...
       'fontweight', 'bold', ...
       'fontsize', 14);
  
  file_txt = uicontrol('Style', 'Edit', ...
		       'Units', 'Normalized',...
		       'Position', [0.02 .5 .65 .045], ...
		       'callback', 'uiresume');
  name_txt = uicontrol('Style', 'Edit', ...
		       'Units', 'Normalized',...
		       'Position', [0.02 .4 .65 .045], ...
		       'callback', 'uiresume');
  file_bt = uicontrol('Style', 'Push', ...
		      'Units', 'Normalized',...
		      'Position', [0.8 .5 .12 .045], ...
		      'String', 'choose file', ...
		      'BackgroundColor', [1 1 1], ...
		      'callback', 'uiresume');
  uicontrol('Style', 'Text', ...
	    'Units', 'Normalized',...
	    'Position', [0.8 .4 .12 .045], ...
	    'String', 'name');
  next_bt = uicontrol('Style', 'Push', ...
		      'Units', 'Normalized',...
		      'Position', [0.6 .0 .05 .045], ...
		      'String', 'OK', ...
		      'enable', 'off', ...
		      'callback', 'uiresume');
  
  while 1
    uiwait
    switch gco
     case file_txt
      file = get(file_txt, 'String');
     
     case file_bt
      [file, path] = uigetfile('*.*', 'Select the data file');
      [p, name] = fileparts(file);
      cd(path)
      set(file_txt, 'String', file);
      set(name_txt, 'String', name)
      set(next_bt, 'enable', 'on');
      
     case next_bt
      file = get(file_txt, 'String');
      name = get(name_txt, 'String');
      close(gcf)
      return
      
    end
  end
  
