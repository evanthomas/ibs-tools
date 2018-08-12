%
% The export to excel gui
%

% Set up initial default values
try, sheetname; catch, sheetname='forsman.xls'; end
try, start_col; catch, start_col=1; end
try, start_row; catch, start_row=1; end


  % Poistion [left bottom width height]
  try, close(1), end
  figure(1)
  set(gcf,'Name','Export to excel', ...
	  'NumberTitle','off', ...
	  'toolbar', 'none')
  % Suppress axes
  set(gca, 'TickLength', [0 0], 'position', [0 0 1 1])

  % Location of gadget columns
  col1      = 0.1;
  col2      = 0.45;
  lpos      = 0.6;
  ldelta    = 0.05;
  gheight   = 0.045;
  
  h = text(col1, lpos+0.02, 'spreadsheet name: ');
  set(h, 'fontweight', 'bold');
  sheet_bt = uicontrol('style', 'edit', ...
		       'Units', 'Normalized', ...
		       'Position', [col2 lpos .25 gheight], ...
		       'String', sheetname);
  lpos = lpos - ldelta;
  
  h = text(col1, lpos+0.02, 'start column: ');
  set(h, 'fontweight', 'bold');
  col_bt = uicontrol('style', 'edit', ...
		     'Units', 'Normalized', ...
		     'Position', [col2 lpos .07 gheight], ...
		     'String', num2str(start_col));
  lpos = lpos - ldelta;
  
  h = text(col1, lpos+0.02, 'start row: ');
  set(h, 'fontweight', 'bold');
  row_bt = uicontrol('style', 'edit', ...
		     'Units', 'Normalized', ...
		     'Position', [col2 lpos .07 gheight], ...
		     'String', num2str(start_row));
  lpos = lpos - ldelta;

  if start_col==1, v=1; else v=0; end
  h = text(col1, lpos+0.02, 'export variable descriptions: ');
  set(h, 'fontweight', 'bold');
  desc_bt = uicontrol('style', 'checkbox', ...
		      'Units', 'Normalized', ...
		      'Position', [col2 lpos .03 gheight], ...
		      'Value', v);
  lpos = lpos - ldelta;
  
  cancel_bt = uicontrol('Style', 'Push', ...
		       'Units', 'Normalized',...
		       'Position', [0.3 .0 .08 .045], ...
		       'String', 'cancel', ...
		       'callback', 'uiresume');
  export_bt = uicontrol('Style', 'Push', ...
			'Units', 'Normalized',...
			'Position', [0.4 .0 .08 .045], ...
			'String', 'export', ...
			'callback', 'uiresume');

  
  uiwait
                
  if gco==cancel_bt, close(1), return, end
  
  sheetname = get(sheet_bt, 'string');
  exptdesc  = get(desc_bt, 'Value');
  start_row = str2num(get(row_bt, 'String'));
  start_col = str2num(get(col_bt, 'String'));
  
  xl = ddeinit('excel', sheetname);
  
  if ~xl
    disp('Unable to make connection to excel.')
    disp(['Note: Excel needs to be running with the named spreadsheet' ...
	  ' loaded'])
    close(1)
    return
  end
  
  % suppress a warning that I have no control over
  % warning off MATLAB:DeprecatedLogicalAPI
  
  i = 1;
  if exptdesc
    item = sprintf('r%dc%d:r%dc%d', ...
		   start_row+i-1, start_col, start_row+i-1, start_col);
    ddepoke(xl, item, 'Experiment name (pulse no)');
    i = i + 1;
    
    item = sprintf('r%dc%d:r%dc%d', ...
		   start_row+i-1, start_col, start_row+i-1, start_col);
    ddepoke(xl, item, 'Mean press');
    i = i + 1;
    
    item = sprintf('r%dc%d:r%dc%d', ...
	   start_row+i-1, start_col, start_row+i-1, start_col);
    ddepoke(xl, item, 'PD at start');
    i = i + 1;
    
    item = sprintf('r%dc%d:r%dc%d', ...
	   start_row+i-1, start_col, start_row+i-1, start_col);
    ddepoke(xl, item, 'PD at end');
    i = i + 1;
    
    item = sprintf('r%dc%d:r%dc%d', ...
	   start_row+i-1, start_col, start_row+i-1, start_col);
    ddepoke(xl, item, 'Maximum at end');
    i = i + 1;
    
    item = sprintf('r%dc%d:r%dc%d', ...
	   start_row+i-1, start_col, start_row+i-1, start_col);
    ddepoke(xl, item, 'PD rebound');
    i = i + 1;
    
    item = sprintf('r%dc%d:r%dc%d', ...
	   start_row+i-1, start_col, start_row+i-1, start_col);
    ddepoke(xl, item, 'Rate of rise');
    i = i + 1;
    
    item = sprintf('r%dc%d:r%dc%d', ...
	   start_row+i-1, start_col, start_row+i-1, start_col);
    ddepoke(xl, item, 'Decay halflife');
    i = i + 1;
    
    start_col = start_col + 1;
  end
  
  
  trigg=0;
  for i=round(ups/4),
    trigg=trigg+1;
    i=1;
    item = sprintf('r%dc%d:r%dc%d', ...
		   start_row+i-1, start_col, start_row+i-1, start_col);
    ddepoke(xl, item, sprintf('%s (%d)', namn, trigg));
    i = i + 1;

    item = sprintf('r%dc%d:r%dc%d', ...
		   start_row+i-1, start_col, start_row+i-1, start_col);
    ddepoke(xl, item, ResF(trigg,PRESSURE));
    i = i + 1;
    
    if ~isnan(ResF(trigg,LOW))
      item = sprintf('r%dc%d:r%dc%d', ...
	     start_row+i-1, start_col, start_row+i-1, start_col);
      ddepoke(xl, item, round(10*ResF(trigg,LOW))/10);
    end
    i = i + 1;
    
    if ~isnan(ResF(trigg,HIGH2))
      item = sprintf('r%dc%d:r%dc%d', ...
	     start_row+i-1, start_col, start_row+i-1, start_col);
      ddepoke(xl, item, ResF(trigg,HIGH2));
    end
      i = i + 1;
    
    if ~isnan(ResF(trigg,HIGH))
      item = sprintf('r%dc%d:r%dc%d', ...
	     start_row+i-1, start_col, start_row+i-1, start_col);
      ddepoke(xl, item, round(10*ResF(trigg,HIGH))/10);
    end
    i = i + 1;
    
    if ~isnan(ResF(trigg, PD_REBOUND_X))
      item = sprintf('r%dc%d:r%dc%d', ...
	     start_row+i-1, start_col, start_row+i-1, start_col);
      ddepoke(xl, item, ResF(trigg,PD_REBOUND));
    end
    i = i + 1;
    
    if ~isnan(ResF(trigg,RISE_X_START))
      item = sprintf('r%dc%d:r%dc%d', ...
	     start_row+i-1, start_col, start_row+i-1, start_col);
      ddepoke(xl, item, ResF(trigg, SLOPE));
    end
    i = i + 1;
    
    if ~isnan(ResF(trigg, T_12))
      item = sprintf('r%dc%d:r%dc%d', ...
	     start_row+i-1, start_col, start_row+i-1, start_col);
      ddepoke(xl, item, ResF(trigg,T_12));
    end
    i = i + 1;
    
    start_col = start_col + 1;
  end    

  ddeterm(xl);
  
  % warning on MATLAB:DeprecatedLogicalAPI
  
  close(1)
