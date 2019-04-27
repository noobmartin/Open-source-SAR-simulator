function plot_simulation_data()
  delimiter = '\t';

  % Open up and parse dimensions file.
  fp = fopen('../output/dimensions.dat');
  
  if(fp == -1)
      disp(strcat('Could not open dimensions.dat for reading - exiting.'));
      return;
  end
  
  while( ~feof(fp) )
      % Read file descriptor and size.
      name = fgetl(fp);
      rows = str2num(fgetl(fp));
      cols = str2num(fgetl(fp));
      
      % Ignore empty matrices.
      if(rows < 1)
          continue;
      end
      
      if(cols < 1)
          continue;
      end
      
      % Read file as a matrix.
      fname = strcat("../output/"+name, '.dat');
      rimg = dlmread(fname, delimiter);
            
      if(cols == 1)
          img = zeros(1,rows);
          
          rvect = 1:rows;
          
          img(1,rvect) = rimg(2*rvect-1)+rimg(2*rvect)*1i;
                    
          figure('Name', name, 'NumberTitle', 'Off');
          title(name);
          plot(abs(img));
          xlabel('Time samples');
          ylabel('Amplitude');

      elseif(cols >= 2)
          r = rows;
          c = cols;
          cols = r;
          rows = c;
          img = zeros(rows, cols);
          
          cvect = 1:cols;
          rvect = 1:rows;
          
          img(rvect, cvect) = rimg(rvect,2*cvect-1)+rimg(rvect,2*cvect)*1i;
          
          figure('Name', name, 'NumberTitle', 'Off');
          imagesc(cvect,rvect,abs(img));
          title(name);
          xlabel('Range (m)');
          ylabel('Azimuth (m)');
          
      end
      
  end
  
  fclose(fp);