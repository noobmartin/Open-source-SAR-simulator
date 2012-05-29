function plot_simulation_data()
  dimensions = dlmread('dimensions.dat');
  
  % Extract parameters from metadata file.
  chirp_length = dimensions(1);
  rows = dimensions(2);
  cols = dimensions(3);
  distance = dimensions(4);
  
  % Read all raw data from file. Note that they contain both real and
  % imaginary parts.
  rchirp = dlmread('chirp.dat', '\t');
  rchirpfft = dlmread('chirpfft.dat', '\t');
  rcompressed = dlmread('compressed.dat', '\t');
  rmatched = dlmread('matched.dat', '\t');
  rmatchedfft = dlmread('matchedfft.dat', '\t');
  rpcimage = dlmread('pulse_compressed_image.dat', '\t');
  rradarimage = dlmread('radar_image.dat', '\t');
  rsarimage = dlmread('sar_image.dat', '\t');
  rsarfft = dlmread('sar_fft.dat', '\t');
  rscene = dlmread('scene_with_waveform.dat', '\t');
  rundistortedimage = dlmread('undistorted_radar_image.dat', '\t');
  runfilteredimage = dlmread('unfiltered_radar_image.dat', '\t');
  
  chirp = zeros(chirp_length);
  chirpfft = zeros(chirp_length);
  compressed = zeros(chirp_length);
  matched = zeros(chirp_length);
  matchedfft = zeros(chirp_length);
  
  chirp = rchirp(:,1)+rchirp(:,2)*1i;
  chirpfft = rchirpfft(:,1)+rchirpfft(:,2)*1i;
  compressed = rcompressed(:,1)+rcompressed(:,2)*1i;
  matched = rmatched(:,1)+rmatched(:,2)*1i;
  matchedfft = rmatchedfft(:,1)+rmatchedfft(:,2)*1i;
  
  pcimage = zeros(rows,cols);
  radarimage = zeros(rows,cols);
  sarimage = zeros(rows,cols);
  sarfft = zeros(rows,cols);
  scene = zeros(rows,cols);
  undistortedimage = zeros(rows,cols);
  unfilteredimage = zeros(rows,cols);
  
  for j=1:rows,
     for k=1:cols,
         pcimage(j,k) = rpcimage(j,2*k-1)+rpcimage(j,2*k)*1i;
         radarimage(j,k) = rradarimage(j,2*k-1)+rradarimage(j,2*k)*1i;
         sarimage(j,k) = rsarimage(j,2*k-1)+rsarimage(j,2*k)*1i;
         sarfft(j,k) = rsarfft(j,2*k-1)+rsarfft(j,2*k)*1i;
         scene(j,k) = rscene(j,2*k-1)+rscene(j,2*k)*1i;
%          undistortedimage = rundistortedimage(j,2*k-1)+rundistortedimage(j,2*k)*1i;
%          unfilteredimage = runfilteredimage(j,2*k-1)+runfilteredimage(j,2*k)*1i;
     end
  end
  
  figure('Name', 'Chirp', 'NumberTitle', 'Off')
  plot(abs(chirp))
  
  figure('Name', 'Chirp FFT', 'NumberTitle', 'Off')
  plot(abs(chirpfft))
  
  figure('Name', 'Matched chirp', 'NumberTitle', 'Off')
  plot(abs(matched))
  
  figure('Name', 'Matched chirp FFT', 'NumberTitle', 'Off')
  plot(abs(matchedfft))
  
  figure('Name', 'Pulse-compressed waveform', 'NumberTitle', 'Off')
  plot(abs(compressed))

  column = pcimage(round(rows/2),:);
  figure('Name', 'Middle column', 'NumberTitle', 'Off');
  plot([1:rows], abs(column));
  
  distance_step = distance/chirp_length;
  end_row_distance = distance_step*rows;
  end_col_distance = distance_step*cols;
  
  x = [0:distance_step:end_row_distance-distance_step];
  y = [0:distance_step:end_col_distance-distance_step];
  
  figure('Name', 'Scene with waveform', 'NumberTitle', 'Off');
  imagesc(x,y,abs(scene))
  title('Scene')
  xlabel('Range (m)')
  ylabel('Azimuth (m)')
  
%   figure('Name', 'Undistorted radar image', 'NumberTitle', 'Off');
%   contour(x,y,abs(undistortedimage),100)
%   title('Undistorted radar image')
%   xlabel('Range (m)')
%   ylabel('Azimuth (m)')
  
%   figure('Name', 'Unfiltered radar image', 'NumberTitle', 'Off');
%   contour(x,y,abs(unfilteredimage),100)
%   title('Unfiltered radar image')
%   xlabel('Range (m)')
%   ylabel('Azimuth (m)')
  
  figure('Name', '2D FFT of radar image', 'NumberTitle', 'Off');
  imagesc(abs(fftshift(fft2(radarimage))));
  title('2D FFT of radar image');

  figure('Name', 'RFI filtered radar image', 'NumberTitle', 'Off');
  imagesc(x,y,abs(radarimage))
  title('Filtered radar image')
  xlabel('Range (m)')
  ylabel('Azimuth (m)')
  
  figure('Name', 'Pulse-compressed radar image', 'NumberTitle', 'Off');
  imagesc(x,y,abs(pcimage))
  title('Pulse-compressed radar image')
  xlabel('Range (m)')
  ylabel('Azimuth (m)')
  
  figure('Name', 'SAR image', 'NumberTitle', 'Off');
  imagesc(x,y,abs(sarimage))
  title('SAR image')
  xlabel('Range (m)')
  ylabel('Azimuth (m)')
  
  figure('Name', 'SAR FFT', 'NumberTitle', 'Off');
  imagesc(abs(sarfft))
  title('SAR image FFT')
  xlabel('Range frequency')
  ylabel('Azimuth frequency')