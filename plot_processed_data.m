function plot_processed_data()

    distance = 10;
    chirp_length = 5;
    
    rows = 512;
    cols = 317;
    
    cradar_image = dlmread('radar_image.dat', '\t');
    csar_image = dlmread('sar_image.dat', '\t');
    csar_fft = dlmread('sar_fft.dat', '\t');
    
    c = [1:2:cols-2];
    
    radar_image = cradar_image(:, c);
    sar_image = csar_image(:, c);
    sar_fft = csar_fft(:, c);
    
    distance_step = distance/chirp_length;
    end_row_distance = distance_step*rows;
    end_col_distance = distance_step*cols;
    
    x = [0:distance_step:end_row_distance-distance_step];
    y = [0:distance_step:end_col_distance-distance_step];
    
    figure('Name', 'Radar image', 'NumberTitle', 'Off');
    imagesc(x,y,abs(radar_image))
    title('Radar image')
    xlabel('Range (m)')
    ylabel('Azimuth (m)')
    
    figure('Name', 'SAR image', 'NumberTitle', 'Off');
    imagesc(x,y,abs(sar_image))
    title('SAR image')
    xlabel('Range (m)')
    ylabel('Azimuth (m)')
    
    figure('Name', 'SAR FFT image', 'NumberTitle', 'Off');
    imagesc(x,y,abs(sar_fft))
    title('Radar image')
    xlabel('Range (m)')
    ylabel('Azimuth (m)')