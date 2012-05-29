function [processed_chirp_data, processed_matched_chirp, processed_pulse_compressed, processed_scene_data, processed_undistorted_image_data, processed_unfiltered_image_data, processed_image_data, processed_compressed_data, processed_sar_data, processed_sar_fft] = plot_simulation_data()
    
    stacksize('max');
    
    data_type = input("Is the data real or complex valued (r/c)? ", "string");
    
    dimensions = read("dimensions.dat", 4, 1);
    chirp_length = dimensions(1);
    radar_rows = dimensions(2);
    radar_cols = dimensions(3);
    signal_distance = dimensions(4);
    
    distance_step = signal_distance/chirp_length;
    end_row_distance = distance_step*radar_rows;
    end_col_distance = distance_step*radar_cols;
    
    printf("Chirp length: %i\n", chirp_length);
    printf("Radar Rows: %i Cols: %i\n", radar_rows, radar_cols);
    
    if data_type == "c"
        chirp_data = read("chirp.dat", chirp_length, 2);
        matched_chirp_data = read("matched.dat", chirp_length, 2);
        pulse_compressed_data = read("compressed.dat", chirp_length, 2);
        scene_with_waveform_data = read("scene_with_waveform.dat", radar_cols, 2*radar_rows);
        undistorted_radar_image_data = read("undistorted_radar_image.dat", radar_cols, 2*radar_rows);
        unfiltered_radar_image_data = read("unfiltered_radar_image.dat", radar_cols, 2*radar_rows);
        radar_image_data = read("radar_image.dat", radar_cols, 2*radar_rows);
        compressed_image_data = read("pulse_compressed_image.dat", radar_cols, 2*radar_rows);
        sar_image_data = read("sar_image.dat", radar_cols, 2*radar_rows);
        sar_fft_data = read("sar_fft.dat", radar_cols, 2*radar_rows);
        
        chirp_fft = read("chirpfft.dat", chirp_length, 2);
        matched_fft = read("matchedfft.dat", chirp_length, 2);
        
        processed_chirp_data = chirp_data(:,1)+%i*chirp_data(:,2);
        processed_matched_chirp = matched_chirp_data(:,1)+%i*matched_chirp_data(:,2);
        processed_pulse_compressed = pulse_compressed_data(:,1)+%i*pulse_compressed_data(:,2);
        processed_scene_data = zeros(radar_rows, radar_cols);
        processed_undistorted_image_data = zeros(radar_rows, radar_cols);
        processed_unfiltered_image_data = zeros(radar_rows, radar_cols);
        processed_image_data = zeros(radar_rows, radar_cols);
        processed_compressed_data = zeros(radar_rows, radar_cols);
        processed_sar_data = zeros(radar_rows, radar_cols);
        processed_sar_fft = zeros(radar_rows, radar_cols);
        
        complex_chirp_fft = chirp_fft(:,1)+%i*chirp_fft(:,2);
        complex_matched_fft = matched_fft(:,1)+%i*matched_fft(:,2);
        
        for col = [1:radar_cols]
            for row = [1:radar_rows]
                processed_scene_data(row,col) = scene_with_waveform_data(col, 2*row-1) + %i*scene_with_waveform_data(col, 2*row);
                processed_undistorted_image_data(row,col) =  undistorted_radar_image_data(col,2*row-1) + %i*undistorted_radar_image_data(col, 2*row);
                processed_unfiltered_image_data(row,col) =  unfiltered_radar_image_data(col,2*row-1) + %i*unfiltered_radar_image_data(col, 2*row);
                processed_image_data(row,col) =  radar_image_data(col,2*row-1) + %i*radar_image_data(col, 2*row);
                processed_compressed_data(row,col) = compressed_image_data(col, 2*row-1) + %i*compressed_image_data(col, 2*row);
                processed_sar_data(row,col) = sar_image_data(col,2*row-1) + %i*radar_image_data(col, 2*row);
                processed_sar_fft(row,col) = sar_fft_data(col, 2*row-1) + %i*sar_fft_data(col, 2*row);
            end
        end
    elseif data_type == "r"
        processed_chirp_data = read("chirp.dat", chirp_length, 1);
        processed_matched_chirp = read("matched.dat", chirp_length, 1);
        processed_pulse_compressed = read("compressed.dat", chirp_length, 1);
        processed_scene_data = read("scene_with_waveform.dat", radar_cols, radar_rows);
        processed_undistorted_image_data = read("undistorted_radar_image.dat", radar_cols, radar_rows);
        processed_unfiltered_image_data = read("unfiltered_radar_image.dat", radar_cols, radar_rows);
        processed_image_data = read("radar_image.dat", radar_cols, radar_rows);
        processed_compressed_data = read("pulse_compressed_image.dat", radar_cols, radar_rows);
        processed_sar_data = read("sar_image.dat", radar_cols, radar_rows);
        processed_sar_fft = read("sar_fft.dat", radar_cols, radar_rows);
        
        complex_chirp_fft = read("chirpfft.dat", chirp_length,1);
        complex_matched_fft = read("matchedfft.dat", chirp_length, 1);
    else
        printf("Invalid data type.");
        return
    end
    
    da = gda();
    da.x_label.text = "Time";
    da.y_label.text = "Amplitude";
    figure("Figure_name", "Chirp signal");
    plot(processed_chirp_data);
    
    da = gda();
    da.x_label.text = "Frequency";
    da.y_label.text = "Amplitude";
    figure("Figure_name", "Chirp FFT");
    plot(abs(complex_chirp_fft));
    
    da = gda();
    da.x_label.text = "Time";
    da.y_label.text = "Amplitude";
    figure("Figure_name", "Matched Chirp signal");
    plot(processed_matched_chirp);
    
    da = gda();
    da.x_label.text = "Frequency";
    da.y_label.text = "Amplitude";
    figure("Figure_name", "Matched Chirp FFT");
    plot(abs(complex_matched_fft));
    
    da = gda();
    da.x_label.text = "Time";
    da.y_label.text = "Amplitude";
    figure("Figure_name", "Pulse-compressed signal");
    plot(abs(processed_pulse_compressed));
    
    [x,y] = meshgrid([0:distance_step:end_row_distance-distance_step], [0:distance_step:end_col_distance-distance_step]);
    
    da = gda();
    da.y_label.text = "Range (m)";
    da.x_label.text = "Crossrange (m)";
    figure("Figure_name", "Scene with waveform");
    //plot3d([0:distance_step:end_row_distance-distance_step], [0:distance_step:end_col_distance-distance_step], 10000*abs(processed_scene_data));
    mesh(x,y,10000*abs(processed_scene_data));
//    contour([0:distance_step:end_row_distance-distance_step],[0:distance_step:end_col_distance-distance_step],abs(processed_scene_data),10);
    
    da = gda();
    da.y_label.text = "Range (m)";
    da.x_label.text = "Crossrange (m)";
    figure("Figure_name", "Undistorted radar image");
    //plot3d([0:distance_step:end_row_distance-distance_step], [0:distance_step:end_col_distance-distance_step], 10000*abs(processed_image_data));
    mesh(x,y,10000*abs(processed_undistorted_image_data));
//    contour([0:distance_step:end_row_distance-distance_step],[0:distance_step:end_col_distance-distance_step],abs(processed_undistorted_image_data),10);
    
    da = gda();
    da.y_label.text = "Range (m)";
    da.x_label.text = "Crossrange (m)";
    figure("Figure_name", "Unfiltered radar image");
    //plot3d([0:distance_step:end_row_distance-distance_step], [0:distance_step:end_col_distance-distance_step], 10000*abs(processed_image_data));
    mesh(x,y,10000*abs(processed_unfiltered_image_data));
    
    da = gda();
    da.y_label.text = "Range (m)";
    da.x_label.text = "Crossrange (m)";
    figure("Figure_name", "Radar image");
    //plot3d([0:distance_step:end_row_distance-distance_step], [0:distance_step:end_col_distance-distance_step], 10000*abs(processed_image_data));
    mesh(x,y,10000*abs(processed_image_data));
    
    da = gda();
    da.y_label.text = "Range (m)";
    da.x_label.text = "Crossrange (m)";
    figure("Figure_name", "Pulse-compressed Radar image");
    //plot3d([0:distance_step:end_row_distance-distance_step], [0:distance_step:end_col_distance-distance_step], 10000*abs(processed_compressed_data));
    mesh(x,y,10000*abs(processed_compressed_data));
    
    da = gda();
    da.y_label.text = "Range (m)";
    da.x_label.text = "Crossrange (m)";
    figure("Figure_name", "SAR image");
    //plot3d([0:distance_step:end_row_distance-distance_step], [0:distance_step:end_col_distance-distance_step], 10000*abs(processed_sar_data));
    mesh(x,y,10000*abs(processed_sar_data));
    
    da = gda();
    da.y_label.text = "Range frequency";
    da.x_label.text = "Crossrange frequency";
    figure("Figure_name", "SAR image FFT");
    mesh(abs(processed_sar_fft'));
    
endfunction
