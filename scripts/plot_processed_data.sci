function plot_processed_data()
    
    //stacksize('max');
    
    //dimensions = read("dimensions.dat", 4, 1);
    //radar_rows = dimensions(2);
    //radar_cols = dimensions(3);
    
    radar_rows = 512;
    radar_cols = 315;

    printf("Radar Rows: %i Cols: %i\n", radar_rows, radar_cols);

        processed_image_data = read("radar_image.dat", radar_cols, radar_rows);
        //processed_compressed_data = read("pulse_compressed_image.dat", radar_cols, radar_rows);
        processed_sar_data = read("sar_image.dat", radar_cols, radar_rows);
        processed_sar_fft = read("sar_fft.dat", radar_cols, radar_rows);
    
    da = gda();
    da.x_label.text = "Range (m)";
    da.y_label.text = "Crossrange (m)";
    figure("Figure_name", "Radar image");
    //plot3d([0:distance_step:end_row_distance-distance_step], [0:distance_step:end_col_distance-distance_step], 10000*abs(processed_image_data));
    mesh(10000*abs(processed_image_data));
    
    //da = gda();
    //da.x_label.text = "Range (m)";
    //da.y_label.text = "Crossrange (m)";
    //figure("Figure_name", "Pulse-compressed Radar image");
    //plot3d([0:distance_step:end_row_distance-distance_step], [0:distance_step:end_col_distance-distance_step], 10000*abs(processed_compressed_data));
    //mesh(10000*abs(processed_compressed_data));
    
    da = gda();
    da.x_label.text = "Range (m)";
    da.y_label.text = "Crossrange (m)";
    figure("Figure_name", "SAR image");
    //plot3d([0:distance_step:end_row_distance-distance_step], [0:distance_step:end_col_distance-distance_step], 10000*abs(processed_sar_data));
    mesh(10000*abs(processed_sar_data));
    
    da = gda();
    da.x_label.text = "Range frequency";
    da.y_label.text = "Crossrange frequency";
    figure("Figure_name", "SAR image FFT");
    mesh(abs(processed_sar_fft'));
    
endfunction
