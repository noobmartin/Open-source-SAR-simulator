/*
 * Author: Alexander Rajula
 * Contact: alexander@rajula.org
 *
 * This code is ugly and really ought to be re-written.
 */

#include "sar_simulator.h"

int write_real_data(data_arrays* data, radar_variables* variables){
  char fmode = 0;
  int ret = 0;

  do{
    fmode = getchar();
  }while(fmode != '\n');

  printf("Would you like to write data in human-readable or binary format (h/b): ");
  do{
    ret = scanf("%c", &fmode);
  }while((fmode != 'h') && (fmode != 'b'));

  FILE* dimensions = fopen("dimensions.dat", "w");
  if(dimensions == NULL){
    printf("Could not open dimensions.dat for writing - exiting.\n");
    return -1;
  }
  FILE* chirpf = fopen("chirp.dat", "wb");
  if(chirpf == NULL){
    printf("Could not open chirp.dat for writing - exiting.\n");
    return -1;
  }
  FILE* matchedf = fopen("matched.dat", "wb");
  if(matchedf == NULL){
    printf("Could not open matched.dat for writing - exiting.\n");
    return -1;
  }
  FILE* chirpfftf = fopen("chirpfft.dat", "wb");
  if(chirpfftf == NULL){
    printf("Could not open chirpfft.dat for writing - exiting.\n");
    return -1;
  }
  FILE* matchedfftf = fopen("matchedfft.dat", "wb");
  if(matchedfftf == NULL){
    printf("Could not open matchedfft.dat for writing - exiting.\n");
    return -1;
  }
  FILE* compressedf = fopen("compressed.dat", "wb");
  if(compressedf == NULL){
    printf("Could not open compressed.dat for writing - exiting.\n");
    return -1;
  }
  FILE* scene_with_waveformf = fopen("scene_with_waveform.dat", "wb");
  if(scene_with_waveformf == NULL){
    printf("Could not open scene_with_waveform.dat for writing - exiting.\n");
    return -1;
  }
  FILE* undistorted_radar_imagef = fopen("undistorted_radar_image.dat", "wb");
  if(undistorted_radar_imagef == NULL){
    printf("Could not open undistorted_radar_image.dat for writing - exiting.\n");
    return -1;
  }
  FILE* radar_imagef = fopen("radar_image.dat", "wb");
  if(radar_imagef == NULL){
    printf("Could not open radar_image.dat for writing - exiting.\n");
    return -1;
  }
  FILE* unfiltered_radar_imagef = fopen("unfiltered_radar_image.dat", "wb");
  if(unfiltered_radar_imagef == NULL){
    printf("Could not open unfiltered_radar_image.dat for writing - exiting.\n");
    return -1;
  }
  FILE* pulse_compressedf = fopen("pulse_compressed_image.dat", "wb");
  if(pulse_compressedf == NULL){
   printf("Could not open pulse_compressed_image.dat for writing - exiting.\n");
   return -1;
  }
  FILE* sar_imagef = fopen("sar_image.dat", "wb");
  if(sar_imagef == NULL){
    printf("Could not open sar_image.dat for writing - exiting.\n");
    return -1;
  }
  FILE* sar_fftf = fopen("sar_fft.dat", "wb");
  if(sar_fftf == NULL){
    printf("Could not open sar_fft.dat for writing - exiting.\n");
    return -1;
  }

  fprintf(dimensions, "%u\n%u\n%u\n%f\n", variables->chirp_length, variables->nrows, variables->ncols, variables->signal_distance);
  
  if(fmode == 'b'){
  	if(variables->mode != 'p'){
	    if(data->chirp_signal)
		ret = fwrite(data->chirp_signal, 1, variables->chirp_length*sizeof(complex double), chirpf);
	    if(data->matched_chirp)
		ret = fwrite(data->matched_chirp, 1, variables->chirp_length*sizeof(complex double), matchedf);
	    if(data->chirp_fft)
	    	ret = fwrite(data->chirp_fft, 1, variables->chirp_length*sizeof(complex double), chirpfftf);
	    if(data->matched_fft)
	    	ret = fwrite(data->matched_fft, 1, variables->chirp_length*sizeof(complex double), matchedfftf);
	    if(data->pulse_compressed_waveform)
	    	ret = fwrite(data->pulse_compressed_waveform, 1, variables->chirp_length*sizeof(complex double), compressedf);
	    if(data->scene_with_waveform)
	    	ret = fwrite(data->scene_with_waveform, 1, variables->nrows*variables->ncols*sizeof(complex double), scene_with_waveformf);
	}
	    if(data->radar_image)
	    	ret = fwrite(data->radar_image, 1, variables->nrows*variables->ncols*sizeof(complex double), radar_imagef);
	    if(data->undistorted_radar_image)
	    	ret = fwrite(data->undistorted_radar_image, 1, variables->nrows*variables->ncols*sizeof(complex double), undistorted_radar_imagef);
	    if(data->unfiltered_radar_image)
	    	ret = fwrite(data->unfiltered_radar_image, 1, variables->nrows*variables->ncols*sizeof(complex double), unfiltered_radar_imagef);
	    if(data->pulse_compressed_radar_image)
	   	ret = fwrite(data->pulse_compressed_radar_image, 1, variables->nrows*variables->ncols*sizeof(complex double), pulse_compressedf);
	    if(data->sar_image)
	    	ret = fwrite(data->sar_image, 1, variables->nrows*variables->ncols*sizeof(complex double), sar_imagef);
	    if(data->sar_fft)
	    	ret = fwrite(data->sar_fft, 1, variables->nrows*variables->ncols*sizeof(double complex), sar_fftf);
  }
  else{
	  int i,j;
  	  if(variables->mode != 'p'){
	  	if(data->chirp_signal){
		  for(i = 0; i < variables->chirp_length; i++){
		    fprintf(chirpf, "%.5f\n", creal(data->chirp_signal[i]));
		  }
		}

		if(data->matched_chirp){
		  for(i = 0; i < variables->chirp_length; i++){
		    fprintf(matchedf, "%.5f\n", creal(data->matched_chirp[i]));
		  }
		}

		if(data->pulse_compressed_waveform){
		  for(i = 0; i < variables->chirp_length; i++){
		    fprintf(compressedf, "%.5f\n", creal(data->pulse_compressed_waveform[i]));
		  }
		}

		if(data->chirp_fft){
		  for(i = 0; i < variables->chirp_length; i++){
		    fprintf(chirpfftf, "%.5f\n", creal(data->chirp_fft[i]));
		  }
		}

		if(data->matched_fft){
		  for(i = 0; i < variables->chirp_length; i++){
		    fprintf(matchedfftf, "%.5f\n", creal(data->matched_fft[i]));
		  }
		}

		if(data->scene_with_waveform){
		  for(i = 0; i < variables->ncols; i++){
		    for(j = 0; j < variables->nrows; j++){
		      fprintf(scene_with_waveformf, "%.5f\n", creal(data->scene_with_waveform[i*variables->nrows+j]));
		    }
		    fprintf(scene_with_waveformf, "\n");
	  	  }
		}
	  }

	if(data->undistorted_radar_image){
	  for(i = 0; i < variables->ncols; i++){
	    for(j = 0; j < variables->nrows; j++){
	      fprintf(undistorted_radar_imagef, "%.5f\t", creal(data->undistorted_radar_image[i*variables->nrows+j]));
	    }
	    fprintf(undistorted_radar_imagef, "\n");
	  }
	}

	if(data->radar_image){
	  for(i = 0; i < variables->ncols; i++){
	    for(j = 0; j < variables->nrows; j++){
	      fprintf(radar_imagef, "%.5f\t", creal(data->radar_image[i*variables->nrows+j]));
	    }
	    fprintf(radar_imagef, "\n");
	  }
	}

	if(data->unfiltered_radar_image){
	  for(i = 0; i < variables->ncols; i++){
	    for(j = 0; j < variables->nrows; j++){
	      fprintf(unfiltered_radar_imagef, "%.5f\t", creal(data->unfiltered_radar_image[i*variables->nrows+j]));
	    }
	    fprintf(unfiltered_radar_imagef, "\n");
	  }
	}

	if(data->pulse_compressed_radar_image){
	  for(i = 0; i < variables->ncols; i++){
	    for(j = 0; j < variables->nrows; j++){
	      fprintf(pulse_compressedf, "%.5f\t", creal(data->pulse_compressed_radar_image[i*variables->nrows+j]));
	    }
	    fprintf(pulse_compressedf, "\n");
	  }
	}

	if(data->sar_image){
	  for(i = 0; i < variables->ncols; i++){
	    for(j = 0; j < variables->nrows; j++){
	      fprintf(sar_imagef, "%.5f\t", creal(data->sar_image[i*variables->nrows+j]));
	    }
	    fprintf(sar_imagef, "\n");
	  }
	}

	if(data->sar_fft){
	  for(i = 0; i < variables->ncols; i++){
	    for(j = 0; j < variables->nrows; j++){
	      fprintf(sar_fftf, "%.5f\t", creal(data->sar_fft[i*variables->nrows+j]));
	    }
	    fprintf(sar_fftf, "\n");
  	  }
	}
  }

  fclose(sar_fftf);
  fclose(pulse_compressedf);
  fclose(dimensions);
  fclose(chirpf);
  fclose(matchedf);
  fclose(chirpfftf);
  fclose(matchedfftf);
  fclose(compressedf);
  fclose(scene_with_waveformf);
  fclose(radar_imagef);

}

int write_complex_data(data_arrays* data, radar_variables* variables){
  char fmode = 0;
  int ret = 0;

  do{
    fmode = getchar();
  }while(fmode != '\n');

  printf("Would you like to write data in human-readable or binary format (h/b): ");
  do{
    ret = scanf("%c", &fmode);
  }while((fmode != 'h') && (fmode != 'b'));

  FILE* dimensions = fopen("dimensions.dat", "w");
  if(dimensions == NULL){
    printf("Could not open dimensions.dat for writing - exiting.\n");
    return -1;
  }
  FILE* chirpf = fopen("chirp.dat", "wb");
  if(chirpf == NULL){
    printf("Could not open chirp.dat for writing - exiting.\n");
    return -1;
  }
  FILE* matchedf = fopen("matched.dat", "wb");
  if(matchedf == NULL){
    printf("Could not open matched.dat for writing - exiting.\n");
    return -1;
  }
  FILE* chirpfftf = fopen("chirpfft.dat", "wb");
  if(chirpfftf == NULL){
    printf("Could not open chirpfft.dat for writing - exiting.\n");
    return -1;
  }
  FILE* matchedfftf = fopen("matchedfft.dat", "wb");
  if(matchedfftf == NULL){
    printf("Could not open matchedfft.dat for writing - exiting.\n");
    return -1;
  }
  FILE* compressedf = fopen("compressed.dat", "wb");
  if(compressedf == NULL){
    printf("Could not open compressed.dat for writing - exiting.\n");
    return -1;
  }
  FILE* scene_with_waveformf = fopen("scene_with_waveform.dat", "wb");
  if(scene_with_waveformf == NULL){
    printf("Could not open scene_with_waveform.dat for writing - exiting.\n");
    return -1;
  }
  FILE* undistorted_radar_imagef = fopen("undistorted_radar_image.dat", "wb");
  if(undistorted_radar_imagef == NULL){
    printf("Could not open undistorted_radar_image.dat for writing - exiting.\n");
    return -1;
  }
  FILE* radar_imagef = fopen("radar_image.dat", "wb");
  if(radar_imagef == NULL){
    printf("Could not open radar_image.dat for writing - exiting.\n");
    return -1;
  }
  FILE* unfiltered_radar_imagef = fopen("unfiltered_radar_image.dat", "wb");
  if(unfiltered_radar_imagef == NULL){
    printf("Could not open unfiltered_radar_image.dat for writing - exiting.\n");
    return -1;
  }
  FILE* pulse_compressedf = fopen("pulse_compressed_image.dat", "wb");
  if(pulse_compressedf == NULL){
   printf("Could not open pulse_compressed_image.dat for writing - exiting.\n");
   return -1;
  }
  FILE* sar_imagef = fopen("sar_image.dat", "wb");
  if(sar_imagef == NULL){
    printf("Could not open sar_image.dat for writing - exiting.\n");
    return -1;
  }
  FILE* sar_fftf = fopen("sar_fft.dat", "wb");
  if(sar_fftf == NULL){
    printf("Could not open sar_fft.dat for writing - exiting.\n");
    return -1;
  }

  fprintf(dimensions, "%u\n%u\n%u\n%f\n", variables->chirp_length, variables->nrows, variables->ncols, variables->signal_distance);
  
  if(fmode == 'b'){
  	if(data->chirp_signal)
	    ret = fwrite(data->chirp_signal, 1, variables->chirp_length*sizeof(complex double), chirpf);
	if(data->matched_chirp)
	    ret = fwrite(data->matched_chirp, 1, variables->chirp_length*sizeof(complex double), matchedf);
	if(data->chirp_fft)
	    ret = fwrite(data->chirp_fft, 1, variables->chirp_length*sizeof(complex double), chirpfftf);
	if(data->matched_fft)
	    ret = fwrite(data->matched_fft, 1, variables->chirp_length*sizeof(complex double), matchedfftf);
	if(data->pulse_compressed_waveform)
	    ret = fwrite(data->pulse_compressed_waveform, 1, variables->chirp_length*sizeof(complex double), compressedf);
	if(data->scene_with_waveform)
	    ret = fwrite(data->scene_with_waveform, 1, variables->nrows*variables->ncols*sizeof(complex double), scene_with_waveformf);
	if(data->undistorted_radar_image)
	    ret = fwrite(data->undistorted_radar_image, 1, variables->nrows*variables->ncols*sizeof(complex double), undistorted_radar_imagef);
	if(data->radar_image)
	    ret = fwrite(data->radar_image, 1, variables->nrows*variables->ncols*sizeof(complex double), radar_imagef);
	if(data->unfiltered_radar_image)
	    ret = fwrite(data->unfiltered_radar_image, 1, variables->nrows*variables->ncols*sizeof(complex double), unfiltered_radar_imagef);
	if(data->pulse_compressed_radar_image)
	    ret = fwrite(data->pulse_compressed_radar_image, 1, variables->nrows*variables->ncols*sizeof(complex double), pulse_compressedf);
	if(data->sar_image)
	    ret = fwrite(data->sar_image, 1, variables->nrows*variables->ncols*sizeof(complex double), sar_imagef);
	if(data->sar_fft)
	    ret = fwrite(data->sar_fft, 1, variables->nrows*variables->ncols*sizeof(double complex), sar_fftf);
  }
  else{
	  int i,j;
	if(data->chirp_signal){
	  for(i = 0; i < variables->chirp_length; i++){
	    fprintf(chirpf, "%.5f\t", creal(data->chirp_signal[i]));
	    fprintf(chirpf, "%.5f\n", cimag(data->chirp_signal[i]));
	  }
	}

	if(data->matched_chirp){
	  for(i = 0; i < variables->chirp_length; i++){
	    fprintf(matchedf, "%.5f\t", creal(data->matched_chirp[i]));
	    fprintf(matchedf, "%.5f\n", cimag(data->matched_chirp[i]));
	  }
	}

	if(data->pulse_compressed_waveform){
	  for(i = 0; i < variables->chirp_length; i++){
	    fprintf(compressedf, "%.5f\t", creal(data->pulse_compressed_waveform[i]));
	    fprintf(compressedf, "%.5f\n", cimag(data->pulse_compressed_waveform[i]));
	  }
	}

	if(data->chirp_fft){
	  for(i = 0; i < variables->chirp_length; i++){
	    fprintf(chirpfftf, "%.5f\t", creal(data->chirp_fft[i]));
	    fprintf(chirpfftf, "%.5f\n", cimag(data->chirp_fft[i]));
	  }
	}

	if(data->matched_fft){
	  for(i = 0; i < variables->chirp_length; i++){
	    fprintf(matchedfftf, "%.5f\t", creal(data->matched_fft[i]));
	    fprintf(matchedfftf, "%.5f\n", cimag(data->matched_fft[i]));
	  }
	}

	if(data->scene_with_waveform){
	  for(i = 0; i < variables->ncols; i++){
	    for(j = 0; j < variables->nrows; j++){
	      fprintf(scene_with_waveformf, "%.5f\t", creal(data->scene_with_waveform[i*variables->nrows+j]));
	      fprintf(scene_with_waveformf, "%.5f\t", cimag(data->scene_with_waveform[i*variables->nrows+j]));
	    }
	    fprintf(scene_with_waveformf, "\n");
	  }
	}

	if(data->undistorted_radar_image){
	  for(i = 0; i < variables->ncols; i++){
	    for(j = 0; j < variables->nrows; j++){
	      fprintf(undistorted_radar_imagef, "%.5f\t", creal(data->undistorted_radar_image[i*variables->nrows+j]));
	      fprintf(undistorted_radar_imagef, "%.5f\t", cimag(data->undistorted_radar_image[i*variables->nrows+j]));
	    }
	    fprintf(undistorted_radar_imagef, "\n");
	  }
	}

	if(data->radar_image){
	  for(i = 0; i < variables->ncols; i++){
	    for(j = 0; j < variables->nrows; j++){
	      fprintf(radar_imagef, "%.5f\t", creal(data->radar_image[i*variables->nrows+j]));
	      fprintf(radar_imagef, "%.5f\t", cimag(data->radar_image[i*variables->nrows+j]));
	    }
	    fprintf(radar_imagef, "\n");
	  }
	}

	if(data->unfiltered_radar_image){
	  for(i = 0; i < variables->ncols; i++){
	    for(j = 0; j < variables->nrows; j++){
	      fprintf(unfiltered_radar_imagef, "%.5f\t", creal(data->unfiltered_radar_image[i*variables->nrows+j]));
	      fprintf(unfiltered_radar_imagef, "%.5f\t", cimag(data->unfiltered_radar_image[i*variables->nrows+j]));
	    }
	    fprintf(unfiltered_radar_imagef, "\n");
	  }
	}

	if(data->pulse_compressed_radar_image){
	  for(i = 0; i < variables->ncols; i++){
	    for(j = 0; j < variables->nrows; j++){
	      fprintf(pulse_compressedf, "%.5f\t", creal(data->pulse_compressed_radar_image[i*variables->nrows+j]));
	      fprintf(pulse_compressedf, "%.5f\t", cimag(data->pulse_compressed_radar_image[i*variables->nrows+j]));
	    }
	    fprintf(pulse_compressedf, "\n");
	  }
	}

	if(data->sar_image){
	  for(i = 0; i < variables->ncols; i++){
	    for(j = 0; j < variables->nrows; j++){
	      fprintf(sar_imagef, "%.5f\t", creal(data->sar_image[i*variables->nrows+j]));
	      fprintf(sar_imagef, "%.5f\t", cimag(data->sar_image[i*variables->nrows+j]));
	    }
	    fprintf(sar_imagef, "\n");
	  }
	}

	if(data->sar_fft){
	  for(i = 0; i < variables->ncols; i++){
	    for(j = 0; j < variables->nrows; j++){
	      fprintf(sar_fftf, "%.5f\t", creal(data->sar_fft[i*variables->nrows+j]));
	      fprintf(sar_fftf, "%.5f\t", cimag(data->sar_fft[i*variables->nrows+j]));
	    }
	    fprintf(sar_fftf, "\n");
  	  }
	}
  }

  fclose(sar_fftf);
  fclose(pulse_compressedf);
  fclose(dimensions);
  fclose(chirpf);
  fclose(matchedf);
  fclose(chirpfftf);
  fclose(matchedfftf);
  fclose(compressedf);
  fclose(scene_with_waveformf);
  fclose(radar_imagef);
  fclose(sar_imagef);
}

int read_radar_file(data_arrays* data, radar_variables* variables){
  FILE* fp = fopen(variables->radar_data_filename, "r");
  if(fp == NULL)
    return -1;

  //fseek(fp, 0L, SEEK_END);
  //unsigned int file_size = ftell(fp);
  //fseek(fp, 0L, SEEK_SET);
  
  radar_metadata meta;
  
  int ret;
  printf("Radar rows: ");
  ret = scanf("%u", &meta.rows);

  printf("Radar cols: ");
  ret = scanf("%u", &meta.cols);

  meta.real_or_complex = 'r';

  data->radar_image = malloc(meta.rows*meta.cols*sizeof(complex double));

  //FILE* mp = fopen("radar_metadata", "r");
  //ret = fread(&meta, sizeof(radar_metadata), 1, mp);

  //fclose(mp);


  int i,j;
  if(meta.real_or_complex == 'r'){
      for(j = 0; j < meta.cols; j++){
    for(i = 0; i < meta.rows; i++){
        //ret = fread(data->radar_image + i*sizeof(complex double), sizeof(double), 1, fp);
        ret = fscanf(fp, "%lf", (double*)(  &data->radar_image[j*meta.rows+i]  ));
      }
    }
  }
  else if(meta.real_or_complex == 'c'){
      for(j = 0; j < meta.rows; j++){
    for(i = 0; i < meta.rows*meta.cols; i++){
      //ret = fread(data->radar_image + i*sizeof(complex double), sizeof(complex double), 1, fp);
        ret = fscanf(fp, "%lf", (double*)(  &data->radar_image[i*meta.rows+j]  ));
      }
    }
  }
  else{
    printf("Invalid data mode, should be real or complex, read %c.\n", meta.real_or_complex);
    fclose(fp);
    return -1;
  }

  variables->nrows = meta.rows;
  variables->ncols = meta.cols;

  fclose(fp);
  return 1;
}
