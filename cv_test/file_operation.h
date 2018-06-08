#ifndef _FILE_OPERATION_H_
#define _FILE_OPERATION_H_
#include <stdio.h>
#include <stdint.h>


#ifdef  __cplusplus  
extern "C" {
#endif 
void read_data_from_file_f32(FILE *fp, float *data, int m, int n);
void read_data_from_file_u8(FILE *fp, uint8_t *data, int m, int n);
void write_data_to_file_f32(FILE *fp, float *data, int m, int n);
void write_data_to_file_u8(FILE *fp, uint8_t *data, int m, int n);
void write_data_to_file_s16(FILE *fp, int16_t *data, int m, int n);

#ifdef  __cplusplus  
}
#endif  


/*
void read_data_from_file_float(FILE *fp, float *data, int m, int n);
void read_data_from_file_uint8(FILE *fp, uint8_t *data, int m, int n);
void write_data_to_file_float(FILE *fp, float *data, int m, int n);
void write_data_to_file_uint8(FILE *fp, uint8_t *data, int m, int n);
void write_data_to_file_int16(FILE *fp, int16_t *data, int m, int n);
*/





void test_file_read_float();
void test_read_file_u8();

#endif
