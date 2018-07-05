#include "iip_wav.h"

WAV* read_WAV(char*file_path)
{
 WAV * wav;
 FILE * f;
 uint32_t read = 0;
 char buffer4[4];
 char buffer2[2];
 UINT num_of_sample;
 //파일 유무확인
 f = fopen(file_path,"r");
 if(f == NULL)
 {
   printf("Failed to open '%s'\n",file_path);
   return NULL;
 }
   
 //WAV 읽기
 wav = (WAV*)malloc(sizeof(WAV)); 
  
 fread(wav->riff_id,sizeof(wav->riff_id),1,f); 
 
 fread(buffer4,sizeof(buffer4),1,f);
 //convert little endial to big endian 4 bytes int;
 wav->riff_size = buffer[0]|(buffer4[1]<<8)|(buffer4[2]<<16)|(buffer4[3]<<24);
 
 fread(wav->wave_id,sizeof(wav->wave_id),1,f);
 
 fread(wav->fmt_id,sizeof(fmt_id),1,f);

 fread(buffer4,sizeof(buffer4),1,f);
 //convert little endial to big endian 4 bytes int;
 wav->fmt_size = buffer[0]|(buffer4[1]<<8)|(buffer4[2]<<16)|(buffer4[3]<<24);
 
 fread(buffer2,sizeof(buffer2),1,f);
 //convert little endial to big endian 2 bytes int;
 wav->fmt_type = buffer2[0] | (buffer2[1]<<8);
  
 fread(buffer2,sizeof(buffer2),1,f);
 //convert little endial to big endian 2 bytes int;
 wav->channels = buffer2[0] | (buffer2[1]<<8);
 
 fread(buffer4,sizeof(buffer4),1,f);
 //convert little endial to big endian 4 bytes int;
 wav->sample_rate = buffer[0]|(buffer4[1]<<8)|(buffer4[2]<<16)|(buffer4[3]<<24);

 fread(buffer4,sizeof(buffer4),1,f);
 //convert little endial to big endian 4 bytes int;
 wav->byte_rate = buffer[0]|(buffer4[1]<<8)|(buffer4[2]<<16)|(buffer4[3]<<24);

 fread(buffer2,sizeof(buffer2),1,f);
 //convert little endial to big endian 2 bytes int;
 wav->block_align = buffer2[0] | (buffer2[1]<<8);
 
 fread(buffer2,sizeof(buffer2),1,f);
 //convert little endial to big endian 2 bytes int;
 wav->bit_per_sample = buffer2[0] | (buffer2[1]<<8);

 fread(wav->data_id,sizeof(wav->data_id),1,f);

 fread(buffer4,sizeof(buffer4),1,f);
 //convert little endial to big endian 4 bytes int;
 wav->data_size = buffer[0]|(buffer4[1]<<8)|(buffer4[2]<<16)|(buffer4[3]<<24);

 //버퍼 만들고 읽기
 num_of_sample = (8 * wav->data_size) / (wav->channels * wav->bit_per_sample);
  
 
 wav->buffer.buf = (short*)malloc(sizeof(short)*wav->num_of_sample);

 free(f);
 return wav;
}

void write_WAV(WAV*){}

MAT* WAV2MAT(WAV*){}

void free_WAV(WAV*wav)
{
  free(wav->buffer.buf);
	free(wav);
}

void print_WAV(WAV*wav)
{
  printf("----------WAV HEADER INFOMATIONi----------\n");
  printf("riff_id :  %s\n",wav->riff_id);
  printf("riff_size : %u \n",wav->riff_size );
  printf("wave_id : %s\n",wav->wave_id  );
  printf("fmt_id : %s\n",wav->fmt_id  );
  printf("fmt_size : %u\n",wav->fmt_size  );
  printf("fmt_type : %u\n",wav->fmt_type  );
  printf("channels : %u\n",wav->channels  );
  printf("sample_rate %u: \n",wav->sample_rate  );
  printf("byte_rate : %u\n",wav->byte_rate  );
  printf("block_align : %u\n",wav->block_align  );
  printf("bit_per_sample : $u\n",wav->bit_per_sample  );
  printf("data_id : %s\n",wav->data_id  );
  printf("data_size : $u\n",wav->data_size  );



}
