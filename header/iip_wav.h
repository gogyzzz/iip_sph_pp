#include "iip_type.h"

typedef struct WAV_BUF {
  UINT buf_size;
  short channels;
  // buf[buf_size * channels]
  short* buf;

} WAV_BUF;

typedef struct WAV {
  // 4bytes fixed size header infomation -> uint32_t

  char riff_id[4];     // riff string
  uint32_t riff_size;  // overall size of file in bytes
  char wave_id[4];     // wave string
  char fmt_id[4];      // fmt string with trailing null char
  uint32_t fmt_size;   // length of the format data;
  short
      fmt_type;  // format type 1-PCM 3-IEEE float 6- 8bit A law, 7- 8bit ,u law
  short channels;        // no of channel
  uint32_t sample_rate;  // SampleRate(blocks per second)
  uint32_t byte_rate;  // ByteRate = SampleRate * NumChannels * BitsPerSample/8
  short block_align;   // NumChannels * BitsPerSample/8
  short bit_per_sample;  // bits per sample, 8 - 8bits, 16-16bits etc
  char data_id[4];       // DATA string or FLLR string
  uint32_t data_size;    // NumSamples * NumChannels * BitsPerSample/8 - size of
                       // the nex chunk that will be read

  WAV_BUF buffer;

} WAV;

WAV* read_WAV(char* file_path);

void write_WAV(WAV* wav, char* file_path);

MAT* WAV2MAT(WAV* wav);

MAT* WAV_BUF2MAT(WAV_BUF* buf);
WAV* MAT2WAV(MAT* mat, UINT sample_rate);

void free_WAV(WAV* wav);

// Display header of WAV
void print_WAV(WAV* wav);
