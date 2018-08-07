/*
 * ===========================================================
 *           Copyright (c) 2018, __IIPLAB__
 *                All rights reserved.
 * 
 * This Source Code Form is subject to the terms of
 * the Mozilla Public License, v. 2.0. 
 * If a copy of the MPL was not distributed with this file,
 *  You can obtain one at http://mozilla.org/MPL/2.0/.
 * ===========================================================
 */
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

/*read wav from file_path and return allocated WAV struct*/
WAV* read_wav(char* file_path);

/*write wav file from WAV struct into file_path */
void write_wav(WAV* wav, char* file_path);

/*return MAT struct with wav's buffer*/
MAT* wav2mat(WAV* wav);

/*return MAT struct with wav_buf*/
MAT* wav_buf2MAT(WAV_BUF* buf);

/*return WAV struct mat as wav_buf with given sample_rate*/
WAV* mat2wav(MAT* mat, UINT sample_rate);

/*free WAV strcut*/
void free_wav(WAV* wav);

/* Display header info of WAV struct */
void print_wav(WAV* wav);
