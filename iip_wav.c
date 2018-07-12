#include "iip_wav.h"

WAV *read_WAV(char *file_path)
{
	WAV *wav;
	FILE *f;
	uint32_t read = 0;
	unsigned char buffer4[4];
	unsigned char buffer2[2];
	UINT num_of_sample;
	UINT size_of_each_sample;
	UINT byte_in_each_channel;

#if DEBUG
	printf("%s\n", __func__);
#endif
	//파일 유무확인
	if (!(f = fopen(file_path, "rb")))
	{
		printf("Failed to open '%s'\nAborting process\n", file_path);
		exit(-1);
		return NULL;
	}

	//WAV 읽기
	wav = (WAV *)malloc(sizeof(WAV));

	fread(wav->riff_id, sizeof(wav->riff_id), 1, f);

	fread(buffer4, sizeof(buffer4), 1, f);
	//convert little endial to big endian 4 bytes int;
	wav->riff_size = buffer4[0] | (buffer4[1] << 8) | (buffer4[2] << 16) | (buffer4[3] << 24);
	// fread(wav->riff_size,sizeof(wav->riff_size),1,f);

	fread(wav->wave_id, sizeof(wav->wave_id), 1, f);

	fread(wav->fmt_id, sizeof(wav->fmt_id), 1, f);

	fread(buffer4, sizeof(buffer4), 1, f);
	//convert little endial to big endian 4 bytes int;
	wav->fmt_size = buffer4[0] | (buffer4[1] << 8) | (buffer4[2] << 16) | (buffer4[3] << 24);

	fread(buffer2, sizeof(buffer2), 1, f);
	//convert little endial to big endian 2 bytes int;
	wav->fmt_type = buffer2[0] | (buffer2[1] << 8);

	fread(buffer2, sizeof(buffer2), 1, f);
	//convert little endial to big endian 2 bytes int;
	wav->channels = buffer2[0] | (buffer2[1] << 8);

	fread(buffer4, sizeof(buffer4), 1, f);
	//convert little endial to big endian 4 bytes int;
	wav->sample_rate = buffer4[0] | (buffer4[1] << 8) | (buffer4[2] << 16) | (buffer4[3] << 24);

	fread(buffer4, sizeof(buffer4), 1, f);
	//convert little endial to big endian 4 bytes int;
	wav->byte_rate = buffer4[0] | (buffer4[1] << 8) | (buffer4[2] << 16) | (buffer4[3] << 24);

	fread(buffer2, sizeof(buffer2), 1, f);
	//convert little endial to big endian 2 bytes int;
	wav->block_align = buffer2[0] | (buffer2[1] << 8);

	fread(buffer2, sizeof(buffer2), 1, f);
	//convert little endial to big endian 2 bytes int;
	wav->bit_per_sample = buffer2[0] | (buffer2[1] << 8);

	fread(wav->data_id, sizeof(wav->data_id), 1, f);

	fread(buffer4, sizeof(buffer4), 1, f);
	//convert little endial to big endian 4 bytes int;
	wav->data_size = buffer4[0] | (buffer4[1] << 8) | (buffer4[2] << 16) | (buffer4[3] << 24);

	//버퍼 만들고 읽기
	num_of_sample = (8 * wav->data_size) / (wav->channels * wav->bit_per_sample);
	size_of_each_sample = (wav->channels * wav->bit_per_sample) / 8;

	byte_in_each_channel = (size_of_each_sample / wav->channels);

	printf("num_of_sample       : %u \n", num_of_sample);
	printf("size_of_each_sample : %u \n", size_of_each_sample);
	printf("duartion in seconds : %u \n", wav->riff_size / wav->byte_rate);
	printf("byte_in_each_channel : %u \n", byte_in_each_channel);
	wav->buffer.buf = (short *)malloc(wav->data_size);

	//printf("Buffer read : %lu\n",fread(wav->buffer.buf,sizeof(short),wav->data_size/sizeof(short),f));
	printf("Buffer read : %lu\n", fread(wav->buffer.buf, sizeof(short), wav->channels * num_of_sample, f));
	wav->buffer.channels = wav->channels;
	wav->buffer.buf_size = num_of_sample;
	fclose(f);
	return wav;
}

void write_WAV(WAV *wav, char *file_path)
{
	FILE *f;

#if DEBUG
	printf("%s\n", __func__);
#endif

	f = fopen(file_path, "wb");
	if (f == NULL)
	{
		printf("Failed to Open : %s\n", file_path);
		return;
	}

	fwrite(wav->riff_id, sizeof(char), 4, f);
	fwrite(&(wav->riff_size), sizeof(uint32_t), 1, f);
	fwrite((wav->wave_id), sizeof(char), 4, f);
	fwrite((wav->fmt_id), sizeof(char), 4, f);
	fwrite(&(wav->fmt_size), sizeof(uint32_t), 1, f);
	fwrite(&(wav->fmt_type), sizeof(short), 1, f);
	fwrite(&(wav->channels), sizeof(short), 1, f);
	fwrite(&(wav->sample_rate), sizeof(uint32_t), 1, f);
	fwrite(&(wav->byte_rate), sizeof(uint32_t), 1, f);
	fwrite(&(wav->block_align), sizeof(short), 1, f);
	fwrite(&(wav->bit_per_sample), sizeof(short), 1, f);
	fwrite(wav->data_id, sizeof(char), 4, f);
	fwrite(&(wav->data_size), sizeof(uint32_t), 1, f);

	printf("WRITE BUFFER : %lu\n", fwrite((wav->buffer.buf), sizeof(short), wav->data_size / (wav->bit_per_sample / 8), f));

	fclose(f);
}
/*
MAT* WAV2MAT(WAV*){}
*/
void free_WAV(WAV *wav)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	free(wav->buffer.buf);
	free(wav);
}

void print_WAV(WAV *wav)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	UINT t;
	printf("----------WAV HEADER INFOMATION----------\n");
	printf("riff_id        : %c%c%c%c\n", wav->riff_id[0], wav->riff_id[1], wav->riff_id[2], wav->riff_id[3]);
	printf("riff_size      : %u \n", wav->riff_size);
	printf("wave_id        : %c%c%c%c\n", wav->wave_id[0], wav->wave_id[1], wav->wave_id[2], wav->wave_id[3]);
	printf("fmt_id         : %c%c%c%c\n", wav->fmt_id[0], wav->fmt_id[1], wav->fmt_id[2], wav->fmt_id[3]);
	printf("fmt_size       : %u\n", wav->fmt_size);
	t = wav->fmt_type;
	switch (t)
	{
	case 1:
		printf("fmt_type       : %u - PCM\n", wav->fmt_type);
		break;
	case 3:
		printf("fmt_type       : %u - IEEE float\n", wav->fmt_type);
		break;
	case 6:
		printf("fmt_type       : %u - 8bit A law\n", wav->fmt_type);
		break;
	case 7:
		printf("fmt_type       : %u - 8 bit U law\n", wav->fmt_type);
		break;
	}
	printf("channels       : %u\n", wav->channels);
	printf("sample_rate    : %u \n", wav->sample_rate);
	printf("byte_rate      : %u\n", wav->byte_rate);
	printf("block_align    : %u\n", wav->block_align);
	printf("bit_per_sample : %u\n", wav->bit_per_sample);
	printf("data_id        : %c%c%c%c\n", wav->data_id[0], wav->data_id[1], wav->data_id[2], wav->data_id[3]);
	printf("data_size      : %u\n", wav->data_size);
}

MAT *WAV_BUF2MAT(WAV_BUF *buf)
{
	MAT *mat;
	ITER i, j;
#if DEBUG
	printf("%s\n", __func__);
#endif
	mat = (MAT *)malloc(sizeof(MAT));
	mat->ndim = 1;
	mat->d0 = buf->buf_size;
	mat->d1 = buf->channels;
	mat->d2 = 1;
	mat->data = (DTYPE *)malloc(sizeof(DTYPE) * mat->d0 * mat->d1);

#pragma omp parallel for shared(mat, buf) private(i, j)
	for (i = 0; i < mat->d1; i++)
	{
		for (j = 0; j < mat->d0; j++)
		{
			mat->data[i * mat->d0 + j] = (DTYPE)(buf->buf[i + mat->d1 * j]);
		}
	}
	return mat;
}

MAT *WAV2MAT(WAV *wav)
{
	return WAV_BUF2MAT(&(wav->buffer));
}

WAV *MAT2WAV(MAT *mat, UINT sample_rate)
{
	ITER i;
	UINT d0, d1;
	WAV *wav;

	d0 = mat->d0;
	d1 = mat->d1;

	wav = (WAV *)malloc(sizeof(WAV));
	wav->riff_id[0] = 'R';
	wav->riff_id[1] = 'I';
	wav->riff_id[2] = 'F';
	wav->riff_id[3] = 'F';

	wav->wave_id[0] = 'W';
	wav->wave_id[1] = 'A';
	wav->wave_id[2] = 'V';
	wav->wave_id[3] = 'E';

	wav->fmt_id[0] = 'f';
	wav->fmt_id[1] = 'm';
	wav->fmt_id[2] = 't';
	wav->fmt_id[3] = ' ';

	//short 16bit ->2 bytes
	wav->fmt_size = 16;

	//1- PCM
	wav->fmt_type = 1;

	wav->channels = d1;

	//have to be Optional
	wav->sample_rate = sample_rate;

	//smaple_rate * channels * fmt_size / 8
	wav->byte_rate = wav->sample_rate * wav->channels * wav->fmt_size / 8;

	//bit per sample, 8 or 16. not sure
	wav->bit_per_sample = 16;

	//bit_per_sample*channels /8
	wav->block_align = wav->bit_per_sample * wav->channels / 8;
	wav->data_id[0] = 'd';
	wav->data_id[1] = 'a';
	wav->data_id[2] = 't';
	wav->data_id[3] = 'a';

	//Number of Samples * Number of Channels * Bit_per_sample / 8
	wav->data_size = mat->d0 * mat->d1 * wav->bit_per_sample / 8;

	wav->riff_size = wav->data_size + 44;

	wav->buffer.buf_size = d0;
	wav->buffer.channels = wav->channels;

	/*
	 * mat  1 1 1 1 1 
	 *      2 2 2 2 2
	 *      3 3 3 3 3
	 *
	 *buf   1 2 3 1 2 3 1 2 3 1 2 3 1 2 3  
	 *
	 * */
	wav->buffer.buf = (short *)malloc(sizeof(short) * d0 * d1);
#pragma omp parallel for shared(wav, mat) private(i)
	for (i = 0; i < d0 * d1; i++)
	{
		wav->buffer.buf[i] = (short)(mat->data[i % d1 * d0 + i / d1]);
	}

	return wav;
}
