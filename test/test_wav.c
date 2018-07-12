#include "mother.h"

int main(int argc, char *argv[])
{
	WAV *wav, *m_wav;
	MAT *mat;
	if (argc == 1)
	{
		printf("ERROR : You need to give wav file as an argument\n");
		return -1;
	}
	wav = read_WAV(argv[1]);
	print_WAV(wav);
	write_WAV(wav, "w2w.wav");
	free_WAV(wav);
	wav = read_WAV("w2w.wav");
	print_WAV(wav);

	mat = WAV2MAT(wav);

	m_wav = MAT2WAV(mat, 48000);
	print_WAV(m_wav);
	write_WAV(m_wav, "m2w.wav");

	free_WAV(m_wav);
	free_MAT(mat);
	free_WAV(wav);
	return 0;
}
