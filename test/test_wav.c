#include "mother.h"

int main(int argc, char*argv[])
{
  WAV * wav;
  MAT * mat;
	if(argc == 1)
		{
		printf("ERROR : You need to give wav file as an argument\n");
		return -1;
		}
	wav = read_WAV(argv[1]);
  print_WAV(wav);
	write_WAV(wav,"temp.wav");
	free_WAV(wav);
  read_WAV("temp.wav");
	print_WAV(wav);

  mat = WAV2MAT(wav);
	free_MAT(mat);

	free_WAV(wav);
return 0;
}
