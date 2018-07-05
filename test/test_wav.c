#include "mother.h"

int main()
{
  WAV * wav;
  wav = read_WAV("wav.wav");
  
  print_WAV(wav);
  free_WAV(wav);
return 0;
}
