#include "mother.h"

int main(int argc, char *argv[]) {
  WAV *wav, *m_wav;
  MAT *mat;
  if (argc == 1) {
    printf("ERROR : You need to give wav file as an argument\n");
    return -1;
  }
  wav = read_wav(argv[1]);
  print_wav(wav);
  mat = wav2mat(wav);

  scal(10., mat);

  m_wav = mat2wav(mat, wav->sample_rate);
  print_wav(m_wav);
  write_wav(m_wav, "m2w.wav");

  free_wav(m_wav);
  free_mat(mat);
  free_wav(wav);
  return 0;
}
