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
  write_wav(wav, "w2w.wav");
  free_wav(wav);
  wav = read_wav("w2w.wav");
  print_wav(wav);

  mat = wav2mat(wav);

  m_wav = mat2wav(mat, 48000);
  print_wav(m_wav);
  write_wav(m_wav, "m2w.wav");

  free_wav(m_wav);
  free_mat(mat);
  free_wav(wav);
  return 0;
}
