#include <stdio.h>

int main(int argc, char *argv[])
{
  if ( argc != 3 ) {
    printf("# usage: prog input.txt output.txt\n");
    return __LINE__;
  }
  FILE* fin = fopen(argv[1], "r");
  if ( fin == NULL ) {
    printf("[error] can not open %s\n", argv[1]);
    return __LINE__;
  }
  FILE* out = fopen(argv[2], "w");
  double time, error;
  char unit[16];
  while ( fscanf(fin, "%*[^:]:%lf%s%*[^\n]", &time, unit) == 2 ) {
    fscanf(fin, "%*[^:]:%lf%*[^\n]", &error);
    if ( strcmp(unit, "msec") == 0 )
      time /= 1000;
    fprintf(out, "%lf %e\n", time, error);
  }
  printf("[info] done\n");
  return 0;
}
