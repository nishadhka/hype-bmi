#include <stdio.h>
extern "C" int get_num_basins();
extern "C" int initialize(char*, int*);

int main()
{
   int iens = 0;
   char path[200];
   path[0] = '.';
   path[1] = '\0';
   initialize(path, &iens);
   int n = get_num_basins();
   printf("The number of subbasins is %1d", n);
}