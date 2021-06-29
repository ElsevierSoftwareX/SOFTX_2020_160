/*
# Copyright (C) 2019 Igor Yakushin, Sergey Klimenko
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[])
{
  int n;
  char fn[1024];

  if(argc==3)
    {
      strcpy(fn,argv[1]);
      n=atoi(argv[2]);
      //      printf("n=%d gps_start=%d gps_end=%d\n",fn,n);
    }
  else
    {
      printf("Usage: extract frame_file_name number_of_samples\n");
      return 1;
    }

  float *a=(float*)malloc(n*sizeof(float));
  FILE *f=fopen(fn,"rb");
  fread(a,sizeof(float),n,f);
  fclose(f);

  for(int i=0;i<n;i++)
    {
      printf("%d\n",int(a[i]));
    }
  free(a);
  return 0;
}
