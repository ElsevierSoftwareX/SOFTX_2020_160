/*
# Copyright (C) 2019 Marco Drago, Igor Yakushin, Sergey Klimenko
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


void mergeTREE()
{
  char list[512], s[512], out[512], nametree[512];
  cin>>list>>out>>nametree;

  TChain tree(nametree);

  FILE *in;
  
  if( (in=fopen(list,"r"))!=NULL )
    { 
      while(fgets(s,512,in) != NULL)
	{
	  for(int i=strlen(s)-1; i<512; i++) s[i]='\0';
	  tree.Add(s);
	}
      fclose(in);
    }
  
  sprintf(s,"merge/%s.root",out);
  tree.Merge(s);
}
