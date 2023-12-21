/*
 * 
 * This toy utility converts a binary data file into an equivalentascii file
 *
 * -- Each row correspond to a particle; the meaning of the columsn is as follows:
 * column 1     : the unique ID of the particle
 * column 2     : the type of the particle (see below)
 * column 3     : the ID of the fof to which the particle belongs
 * column 4     : the ID of the subhalo to which the particle belongs
 * column 5,6,7 : the spatial position of the particle
 * column 8     : the kinetic energy of the particle
 * column 9     : the potential energy of the particle 
 *
 * -- particle types:
 * 0     : gas particles
 * 1,2,3 : dark matter (there will be only type 1
 * 4     : stars
 * 5     : black holes
 *
 * -- meaning of the subhalo IDs:
 *
 * positive values, i.e. [1, ...] ==> a gravitationally bounded subhalo, i.e. "a galaxy"
 * 0 (zero)                       ==> particles gravitationally bound to the FOF that are not
 *                                    bound to a subhalo ("diffuse" particles)
 * -1                             ==> "fuzzy" particles: they have been included by the fof
 *                                    algorithm because of gemoetrical reasons but are not
 *                                    bound to the FOF (nor to any subhalos)
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../include/read_fof_snapshot.h"

file_info read_data_from_file(char* filename){
  data_t *data;
  num_t   N;
  int     id_size;
  int     float_size;
  size_t  ret;
  //aaa 
  FILE *file = fopen( filename, "r" );

  if( file == NULL ) {
    printf("unable to open file %s\n", filename ); exit(2); }
  // get the number of particles in the file
  ret = fread ( &N, sizeof(num_t), 1, file );

  // get the size of IDs
  ret = fread ( &id_size, sizeof(int), 1, file );

  // get the size of floats
  ret = fread ( &float_size, sizeof(int), 1, file );

  // N.B. obviously how the data_t and list_t structures are made up should depend on the
  // size of IDs and floats, that we just got from the file
  // For the sake of simplicity here we have defined it through #define at compile time.  
  
  // get the data
  data = (data_t*)malloc( sizeof(data_t) * N );

  ret = fread( data, sizeof(data_t), N, file );

  file_info f = {N, id_size, float_size, data};

  return f; 
}
