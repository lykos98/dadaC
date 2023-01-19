typedef struct 
{
   unsigned long int idx;
   float d; 
} datapoint;


#if defined(LONG_IDS)
typedef unsigned long long PID_t;
#else
typedef unsigned int       PID_t;
#endif

#if defined(DOUBLEPRECISION)
typedef double float_out;
#else
typedef float float_out;
#endif

typedef unsigned long long num_t;
typedef int                fgid_t;

typedef struct { PID_t pid; int fofid; fgid_t gid; } list_t;
typedef struct { list_t idx; int type; float_out pos[3], ekin, pot; } data_t;
typedef struct { int N; int id_size; int float_size; data_t * data} file_info;


