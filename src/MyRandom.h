#ifndef MyRandomF
#define MyRandomF

#ifdef __cplusplus
extern "C" {
#endif
void my_rand_init();
void *my_rand_get_data();
void my_rand_free_data(void *data);
double my_rand_unif_unit(void *data);
int my_rand_bern(void *data, double prob);
int my_rand_unif_range(void *data, int max);

#ifdef __cplusplus
}
#endif

#endif
