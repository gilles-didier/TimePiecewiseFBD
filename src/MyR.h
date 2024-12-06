#ifndef MyRF
#define MyRF

#ifdef __cplusplus
extern "C" {
#endif

void error(const char * format, ...);
void warning(const char * format, ...);
void Rprintf(const char * format, ...);
double lgammafn(double x);
double unif_rand();
int bernoulli(double prob);
int unif_int_rand(int max);

#ifdef __cplusplus
}
#endif

#endif
