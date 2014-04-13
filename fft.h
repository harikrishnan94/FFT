#ifndef _FFT_H
#define _FFT_H

#define VERSION "FFT-1.0.0"

#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpc.h>

typedef mpc_ptr Sequence;
typedef mpc_t Element;

void fft_init(size_t, mpfr_prec_t);
void print(Sequence, size_t);
void bitr(Sequence, size_t);
void fft(Sequence, size_t);
void fft_free(size_t);
void seq_free(Sequence, size_t);

#endif
