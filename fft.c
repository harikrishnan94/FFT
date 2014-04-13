#include "fft.h"
#include <stdlib.h>
#include <math.h>

mpfr_t ZERO, tmp;
Sequence twid_fact;
size_t LEN = 2;
mpfr_prec_t precision = 53;
mpc_rnd_t RND = MPC_RNDDU;
Element min2pii; //(0, -1 * 2 * PI / N);
Element temp;
Sequence new_seq;

#define swap(a,b) mpc_swap(a,b)

void seq_free(Sequence seq, size_t N)
{
	while (N--)
		mpc_clear(seq + N);
	free(seq);
}

void fft_init(size_t N, mpfr_prec_t prec)
{
	if (prec) precision = prec;
	if (N) LEN = N;
	twid_fact = (Sequence) calloc(LEN, sizeof(mpc_t));
	mpfr_init_set_d(ZERO, 0.0, MPFR_RNDA);
	mpfr_init2(tmp, precision);
	mpfr_const_pi(tmp, MPFR_RNDA);
	mpfr_mul_si(tmp, tmp, -2, MPFR_RNDA);
	mpfr_div_ui(tmp, tmp, N, MPFR_RNDA);
	mpc_init2(min2pii, precision);
	mpc_set_fr_fr(min2pii, ZERO, tmp, RND);
	mpc_init2(temp, precision);
	new_seq = (Sequence) calloc(N, sizeof(mpc_t));
	size_t n = N / 2;
	while (n--)
	{
		mpc_init2(twid_fact + n, precision);
		mpc_mul_ui(twid_fact + n, min2pii, n, RND);
		mpc_exp(twid_fact + n, twid_fact + n, RND);
	}
}

void add_diff(Sequence seq, size_t N)
{
	Sequence hseq = seq + N / 2, end = seq + N;
	while (hseq < end)
	{
		mpc_set(temp, seq, RND); //temp = *seq
		mpc_add(seq, seq, hseq, RND); //*seq += *hseq
		mpc_sub(hseq, temp, hseq, RND); //*hseq = temp - *hseq
		seq++, hseq++;
	}
}

void mul_twid(Sequence seq, size_t N)
{
	unsigned int mul_fact = LEN / (N + N);
	Sequence end = seq + N;
	for (unsigned i = 0; seq < end; i += mul_fact)
	{
		//(*seq++) *= twid_fact[i];
		mpc_mul(seq, seq, twid_fact + i, RND);
		seq++;
	}
}

void print(Sequence seq, size_t N)
{
	size_t n = 0;
	while (n++ < N)
	{
		mpc_out_str(stdout, 10, 16, seq, RND); //cout << *seq++ << endl;
		seq++;
		putchar('\n');
	}
}

unsigned reverse(unsigned int v)
{
// swap odd and even bits
	v = ((v >> 1) & 0x55555555) | ((v & 0x55555555) << 1);
// swap consecutive pairs
	v = ((v >> 2) & 0x33333333) | ((v & 0x33333333) << 2);
// swap nibbles ...
	v = ((v >> 4) & 0x0F0F0F0F) | ((v & 0x0F0F0F0F) << 4);
// swap bytes
	v = ((v >> 8) & 0x00FF00FF) | ((v & 0x00FF00FF) << 8);
// swap 2-byte long pairs
	v = (v >> 16) | (v << 16);
	return (v);
}

void bitr(Sequence seq, size_t N)
{
	unsigned n = ceil(log((double) (N)) / log(2.0));
	for (unsigned i = N - 1; i + 1; i--)
	{
		mpc_init2(new_seq + i, precision);
		mpc_set(new_seq + i, seq + (reverse(i) >> (32 - n)), RND);
	}
	for (unsigned i = N - 1; i + 1; i--)
		mpc_set(seq + i, new_seq + i, RND);
}

void fft(Sequence seq, size_t N)
{
	if (N == 2)
	{
		mpc_set(temp, seq, RND); //temp = seq[0]
		mpc_add(seq, seq, seq + 1, RND); //seq[0] += seq[1]
		mpc_sub(seq + 1, temp, seq + 1, RND); //seq[1] = temp - seq[1]
		return;
	}
	add_diff(seq, N);
	fft(seq, N / 2);
	mul_twid(seq + N / 2, N / 2);
	fft(seq + N / 2, N / 2);
}

void fft_free(size_t N)
{
	mpfr_clear(ZERO);
	mpfr_clear(tmp);
	mpc_clear(temp);
	mpc_clear(min2pii);
	seq_free(twid_fact, N / 2);
	seq_free(new_seq, N);
}
