/*
 * I am writing this thing in C
 * because I need to do tens of billions of things.
 *
 * This C code is very hard-coded style.
 * If something more flexible is needed,
 * use Python instead.
 */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#define NPHEN 40
#define NSNP 30000
#define LINESIZE 1000

double get_mean(double *phen, int *snp, int marker)
{
  double total = 0.0;
  int count = 0;
  int i;
  for (i=0; i<NPHEN; i++)
  {
    if (snp[i] == marker)
    {
      total += phen[i];
      count += 1;
    }
  }
  return total / count;
}

int snp_is_valid(int *snp)
{
  /*
   * A snp is said to be valid if it has a zero and a one.
   */
  int i;
  int has[2] = {0, 0};
  for (i=0; i<NPHEN; i++)
  {
    has[snp[i]] = 1;
  }
  return (has[0] && has[1]);
}

int get_additive_parameters(double *params, double *phen, int *snpa, int *snpb, int verbose)
{
  /*
   * The additive parameters are mu, x1, and x2.
   * These will be passed to the caller through the params argument.
   */
  int nsnp1 = 0;
  int nsnp2 = 0;
  int nboth = 0;
  double a11, a12, a13;
  double a21, a22, a23;
  double a31, a32, a33;
  double b11, b12, b13;
  double b21, b22, b23;
  double b31, b32, b33;
  double determinant;
  double phentotal, phen_snpa, phen_snpb;
  int i;
  phentotal = 0.0;
  phen_snpa = 0.0;
  phen_snpb = 0.0;
  for (i=0; i<NPHEN; i++)
  {
    nsnp1 += snpa[i];
    nsnp2 += snpb[i];
    nboth += snpa[i] * snpb[i];
    phentotal += phen[i];
    phen_snpa += phen[i]*snpa[i];
    phen_snpb += phen[i]*snpb[i];
  }
  /* Calculate X'X the fast way */
  a11 = NPHEN; a12 = nsnp1; a13 = nsnp2;
  a21 = nsnp1; a22 = nsnp1; a23 = nboth;
  a31 = nsnp2; a32 = nboth; a33 = nsnp2;
  /* Calculate two parts of the inverse of X'X */
  b11 = a33*a22 - a32*a23; b12 = a32*a13 - a33*a12; b13 = a23*a12 - a22*a13;
  b21 = a31*a23 - a33*a21; b22 = a33*a11 - a31*a13; b23 = a21*a13 - a23*a11;
  b31 = a32*a21 - a31*a22; b32 = a31*a12 - a32*a11; b33 = a22*a11 - a21*a12;
  determinant = a11*(a33*a22 - a32*a23) - a21*(a33*a12 - a32*a13) + a31*(a23*a12 - a22*a13);
  /* The determinant should not be zero */
  if (determinant == 0.0) return -1;
  if (verbose)
  {
    printf("%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", a11, a12, a13, a21, a22, a23, a31, a32, a33);
    printf("%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", b11, b12, b13, b21, b22, b23, b31, b32, b33);
  }
  /* Calculate the parameters */
  params[0] = (b11*phentotal + b12*phen_snpa + b13*phen_snpb) / determinant;
  params[1] = (b21*phentotal + b22*phen_snpa + b23*phen_snpb) / determinant;
  params[2] = (b31*phentotal + b32*phen_snpa + b33*phen_snpb) / determinant;
  return 0;
}

double get_additive_squared_error(double *params, double *phen, int *snpa, int *snpb)
{
  /*
   * The additive parameters are mu, x1, and x2.
   */
  double squared_error = 0;
  double estimate;
  double d;
  int i;
  for (i=0; i<NPHEN; i++)
  {
    estimate = params[0] + params[1] * snpa[i] + params[2] * snpb[i];
    d = phen[i] - estimate;
    squared_error += d*d;
  }
  return squared_error;
}

double get_additive_snp_error(double *phen, int *snpa, int *snpb)
{
  /*
   * This counts on the snps not being degenerate.
   * This would include snps that are all ones or zeros.
   */
  double params[3];
  get_additive_parameters(params, phen, snpa, snpb, 0);
  return get_additive_squared_error(params, phen, snpa, snpb);
}

double get_interacting_snp_error(double *phen, int *snpa, int *snpb)
{
  int count[4] = {0, 0, 0, 0};
  double sum[4] = {0.0, 0.0, 0.0, 0.0};
  double sumsq[4] = {0.0, 0.0, 0.0, 0.0};
  int i;
  double sum_squared_error = 0.0;
  double ph;
  int index;
  for (i=0; i<NPHEN; i++)
  {
    ph = phen[i];
    index = (snpa[i]<<1) + snpb[i];
    count[index]++;
    sum[index] += ph;
    sumsq[index] += ph*ph;
  }
  for (i=0; i<4; i++)
  {
    if (count[i])
    {
      sum_squared_error += sumsq[i] - (sum[i]*sum[i])/count[i];
    }
  }
  return sum_squared_error;
}

int do_additive_snp_analysis(double *phen, int *snp)
{
  int i, j, k;
  int any_squared_error = 0;
  double best_squared_error = -1;
  double se;
  int bins[4];
  int ncombos;
  int index;
  int *snpa;
  int *snpb;
  for (i=0; i<NSNP; i++) {
    snpa = snp + i*NPHEN;
    for (j=i+1; j<NSNP; j++) {
      snpb = snp + j*NPHEN;
      /* Check the compatibility of the snps. */
      for (k=0; k<4; k++) bins[k] = 0;
      for (k=0; k<NPHEN; k++)
      {
        index = (snpa[k]<<1) + snpb[k];
        bins[index] += 1;
      }
      ncombos = 0;
      for (k=0; k<4; k++)
      {
        if (bins[k])
        {
          ncombos++;
        }
      }
      /* 
       * If all four combos are present then use additive analysis.
       * If only three combos are present then additive analysis would
       * give the same result as interacting analysis.
       * For the three combo case, the additive analysis would break on the
       * singular matrix, while the interacting analysis would work fine.
       * If fewer than three combos are present then do not analyze the pair.
       */
      if (ncombos < 3) continue;
      if (ncombos < 4)
      {
        se = get_interacting_snp_error(phen, snpa, snpb);
      }
      else
      {
        se = get_additive_snp_error(phen, snpa, snpb);
      }
      /* Report a new best score if one was found. */
      if ((!any_squared_error) || (se < best_squared_error))
      {
        printf("SNP pair {%d, %d} has squared error %lf\n", i, j, se);
        any_squared_error = 1;
        best_squared_error = se;
      }
    }
  }
  return 0;
}

int do_interacting_snp_analysis(double *phen, int *snp)
{
  int i, j;
  int any_squared_error = 0;
  double best_squared_error = -1;
  double se;
  for (i=0; i<NSNP; i++) {
    for (j=i+1; j<NSNP; j++) {
      se = get_interacting_snp_error(phen, snp + i*NPHEN, snp + j*NPHEN);
      if ((!any_squared_error) || (se < best_squared_error))
      {
        printf("SNP pair {%d, %d} has squared error %lf\n", i, j, se);
        any_squared_error = 1;
        best_squared_error = se;
      }
    }
  }
}

int do_single_snp_analysis(double *phen, int *snp)
{
  int i, j;
  int any_found = 0;
  double best_squared_error = -1;
  int best_snp_index;
  double squared_error;
  double mu_0;
  double mu_1;
  double mu;
  double d;
  int *snpa;
  int marker;
  for (i=0; i<NSNP; i++)
  {
    snpa = snp + i*NPHEN;
    if (!snp_is_valid(snpa)) continue;
    mu_0 = get_mean(phen, snpa, 0);
    mu_1 = get_mean(phen, snpa, 1);
    squared_error = 0.0;
    for (j=0; j<NPHEN; j++)
    {
      mu = (snpa[j] == 0) ? mu_0 : mu_1;
      d = phen[j] - mu;
      squared_error += d*d;
    }
    if ((!any_found) || (squared_error < best_squared_error))
    {
      any_found = 1;
      best_squared_error = squared_error;
      printf("best snp index: %d\n", i);
      printf("best squared_error: %lf\n", squared_error);
      printf("\n");
    }
  }
}

double get_balancedness(int *snp)
{
  int n[2] = {0, 0};
  int i;
  double mu;
  for (i=0; i<NPHEN; i++)
  {
    n[snp[i]]++;
  }
  mu = n[1] / (float) NPHEN;
  return n[0]*mu*mu + n[1]*(1-mu)*(1-mu);
}

int do_single_snp_t_analysis(double *phen, int *snp)
{
  int i, j;
  int any_found = 0;
  double squared_error;
  double best_t;
  double mu_0;
  double mu_1;
  double mu;
  double d;
  double t;
  double balancedness;
  double effect_size;
  double sumphena;
  double sumphenb;
  double sumphen;
  double sumphen2;
  double tscorehat_num;
  double tscorehat_denom;
  int n0;
  int n1;
  double monofoo;
  int *snpa;
  int marker;
  for (i=0; i<NPHEN; i++)
  {
    sumphen += phen[i];
    sumphen2 += phen[i]*phen[i];
  }
  for (i=0; i<NSNP; i++)
  {
    snpa = snp + i*NPHEN;
    if (!snp_is_valid(snpa)) continue;
    /* get the balancedness */
    balancedness = get_balancedness(snpa);
    /* get the effect size */
    mu_0 = get_mean(phen, snpa, 0);
    mu_1 = get_mean(phen, snpa, 1);
    effect_size = mu_1 - mu_0;
    /* get the sum of squares of errors */
    squared_error = 0.0;
    for (j=0; j<NPHEN; j++)
    {
      mu = (snpa[j] == 0) ? mu_0 : mu_1;
      d = phen[j] - mu;
      squared_error += d*d;
    }
    /*
     * Combine these to get the statistic of interest.
     * This isn't really the t statistic.
     */
    t = squared_error;
    t /= effect_size*effect_size*balancedness;
    /* compare the statistic to that of other models */
    if ((!any_found) || (t < best_t))
    {
      any_found = 1;
      best_t = t;
      printf("snp index: %d\n", i);
      printf("t-like statistic: %lf\n", t);
      printf("squared error: %lf\n", squared_error);
      /* get some temp stats */
      sumphena = 0.0;
      sumphenb = 0.0;
      n0 = 0;
      n1 = 0;
      for (j=0; j<NPHEN; j++)
      {
        if (snpa[j])
        {
          n1++;
          sumphenb += phen[j];
        }
        else
        {
          n0++;
          sumphena += phen[j];
        }
      }
      monofoo = (sumphena*sumphena)/n0 + (sumphenb*sumphenb)/n1;
      tscorehat_num = sumphen2 - monofoo;
      tscorehat_denom = sumphen*sumphen + monofoo;
      printf("numerator: %lf\n", tscorehat_num);
      printf("denominator: %lf\n", tscorehat_denom);
      printf("ratio: %lf\n", tscorehat_num / tscorehat_denom);
      printf("\n");
    }
  }
}

int get_data_from_file(double *phen, int *snp, const char *filename_in)
{
  /*
   * @param phen: fill this array of phenotypes
   * @param snp: fill this array of snps
   * @param filename_in: read this filename
   * @return: nonzero if error
   */
  int i, j;
  char line[LINESIZE+1];
  int marker;
  FILE *fin;
  char c;
  int nread;
  fin = fopen(filename_in, "rt");
  for (i=0; i<NPHEN; i++)
  {
    nread = fscanf(fin, "%lf", phen + i);
  }
  for (i=0; i<NSNP; i++)
  {
    nread = fscanf(fin, "%s", line);
    for (j=0; j<NPHEN; j++)
    {
      c = line[j];
      if (c == '0')
      {
        marker = 0;
      }
      else if (c == '1')
      {
        marker = 1;
      }
      else
      {
        return -1;
      }
      snp[i*NPHEN + j] = marker;
    }
  }
  fclose(fin);
  return 0;
}

int test_additive_parameters(void)
{
  /*
   * Call this directly from main().
   */
  srand(time(NULL));
  int snpa[NPHEN];
  int snpb[NPHEN];
  double phen[NPHEN];
  int i;
  double mu = 1.0;
  double alpha = 2.0;
  double beta = 3.0;
  double error;
  double params[3];
  int result;
  for (i=0; i<NPHEN; i++)
  {
    snpa[i] = rand() % 2;
    snpb[i] = rand() % 2;
    error = (double) rand();
    error /= RAND_MAX;
    error -= 0.5;
    error /= 5;
    phen[i] = mu + alpha * snpa[i] + beta * snpb[i] + error;
  }
  result = get_additive_parameters(params, phen, snpa, snpb, 1);
  if (result) return result;
  printf("estimated mu: %lf\n", params[0]);
  printf("estimated alpha: %lf\n", params[1]);
  printf("estimated beta: %lf\n", params[2]);
  return 0;
}

int do_analysis(void)
{
  /*
   * Call this directly from main().
   */
  double phen[NPHEN];
  int snp[NPHEN*NSNP];
  const char *filename_in = "/home/argriffi/phenotype-challenge/testphen-simple.txt";
  int result;
  printf("reading the data...\n");
  result = get_data_from_file(phen, snp, filename_in);
  if (result) return result;
  printf("doing the analysis...\n");
  /* return do_additive_snp_analysis(phen, snp); */
  return do_single_snp_t_analysis(phen, snp);
}

int main(int argc, char *argv[])
{
  return do_analysis();
}
