/* ******************************************************************
NGA Subduction Model for Bozorgnia, Campbell, Gegor, Kuehn

magnitude scaling is modeled as logistic hinge function
magnitude break point is different for slab/interface, and is different for different subduction zones
parameter delta (controls the smoothness of the transition) is fixed
slope for magnitude scaling below the break point is different for slab/interface

depth scaling is modeled only for slab events
depth scaling is modeled as a logisic hinge function
scaling beyond the break point is fixed to zero
parameter delta is fixed

geometrical spreading is magnitude and slab dependent

overall constant, anelastic attenuation and linear Vs30-scaling vary by region
prior for regional standard deviation is exponential distribution
non-centered parameterisation

includes arc-crossing model for Japan, Central America, South America
other regions only have one attenuation coefficient

nonlinear site amplification is based on Campbell and Bozorgina 2014
psa for sie amplification includes event term

robust regression, using Student-t distribution for the record likelihood
different value of nu for different data base regions


region assignment:
1 Alaska (no FX)
2 Cascadia (no FX)
3 Central America & Mexico (FX)
4 Japan (FX)
5 New Zealand (no FX)
6 South America (FX)
7 Taiwan (FX)

find delta mb and delta zb

depth scaling for both if/slab

read in parameters for PGA to calculate median PGA values for nonlinear site amplification

estimate nft terms -- important: nft has M-6 in he exponent


no mb adjustment

centered on different M/R values

regional parameters are correlated

Cascadia is split into two regions for the constant term for Intraslab

9 regional coefficients:
  1 par1 constant interface
  2 par1a constant intraslab
  3 par7 Vs30-scaling
  4 par6x1 crossing arc r1
  5 par6x2 crossing arc r2
  6 par6x3 crossing arc r3
  7 par61 not crossing arc r1
  8 par62 not crossing arc r2
  9 par63 not crossing arc r3

apply magnitude break adjustment for if only for Japan and South America

****************************************************************** */

functions {
  real logistic_hinge(real x, real x0, real a, real b0, real b1, real delta);

  real logistic_hinge(real x, real x0, real a, real b0, real b1, real delta) { 
  real xdiff = x - x0;
  return a + b0 * xdiff + (b1 - b0) * delta * log1p_exp(xdiff / delta);
}
}

data {
  int<lower=1> N;  // number of records
  int<lower=1> NEQ;  // number of events
  int<lower=1> NREG;  // number of regions

  vector[N] R;
  vector[N] R1;
  vector[N] R2;
  vector[N] R3;
  vector[N] M;
  vector[N] VS;
  vector[N] HD;
  vector[N] FS;
  vector[N] mbreak_slab;
  vector[N] mbreak_if;
  vector[N] FX;

  vector[N] Y;  // PSA/FAS values

  int<lower=1,upper=NEQ> eq[N];
  int<lower=1,upper=NREG> reg[N]; // region index
  int<lower=1,upper=NREG+1> reg2[N]; // region index2 - for Cascadia small events

  real k1; // parameter for nonlinear-sie amplification
  real k2; // parameter for nonlinear-sie amplification

  real zbreak_if;
  real zbreak_slab;

  real dmb;

  // read in mean and standard deviations of prior distributions
  real mu_prior_1;
  real mu_prior_2;
  real mu_prior_2a;
  real mu_prior_3;
  real mu_prior_4;
  real mu_prior_4a;
  real mu_prior_5;
  real mu_prior_6;
  real mu_prior_7;
  real mu_prior_8;
  real mu_prior_9;
  real mu_prior_6xc;

  real<lower=0> sd_prior_1;
  real<lower=0> sd_prior_2;
  real<lower=0> sd_prior_2a;
  real<lower=0> sd_prior_3;
  real<lower=0> sd_prior_4;
  real<lower=0> sd_prior_4a;
  real<lower=0> sd_prior_5;
  real<lower=0> sd_prior_6;
  real<lower=0> sd_prior_7;
  real<lower=0> sd_prior_8;
  real<lower=0> sd_prior_9;
  real<lower=0> sd_prior_6xc;

  real<lower=0> sd_exp_1;
  real<lower=0> sd_exp_6;
  real<lower=0> sd_exp_7;

  real mu_nft_1; // parameter for near-fault term (h(m) = 10^(nft_1 + nft_2 * m))
  real mu_nft_2;

  real<lower=0> sd_nft_1;
  real<lower=0> sd_nft_2;

  // parameters for pga
  real par1_pga;
  real par1a_pga;
  real par2_pga;
  real par2a_pga;
  real par3_pga;
  real par4_pga;
  real par4a_pga;
  real par5_pga;
  real par6_pga;
  real par6b_pga;
  real par7_pga;
  real par9_pga;
  real par9a_pga;
  real par6xc_pga;

  real dzb_if_pga;
  real dzb_slab_pga;

  real nft_1_pga;
  real nft_2_pga;

  vector[NREG] par1_reg_pga; // regional constant
  vector[NREG+1] par1a_reg_pga; // regional constant
  vector[NREG] par7_reg_pga; // regional site scaling

  vector[NREG] par6x1_reg_pga; // regional attenuation
  vector[NREG] par6x2_reg_pga; // regional attenuation
  vector[NREG] par6x3_reg_pga; // regional attenuation
  vector[NREG] par61_reg_pga; // regional attenuation
  vector[NREG] par62_reg_pga; // regional attenuation
  vector[NREG] par63_reg_pga; // regional attenuation

  vector[NEQ] eqterm_pga;
  vector<lower=0>[NREG] nu;
}

transformed data {
  real<lower=0> delta_m;
  real<lower=0> delta_z;
  real refmb_lh;
  real ref_r;
  real par10;
  real ref_z_if;
  real ref_z_slab;

  real VSrock;
  real c;
  real n;
  real lnVSrock_k1;
  real k2n;

  vector[N] lnVSk1;
  vector[N] cVSk1n;
  vector[N] pga_rock;

  vector[N] dmb_rec;

  real k1_pga;
  real k2_pga;

  delta_m = 0.1;
  delta_z = 1;
  refmb_lh = 6.;
  ref_r = 100;
  ref_z_if = 15;
  ref_z_slab = 50;

  par10 = 0;

  VSrock = 1100;
  c = 1.88;
  n = 1.18;
  k2n = k2 * n;

  k1_pga = 865.;
  k2_pga = -1.186;

  lnVSrock_k1 = log(VSrock/k1);

  for(i in 1:N) {
    lnVSk1[i] = log(VS[i]/k1);
    cVSk1n[i] = c * (VS[i]/k1)^n;
  }


  // Rock PGA
  for(i in 1:N) {
    real mu1;
    real mu1a;
    real mu2;
    real mu3;
    real mu3a;
    real mu4;
    real mu5;
    real mu_geom;
    real mu_geoma;
    real siterock;

    // base model
    mu1 = (1 - FS[i]) * (par1_reg_pga[reg[i]] + logistic_hinge(M[i],mbreak_if[i],par4_pga * (mbreak_if[i] - refmb_lh),par4_pga,par5_pga,delta_m));
    mu1a = FS[i] * (par1a_reg_pga[reg2[i]] + logistic_hinge(M[i],mbreak_slab[i],par4a_pga * (mbreak_slab[i] - refmb_lh),par4a_pga,par5_pga,delta_m));

    mu4 = (1 - FS[i]) * logistic_hinge(HD[i],zbreak_if + dzb_if_pga,par9_pga * (zbreak_if + dzb_if_pga - ref_z_if),par9_pga,par10,delta_z);
    mu5 = FS[i] * logistic_hinge(HD[i],zbreak_slab + dzb_slab_pga,par9a_pga * (zbreak_slab + dzb_slab_pga - ref_z_slab),par9a_pga,par10,delta_z);

    mu_geom = (1 - FS[i]) * (par2_pga + par3_pga * M[i]) * log(R[i] + 10^(nft_1_pga + nft_2_pga * (M[i] - 6)));
    mu_geoma = FS[i] * (par2a_pga + par3_pga * M[i]) * log(R[i] + 10^(nft_1_pga + nft_2_pga * (M[i] - 6)));

    // anelastic attenuation
    mu2 = FX[i] * (par6xc_pga + par6x1_reg_pga[reg[i]] * R1[i] + par6x2_reg_pga[reg[i]] * R2[i] + par6x3_reg_pga[reg[i]] * R3[i]) + (1 - FX[i]) * (par61_reg_pga[reg[i]] * R1[i] + par62_reg_pga[reg[i]] * R2[i] + par63_reg_pga[reg[i]] * R3[i]);

    // site amplification
    siterock = (par7_reg_pga[reg[i]] + k2_pga * n) * log(VSrock/k1_pga);
    pga_rock[i] = exp(mu1 + mu1a + mu2 + mu4 + mu5 + siterock + mu_geom + mu_geoma + eqterm_pga[eq[i]]);


    // magnitude break adjustment
    if(reg[i] == 4 || reg[i] == 6)
      dmb_rec[i] = dmb;
    else
      dmb_rec[i] = 0;

  }
}

parameters {
  real par1;
  real par1a;
  real par2;
  real par2a;
  real par3;
  real<lower=0> par4;
  real<lower=0> par4a;
  real<lower=0,upper=par4> par5;
  real<upper=0> par6;
  real<upper=0> par6b;
  real par7;
  real par9;
  real par9a;
  real par6xc;
  real par1a_ca;

  real dzb_if;
  real dzb_slab;

  real nft_1;
  real nft_2;

  real<lower=0> sigma_rec;
  real<lower=0> sigma_eq;

  real<lower=0> sigma_par1; // regional standard deviation for constant
  real<lower=0> sigma_par1a; // regional standard deviation for constant
  real<lower=0> sigma_par6; // regional standard deviation for attenuation
  real<lower=0> sigma_par6b; // regional standard deviation for attenuation
  real<lower=0> sigma_par7; // regional standard deviation for site scaling

  vector[NEQ] eqterm;

  cholesky_factor_corr[9] L_reg;

  vector[NREG] z1; // regional site scaling
  vector[NREG] z1a; // regional site scaling
  vector[NREG] z7; // regional site scaling

  vector<upper=-par6b/sigma_par6b>[NREG] z6x1; // regional attenuation
  vector<upper=-par6b/sigma_par6b>[NREG] z6x2; // regional attenuation
  vector<upper=-par6b/sigma_par6b>[NREG] z6x3; // regional attenuation
  vector<upper=-par6b/sigma_par6b>[NREG] z61; // regional attenuation
  vector<upper=-par6/sigma_par6>[NREG] z62; // regional attenuation
  vector<upper=-par6b/sigma_par6b>[NREG] z63; // regional attenuation
}

transformed parameters {
  matrix[9,NREG] deltaR;
  vector[9] sigma_reg;
  matrix[9,NREG] z_reg;

  vector[NREG] par1_reg; // regional constant
  vector[NREG + 1] par1a_reg; // regional constant
  vector[NREG] par7_reg; // regional site scaling

  vector[NREG] par6x1_reg; // regional attenuation
  vector[NREG] par6x2_reg; // regional attenuation
  vector[NREG] par6x3_reg; // regional attenuation
  vector[NREG] par61_reg; // regional attenuation
  vector[NREG] par62_reg; // regional attenuation
  vector[NREG] par63_reg; // regional attenuation

  sigma_reg[1] = sigma_par1;
  sigma_reg[2] = sigma_par1a;
  sigma_reg[3] = sigma_par7;
  sigma_reg[4] = sigma_par6b;
  sigma_reg[5] = sigma_par6b;
  sigma_reg[6] = sigma_par6b;
  sigma_reg[7] = sigma_par6b;
  sigma_reg[8] = sigma_par6;
  sigma_reg[9] = sigma_par6b;

  z_reg[1] = z1';
  z_reg[2] = z1a';
  z_reg[3] = z7';
  z_reg[4] = z6x1';
  z_reg[5] = z6x2';
  z_reg[6] = z6x3';
  z_reg[7] = z61';
  z_reg[8] = z62';
  z_reg[9] = z63';

  deltaR = diag_pre_multiply(sigma_reg, L_reg) * z_reg;

  par1_reg = par1 + deltaR[1]';
  par7_reg = par7 + deltaR[3]';

  par6x1_reg = par6b + deltaR[4]';
  par6x2_reg = par6b + deltaR[5]';
  par6x3_reg = par6b + deltaR[6]';
  par61_reg = par6b + deltaR[7]';
  par62_reg = par6 + deltaR[8]';
  par63_reg = par6b + deltaR[9]';

  for(i in 1:NREG) {
    par1a_reg[i] = par1a + deltaR[2,i];
  }
  par1a_reg[NREG+1] = par1a_ca;

}

model {
  vector[N] mu_rec;
  vector[N] nu_rec;
  vector[N] mu_slab;

  sigma_eq ~ cauchy(0,0.5);
  sigma_rec ~ cauchy(0,0.5);

  sigma_par1 ~ exponential(sd_exp_1);
  sigma_par1a ~ exponential(sd_exp_1);
  sigma_par6 ~ exponential(sd_exp_6);
  sigma_par6b ~ exponential(sd_exp_6);
  sigma_par7 ~ exponential(sd_exp_7);

  par1 ~ normal(mu_prior_1,sd_prior_1);
  par1a ~ normal(mu_prior_1,sd_prior_1);
  par2 ~ normal(mu_prior_2,sd_prior_2);
  par2a ~ normal(mu_prior_2,sd_prior_2);
  par3 ~ normal(mu_prior_3,sd_prior_3);
  par4 ~ normal(mu_prior_4,sd_prior_4);
  par4a ~ normal(mu_prior_4a,sd_prior_4a);
  par5 ~ normal(mu_prior_5,sd_prior_5);
  par6 ~ normal(mu_prior_6,sd_prior_6);
  par6b ~ normal(mu_prior_6,sd_prior_6);
  par7 ~ normal(mu_prior_7,sd_prior_7);
  par9 ~ normal(mu_prior_9,sd_prior_9);
  par9a ~ normal(mu_prior_9,sd_prior_9);
  par10 ~ normal(mu_prior_9,sd_prior_9);
  par6xc ~ normal(mu_prior_6xc,sd_prior_6xc);

  dzb_if ~ normal(0,10);
  dzb_slab ~ normal(0,10);

  nft_1 ~ normal(mu_nft_1,sd_nft_1);
  nft_2 ~ normal(mu_nft_2,sd_nft_2);

  z1 ~ normal(0,1);
  z1a ~ normal(0,1);
  z7 ~ normal(0,1);

  z6x1 ~ normal(0,1);
  z6x2 ~ normal(0,1);
  z6x3 ~ normal(0,1);
  z61 ~ normal(0,1);
  z62 ~ normal(0,1);
  z63 ~ normal(0,1);

  L_reg ~ lkj_corr_cholesky(2);

  eqterm ~ normal(0,sigma_eq);

  for(i in 1:N) {
    real mu1;
    real mu1a;
    real mu2;
    real mu3;
    real mu3a;
    real mu4;
    real mu5;
    real mu_geom;
    real mu_geoma;

    // base model
    mu1 = (1 - FS[i]) * (par1_reg[reg[i]] + logistic_hinge(M[i],mbreak_if[i] + dmb_rec[i],par4 * (mbreak_if[i] + dmb_rec[i] - refmb_lh),par4,par5,delta_m));
    mu1a = FS[i] * (par1a_reg[reg2[i]] + logistic_hinge(M[i],mbreak_slab[i],par4a * (mbreak_slab[i] - refmb_lh),par4a,par5,delta_m));

    mu4 = (1 - FS[i]) * logistic_hinge(HD[i],zbreak_if + dzb_if,par9 * (zbreak_if + dzb_if - ref_z_if),par9,par10,delta_z);
    mu5 = FS[i] * logistic_hinge(HD[i],zbreak_slab + dzb_slab,par9a * (zbreak_slab + dzb_slab - ref_z_slab),par9a,par10,delta_z);

    mu_geom = (1 - FS[i]) * (par2 + par3 * M[i]) * log(R[i] + 10^(nft_1 + nft_2 * (M[i] - 6)));
    mu_geoma = FS[i] * (par2a + par3 * M[i]) * log(R[i] + 10^(nft_1 + nft_2 * (M[i] - 6)));

    // anelastic attenuation
    mu2 = FX[i] * (par6xc + par6x1_reg[reg[i]] * R1[i] + par6x2_reg[reg[i]] * R2[i] + par6x3_reg[reg[i]] * R3[i]) + (1 - FX[i]) * (par61_reg[reg[i]] * R1[i] + par62_reg[reg[i]] * R2[i] + par63_reg[reg[i]] * R3[i]);

    // site amplification
    if(VS[i] < k1) {
      mu3 = par7_reg[reg[i]] * lnVSk1[i] + k2 * (log(pga_rock[i] + cVSk1n[i]) - log(pga_rock[i] + c));
    }
    else {
      mu3 = (par7_reg[reg[i]] + k2n) * lnVSk1[i];
    }

    nu_rec[i] = nu[reg[i]];
    
    mu_rec[i] = mu1 + mu1a + mu2 + mu3 + mu4 + mu5 + mu_geom + mu_geoma + eqterm[eq[i]];
  }
  Y ~ student_t(nu_rec,mu_rec,sigma_rec);
}
