/*
	Genome-wide Efficient Mixed Model Association (GEMMA)
    Copyright (C) 2011  Xiang Zhou

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * @file gemma_lmm.cpp
 * @brief Linear Mixed Model implementation for eQTL/multi-tissue analysis
 * 
 * This file implements the statistical methods for multi-tissue eQTL analysis
 * using linear mixed models as described in the mmQTL paper.
 * 
 * MATHEMATICAL MODEL:
 * The core model (Equation 1-3 from mmqtl_extracted.md) is:
 *   Y_t = X_j*β_{j,t} + ε̂_t
 * where:
 *   - Y_t: normalized gene expression for tissue t
 *   - X_j: normalized genotype dosage for variant j
 *   - β_{j,t}: variant effect size (estimated in this code)
 *   - ε̂_t ~ N(0, K*σ²_g + σ²_e*I): residual with genetic relatedness K
 * 
 * KEY TRANSFORMATIONS:
 * The covariance matrix V = K*σ²_g + σ²_e*I = K*λ*σ²_e + σ²_e*I is decomposed as:
 *   V = U*(D*λ + I)*σ²_e*U^T
 * where K = U*D*U^T (eigendecomposition), λ = σ²_g/σ²_e (variance ratio)
 * 
 * This allows transformation to:
 *   U^T*Y_t = U^T*X_j*β_{j,t} + U^T*ε̂_t
 * making the covariance diagonal: (D*λ + I)*σ²_e
 * 
 * ESTIMATOR (Equation 4):
 *   β̂_{j,t} = (X_j^T * V^{-1} * X_j)^{-1} * (X_j^T * V^{-1} * Y_t)
 * 
 * ANALYSIS MODES:
 *   1. Wald test (REMLE)
 *   2. Likelihood ratio test (MLE)
 *   3. Score test
 *   4. All three tests combined
 */



#include <iostream>
#include <fstream>
#include <sstream>

#include <iomanip>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h> 
#include <bitset>
#include <cstring>

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"

#include "gsl/gsl_cdf.h"
#include "gsl/gsl_roots.h"
#include "gsl/gsl_min.h"
#include "gsl/gsl_integration.h"

#include "gemma_io.h"
#include "gemma_lapack.h"
#include "gemma_gzstream.h"

#ifdef FORCE_FLOAT
#include "lmm_float.h"
#else
#include "gemma_lmm.h"
#endif


using namespace std;


/**
 * @brief Copy parameters from PARAM structure to LMM object
 * @param cPar Reference to PARAM object containing analysis configuration
 * 
 * Initializes the LMM object with analysis settings including:
 * - File paths for input/output
 * - Analysis mode (Wald/LRT/Score tests)
 * - Sample and SNP filtering indicators
 * - Search range for variance ratio parameter λ = σ²_g/σ²_e
 */
void LMM::CopyFromParam (PARAM &cPar) 
{
	a_mode=cPar.a_mode;
	d_pace=cPar.d_pace;
	
	file_bfile=cPar.file_bfile;
	file_geno=cPar.file_geno;
	file_out=cPar.file_out;
	file_gene=cPar.file_gene;
	
	l_min=cPar.l_min;
	l_max=cPar.l_max;
	n_region=cPar.n_region;	
	l_mle_null=cPar.l_mle_null;
	logl_mle_H0=cPar.logl_mle_H0;
	
	time_UtX=0.0;
	time_opt=0.0;
	
	ni_total=cPar.ni_total;
	ns_total=cPar.ns_total;
	ni_test=cPar.ni_test;
	ns_test=cPar.ns_test;
	n_cvt=cPar.n_cvt;
	
	ng_total=cPar.ng_total;
	ng_test=0;
	
	indicator_idv=cPar.indicator_idv;	
	indicator_snp=cPar.indicator_snp;	
	snpInfo=cPar.snpInfo;
	
	return;
}

/**
 * @brief Copy results back to PARAM structure
 * @param cPar Reference to PARAM object to receive results
 * 
 * Updates PARAM with:
 * - Timing information for matrix operations and optimization
 * - Number of genes/SNPs tested
 */
void LMM::CopyToParam (PARAM &cPar) 
{
	cPar.time_UtX=time_UtX;
	cPar.time_opt=time_opt;	
	
	cPar.ng_test=ng_test;
	
	return;
}


/**
 * @brief Write association test results to output file
 * 
 * Output format depends on analysis mode:
 * - Mode 1: beta, SE, lambda_remle, p_wald
 * - Mode 2: lambda_mle, p_lrt  
 * - Mode 3: beta, SE, p_score
 * - Mode 4: All statistics combined
 * 
 * For SNP analysis: includes chr, rs, position, alleles, MAF
 * For gene analysis: includes gene ID only
 */
void LMM::WriteFiles () 
{
	string file_str;
	file_str="./output/"+file_out;
	file_str+=".assoc.txt";

	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}

	if (!file_gene.empty()) {
		outfile<<"geneID"<<"\t";
		
		if (a_mode==1) {
			outfile<<"beta"<<"\t"<<"se"<<"\t"<<"l_remle"<<"\t"<<"p_wald"<<endl;
		} else if (a_mode==2) {
			outfile<<"l_mle"<<"\t"<<"p_lrt"<<endl;
		} else if (a_mode==3) {
			outfile<<"beta"<<"\t"<<"se"<<"\t"<<"p_score"<<endl;
		} else if (a_mode==4) {
			outfile<<"beta"<<"\t"<<"se"<<"\t"<<"l_remle"<<"\t"<<"l_mle"<<"\t"<<"p_wald"<<"\t"<<"p_lrt"<<"\t"<<"p_score"<<endl;
		} else {}
				
		for (vector<SUMSTAT>::size_type t=0; t<sumStat.size(); ++t) {	
			outfile<<snpInfo[t].rs_number<<"\t";
			
			if (a_mode==1) {
				outfile<<scientific<<setprecision(6)<<sumStat[t].beta<<"\t"<<sumStat[t].se<<"\t"<<sumStat[t].lambda_remle<<"\t"<<sumStat[t].p_wald <<endl;
			} else if (a_mode==2) {
				outfile<<scientific<<setprecision(6)<<sumStat[t].lambda_mle<<"\t"<<sumStat[t].p_lrt<<endl;
			} else if (a_mode==3) {
				outfile<<scientific<<setprecision(6)<<sumStat[t].beta<<"\t"<<sumStat[t].se<<"\t"<<sumStat[t].p_score<<endl;
			} else if (a_mode==4) {
				outfile<<scientific<<setprecision(6)<<sumStat[t].beta<<"\t"<<sumStat[t].se<<"\t"<<sumStat[t].lambda_remle<<"\t"<<sumStat[t].lambda_mle<<"\t"<<sumStat[t].p_wald <<"\t"<<sumStat[t].p_lrt<<"\t"<<sumStat[t].p_score<<endl;
			} else {}
		}	
	}  else {
		outfile<<"chr"<<"\t"<<"rs"<<"\t"<<"ps"<<"\t"<<"n_miss"<<"\t"<<"allele1"<<"\t"<<"allele0"<<"\t"<<"af"<<"\t";
		
		if (a_mode==1) {
			outfile<<"beta"<<"\t"<<"se"<<"\t"<<"l_remle"<<"\t"<<"p_wald"<<endl;
		} else if (a_mode==2) {
			outfile<<"l_mle"<<"\t"<<"p_lrt"<<endl;
		} else if (a_mode==3) {
			outfile<<"beta"<<"\t"<<"se"<<"\t"<<"p_score"<<endl;
		} else if (a_mode==4) {
			outfile<<"beta"<<"\t"<<"se"<<"\t"<<"l_remle"<<"\t"<<"l_mle"<<"\t"<<"p_wald"<<"\t"<<"p_lrt"<<"\t"<<"p_score"<<endl;
		} else {}
		
		size_t t=0;
		for (size_t i=0; i<snpInfo.size(); ++i) {
			if (indicator_snp[i]==0) {continue;}
			
			outfile<<snpInfo[i].chr<<"\t"<<snpInfo[i].rs_number<<"\t"<<snpInfo[i].base_position<<"\t"<<snpInfo[i].n_miss<<"\t"<<snpInfo[i].a_minor<<"\t"<<snpInfo[i].a_major<<"\t"<<fixed<<setprecision(3)<<snpInfo[i].maf<<"\t";
			
			if (a_mode==1) {
				outfile<<scientific<<setprecision(6)<<sumStat[t].beta<<"\t"<<sumStat[t].se<<"\t"<<sumStat[t].lambda_remle<<"\t"<<sumStat[t].p_wald <<endl;
			} else if (a_mode==2) {
				outfile<<scientific<<setprecision(6)<<sumStat[t].lambda_mle<<"\t"<<sumStat[t].p_lrt<<endl;
			} else if (a_mode==3) {
				outfile<<scientific<<setprecision(6)<<sumStat[t].beta<<"\t"<<sumStat[t].se<<"\t"<<sumStat[t].p_score<<endl;
			} else if (a_mode==4) {
				outfile<<scientific<<setprecision(6)<<sumStat[t].beta<<"\t"<<sumStat[t].se<<"\t"<<sumStat[t].lambda_remle<<"\t"<<sumStat[t].lambda_mle<<"\t"<<sumStat[t].p_wald <<"\t"<<sumStat[t].p_lrt<<"\t"<<sumStat[t].p_score<<endl;
			} else {}
			t++;
		}
	}
	
		
	outfile.close();
	outfile.clear();
	return;
}







/**
 * @brief Map two indices (a,b) to a single index for symmetric matrix storage
 * @param a First index (1 to n_cvt+2)
 * @param b Second index (1 to n_cvt+2)
 * @param n_cvt Number of covariates
 * @return Single index for compact storage of upper triangular matrix
 * 
 * This function enables efficient storage of symmetric matrices by storing
 * only the upper triangle. Used for storing products like U^T*W*W^T*U where
 * W includes covariates and test variants.
 * 
 * Index mapping:
 *   a=1,b=1 → 0
 *   a=1,b=2 → 1  
 *   a=2,b=2 → n_cvt+1
 *   etc.
 */
size_t GetabIndex (const size_t a, const size_t b, const size_t n_cvt) {
	if (a>n_cvt+2 || b>n_cvt+2 || a<=0 || b<=0) {cout<<"error in GetabIndex."<<endl; return 0;}
	size_t index;
	size_t l, h;
	if (b>a) {l=a; h=b;} else {l=b; h=a;}
	
	size_t n=n_cvt+2;
	index=(2*n-l+2)*(l-1)/2+h-l;	
	
	return index;
}

/**
 * @brief Calculate P matrices for association tests with covariates
 * @param n_cvt Number of covariates
 * @param e_mode Error mode (0 or non-zero for different variance structures)
 * @param Hi_eval Vector containing (λ*D + I)^{-1} where D are eigenvalues
 * @param Uab Matrix containing transformed cross-products U^T*[W,X,Y]*[W,X,Y]^T*U
 * @param ab Vector containing original cross-products [W,X,Y]*[W,X,Y]^T (for e_mode!=0)
 * @param Pab Output matrix of adjusted P statistics after accounting for covariates
 * 
 * MATHEMATICAL BACKGROUND:
 * This computes projection matrices that account for covariates in the mixed model.
 * The P matrix represents the variance-adjusted inner products after marginalizing
 * over covariates using the Woodbury matrix identity iteratively.
 * 
 * For the LMM: V^{-1} = H^{-1} = (λ*K + I)^{-1}
 * In eigenspace: U^T*H*U = λ*D + I (diagonal)
 * 
 * The function computes statistics like:
 *   P_{yy} = Y^T * (H^{-1} - H^{-1}*W*(W^T*H^{-1}*W)^{-1}*W^T*H^{-1}) * Y
 * which is Y^T*V^{-1}*Y adjusted for covariates W.
 * 
 * This is essential for calculating test statistics while controlling for
 * confounders (covariates).
 */
void CalcPab (const size_t n_cvt, const size_t e_mode, const gsl_vector *Hi_eval, const gsl_matrix *Uab, const gsl_vector *ab, gsl_matrix *Pab)
{
	size_t index_ab, index_aw, index_bw, index_ww;
	double p_ab;
	double ps_ab, ps_aw, ps_bw, ps_ww;
	
	for (size_t p=0; p<=n_cvt+1; ++p) {
		for (size_t a=p+1; a<=n_cvt+2; ++a) {
			for (size_t b=a; b<=n_cvt+2; ++b) {
				index_ab=GetabIndex (a, b, n_cvt);
				if (p==0) {			
					gsl_vector_const_view Uab_col=gsl_matrix_const_column (Uab, index_ab);
					gsl_blas_ddot (Hi_eval, &Uab_col.vector, &p_ab);
					if (e_mode!=0) {p_ab=gsl_vector_get (ab, index_ab)-p_ab;}
					gsl_matrix_set (Pab, 0, index_ab, p_ab);
				}
				else {
					index_aw=GetabIndex (a, p, n_cvt);
					index_bw=GetabIndex (b, p, n_cvt);
					index_ww=GetabIndex (p, p, n_cvt);
					
					ps_ab=gsl_matrix_get (Pab, p-1, index_ab);
					ps_aw=gsl_matrix_get (Pab, p-1, index_aw);
					ps_bw=gsl_matrix_get (Pab, p-1, index_bw);
					ps_ww=gsl_matrix_get (Pab, p-1, index_ww);
					
					p_ab=ps_ab-ps_aw*ps_bw/ps_ww;
					gsl_matrix_set (Pab, p, index_ab, p_ab);
				}
			}
		}
	}
	return;
}

/**
 * @brief Calculate PPab matrices (second-order derivatives w.r.t. λ)
 * @param n_cvt Number of covariates
 * @param e_mode Error mode
 * @param HiHi_eval Vector containing (λ*D + I)^{-2}
 * @param Uab Transformed cross-products
 * @param ab Original cross-products
 * @param Pab First-order P matrices
 * @param PPab Output: second derivatives of P matrices w.r.t. λ
 * 
 * MATHEMATICAL PURPOSE:
 * Computes d²/dλ² of the projection matrices. These second derivatives are
 * needed for:
 * 1. Newton-Raphson optimization of λ (requires Hessian)
 * 2. Standard error estimation via Fisher information
 * 
 * The variance ratio λ = σ²_g/σ²_e is estimated by maximizing the (restricted)
 * likelihood. This requires computing derivatives of log-likelihood w.r.t. λ,
 * which in turn requires derivatives of quadratic forms like Y^T*P*Y.
 */
void CalcPPab (const size_t n_cvt, const size_t e_mode, const gsl_vector *HiHi_eval, const gsl_matrix *Uab, const gsl_vector *ab, const gsl_matrix *Pab, gsl_matrix *PPab)
{
	size_t index_ab, index_aw, index_bw, index_ww;
	double p2_ab;
	double ps2_ab, ps_aw, ps_bw, ps_ww, ps2_aw, ps2_bw, ps2_ww;
	
	for (size_t p=0; p<=n_cvt+1; ++p) {
		for (size_t a=p+1; a<=n_cvt+2; ++a) {
			for (size_t b=a; b<=n_cvt+2; ++b) {
				index_ab=GetabIndex (a, b, n_cvt);
				if (p==0) {					
					gsl_vector_const_view Uab_col=gsl_matrix_const_column (Uab, index_ab);
					gsl_blas_ddot (HiHi_eval, &Uab_col.vector, &p2_ab);
					if (e_mode!=0) {p2_ab=p2_ab-gsl_vector_get (ab, index_ab)+2.0*gsl_matrix_get (Pab, 0, index_ab);}
					gsl_matrix_set (PPab, 0, index_ab, p2_ab);
				}
				else {
					index_aw=GetabIndex (a, p, n_cvt);
					index_bw=GetabIndex (b, p, n_cvt);
					index_ww=GetabIndex (p, p, n_cvt);
					
					ps2_ab=gsl_matrix_get (PPab, p-1, index_ab);
					ps_aw=gsl_matrix_get (Pab, p-1, index_aw);
					ps_bw=gsl_matrix_get (Pab, p-1, index_bw);
					ps_ww=gsl_matrix_get (Pab, p-1, index_ww);
					ps2_aw=gsl_matrix_get (PPab, p-1, index_aw);
					ps2_bw=gsl_matrix_get (PPab, p-1, index_bw);
					ps2_ww=gsl_matrix_get (PPab, p-1, index_ww);
					
					p2_ab=ps2_ab+ps_aw*ps_bw*ps2_ww/(ps_ww*ps_ww);
					p2_ab-=(ps_aw*ps2_bw+ps_bw*ps2_aw)/ps_ww;
					gsl_matrix_set (PPab, p, index_ab, p2_ab);
					
				}
			}
		}
	}
	return;
}

/**
 * @brief Calculate PPPab matrices (third-order derivatives w.r.t. λ)
 * @param n_cvt Number of covariates
 * @param e_mode Error mode
 * @param HiHiHi_eval Vector containing (λ*D + I)^{-3}
 * @param Uab Transformed cross-products
 * @param ab Original cross-products
 * @param Pab First-order P matrices
 * @param PPab Second-order P matrices
 * @param PPPab Output: third derivatives of P matrices w.r.t. λ
 * 
 * MATHEMATICAL PURPOSE:
 * Third derivatives are used for more accurate optimization and checking
 * convergence properties. They help ensure the likelihood maximum is not
 * a saddle point and can improve numerical stability in difficult cases.
 */
void CalcPPPab (const size_t n_cvt, const size_t e_mode, const gsl_vector *HiHiHi_eval, const gsl_matrix *Uab, const gsl_vector *ab, const gsl_matrix *Pab, const gsl_matrix *PPab, gsl_matrix *PPPab)
{
	size_t index_ab, index_aw, index_bw, index_ww;
	double p3_ab;
	double ps3_ab, ps_aw, ps_bw, ps_ww, ps2_aw, ps2_bw, ps2_ww, ps3_aw, ps3_bw, ps3_ww;
	
	for (size_t p=0; p<=n_cvt+1; ++p) {
		for (size_t a=p+1; a<=n_cvt+2; ++a) {
			for (size_t b=a; b<=n_cvt+2; ++b) {
				index_ab=GetabIndex (a, b, n_cvt);
				if (p==0) {					
					gsl_vector_const_view Uab_col=gsl_matrix_const_column (Uab, index_ab);
					gsl_blas_ddot (HiHiHi_eval, &Uab_col.vector, &p3_ab);
					if (e_mode!=0) {p3_ab=gsl_vector_get (ab, index_ab)-p3_ab+3.0*gsl_matrix_get (PPab, 0, index_ab)-3.0*gsl_matrix_get (Pab, 0, index_ab);}
					gsl_matrix_set (PPPab, 0, index_ab, p3_ab);
				}
				else {
					index_aw=GetabIndex (a, p, n_cvt);
					index_bw=GetabIndex (b, p, n_cvt);
					index_ww=GetabIndex (p, p, n_cvt);
					
					ps3_ab=gsl_matrix_get (PPPab, p-1, index_ab);
					ps_aw=gsl_matrix_get (Pab, p-1, index_aw);
					ps_bw=gsl_matrix_get (Pab, p-1, index_bw);
					ps_ww=gsl_matrix_get (Pab, p-1, index_ww);
					ps2_aw=gsl_matrix_get (PPab, p-1, index_aw);
					ps2_bw=gsl_matrix_get (PPab, p-1, index_bw);
					ps2_ww=gsl_matrix_get (PPab, p-1, index_ww);
					ps3_aw=gsl_matrix_get (PPPab, p-1, index_aw);
					ps3_bw=gsl_matrix_get (PPPab, p-1, index_bw);
					ps3_ww=gsl_matrix_get (PPPab, p-1, index_ww);
					
					p3_ab=ps3_ab-ps_aw*ps_bw*ps2_ww*ps2_ww/(ps_ww*ps_ww*ps_ww);
					p3_ab-=(ps_aw*ps3_bw+ps_bw*ps3_aw+ps2_aw*ps2_bw)/ps_ww;
					p3_ab+=(ps_aw*ps2_bw*ps2_ww+ps_bw*ps2_aw*ps2_ww+ps_aw*ps_bw*ps3_ww)/(ps_ww*ps_ww);
					
					gsl_matrix_set (PPPab, p, index_ab, p3_ab);
				}
			}
		}
	}
	return;
}


/**
 * @brief Log-likelihood function for variance ratio λ
 * @param l Current value of λ = σ²_g/σ²_e (variance ratio parameter)
 * @param params Pointer to FUNC_PARAM structure with data and settings
 * @return Log-likelihood value at λ=l
 * 
 * MATHEMATICAL DERIVATION:
 * For the linear mixed model Y = Xβ + ε where ε ~ N(0, V),
 * V = λ*σ²_e*K + σ²_e*I = σ²_e*(λ*K + I)
 * 
 * The log-likelihood is:
 *   ℓ(λ) = -n/2*log(2π) - 1/2*log|V| - 1/2*Y^T*V^{-1}*Y
 * 
 * After eigendecomposition V = U*(λ*D + I)*σ²_e*U^T and profiling out σ²_e:
 *   ℓ(λ) ∝ -1/2*log|λ*D + I| - n/2*log(Y^T*(λ*D + I)^{-1}*Y)
 * 
 * This function evaluates this profile log-likelihood for optimization.
 * The MLE of λ maximizes this function (a_mode=2).
 */
double LogL_f (double l, void *params)
{
	FUNC_PARAM *p=(FUNC_PARAM *) params;
	size_t n_cvt=p->n_cvt;
	size_t ni_test=p->ni_test;	
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	size_t nc_total;
	if (p->calc_null==true) {nc_total=n_cvt;} else {nc_total=n_cvt+1;}
	
	double f=0.0, logdet_h=0.0, d;
	size_t index_yy;
	
	gsl_matrix *Pab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_vector *Hi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *v_temp=gsl_vector_alloc((p->eval)->size);
				
	gsl_vector_memcpy (v_temp, p->eval);
	gsl_vector_scale (v_temp, l);
	if (p->e_mode==0) {gsl_vector_set_all (Hi_eval, 1.0);} else {gsl_vector_memcpy (Hi_eval, v_temp);}
	gsl_vector_add_constant (v_temp, 1.0);
	gsl_vector_div (Hi_eval, v_temp);	
	
	for (size_t i=0; i<(p->eval)->size; ++i) {
		d=gsl_vector_get (v_temp, i);
		logdet_h+=log(fabs(d));
	}	
	
	CalcPab (n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, Pab);	
	
	double c=0.5*(double)ni_test*(log((double)ni_test)-log(2*M_PI)-1.0);
	
	index_yy=GetabIndex (n_cvt+2, n_cvt+2, n_cvt);	
	double P_yy=gsl_matrix_get (Pab, nc_total, index_yy);
	f=c-0.5*logdet_h-0.5*(double)ni_test*log(P_yy);
	
	gsl_matrix_free (Pab);
	gsl_vector_free (Hi_eval);
	gsl_vector_free (v_temp);
	return f;
}

 
/**
 * @brief First derivative of log-likelihood w.r.t. λ
 * @param l Current value of λ
 * @param params Pointer to FUNC_PARAM structure
 * @return dℓ/dλ evaluated at λ=l
 * 
 * MATHEMATICAL DERIVATION:
 * The first derivative is:
 *   dℓ/dλ = -1/2 * tr(H^{-1}*K) + 1/2 * Y^T*H^{-1}*K*H^{-1}*Y / (Y^T*P*Y)
 * 
 * where H = λ*K + I, P is the covariate-adjusted projection.
 * 
 * This gradient is used for:
 * 1. Finding roots (MLE/REMLE where dℓ/dλ = 0)
 * 2. Gradient-based optimization
 * 3. Checking convergence
 */
double LogL_dev1 (double l, void *params)
{
	FUNC_PARAM *p=(FUNC_PARAM *) params;	
	size_t n_cvt=p->n_cvt;
	size_t ni_test=p->ni_test;	
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	size_t nc_total;
	if (p->calc_null==true) {nc_total=n_cvt;} else {nc_total=n_cvt+1;}
	
	double dev1=0.0, trace_Hi=0.0;
	size_t index_yy;
	
	gsl_matrix *Pab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_matrix *PPab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_vector *Hi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *HiHi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *v_temp=gsl_vector_alloc((p->eval)->size);
	
	gsl_vector_memcpy (v_temp, p->eval);
	gsl_vector_scale (v_temp, l);
	if (p->e_mode==0) {gsl_vector_set_all (Hi_eval, 1.0);} else {gsl_vector_memcpy (Hi_eval, v_temp);}
	gsl_vector_add_constant (v_temp, 1.0);
	gsl_vector_div (Hi_eval, v_temp);
	
	gsl_vector_memcpy (HiHi_eval, Hi_eval);
	gsl_vector_mul (HiHi_eval, Hi_eval);	
	
	gsl_vector_set_all (v_temp, 1.0);
	gsl_blas_ddot (Hi_eval, v_temp, &trace_Hi);
		
	if (p->e_mode!=0) {trace_Hi=(double)ni_test-trace_Hi;}
	
	CalcPab (n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, Pab);	
	CalcPPab (n_cvt, p->e_mode, HiHi_eval, p->Uab, p->ab, Pab, PPab);	
	
	double trace_HiK=((double)ni_test-trace_Hi)/l;	
	
	index_yy=GetabIndex (n_cvt+2, n_cvt+2, n_cvt);
	
	double P_yy=gsl_matrix_get (Pab, nc_total, index_yy);
	double PP_yy=gsl_matrix_get (PPab, nc_total, index_yy);
	double yPKPy=(P_yy-PP_yy)/l;	
	dev1=-0.5*trace_HiK+0.5*(double)ni_test*yPKPy/P_yy;
			
	gsl_matrix_free (Pab);
	gsl_matrix_free (PPab);
	gsl_vector_free (Hi_eval);
	gsl_vector_free (HiHi_eval);
	gsl_vector_free (v_temp);	
	
	return dev1;
}
	
/**
 * @brief Second derivative of log-likelihood w.r.t. λ (Hessian)
 * @param l Current value of λ
 * @param params Pointer to FUNC_PARAM structure
 * @return d²ℓ/dλ² evaluated at λ=l
 * 
 * MATHEMATICAL PURPOSE:
 * The Hessian (second derivative) is used for:
 * 1. Newton-Raphson optimization: λ_new = λ_old - (d²ℓ/dλ²)^{-1} * (dℓ/dλ)
 * 2. Fisher information: I(λ) = -E[d²ℓ/dλ²]
 * 3. Standard error: SE(λ) = sqrt(-1/(d²ℓ/dλ²))
 * 4. Checking if critical point is maximum (need d²ℓ/dλ² < 0)
 * 
 * A negative Hessian at the MLE confirms we have a maximum, not a minimum.
 */
double LogL_dev2 (double l, void *params)
{
	FUNC_PARAM *p=(FUNC_PARAM *) params;	
	size_t n_cvt=p->n_cvt;
	size_t ni_test=p->ni_test;	
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	size_t nc_total;
	if (p->calc_null==true) {nc_total=n_cvt;} else {nc_total=n_cvt+1;}
	
	double dev2=0.0, trace_Hi=0.0, trace_HiHi=0.0;
	size_t index_yy;
	
	gsl_matrix *Pab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_matrix *PPab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_matrix *PPPab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_vector *Hi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *HiHi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *HiHiHi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *v_temp=gsl_vector_alloc((p->eval)->size);
	
	gsl_vector_memcpy (v_temp, p->eval);
	gsl_vector_scale (v_temp, l);
	if (p->e_mode==0) {gsl_vector_set_all (Hi_eval, 1.0);} else {gsl_vector_memcpy (Hi_eval, v_temp);}
	gsl_vector_add_constant (v_temp, 1.0);
	gsl_vector_div (Hi_eval, v_temp);
		
	gsl_vector_memcpy (HiHi_eval, Hi_eval);
	gsl_vector_mul (HiHi_eval, Hi_eval);	
	gsl_vector_memcpy (HiHiHi_eval, HiHi_eval);
	gsl_vector_mul (HiHiHi_eval, Hi_eval);
	
	gsl_vector_set_all (v_temp, 1.0);
	gsl_blas_ddot (Hi_eval, v_temp, &trace_Hi);
	gsl_blas_ddot (HiHi_eval, v_temp, &trace_HiHi);
	
	if (p->e_mode!=0) {		
		trace_Hi=(double)ni_test-trace_Hi;
		trace_HiHi=2*trace_Hi+trace_HiHi-(double)ni_test;
	}
	
	CalcPab (n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, Pab);	
	CalcPPab (n_cvt, p->e_mode, HiHi_eval, p->Uab, p->ab, Pab, PPab);	
	CalcPPPab (n_cvt, p->e_mode, HiHiHi_eval, p->Uab, p->ab, Pab, PPab, PPPab);	
	
	double trace_HiKHiK=((double)ni_test+trace_HiHi-2*trace_Hi)/(l*l);
	
	index_yy=GetabIndex (n_cvt+2, n_cvt+2, n_cvt);
	double P_yy=gsl_matrix_get (Pab, nc_total, index_yy);
	double PP_yy=gsl_matrix_get (PPab, nc_total, index_yy);
	double PPP_yy=gsl_matrix_get (PPPab, nc_total, index_yy);		
		
	double yPKPy=(P_yy-PP_yy)/l;
	double yPKPKPy=(P_yy+PPP_yy-2.0*PP_yy)/(l*l);
		
	dev2=0.5*trace_HiKHiK-0.5*(double)ni_test*(2.0*yPKPKPy*P_yy-yPKPy*yPKPy)/(P_yy*P_yy);
		
	gsl_matrix_free (Pab);
	gsl_matrix_free (PPab);
	gsl_matrix_free (PPPab);
	gsl_vector_free (Hi_eval);
	gsl_vector_free (HiHi_eval);
	gsl_vector_free (HiHiHi_eval);
	gsl_vector_free (v_temp);	
	
	return dev2;
}
	
	
/**
 * @brief Compute both first and second derivatives simultaneously
 * @param l Current value of λ
 * @param params Pointer to FUNC_PARAM structure
 * @param dev1 Output: first derivative dℓ/dλ
 * @param dev2 Output: second derivative d²ℓ/dλ²
 * 
 * EFFICIENCY:
 * Computing both derivatives together is more efficient than calling
 * LogL_dev1 and LogL_dev2 separately because intermediate matrices
 * (Pab, PPab, PPPab) are reused.
 * 
 * This is the preferred function for Newton-Raphson optimization which
 * requires both gradient and Hessian at each iteration.
 */	
void LogL_dev12 (double l, void *params, double *dev1, double *dev2)
{
	FUNC_PARAM *p=(FUNC_PARAM *) params;
	size_t n_cvt=p->n_cvt;
	size_t ni_test=p->ni_test;	
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	size_t nc_total;
	if (p->calc_null==true) {nc_total=n_cvt;} else {nc_total=n_cvt+1;}
	
	double trace_Hi=0.0, trace_HiHi=0.0;
	size_t index_yy;
	
	gsl_matrix *Pab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_matrix *PPab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_matrix *PPPab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_vector *Hi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *HiHi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *HiHiHi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *v_temp=gsl_vector_alloc((p->eval)->size);
	
	gsl_vector_memcpy (v_temp, p->eval);
	gsl_vector_scale (v_temp, l);
	if (p->e_mode==0) {gsl_vector_set_all (Hi_eval, 1.0);} else {gsl_vector_memcpy (Hi_eval, v_temp);}
	gsl_vector_add_constant (v_temp, 1.0);
	gsl_vector_div (Hi_eval, v_temp);
		
	gsl_vector_memcpy (HiHi_eval, Hi_eval);
	gsl_vector_mul (HiHi_eval, Hi_eval);	
	gsl_vector_memcpy (HiHiHi_eval, HiHi_eval);
	gsl_vector_mul (HiHiHi_eval, Hi_eval);
	
	gsl_vector_set_all (v_temp, 1.0);
	gsl_blas_ddot (Hi_eval, v_temp, &trace_Hi);
	gsl_blas_ddot (HiHi_eval, v_temp, &trace_HiHi);
	
	if (p->e_mode!=0) {	
		trace_Hi=(double)ni_test-trace_Hi;
		trace_HiHi=2*trace_Hi+trace_HiHi-(double)ni_test;
	}
	
	CalcPab (n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, Pab);	
	CalcPPab (n_cvt, p->e_mode, HiHi_eval, p->Uab, p->ab, Pab, PPab);	
	CalcPPPab (n_cvt, p->e_mode, HiHiHi_eval, p->Uab, p->ab, Pab, PPab, PPPab);	
	
	double trace_HiK=((double)ni_test-trace_Hi)/l;
	double trace_HiKHiK=((double)ni_test+trace_HiHi-2*trace_Hi)/(l*l);
	
	index_yy=GetabIndex (n_cvt+2, n_cvt+2, n_cvt);
	
	double P_yy=gsl_matrix_get (Pab, nc_total, index_yy);
	double PP_yy=gsl_matrix_get (PPab, nc_total, index_yy);
	double PPP_yy=gsl_matrix_get (PPPab, nc_total, index_yy);		
		
	double yPKPy=(P_yy-PP_yy)/l;	
	double yPKPKPy=(P_yy+PPP_yy-2.0*PP_yy)/(l*l);
		
	*dev1=-0.5*trace_HiK+0.5*(double)ni_test*yPKPy/P_yy;
	*dev2=0.5*trace_HiKHiK-0.5*(double)ni_test*(2.0*yPKPKPy*P_yy-yPKPy*yPKPy)/(P_yy*P_yy);
			
	gsl_matrix_free (Pab);
	gsl_matrix_free (PPab);
	gsl_matrix_free (PPPab);
	gsl_vector_free (Hi_eval);
	gsl_vector_free (HiHi_eval);
	gsl_vector_free (HiHiHi_eval);
	gsl_vector_free (v_temp);	
	
	return;
}


/**
 * @brief Restricted log-likelihood (REML) function
 * @param l Current value of λ
 * @param params Pointer to FUNC_PARAM structure
 * @return Restricted log-likelihood value at λ=l
 * 
 * MATHEMATICAL BACKGROUND:
 * REML differs from ML by integrating out fixed effects (β), treating them
 * as nuisance parameters. This reduces bias in variance component estimation.
 * 
 * The restricted likelihood is proportional to:
 *   ℓ_R(λ) = ℓ(λ) - 1/2*log|X^T*V^{-1}*X|
 * 
 * where the second term accounts for uncertainty in β.
 * 
 * REML is preferred for variance estimation (a_mode=1) because:
 * 1. Unbiased for variance components in balanced designs
 * 2. Accounts for degrees of freedom lost estimating fixed effects
 * 3. More conservative in small samples
 * 
 * The REMLE (Restricted Maximum Likelihood Estimator) maximizes this function.
 */
double LogRL_f (double l, void *params)
{
	FUNC_PARAM *p=(FUNC_PARAM *) params;	
	size_t n_cvt=p->n_cvt;
	size_t ni_test=p->ni_test;	
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	double df;
	size_t nc_total;
	if (p->calc_null==true) {nc_total=n_cvt; df=(double)ni_test-(double)n_cvt; }
	else {nc_total=n_cvt+1; df=(double)ni_test-(double)n_cvt-1.0;}
	
	double f=0.0, logdet_h=0.0, logdet_hiw=0.0, d;
	size_t index_ww;
	
	gsl_matrix *Pab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_matrix *Iab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_vector *Hi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *v_temp=gsl_vector_alloc((p->eval)->size);
	
	gsl_vector_memcpy (v_temp, p->eval);
	gsl_vector_scale (v_temp, l);
	if (p->e_mode==0) {gsl_vector_set_all (Hi_eval, 1.0);} else {gsl_vector_memcpy (Hi_eval, v_temp);}
	gsl_vector_add_constant (v_temp, 1.0);
	gsl_vector_div (Hi_eval, v_temp);	
	
	for (size_t i=0; i<(p->eval)->size; ++i) {
		d=gsl_vector_get (v_temp, i);
		logdet_h+=log(fabs(d));
	}
	
	CalcPab (n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, Pab);	
	gsl_vector_set_all (v_temp, 1.0);
	CalcPab (n_cvt, p->e_mode, v_temp, p->Uab, p->ab, Iab);	
	
	//calculate |WHiW|-|WW|
	logdet_hiw=0.0;
	for (size_t i=0; i<nc_total; ++i) {
		index_ww=GetabIndex (i+1, i+1, n_cvt);
		d=gsl_matrix_get (Pab, i, index_ww);
		logdet_hiw+=log(d);
		d=gsl_matrix_get (Iab, i, index_ww);
		logdet_hiw-=log(d);
	}
	index_ww=GetabIndex (n_cvt+2, n_cvt+2, n_cvt);	
	double P_yy=gsl_matrix_get (Pab, nc_total, index_ww);
	
	double c=0.5*df*(log(df)-log(2*M_PI)-1.0);		
	f=c-0.5*logdet_h-0.5*logdet_hiw-0.5*df*log(P_yy);
		
	gsl_matrix_free (Pab);
	gsl_matrix_free (Iab);
	gsl_vector_free (Hi_eval);
	gsl_vector_free (v_temp);
	return f;
}


/**
 * @brief First derivative of restricted log-likelihood
 * @param l Current value of λ
 * @param params Pointer to FUNC_PARAM structure
 * @return dℓ_R/dλ evaluated at λ=l
 * 
 * Used to find REMLE where dℓ_R/dλ = 0. The REMLE is preferred over MLE
 * for Wald tests (a_mode=1) as it provides better variance estimates.
 */
double LogRL_dev1 (double l, void *params)
{
	FUNC_PARAM *p=(FUNC_PARAM *) params;	
	size_t n_cvt=p->n_cvt;
	size_t ni_test=p->ni_test;	
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	double df;
	size_t nc_total;
	if (p->calc_null==true) {nc_total=n_cvt; df=(double)ni_test-(double)n_cvt; }
	else {nc_total=n_cvt+1; df=(double)ni_test-(double)n_cvt-1.0;}
	
	double dev1=0.0, trace_Hi=0.0;
	size_t index_ww;
	
	gsl_matrix *Pab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_matrix *PPab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_vector *Hi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *HiHi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *v_temp=gsl_vector_alloc((p->eval)->size);
	
	gsl_vector_memcpy (v_temp, p->eval);
	gsl_vector_scale (v_temp, l);
	if (p->e_mode==0) {gsl_vector_set_all (Hi_eval, 1.0);} else {gsl_vector_memcpy (Hi_eval, v_temp);}
	gsl_vector_add_constant (v_temp, 1.0);
	gsl_vector_div (Hi_eval, v_temp);
	
	gsl_vector_memcpy (HiHi_eval, Hi_eval);
	gsl_vector_mul (HiHi_eval, Hi_eval);	
	
	gsl_vector_set_all (v_temp, 1.0);
	gsl_blas_ddot (Hi_eval, v_temp, &trace_Hi);
	
	if (p->e_mode!=0) {	
		trace_Hi=(double)ni_test-trace_Hi;
	}
	
	CalcPab (n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, Pab);	
	CalcPPab (n_cvt, p->e_mode, HiHi_eval, p->Uab, p->ab, Pab, PPab);	
	
	//calculate tracePK and trace PKPK
	double trace_P=trace_Hi;
	double ps_ww, ps2_ww;
	for (size_t i=0; i<nc_total; ++i) {
		index_ww=GetabIndex (i+1, i+1, n_cvt);
		ps_ww=gsl_matrix_get (Pab, i, index_ww);
		ps2_ww=gsl_matrix_get (PPab, i, index_ww);
		trace_P-=ps2_ww/ps_ww;
	}
	double trace_PK=(df-trace_P)/l;
	
	//calculate yPKPy, yPKPKPy
	index_ww=GetabIndex (n_cvt+2, n_cvt+2, n_cvt);
	double P_yy=gsl_matrix_get (Pab, nc_total, index_ww);
	double PP_yy=gsl_matrix_get (PPab, nc_total, index_ww);		
	double yPKPy=(P_yy-PP_yy)/l;	
	
	dev1=-0.5*trace_PK+0.5*df*yPKPy/P_yy;	
			
	gsl_matrix_free (Pab);
	gsl_matrix_free (PPab);
	gsl_vector_free (Hi_eval);
	gsl_vector_free (HiHi_eval);
	gsl_vector_free (v_temp);	
	
	return dev1;
}



/**
 * @brief Second derivative of restricted log-likelihood (REML Hessian)
 * @param l Current value of λ
 * @param params Pointer to FUNC_PARAM structure
 * @return d²ℓ_R/dλ² evaluated at λ=l
 * 
 * Used for Newton-Raphson optimization of REMLE and computing standard
 * errors of variance component estimates.
 */
double LogRL_dev2 (double l, void *params)
{
	FUNC_PARAM *p=(FUNC_PARAM *) params;	
	size_t n_cvt=p->n_cvt;
	size_t ni_test=p->ni_test;	
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	double df;
	size_t nc_total;
	if (p->calc_null==true) {nc_total=n_cvt; df=(double)ni_test-(double)n_cvt; }
	else {nc_total=n_cvt+1; df=(double)ni_test-(double)n_cvt-1.0;}
	
	double dev2=0.0, trace_Hi=0.0, trace_HiHi=0.0;
	size_t index_ww;
	
	gsl_matrix *Pab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_matrix *PPab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_matrix *PPPab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_vector *Hi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *HiHi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *HiHiHi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *v_temp=gsl_vector_alloc((p->eval)->size);
	
	gsl_vector_memcpy (v_temp, p->eval);
	gsl_vector_scale (v_temp, l);
	if (p->e_mode==0) {gsl_vector_set_all (Hi_eval, 1.0);} else {gsl_vector_memcpy (Hi_eval, v_temp);}
	gsl_vector_add_constant (v_temp, 1.0);
	gsl_vector_div (Hi_eval, v_temp);
		
	gsl_vector_memcpy (HiHi_eval, Hi_eval);
	gsl_vector_mul (HiHi_eval, Hi_eval);	
	gsl_vector_memcpy (HiHiHi_eval, HiHi_eval);
	gsl_vector_mul (HiHiHi_eval, Hi_eval);
	
	gsl_vector_set_all (v_temp, 1.0);
	gsl_blas_ddot (Hi_eval, v_temp, &trace_Hi);
	gsl_blas_ddot (HiHi_eval, v_temp, &trace_HiHi);
	
	if (p->e_mode!=0) {	
		trace_Hi=(double)ni_test-trace_Hi;
		trace_HiHi=2*trace_Hi+trace_HiHi-(double)ni_test;
	}
	
	CalcPab (n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, Pab);	
	CalcPPab (n_cvt, p->e_mode, HiHi_eval, p->Uab, p->ab, Pab, PPab);	
	CalcPPPab (n_cvt, p->e_mode, HiHiHi_eval, p->Uab, p->ab, Pab, PPab, PPPab);	
	
	//calculate tracePK and trace PKPK
	double trace_P=trace_Hi, trace_PP=trace_HiHi;
	double ps_ww, ps2_ww, ps3_ww;
	for (size_t i=0; i<nc_total; ++i) {
		index_ww=GetabIndex (i+1, i+1, n_cvt);
		ps_ww=gsl_matrix_get (Pab, i, index_ww);
		ps2_ww=gsl_matrix_get (PPab, i, index_ww);
		ps3_ww=gsl_matrix_get (PPPab, i, index_ww);
		trace_P-=ps2_ww/ps_ww;
		trace_PP+=ps2_ww*ps2_ww/(ps_ww*ps_ww)-2.0*ps3_ww/ps_ww;
	}
	double trace_PKPK=(df+trace_PP-2.0*trace_P)/(l*l);
	
	//calculate yPKPy, yPKPKPy
	index_ww=GetabIndex (n_cvt+2, n_cvt+2, n_cvt);
	double P_yy=gsl_matrix_get (Pab, nc_total, index_ww);
	double PP_yy=gsl_matrix_get (PPab, nc_total, index_ww);
	double PPP_yy=gsl_matrix_get (PPPab, nc_total, index_ww);				
	double yPKPy=(P_yy-PP_yy)/l;	
	double yPKPKPy=(P_yy+PPP_yy-2.0*PP_yy)/(l*l);
	
	dev2=0.5*trace_PKPK-0.5*df*(2.0*yPKPKPy*P_yy-yPKPy*yPKPy)/(P_yy*P_yy);
	
	gsl_matrix_free (Pab);
	gsl_matrix_free (PPab);
	gsl_matrix_free (PPPab);
	gsl_vector_free (Hi_eval);
	gsl_vector_free (HiHi_eval);
	gsl_vector_free (HiHiHi_eval);
	gsl_vector_free (v_temp);	
	
	return dev2;
}
	


/**
 * @brief Compute both REML derivatives simultaneously
 * @param l Current value of λ
 * @param params Pointer to FUNC_PARAM structure  
 * @param dev1 Output: first derivative dℓ_R/dλ
 * @param dev2 Output: second derivative d²ℓ_R/dλ²
 * 
 * Efficient computation of both derivatives for Newton-Raphson optimization
 * of the restricted likelihood (used in REMLE estimation).
 */
void LogRL_dev12 (double l, void *params, double *dev1, double *dev2)
{
	FUNC_PARAM *p=(FUNC_PARAM *) params;	
	size_t n_cvt=p->n_cvt;
	size_t ni_test=p->ni_test;	
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	double df;
	size_t nc_total;
	if (p->calc_null==true) {nc_total=n_cvt; df=(double)ni_test-(double)n_cvt; }
	else {nc_total=n_cvt+1; df=(double)ni_test-(double)n_cvt-1.0;}
	
	double trace_Hi=0.0, trace_HiHi=0.0;
	size_t index_ww;
	
	gsl_matrix *Pab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_matrix *PPab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_matrix *PPPab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_vector *Hi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *HiHi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *HiHiHi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *v_temp=gsl_vector_alloc((p->eval)->size);
	
	gsl_vector_memcpy (v_temp, p->eval);
	gsl_vector_scale (v_temp, l);
	if (p->e_mode==0) {gsl_vector_set_all (Hi_eval, 1.0);} else {gsl_vector_memcpy (Hi_eval, v_temp);}
	gsl_vector_add_constant (v_temp, 1.0);
	gsl_vector_div (Hi_eval, v_temp);
		
	gsl_vector_memcpy (HiHi_eval, Hi_eval);
	gsl_vector_mul (HiHi_eval, Hi_eval);	
	gsl_vector_memcpy (HiHiHi_eval, HiHi_eval);
	gsl_vector_mul (HiHiHi_eval, Hi_eval);
	
	gsl_vector_set_all (v_temp, 1.0);
	gsl_blas_ddot (Hi_eval, v_temp, &trace_Hi);
	gsl_blas_ddot (HiHi_eval, v_temp, &trace_HiHi);
	
	if (p->e_mode!=0) {	
		trace_Hi=(double)ni_test-trace_Hi;
		trace_HiHi=2*trace_Hi+trace_HiHi-(double)ni_test;
	}
	
	CalcPab (n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, Pab);	
	CalcPPab (n_cvt, p->e_mode, HiHi_eval, p->Uab, p->ab, Pab, PPab);	
	CalcPPPab (n_cvt, p->e_mode, HiHiHi_eval, p->Uab, p->ab, Pab, PPab, PPPab);	
		
	//calculate tracePK and trace PKPK
	double trace_P=trace_Hi, trace_PP=trace_HiHi;
	double ps_ww, ps2_ww, ps3_ww;
	for (size_t i=0; i<nc_total; ++i) {
		index_ww=GetabIndex (i+1, i+1, n_cvt);
		ps_ww=gsl_matrix_get (Pab, i, index_ww);
		ps2_ww=gsl_matrix_get (PPab, i, index_ww);
		ps3_ww=gsl_matrix_get (PPPab, i, index_ww);
		trace_P-=ps2_ww/ps_ww;
		trace_PP+=ps2_ww*ps2_ww/(ps_ww*ps_ww)-2.0*ps3_ww/ps_ww;
	}
	double trace_PK=(df-trace_P)/l;
	double trace_PKPK=(df+trace_PP-2.0*trace_P)/(l*l);
	
	//calculate yPKPy, yPKPKPy
	index_ww=GetabIndex (n_cvt+2, n_cvt+2, n_cvt);
	double P_yy=gsl_matrix_get (Pab, nc_total, index_ww);
	double PP_yy=gsl_matrix_get (PPab, nc_total, index_ww);
	double PPP_yy=gsl_matrix_get (PPPab, nc_total, index_ww);				
	double yPKPy=(P_yy-PP_yy)/l;	
	double yPKPKPy=(P_yy+PPP_yy-2.0*PP_yy)/(l*l);
	
	*dev1=-0.5*trace_PK+0.5*df*yPKPy/P_yy;
	*dev2=0.5*trace_PKPK-0.5*df*(2.0*yPKPKPy*P_yy-yPKPy*yPKPy)/(P_yy*P_yy);
	
	gsl_matrix_free (Pab);
	gsl_matrix_free (PPab);
	gsl_matrix_free (PPPab);
	gsl_vector_free (Hi_eval);
	gsl_vector_free (HiHi_eval);
	gsl_vector_free (HiHiHi_eval);
	gsl_vector_free (v_temp);	
	
	return ;
}
	







/**
 * @brief Calculate Wald test statistics using REMLE
 * @param l Variance ratio λ (typically REMLE estimate)
 * @param params FUNC_PARAM structure with transformed data
 * @param beta Output: effect size estimate β̂_{j,t} (Equation 4 from model)
 * @param se Output: standard error of β̂
 * @param p_wald Output: Wald test p-value
 * 
 * MATHEMATICAL IMPLEMENTATION:
 * Implements the estimator from Equation 4:
 *   β̂ = (X^T * V^{-1} * X)^{-1} * (X^T * V^{-1} * Y)
 * 
 * In transformed space (after eigendecomposition):
 *   β̂ = P_xy / P_xx
 * where P_xy = (U^T*X)^T * (λ*D + I)^{-1} * (U^T*Y)
 * 
 * Standard error:
 *   SE(β̂) = sqrt(σ²_e / (X^T * V^{-1} * X)) = sqrt(1 / (τ * P_xx))
 * where τ = df / P_yy is the residual variance estimate.
 * 
 * Wald statistic:
 *   W = (P_yy - Px_yy) * τ ~ F(1, df)
 * which tests H0: β = 0 vs H1: β ≠ 0
 * 
 * This is the primary test for detecting eQTLs in a_mode=1.
 */
void LMM::CalcRLWald (const double &l, const FUNC_PARAM &params, double &beta, double &se, double &p_wald)
{
	size_t n_cvt=params.n_cvt;
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	int df=(int)ni_test-(int)n_cvt-1;
			
	gsl_matrix *Pab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_vector *Hi_eval=gsl_vector_alloc(params.eval->size);
	gsl_vector *v_temp=gsl_vector_alloc(params.eval->size);
	
	gsl_vector_memcpy (v_temp, params.eval);
	gsl_vector_scale (v_temp, l);
	if (params.e_mode==0) {gsl_vector_set_all (Hi_eval, 1.0);} else {gsl_vector_memcpy (Hi_eval, v_temp);}
	gsl_vector_add_constant (v_temp, 1.0);
	gsl_vector_div (Hi_eval, v_temp);	
	
	CalcPab (n_cvt, params.e_mode, Hi_eval, params.Uab, params.ab, Pab);	
	
	size_t index_yy=GetabIndex (n_cvt+2, n_cvt+2, n_cvt);	
	size_t index_xx=GetabIndex (n_cvt+1, n_cvt+1, n_cvt);
	size_t index_xy=GetabIndex (n_cvt+2, n_cvt+1, n_cvt);
	double P_yy=gsl_matrix_get (Pab, n_cvt, index_yy);
	double P_xx=gsl_matrix_get (Pab, n_cvt, index_xx);
	double P_xy=gsl_matrix_get (Pab, n_cvt, index_xy);	
	double Px_yy=gsl_matrix_get (Pab, n_cvt+1, index_yy);	
	
	beta=P_xy/P_xx;
	double tau=(double)df/Px_yy;
	se=sqrt(1.0/(tau*P_xx));	
	p_wald=gsl_cdf_fdist_Q ((P_yy-Px_yy)*tau, 1.0, df);	
//	p_wald=gsl_cdf_chisq_Q ((P_yy-Px_yy)*tau, 1);	
	
	gsl_matrix_free (Pab);
	gsl_vector_free (Hi_eval);
	gsl_vector_free (v_temp);
	return ;
}

/**
 * @brief Calculate Score test statistics
 * @param l Variance ratio λ (typically MLE from null model)
 * @param params FUNC_PARAM structure with transformed data
 * @param beta Output: effect size estimate
 * @param se Output: standard error
 * @param p_score Output: Score test p-value
 * 
 * MATHEMATICAL BACKGROUND:
 * The Score test evaluates the derivative of the likelihood at the null
 * hypothesis (β=0), without requiring MLE under the alternative.
 * 
 * Advantages over Wald/LRT:
 * 1. Only requires fitting null model (more efficient)
 * 2. More powerful when alternative is far from null
 * 3. Invariant to parameterization
 * 
 * Score statistic:
 *   S = (U^T*X)^T * (λ*D + I)^{-1} * (U^T*Y)  [score/gradient]
 *   I = (U^T*X)^T * (λ*D + I)^{-1} * (U^T*X)  [Fisher information]
 *   
 * Test statistic:
 *   T = S^2 / I = n * P_xy² / (P_yy * P_xx) ~ F(1, df)
 * 
 * The Score test (a_mode=3) is efficient for genome-wide scans where
 * the null model is fit once and reused for all SNPs.
 */
void LMM::CalcRLScore (const double &l, const FUNC_PARAM &params, double &beta, double &se, double &p_score)
{
	size_t n_cvt=params.n_cvt;
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	int df=(int)ni_test-(int)n_cvt-1;
			
	gsl_matrix *Pab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_vector *Hi_eval=gsl_vector_alloc(params.eval->size);
	gsl_vector *v_temp=gsl_vector_alloc(params.eval->size);
	
	gsl_vector_memcpy (v_temp, params.eval);
	gsl_vector_scale (v_temp, l);
	if (params.e_mode==0) {gsl_vector_set_all (Hi_eval, 1.0);} else {gsl_vector_memcpy (Hi_eval, v_temp);}
	gsl_vector_add_constant (v_temp, 1.0);
	gsl_vector_div (Hi_eval, v_temp);	
	
	CalcPab (n_cvt, params.e_mode, Hi_eval, params.Uab, params.ab, Pab);	
	
	size_t index_yy=GetabIndex (n_cvt+2, n_cvt+2, n_cvt);	
	size_t index_xx=GetabIndex (n_cvt+1, n_cvt+1, n_cvt);
	size_t index_xy=GetabIndex (n_cvt+2, n_cvt+1, n_cvt);
	double P_yy=gsl_matrix_get (Pab, n_cvt, index_yy);
	double P_xx=gsl_matrix_get (Pab, n_cvt, index_xx);
	double P_xy=gsl_matrix_get (Pab, n_cvt, index_xy);	
	double Px_yy=gsl_matrix_get (Pab, n_cvt+1, index_yy);	
	
	beta=P_xy/P_xx;
	double tau=(double)df/Px_yy;
	se=sqrt(1.0/(tau*P_xx));	
	
	p_score=gsl_cdf_fdist_Q ((double)ni_test*P_xy*P_xy/(P_yy*P_xx), 1.0, df);
//	p_score=gsl_cdf_chisq_Q ((double)ni_test*P_xy*P_xy/(P_yy*P_xx), 1);	
	
	gsl_matrix_free (Pab);
	gsl_vector_free (Hi_eval);
	gsl_vector_free (v_temp);
	return ;
}







/**
 * @brief Calculate transformed cross-products U^T * [W,X,Y] * [W,X,Y]^T * U
 * @param UtW Matrix U^T*W where W contains covariates (n_test × n_cvt)
 * @param Uty Vector U^T*y where y is the phenotype (n_test × 1)
 * @param Uab Output: matrix storing all pairwise products in compact form
 * 
 * MATHEMATICAL PURPOSE:
 * After eigendecomposition of K = U*D*U^T, we work in transformed space.
 * This function computes all cross-products needed for the LMM:
 * 
 * Products stored (using GetabIndex mapping):
 *   U^T*W_i * (U^T*W_j)^T  for all covariate pairs i,j
 *   U^T*W_i * (U^T*y)^T    for covariate-phenotype products
 *   U^T*y * (U^T*y)^T      for phenotype variance
 * 
 * These are the building blocks for:
 *   - P matrices used in likelihood calculations
 *   - Effect size estimates: β̂ = (X^T*V^{-1}*X)^{-1} * X^T*V^{-1}*Y
 *   - Test statistics
 * 
 * The compact storage exploits symmetry: only upper triangle is stored.
 */
void CalcUab (const gsl_matrix *UtW, const gsl_vector *Uty, gsl_matrix *Uab) 
{
	size_t index_ab;
	size_t n_cvt=UtW->size2;
	
	gsl_vector *u_a=gsl_vector_alloc (Uty->size);
	
	for (size_t a=1; a<=n_cvt+2; ++a) {
		if (a==n_cvt+1) {continue;}
		
		if (a==n_cvt+2) {gsl_vector_memcpy (u_a, Uty);}
		else {
			gsl_vector_const_view UtW_col=gsl_matrix_const_column (UtW, a-1);
			gsl_vector_memcpy (u_a, &UtW_col.vector);
		}
		
		for (size_t b=a; b>=1; --b) {		
			if (b==n_cvt+1) {continue;}
			
			index_ab=GetabIndex (a, b, n_cvt);
			gsl_vector_view Uab_col=gsl_matrix_column (Uab, index_ab);
			
			if (b==n_cvt+2) {gsl_vector_memcpy (&Uab_col.vector, Uty);}
			else {
				gsl_vector_const_view UtW_col=gsl_matrix_const_column (UtW, b-1);
				gsl_vector_memcpy (&Uab_col.vector, &UtW_col.vector);
			}			
			
			gsl_vector_mul(&Uab_col.vector, u_a);
		}
	}
	
	gsl_vector_free (u_a);
	return;
}

/**
 * @brief Extend Uab to include SNP genotype cross-products
 * @param UtW Matrix U^T*W (covariates in transformed space)
 * @param Uty Vector U^T*y (phenotype in transformed space)
 * @param Utx Vector U^T*x (SNP genotype in transformed space)
 * @param Uab Output: extended cross-product matrix including SNP
 * 
 * MATHEMATICAL PURPOSE:
 * This function adds the test SNP (x) to the cross-product calculations.
 * Computes additional products:
 *   U^T*x * (U^T*W_i)^T  for SNP-covariate products
 *   U^T*x * (U^T*x)^T    for SNP variance
 *   U^T*x * (U^T*y)^T    for SNP-phenotype covariance
 * 
 * These are needed to test the alternative hypothesis H1: β ≠ 0
 * by incorporating the SNP into the model alongside covariates.
 * 
 * Called for each SNP in genome-wide scan.
 */
void CalcUab (const gsl_matrix *UtW, const gsl_vector *Uty, const gsl_vector *Utx, gsl_matrix *Uab) 
{	
	size_t index_ab;
	size_t n_cvt=UtW->size2;
	
	for (size_t b=1; b<=n_cvt+2; ++b) {			
		index_ab=GetabIndex (n_cvt+1, b, n_cvt);
		gsl_vector_view Uab_col=gsl_matrix_column (Uab, index_ab);
		
		if (b==n_cvt+2) {gsl_vector_memcpy (&Uab_col.vector, Uty);}
		else if (b==n_cvt+1) {gsl_vector_memcpy (&Uab_col.vector, Utx);}
		else {
			gsl_vector_const_view UtW_col=gsl_matrix_const_column (UtW, b-1);
			gsl_vector_memcpy (&Uab_col.vector, &UtW_col.vector);
		}
		
		gsl_vector_mul(&Uab_col.vector, Utx);
	}
	
	return;
}



void Calcab (const gsl_matrix *W, const gsl_vector *y, gsl_vector *ab) 
{
	size_t index_ab;
	size_t n_cvt=W->size2;
	
	double d;
	gsl_vector *v_a=gsl_vector_alloc (y->size);
	gsl_vector *v_b=gsl_vector_alloc (y->size);
	
	for (size_t a=1; a<=n_cvt+2; ++a) {
		if (a==n_cvt+1) {continue;}
		
		if (a==n_cvt+2) {gsl_vector_memcpy (v_a, y);}
		else {
			gsl_vector_const_view W_col=gsl_matrix_const_column (W, a-1);
			gsl_vector_memcpy (v_a, &W_col.vector);
		}
		
		for (size_t b=a; b>=1; --b) {		
			if (b==n_cvt+1) {continue;}
			
			index_ab=GetabIndex (a, b, n_cvt);
			
			if (b==n_cvt+2) {gsl_vector_memcpy (v_b, y);}
			else {
				gsl_vector_const_view W_col=gsl_matrix_const_column (W, b-1);
				gsl_vector_memcpy (v_b, &W_col.vector);
			}			
			
			gsl_blas_ddot (v_a, v_b, &d);
			gsl_vector_set(ab, index_ab, d);
		}
	}
	
	gsl_vector_free (v_a);
	gsl_vector_free (v_b);
	return;
}


void Calcab (const gsl_matrix *W, const gsl_vector *y, const gsl_vector *x, gsl_vector *ab) 
{	
	size_t index_ab;
	size_t n_cvt=W->size2;
	
	double d;
	gsl_vector *v_b=gsl_vector_alloc (y->size);
	
	for (size_t b=1; b<=n_cvt+2; ++b) {			
		index_ab=GetabIndex (n_cvt+1, b, n_cvt);
		
		if (b==n_cvt+2) {gsl_vector_memcpy (v_b, y);}
		else if (b==n_cvt+1) {gsl_vector_memcpy (v_b, x);}
		else {
			gsl_vector_const_view W_col=gsl_matrix_const_column (W, b-1);
			gsl_vector_memcpy (v_b, &W_col.vector);
		}
		
		gsl_blas_ddot (x, v_b, &d);
		gsl_vector_set(ab, index_ab, d);
	}
	
	gsl_vector_free (v_b);
	
	return;
}




/**
 * @brief Perform eQTL analysis for gene expression data
 * @param U Eigenvectors of kinship matrix K (n_test × n_test)
 * @param eval Eigenvalues of K (n_test × 1)
 * @param UtW Transformed covariates U^T*W
 * @param Utx Transformed genotype U^T*x (for gene expression as "genotype")
 * @param W Original covariate matrix
 * @param x Original gene expression vector (treated as predictor)
 * 
 * ANALYSIS WORKFLOW:
 * For each gene expression phenotype:
 * 
 * 1. TRANSFORM TO EIGENSPACE:
 *    Compute U^T*y for current gene
 * 
 * 2. FIT NULL MODEL (no genotype effect):
 *    Estimate λ_0 maximizing ℓ(λ) or ℓ_R(λ) under H0: β = 0
 *    Only done for LRT (a_mode=2) and Score test (a_mode=3)
 * 
 * 3. FIT ALTERNATIVE MODEL (with genotype):
 *    Estimate λ_1 and test statistics under H1: β ≠ 0
 * 
 * 4. COMPUTE TEST STATISTICS:
 *    - Wald test (a_mode=1): uses REMLE, tests β̂/SE(β̂)
 *    - LRT (a_mode=2): tests 2*(ℓ_1 - ℓ_0) ~ χ²(1)
 *    - Score test (a_mode=3): efficient, uses null model only
 *    - All tests (a_mode=4)
 * 
 * Results stored in sumStat vector for later output.
 * 
 * This function is used when gene expression is the PREDICTOR (not response),
 * which is less common but useful for reverse eQTL analysis.
 */
void LMM::AnalyzeGene (const gsl_matrix *U, const gsl_vector *eval, const gsl_matrix *UtW, const gsl_vector *Utx, const gsl_matrix *W, const gsl_vector *x) 
{
	ifstream infile (file_gene.c_str(), ifstream::in);
	if (!infile) {cout<<"error reading gene expression file:"<<file_gene<<endl; return;}
	
	clock_t time_start=clock();
	
	string line;
	char *ch_ptr;
	
	double lambda_mle=0, lambda_remle=0, beta=0, se=0, p_wald=0, p_lrt=0, p_score=0;
	double logl_H1=0.0, logl_H0=0.0, l_H0;
	int c_phen;
	string rs; //gene id
	double d;
	
	//Calculate basic quantities
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	gsl_vector *y=gsl_vector_alloc (U->size1);
	gsl_vector *Uty=gsl_vector_alloc (U->size2);
	gsl_matrix *Uab=gsl_matrix_alloc (U->size2, n_index);
	gsl_vector *ab=gsl_vector_alloc (n_index);	
		
	//header
	getline(infile, line);
	
	for (size_t t=0; t<ng_total; t++) {
		!safeGetline(infile, line).eof();
		if (t%d_pace==0 || t==ng_total-1) {ProgressBar ("Performing Analysis ", t, ng_total-1);}
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		rs=ch_ptr;
		
		c_phen=0; 
		for (size_t i=0; i<indicator_idv.size(); ++i) {
			ch_ptr=strtok (NULL, " , \t");
			if (indicator_idv[i]==0) {continue;}
			
			d=atof(ch_ptr); 			
			gsl_vector_set(y, c_phen, d);
			
			c_phen++;
		}
		
		time_start=clock();
		gsl_blas_dgemv (CblasTrans, 1.0, U, y, 0.0, Uty);		
		time_UtX+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
	
		//calculate null
		time_start=clock();
		
		gsl_matrix_set_zero (Uab);
		
		CalcUab (UtW, Uty, Uab);
		FUNC_PARAM param0={false, ni_test, n_cvt, eval, Uab, ab, 0};
		
		if (a_mode==2 || a_mode==3 || a_mode==4) {
			CalcLambda('L', param0, l_min, l_max, n_region, l_H0, logl_H0);
		}
		
		//calculate alternative
		CalcUab(UtW, Uty, Utx, Uab);
		FUNC_PARAM param1={false, ni_test, n_cvt, eval, Uab, ab, 0};
		
		//3 is before 1
		if (a_mode==3 || a_mode==4) {
			CalcRLScore (l_H0, param1, beta, se, p_score);
		}
		
		if (a_mode==1 || a_mode==4) {
			CalcLambda ('R', param1, l_min, l_max, n_region, lambda_remle, logl_H1);
			CalcRLWald (lambda_remle, param1, beta, se, p_wald);
		}
		
		if (a_mode==2 || a_mode==4) {
			CalcLambda ('L', param1, l_min, l_max, n_region, lambda_mle, logl_H1);
			p_lrt=gsl_cdf_chisq_Q (2.0*(logl_H1-logl_H0), 1);	
		}
		
		time_opt+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
		
		//store summary data
		SUMSTAT SNPs={beta, se, lambda_remle, lambda_mle, p_wald, p_lrt, p_score};
		sumStat.push_back(SNPs);
    }
	cout<<endl;
	
	gsl_vector_free (y);
	gsl_vector_free (Uty);
	gsl_matrix_free (Uab);
	gsl_vector_free (ab);
	
	infile.close();
	infile.clear();
	
	return;
}




/**
 * @brief Perform genome-wide association scan using BIMBAM format genotypes
 * @param U Eigenvectors of kinship matrix K
 * @param eval Eigenvalues of K
 * @param UtW Transformed covariates U^T*W
 * @param Uty Transformed phenotype U^T*y (fixed for all SNPs)
 * @param W Original covariate matrix
 * @param y Original phenotype vector
 * 
 * GENOME-WIDE eQTL ANALYSIS (Standard Case):
 * Tests association between each SNP and a gene expression phenotype.
 * 
 * WORKFLOW FOR EACH SNP:
 * 
 * 1. READ GENOTYPE:
 *    - BIMBAM format: chr rs pos geno1 geno2 ... genoN
 *    - Genotypes are dosages (0, 1, 2) or probabilities
 *    - Missing values (NA) imputed with mean genotype
 * 
 * 2. QUALITY CONTROL:
 *    - Mean-center genotypes
 *    - Flip alleles if mean > 1 (ensure minor allele coding)
 *    - Skip SNPs based on indicator_snp filter
 * 
 * 3. TRANSFORM GENOTYPE:
 *    Compute U^T*x for current SNP (time tracked in time_UtX)
 * 
 * 4. CALCULATE STATISTICS:
 *    Using pre-computed null model (l_mle_null, logl_mle_H0):
 *    
 *    a) Score test (a_mode=3,4): 
 *       - Most efficient: reuses null model λ
 *       - Computed first to get β estimate
 *    
 *    b) Wald test (a_mode=1,4):
 *       - Estimates λ via REML for this SNP
 *       - Computes β̂ and SE(β̂) (Equation 4 from model)
 *       - Tests H0: β = 0 using t-statistic
 *    
 *    c) LRT (a_mode=2,4):
 *       - Estimates λ via ML for this SNP
 *       - Compares likelihoods: 2*(ℓ_1 - ℓ_0) ~ χ²(1)
 * 
 * 5. STORE RESULTS:
 *    Summary statistics saved to sumStat vector
 * 
 * MULTI-TISSUE CONTEXT:
 * This implements the single-tissue analysis from Equation 1-4.
 * For multi-tissue analysis, run separately for each tissue t,
 * then combine using covariance structure (Equations 5-6).
 */
void LMM::AnalyzeBimbam (const gsl_matrix *U, const gsl_vector *eval, const gsl_matrix *UtW, const gsl_vector *Uty, const gsl_matrix *W, const gsl_vector *y) 
{
	igzstream infile (file_geno.c_str(), igzstream::in);
//	ifstream infile (file_geno.c_str(), ifstream::in);
	if (!infile) {cout<<"error reading genotype file:"<<file_geno<<endl; return;}

	clock_t time_start=clock();
	
	string line;
	char *ch_ptr;
	
	double lambda_mle=0, lambda_remle=0, beta=0, se=0, p_wald=0, p_lrt=0, p_score=0;
	double logl_H1=0.0;
	int n_miss, c_phen;
	double geno, x_mean;
	
	//Calculate basic quantities
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;

	gsl_vector *x=gsl_vector_alloc (U->size1);
	gsl_vector *x_miss=gsl_vector_alloc (U->size1);
	gsl_vector *Utx=gsl_vector_alloc (U->size2);
	gsl_matrix *Uab=gsl_matrix_alloc (U->size2, n_index);
	gsl_vector *ab=gsl_vector_alloc (n_index);	
	
	gsl_matrix_set_zero (Uab);
	CalcUab (UtW, Uty, Uab);
//	if (e_mode!=0) {
//		gsl_vector_set_zero (ab);
//		Calcab (W, y, ab);
//	}	
	
	//start reading genotypes and analyze	
	for (size_t t=0; t<indicator_snp.size(); ++t) {
//		if (t>1) {break;}
		!safeGetline(infile, line).eof();
		if (t%d_pace==0 || t==(ns_total-1)) {ProgressBar ("Reading SNPs  ", t, ns_total-1);}
		if (indicator_snp[t]==0) {continue;}
		
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		ch_ptr=strtok (NULL, " , \t");
		ch_ptr=strtok (NULL, " , \t");		
		
		x_mean=0.0; c_phen=0; n_miss=0;
		gsl_vector_set_zero(x_miss);
		for (size_t i=0; i<ni_total; ++i) {
			ch_ptr=strtok (NULL, " , \t");
			if (indicator_idv[i]==0) {continue;}
			
			if (strcmp(ch_ptr, "NA")==0) {gsl_vector_set(x_miss, c_phen, 0.0); n_miss++;}
			else {
				geno=atof(ch_ptr); 				
				
				gsl_vector_set(x, c_phen, geno); 
				gsl_vector_set(x_miss, c_phen, 1.0); 
				x_mean+=geno;
			}
			c_phen++;
		}	
		
		x_mean/=(double)(ni_test-n_miss);
		
		for (size_t i=0; i<ni_test; ++i) {
			if (gsl_vector_get (x_miss, i)==0) {gsl_vector_set(x, i, x_mean);}
			geno=gsl_vector_get(x, i);
			if (x_mean>1) {
				gsl_vector_set(x, i, 2-geno);
			}
		}
		
		
		//calculate statistics
		time_start=clock();
		gsl_blas_dgemv (CblasTrans, 1.0, U, x, 0.0, Utx);		
		time_UtX+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
		
		CalcUab(UtW, Uty, Utx, Uab);
//		if (e_mode!=0) {
//			Calcab (W, y, x, ab);
//		}
		
		time_start=clock();
		FUNC_PARAM param1={false, ni_test, n_cvt, eval, Uab, ab, 0};
		
		//3 is before 1
		if (a_mode==3 || a_mode==4) {
			CalcRLScore (l_mle_null, param1, beta, se, p_score);
		}
		
		if (a_mode==1 || a_mode==4) {
			CalcLambda ('R', param1, l_min, l_max, n_region, lambda_remle, logl_H1);	
			CalcRLWald (lambda_remle, param1, beta, se, p_wald);
		}
		
		if (a_mode==2 || a_mode==4) {
			CalcLambda ('L', param1, l_min, l_max, n_region, lambda_mle, logl_H1);
			p_lrt=gsl_cdf_chisq_Q (2.0*(logl_H1-logl_mle_H0), 1);	
		}			
		
		if (x_mean>1) {beta*=-1;}
		
		time_opt+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
		
		//store summary data
		SUMSTAT SNPs={beta, se, lambda_remle, lambda_mle, p_wald, p_lrt, p_score};
		sumStat.push_back(SNPs);
    }	
	cout<<endl;
	
	gsl_vector_free (x);
	gsl_vector_free (x_miss);
	gsl_vector_free (Utx);
	gsl_matrix_free (Uab);
	gsl_vector_free (ab);
	
	infile.close();
	infile.clear();
	
	return;
}






/**
 * @brief Perform genome-wide association scan using PLINK binary format
 * @param U Eigenvectors of kinship matrix K
 * @param eval Eigenvalues of K
 * @param UtW Transformed covariates U^T*W
 * @param Uty Transformed phenotype U^T*y
 * @param W Original covariate matrix
 * @param y Original phenotype vector
 * 
 * PLINK BED FORMAT:
 * Binary genotype file with 2 bits per genotype:
 *   00: homozygous major allele (0 copies of minor allele)
 *   01: missing genotype
 *   10: heterozygous (1 copy)
 *   11: homozygous minor allele (2 copies)
 * 
 * File structure:
 *   - 3 magic bytes header
 *   - Genotypes in SNP-major order (all individuals for SNP1, then SNP2, ...)
 *   - Efficiently stores large datasets
 * 
 * GENOTYPE PROCESSING:
 * 1. Read 2-bit encoded genotypes from .bed file
 * 2. Decode: 00→2, 10→1, 11→0, 01→missing
 * 3. Impute missing with mean genotype
 * 4. Flip if mean > 1 (ensure minor allele)
 * 5. Transform: U^T*x
 * 6. Run association tests (same as AnalyzeBimbam)
 * 
 * This is the most memory-efficient format for large GWAS/eQTL studies.
 * Functionally equivalent to AnalyzeBimbam but for binary genotype files.
 */
void LMM::AnalyzePlink (const gsl_matrix *U, const gsl_vector *eval, const gsl_matrix *UtW, const gsl_vector *Uty, const gsl_matrix *W, const gsl_vector *y) 
{
	string file_bed=file_bfile+".bed";
	ifstream infile (file_bed.c_str(), ios::binary);
	if (!infile) {cout<<"error reading bed file:"<<file_bed<<endl; return;}
	
	clock_t time_start=clock();
	
	char ch[1];
	bitset<8> b;	
	
	double lambda_mle=0, lambda_remle=0, beta=0, se=0, p_wald=0, p_lrt=0, p_score=0;
	double logl_H1=0.0;
	int n_bit, n_miss, ci_total, ci_test;
	double geno, x_mean;
		
	//Calculate basic quantities
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;

	gsl_vector *x=gsl_vector_alloc (U->size1);
	gsl_vector *Utx=gsl_vector_alloc (U->size2);
	gsl_matrix *Uab=gsl_matrix_alloc (U->size2, n_index);	
	gsl_vector *ab=gsl_vector_alloc (n_index);	
	
	gsl_matrix_set_zero (Uab);
	CalcUab (UtW, Uty, Uab);
//	if (e_mode!=0) {
//		gsl_vector_set_zero (ab);
//		Calcab (W, y, ab);
//	}
		
	//calculate n_bit and c, the number of bit for each snp
	if (ni_total%4==0) {n_bit=ni_total/4;}
	else {n_bit=ni_total/4+1; }

	//print the first three majic numbers
	for (int i=0; i<3; ++i) {
		infile.read(ch,1);
		b=ch[0];
	}
	
	
	for (vector<SNPINFO>::size_type t=0; t<snpInfo.size(); ++t) {
		if (t%d_pace==0 || t==snpInfo.size()-1) {ProgressBar ("Reading SNPs  ", t, snpInfo.size()-1);}
		if (indicator_snp[t]==0) {continue;}
		
		infile.seekg(t*n_bit+3);		//n_bit, and 3 is the number of magic numbers
		
		//read genotypes
		x_mean=0.0;	n_miss=0; ci_total=0; ci_test=0; 
		for (int i=0; i<n_bit; ++i) {
			infile.read(ch,1);
			b=ch[0];
			for (size_t j=0; j<4; ++j) {                //minor allele homozygous: 2.0; major: 0.0;
				if ((i==(n_bit-1)) && ci_total==(int)ni_total) {break;}
				if (indicator_idv[ci_total]==0) {ci_total++; continue;}

				if (b[2*j]==0) {
					if (b[2*j+1]==0) {gsl_vector_set(x, ci_test, 2); x_mean+=2.0; }
					else {gsl_vector_set(x, ci_test, 1); x_mean+=1.0; }
				}
				else {
					if (b[2*j+1]==1) {gsl_vector_set(x, ci_test, 0); }                                  
					else {gsl_vector_set(x, ci_test, -9); n_miss++; }
				}

				ci_total++;
				ci_test++;
			}
		}
		
		x_mean/=(double)(ni_test-n_miss);
				
		for (size_t i=0; i<ni_test; ++i) {			
			geno=gsl_vector_get(x,i);
			if (geno==-9) {gsl_vector_set(x, i, x_mean); geno=x_mean;}
			if (x_mean>1) {
				gsl_vector_set(x, i, 2-geno);
			}
		}
		
		//calculate statistics
		time_start=clock();
		gsl_blas_dgemv (CblasTrans, 1.0, U, x, 0.0, Utx);
		time_UtX+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
		
		CalcUab(UtW, Uty, Utx, Uab);
//		if (e_mode!=0) {
//			Calcab (W, y, x, ab);
//		}
		
		time_start=clock();
		FUNC_PARAM param1={false, ni_test, n_cvt, eval, Uab, ab, 0};
		
		//3 is before 1, for beta
		if (a_mode==3 || a_mode==4) {
			CalcRLScore (l_mle_null, param1, beta, se, p_score);
		}
		
		if (a_mode==1 || a_mode==4) {
			CalcLambda ('R', param1, l_min, l_max, n_region, lambda_remle, logl_H1);	
			CalcRLWald (lambda_remle, param1, beta, se, p_wald);
		}
		
		if (a_mode==2 || a_mode==4) {
			CalcLambda ('L', param1, l_min, l_max, n_region, lambda_mle, logl_H1);
			p_lrt=gsl_cdf_chisq_Q (2.0*(logl_H1-logl_mle_H0), 1);	
		}		
		
		if (x_mean>1) {beta*=-1;}		
		
		time_opt+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
		
		//store summary data
		SUMSTAT SNPs={beta, se, lambda_remle, lambda_mle, p_wald, p_lrt, p_score};
		sumStat.push_back(SNPs);
    }	
	cout<<endl;
	
	gsl_vector_free (x);
	gsl_vector_free (Utx);
	gsl_matrix_free (Uab);
	gsl_vector_free (ab);
	
	infile.close();
	infile.clear();	
	
	return;
}




/**
 * @brief Calculate log-likelihood ratios for multiple SNPs (batch processing)
 * @param U Eigenvectors of kinship matrix K
 * @param UtX Matrix of transformed genotypes U^T*X for multiple SNPs
 * @param Uty Transformed phenotype U^T*y
 * @param K_eval Eigenvalues of kinship matrix
 * @param l_min Minimum λ for search
 * @param l_max Maximum λ for search
 * @param n_region Number of search intervals
 * @param pos_loglr Output: vector of (position, log-LR) pairs
 * 
 * BATCH LRT ANALYSIS:
 * Performs likelihood ratio tests for multiple SNPs efficiently.
 * 
 * For each SNP column in UtX:
 * 1. Fit null model (no SNP): estimate λ_0, compute ℓ_0
 * 2. Fit alternative model (with SNP): estimate λ_1, compute ℓ_1
 * 3. Calculate LRT statistic: LR = 2*(ℓ_1 - ℓ_0)
 * 
 * The log-LR can be converted to p-value: p = P(χ²_1 > LR)
 * 
 * This is useful for:
 * - Ranking SNPs by strength of association
 * - Fine-mapping: identifying causal variants in a region
 * - Conditional analysis: finding independent signals
 * 
 * More efficient than calling CalcLambda individually for each SNP
 * when multiple SNPs need to be evaluated.
 */
void MatrixCalcLR (const gsl_matrix *U, const gsl_matrix *UtX, const gsl_vector *Uty, const gsl_vector *K_eval, const double l_min, const double l_max, const size_t n_region, vector<pair<size_t, double> > &pos_loglr) 
{
	double logl_H0, logl_H1, log_lr, lambda0, lambda1;
	
	gsl_vector *w=gsl_vector_alloc (Uty->size);
	gsl_matrix *Utw=gsl_matrix_alloc (Uty->size, 1);	
	gsl_matrix *Uab=gsl_matrix_alloc (Uty->size, 6);
	gsl_vector *ab=gsl_vector_alloc (6);	
	
	gsl_vector_set_zero(ab);
	gsl_vector_set_all (w, 1.0);
	gsl_vector_view Utw_col=gsl_matrix_column (Utw, 0);	
	gsl_blas_dgemv (CblasTrans, 1.0, U, w, 0.0, &Utw_col.vector);		
	
	CalcUab (Utw, Uty, Uab) ;	
	FUNC_PARAM param0={true, Uty->size, 1, K_eval, Uab, ab, 0};	
	
	CalcLambda('L', param0, l_min, l_max, n_region, lambda0, logl_H0);
	
	for (size_t i=0; i<UtX->size2; ++i) {
		gsl_vector_const_view UtX_col=gsl_matrix_const_column (UtX, i);
		CalcUab(Utw, Uty, &UtX_col.vector, Uab);
		FUNC_PARAM param1={false, UtX->size1, 1, K_eval, Uab, ab, 0};
		
		CalcLambda ('L', param1, l_min, l_max, n_region, lambda1, logl_H1);
		log_lr=logl_H1-logl_H0;				
		
		pos_loglr.push_back(make_pair(i,log_lr) );
	}
	
	gsl_vector_free (w);
	gsl_matrix_free (Utw);
	gsl_matrix_free (Uab);
	gsl_vector_free (ab);
	
	return;
}



/**
 * @brief Estimate variance ratio parameter λ = σ²_g/σ²_e
 * @param func_name 'L' for MLE (max likelihood), 'R' for REMLE (restricted max likelihood)
 * @param params FUNC_PARAM with data and covariance structure
 * @param l_min Minimum search value for λ
 * @param l_max Maximum search value for λ
 * @param n_region Number of intervals for initial grid search
 * @param lambda Output: estimated λ value
 * @param logf Output: log-likelihood (or log-restricted-likelihood) at λ
 * 
 * PARAMETER INTERPRETATION:
 * λ = σ²_g/σ²_e is the ratio of genetic to environmental variance.
 * 
 * From the model (Equations 1-3):
 *   V = K*σ²_g + σ²_e*I = σ²_e*(λ*K + I)
 * 
 * - λ = 0: no genetic effects (pure environmental variance)
 * - λ → ∞: strong genetic effects dominate
 * - Typical range: [10^-5, 10^5] for robust search
 * 
 * ESTIMATION STRATEGY:
 * 
 * 1. GRID SEARCH PHASE:
 *    Divide [l_min, l_max] into n_region intervals (log-scale)
 *    Find intervals where derivative changes sign (potential maxima)
 * 
 * 2. ROOT FINDING (Brent's method):
 *    For each sign-change interval:
 *    - Use gsl_root_fsolver_brent to find where dℓ/dλ = 0
 *    - Robust bracketing method, no derivatives needed
 * 
 * 3. REFINEMENT (Newton-Raphson):
 *    Polish the root using derivatives:
 *    - λ_new = λ_old - (dℓ/dλ) / (d²ℓ/dλ²)
 *    - Faster convergence, needs gradient and Hessian
 * 
 * 4. BOUNDARY CHECK:
 *    Evaluate likelihood at l_min and l_max
 *    Return boundary value if better than interior maximum
 * 
 * MLE vs REMLE:
 * - MLE: maximizes full likelihood (used in LRT, a_mode=2)
 * - REMLE: accounts for fixed effects uncertainty (used in Wald, a_mode=1)
 * - REMLE generally preferred for variance estimation (less biased)
 * 
 * This is the core optimization for the LMM, analogous to estimating
 * heritability h² = σ²_g/(σ²_g + σ²_e) = λ/(λ+1).
 */
void CalcLambda (const char func_name, FUNC_PARAM &params, const double l_min, const double l_max, const size_t n_region, double &lambda, double &logf)
{
	if (func_name!='R' && func_name!='L' && func_name!='r' && func_name!='l') {cout<<"func_name only takes 'R' or 'L': 'R' for log-restricted likelihood, 'L' for log-likelihood."<<endl; return;}
	
	vector<pair<double, double> > lambda_lh;
	
	//evaluate first order derivates in different intervals
	double lambda_l, lambda_h, lambda_interval=log(l_max/l_min)/(double)n_region;
	double dev1_l, dev1_h, logf_l, logf_h;
	
	for (size_t i=0; i<n_region; ++i) {
		lambda_l=l_min*exp(lambda_interval*i);
		lambda_h=l_min*exp(lambda_interval*(i+1.0));
		
		if (func_name=='R' || func_name=='r') {
			dev1_l=LogRL_dev1 (lambda_l, &params);
			dev1_h=LogRL_dev1 (lambda_h, &params);
		}
		else {
			dev1_l=LogL_dev1 (lambda_l, &params);
			dev1_h=LogL_dev1 (lambda_h, &params);
		}
		
		if (dev1_l*dev1_h<=0) {
			lambda_lh.push_back(make_pair(lambda_l, lambda_h));
		}
	}
	
	//if derivates do not change signs in any interval
	if (lambda_lh.empty()) {
		if (func_name=='R' || func_name=='r') {
			logf_l=LogRL_f (l_min, &params);
			logf_h=LogRL_f (l_max, &params);
		}
		else {
			logf_l=LogL_f (l_min, &params);
			logf_h=LogL_f (l_max, &params);
		}
		
		if (logf_l>=logf_h) {lambda=l_min; logf=logf_l;} else {lambda=l_max; logf=logf_h;}
	}
	else {
		//if derivates change signs
		int status;
		int iter=0, max_iter=100;
		double l, l_temp;	
		
		gsl_function F;
		gsl_function_fdf FDF;
		
		F.params=&params;
		FDF.params=&params;
		
		if (func_name=='R' || func_name=='r') {
			F.function=&LogRL_dev1;
			FDF.f=&LogRL_dev1;
			FDF.df=&LogRL_dev2;
			FDF.fdf=&LogRL_dev12;
		}
		else {
			F.function=&LogL_dev1;
			FDF.f=&LogL_dev1;
			FDF.df=&LogL_dev2;
			FDF.fdf=&LogL_dev12;
		}
		
		const gsl_root_fsolver_type *T_f;
		gsl_root_fsolver *s_f;
		T_f=gsl_root_fsolver_brent;
		s_f=gsl_root_fsolver_alloc (T_f);
		
		const gsl_root_fdfsolver_type *T_fdf;
		gsl_root_fdfsolver *s_fdf;
		T_fdf=gsl_root_fdfsolver_newton;
		s_fdf=gsl_root_fdfsolver_alloc(T_fdf);	
		
		for (vector<double>::size_type i=0; i<lambda_lh.size(); ++i) {
			lambda_l=lambda_lh[i].first; lambda_h=lambda_lh[i].second;
			
			gsl_root_fsolver_set (s_f, &F, lambda_l, lambda_h);
			
			do {
				iter++;
				status=gsl_root_fsolver_iterate (s_f);
				l=gsl_root_fsolver_root (s_f);
				lambda_l=gsl_root_fsolver_x_lower (s_f);
				lambda_h=gsl_root_fsolver_x_upper (s_f);
				status=gsl_root_test_interval (lambda_l, lambda_h, 0, 1e-1);		
			}
			while (status==GSL_CONTINUE && iter<max_iter); 				
			
			iter=0;
			
			gsl_root_fdfsolver_set (s_fdf, &FDF, l);	
			
			do {
				iter++;
				status=gsl_root_fdfsolver_iterate (s_fdf);
				l_temp=l;
				l=gsl_root_fdfsolver_root (s_fdf);
				status=gsl_root_test_delta (l, l_temp, 0, 1e-5);		
			}
			while (status==GSL_CONTINUE && iter<max_iter && l>l_min && l<l_max); 
			
			l=l_temp;
			if (l<l_min) {l=l_min;}
			if (l>l_max) {l=l_max;}
			if (func_name=='R' || func_name=='r') {logf_l=LogRL_f (l, &params);} else {logf_l=LogL_f (l, &params);}			
			
			if (i==0) {logf=logf_l; lambda=l;}
			else if (logf<logf_l) {logf=logf_l; lambda=l;}
			else {}
		}
		gsl_root_fsolver_free (s_f);	
		gsl_root_fdfsolver_free (s_fdf);		
		
		if (func_name=='R' || func_name=='r') {
			logf_l=LogRL_f (l_min, &params);
			logf_h=LogRL_f (l_max, &params);
		}
		else {
			logf_l=LogL_f (l_min, &params);
			logf_h=LogL_f (l_max, &params);
		}
		
		if (logf_l>logf) {lambda=l_min; logf=logf_l;} 
		if (logf_h>logf) {lambda=l_max; logf=logf_h;}
	}
	
	return;
}





//calculate lambda in the null model
void CalcLambda (const char func_name, const gsl_vector *eval, const gsl_matrix *UtW, const gsl_vector *Uty, const double l_min, const double l_max, const size_t n_region, double &lambda, double &logl_H0)
{
	if (func_name!='R' && func_name!='L' && func_name!='r' && func_name!='l') {cout<<"func_name only takes 'R' or 'L': 'R' for log-restricted likelihood, 'L' for log-likelihood."<<endl; return;}

	size_t n_cvt=UtW->size2, ni_test=UtW->size1;
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	gsl_matrix *Uab=gsl_matrix_alloc (ni_test, n_index);	
	gsl_vector *ab=gsl_vector_alloc (n_index);	
	
	gsl_matrix_set_zero (Uab);
	CalcUab (UtW, Uty, Uab);
//	if (e_mode!=0) {
//		gsl_vector_set_zero (ab);
//		Calcab (W, y, ab);
//	}
		
	FUNC_PARAM param0={true, ni_test, n_cvt, eval, Uab, ab, 0};
	
	CalcLambda(func_name, param0, l_min, l_max, n_region, lambda, logl_H0);
	
	gsl_matrix_free(Uab);	
	gsl_vector_free(ab);	
	
	return;
}




/**
 * @brief Estimate λ for null model (wrapper function)
 * @param func_name 'L' for MLE, 'R' for REMLE
 * @param eval Eigenvalues of kinship matrix K
 * @param UtW Transformed covariates U^T*W
 * @param Uty Transformed phenotype U^T*y
 * @param l_min Minimum λ to search
 * @param l_max Maximum λ to search
 * @param n_region Number of search intervals
 * @param lambda Output: estimated λ under null model
 * @param logl_H0 Output: log-likelihood at null model
 * 
 * NULL MODEL ESTIMATION:
 * Fits the model without any test SNP (only covariates):
 *   Y = W*α + ε, where ε ~ N(0, λ*K*σ²_e + σ²_e*I)
 * 
 * This null model is fit ONCE per phenotype, then reused for all SNPs
 * in Score test (a_mode=3) and LRT (a_mode=2).
 * 
/**
 * @brief Calculate proportion of variance explained (PVE) / heritability
 * @param eval Eigenvalues of kinship matrix K
 * @param UtW Transformed covariates
 * @param Uty Transformed phenotype
 * @param lambda Estimated variance ratio λ (typically REMLE)
 * @param trace_G Trace of kinship matrix K (sum of eigenvalues)
 * @param pve Output: proportion of variance explained by genetics
 * @param pve_se Output: standard error of PVE estimate
 * 
 * MATHEMATICAL DERIVATION:
 * PVE = σ²_g / (σ²_g + σ²_e) = (λ·tr(K)) / (λ·tr(K) + 1)
 * where λ = σ²_g/σ²_e is the variance ratio
 * 
 * Standard error computed via delta method from λ uncertainty.
 */
/*
// DUPLICATE FUNCTION - COMMENTED OUT
void CalcLambda (const char func_name, const gsl_vector *eval, const gsl_matrix *UtW, const gsl_vector *Uty, const double l_min, const double l_max, const size_t n_region, double &lambda, double &logl_H0)
{
	size_t n_cvt=UtW->size2, ni_test=UtW->size1;
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	gsl_matrix *Uab=gsl_matrix_alloc (ni_test, n_index);	
	gsl_vector *ab=gsl_vector_alloc (n_index);	
	
	gsl_matrix_set_zero (Uab);
	CalcUab (UtW, Uty, Uab);
	//	if (e_mode!=0) {
	//		gsl_vector_set_zero (ab);
	//		Calcab (W, y, ab);
	//	}
	
	FUNC_PARAM param0={true, ni_test, n_cvt, eval, Uab, ab, 0};
	
	double se=sqrt(-1.0/LogRL_dev2 (lambda, &param0));
	
	pve=trace_G*lambda/(trace_G*lambda+1.0);
	pve_se=trace_G/((trace_G*lambda+1.0)*(trace_G*lambda+1.0))*se;
	
	gsl_matrix_free (Uab);
	gsl_vector_free (ab);	
	return;
}
*/


/**
 * @param eval Eigenvalues of kinship matrix K
 * @param UtW Transformed covariates
 * @param Uty Transformed phenotype
 * @param lambda Estimated variance ratio λ (typically REMLE)
 * @param trace_G Trace of kinship matrix K (sum of eigenvalues)
 * @param pve Output: proportion of variance explained by genetics
 * @param pve_se Output: standard error of PVE estimate
 * 
 * MATHEMATICAL DERIVATION:
 * Total phenotypic variance: σ²_P = σ²_g + σ²_e
 * 
 * PVE (also called narrow-sense heritability h²):
 *   h² = σ²_g / σ²_P = σ²_g / (σ²_g + σ²_e)
 *      = λ*σ²_e / (λ*σ²_e + σ²_e)
 *      = λ / (λ + 1)
 * 
 * With kinship scaling (trace_G adjustment):
 *   PVE = trace(K)*λ / (trace(K)*λ + 1)
 * 
 * This accounts for the scale of the kinship matrix K.
 * 
 * Standard error via delta method:
 *   SE(PVE) = |∂PVE/∂λ| * SE(λ)
 *   where SE(λ) = sqrt(-1 / (d²ℓ/dλ²))  [from Fisher information]
 * 
 * INTERPRETATION:
 * - PVE = 0: genetics explains no variance (all environmental)
 * - PVE = 1: genetics explains all variance  
 * - PVE = 0.5: half of variance is genetic, half environmental
 * 
 * This is a key quantity in quantitative genetics and is related to
 * the concept of "SNP heritability" in GWAS studies.
 */
void CalcPve (const gsl_vector *eval, const gsl_matrix *UtW, const gsl_vector *Uty, const double lambda, const double trace_G, double &pve, double &pve_se)
{
	size_t n_cvt=UtW->size2, ni_test=UtW->size1;
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	gsl_matrix *Uab=gsl_matrix_alloc (ni_test, n_index);	
	gsl_vector *ab=gsl_vector_alloc (n_index);	
	
	gsl_matrix_set_zero (Uab);
	CalcUab (UtW, Uty, Uab);
	//	if (e_mode!=0) {
	//		gsl_vector_set_zero (ab);
	//		Calcab (W, y, ab);
	//	}
	
	FUNC_PARAM param0={true, ni_test, n_cvt, eval, Uab, ab, 0};
	
	double se=sqrt(-1.0/LogRL_dev2 (lambda, &param0));
	
	pve=trace_G*lambda/(trace_G*lambda+1.0);
	pve_se=trace_G/((trace_G*lambda+1.0)*(trace_G*lambda+1.0))*se;
	
	gsl_matrix_free (Uab);
	gsl_vector_free (ab);	
	return;
}

/**
 * @brief Estimate variance components and fixed effects
 * @param eval Eigenvalues of kinship matrix K
 * @param UtW Transformed covariates U^T*W
 * @param Uty Transformed phenotype U^T*y
 * @param lambda Estimated variance ratio λ (REMLE)
 * @param vg Output: genetic variance σ²_g
 * @param ve Output: environmental variance σ²_e
 * @param beta Output: vector of fixed effect estimates (covariate effects)
 * @param se_beta Output: standard errors for beta coefficients
 * 
 * MATHEMATICAL DETAILS:
 * 
 * 1. VARIANCE COMPONENTS:
 *    From REML estimation at λ:
 *    - Residual variance: σ²_e = Y^T*P*Y / df
 *    - Genetic variance: σ²_g = λ * σ²_e
 *    where P is the projection matrix adjusted for covariates and λ
 * 
 * 2. FIXED EFFECTS (β for covariates):
 *    β̂ = (W^T * V^{-1} * W)^{-1} * W^T * V^{-1} * Y
 *    
 *    In transformed space:
 *    - Compute H^{-1}*W where H = λ*D + I
 *    - Solve: (W^T*H^{-1}*W) * β̂ = W^T*H^{-1}*Y
 *    - Use LU decomposition for numerical stability
 * 
 * 3. STANDARD ERRORS:
 *    Var(β̂) = σ²_e * (W^T * V^{-1} * W)^{-1}
 *    SE(β̂_i) = sqrt(Var(β̂)_{ii})
 * 
 * USAGE:
 * This provides complete parameter estimates for the LMM:
 * - Variance components (vg, ve) for understanding genetic architecture
 * - Fixed effects (beta) for covariate adjustments
 * - Standard errors for inference
 * 
 * These are typically reported alongside eQTL results to characterize
 * the overall model fit and partition variance.
 */
void CalcLmmVgVeBeta (const gsl_vector *eval, const gsl_matrix *UtW, const gsl_vector *Uty, const double lambda, double &vg, double &ve, gsl_vector *beta, gsl_vector *se_beta)
{
	size_t n_cvt=UtW->size2, ni_test=UtW->size1;
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	gsl_matrix *Uab=gsl_matrix_alloc (ni_test, n_index);	
	gsl_vector *ab=gsl_vector_alloc (n_index);	
	gsl_matrix *Pab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_vector *Hi_eval=gsl_vector_alloc(eval->size);
	gsl_vector *v_temp=gsl_vector_alloc(eval->size);
	gsl_matrix *HiW=gsl_matrix_alloc(eval->size, UtW->size2);
	gsl_matrix *WHiW=gsl_matrix_alloc(UtW->size2, UtW->size2);
	gsl_vector *WHiy=gsl_vector_alloc(UtW->size2);
	gsl_matrix *Vbeta=gsl_matrix_alloc(UtW->size2, UtW->size2);
	
	gsl_matrix_set_zero (Uab);
	CalcUab (UtW, Uty, Uab);	
	
	gsl_vector_memcpy (v_temp, eval);
	gsl_vector_scale (v_temp, lambda);
	gsl_vector_set_all (Hi_eval, 1.0);
	gsl_vector_add_constant (v_temp, 1.0);
	gsl_vector_div (Hi_eval, v_temp);
	
	//calculate beta
	gsl_matrix_memcpy (HiW, UtW);
	for (size_t i=0; i<UtW->size2; i++) {
		gsl_vector_view HiW_col=gsl_matrix_column(HiW, i);
		gsl_vector_mul(&HiW_col.vector, Hi_eval);
	}
	gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, HiW, UtW, 0.0, WHiW);
	gsl_blas_dgemv (CblasTrans, 1.0, HiW, Uty, 0.0, WHiy);
	
	int sig;
	gsl_permutation * pmt=gsl_permutation_alloc (UtW->size2);
	LUDecomp (WHiW, pmt, &sig);
	LUSolve (WHiW, pmt, WHiy, beta);
	LUInvert (WHiW, pmt, Vbeta);
		
	//calculate vg and ve
	CalcPab (n_cvt, 0, Hi_eval, Uab, ab, Pab);	
	
	size_t index_yy=GetabIndex (n_cvt+2, n_cvt+2, n_cvt);	
	double P_yy=gsl_matrix_get (Pab, n_cvt, index_yy);	
	
	ve=P_yy/(double)(ni_test-n_cvt);
	vg=ve*lambda;
	
	//with ve, calculate se(beta)
	gsl_matrix_scale(Vbeta, ve);
	
	//obtain se_beta
	for (size_t i=0; i<Vbeta->size1; i++) {
		gsl_vector_set (se_beta, i, sqrt(gsl_matrix_get(Vbeta, i, i) ) );
	}
	
	gsl_matrix_free(Uab);
	gsl_matrix_free(Pab);
	gsl_vector_free(ab);
	gsl_vector_free(Hi_eval);
	gsl_vector_free(v_temp);
	gsl_matrix_free(HiW);
	gsl_matrix_free(WHiW);
	gsl_vector_free(WHiy);
	gsl_matrix_free(Vbeta);
	
	gsl_permutation_free(pmt);
	return;
}

/*
 * ============================================================================
 * END OF FILE: gemma_lmm.cpp
 * ============================================================================
 * 
 * SUMMARY OF IMPLEMENTATION:
 * 
 * This file implements the statistical framework for multi-tissue eQTL analysis
 * using linear mixed models (LMMs) as described in the mmQTL paper.
 * 
 * KEY MATHEMATICAL CONCEPTS IMPLEMENTED:
 * 
 * 1. LINEAR MIXED MODEL (Equations 1-3):
 *    Y_t = X_j*β_{j,t} + ε̂_t
 *    where ε̂_t ~ N(0, K*σ²_g + σ²_e*I)
 *    
 *    - Accounts for population structure via kinship K
 *    - Estimates variant effect sizes β controlling for relatedness
 *    - Partitions variance into genetic (σ²_g) and environmental (σ²_e)
 * 
 * 2. EIGENDECOMPOSITION TRICK:
 *    K = U*D*U^T transforms problem to diagonal covariance:
 *    - Original: V = K*σ²_g + σ²_e*I (requires O(n³) matrix inversion)
 *    - Transformed: Ṽ = (D*λ + I)*σ²_e (diagonal, O(n) inversion)
 *    where λ = σ²_g/σ²_e
 * 
 * 3. EFFECT SIZE ESTIMATION (Equation 4):
 *    β̂ = (X^T*V^{-1}*X)^{-1} * X^T*V^{-1}*Y
 *    Computed efficiently in eigenspace as P_xy/P_xx
 * 
 * 4. VARIANCE RATIO ESTIMATION:
 *    λ = σ²_g/σ²_e estimated by maximizing:
 *    - Log-likelihood (MLE) for LRT (a_mode=2)
 *    - Log-restricted-likelihood (REMLE) for Wald (a_mode=1)
 *    Using Newton-Raphson with derivatives up to order 3
 * 
 * 5. THREE HYPOTHESIS TESTS:
 *    a) Wald test: (β̂/SE)² ~ F(1,df) - requires REMLE
 *    b) LRT: 2*(ℓ₁ - ℓ₀) ~ χ²(1) - compares likelihoods
 *    c) Score test: S²/I ~ F(1,df) - most efficient, uses null model only
 * 
 * MULTI-TISSUE EXTENSIONS (from mmqtl_extracted.md):
 * 
 * While this file implements single-tissue analysis, it provides the foundation
 * for multi-tissue meta-analysis:
 * 
 * - Covariance across tissues (Equations 5-6): estimates from this code are
 *   combined using empirical covariance matrix
 * - Fixed/Random effects meta-analysis: effect sizes β̂_{j,t} from each tissue
 *   are combined accounting for correlation structure
 * - Conditional analysis: iteratively identifies independent eQTL signals
 * 
 * COMPUTATIONAL EFFICIENCY:
 * 
 * - Eigendecomposition done ONCE per phenotype
 * - Matrix transformations (U^T*X) tracked separately (time_UtX)
 * - Optimization iterations tracked (time_opt)
 * - Score test reuses null model for all SNPs (most efficient)
 * - Compact storage of symmetric matrices via GetabIndex
 * 
 * FILE FORMATS SUPPORTED:
 * 
 * - PLINK binary (.bed): most memory-efficient for large GWAS
 * - BIMBAM: flexible text format with dosages
 * - Gene expression: custom format for reverse eQTL analysis
 * 
 * RELATIONSHIP TO BROADER mmQTL FRAMEWORK:
 * 
 * This LMM implementation forms the core statistical engine. Other components:
 * - gemma_eigenlib.cpp: eigendecomposition of kinship matrix
 * - PostCal.cpp: posterior probability calculations
 * - caviar_PostCal.cpp: fine-mapping causal variants
 * - conditional_function.cpp: iterative conditional analysis
 * - TopKSNP.cpp: selecting top K independent signals
 * 
 * Together these implement the full mmQTL pipeline for identifying and
 * characterizing multi-tissue genetic effects on gene expression.
 */
