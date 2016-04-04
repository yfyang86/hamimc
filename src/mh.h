#pragma once


#include "comm.h"


namespace biostacs
{
	

	class Likelihood {
	private:
		double SumOfLogLik;
		double SumOfPriorLogLik;
		random_device rd;
		vector<double> PriorParam;
		int SampleSizeNum;  // n: # sample size
		int ParamNum;		// p: # of params
		mat Xval;

		double sumloglik(double pdf(vector<double> param, mat x), vector<double>param) {
			for (int i = 0; i < SampleSizeNum; i++) {
				SumOfLogLik += pdf(param, Xval.row(i));
			}
			return(SumOfLogLik);
		}

		double sumloglik(double pdf(vector<double> param, mat x)) {
			sumloglik(pdf, PriorParam);
		}

	public:

		Likelihood() {
			SumOfLogLik = 0.;
			SumOfPriorLogLik = 0.;
		}
		Likelihood(vector<double> parameters, mat &data) {
			SampleSizeNum = data.n_rows;
			ParamNum = data.n_cols;
			Xval = data;
			PriorParam = parameters;

			SumOfLogLik = 0.;
			SumOfPriorLogLik = 0.;
		}

	};

	


	class MetropolisMCMC :public Likelihood{
	private:
		double SumOfLogLik;
		double SumOfPriorLogLik;
		int iter;
		int maxiter;
		vector<double> StartingParam;
		random_device rd;
		vector<double> PriorParam;
		vector<double> Proposal;
		mat MonteCarloChain;

		int SampleSizeNum;  // n: # sample size
		int ParamNum;		// p: # of params, redesign need to be dropped.
		int PriorParamNum;		// p: # of params
		mat Xval;
		mh_rnd::default_rnd rnd;
		double sumloglik(double pdf(vector<double> param, mat x), vector<double>param) {
			SumOfLogLik = 0.;
			for (int i = 0; i < SampleSizeNum; i++) {
				SumOfLogLik += pdf(param, Xval.row(i));
			}
			return(SumOfLogLik);
		}

	public:
		double priorlog(double priorpdf(vector<double> param), vector<double>param) {
			SumOfPriorLogLik = priorpdf(param);
			return(SumOfPriorLogLik);
		}

		double posterior(double pdf(vector<double> param, mat x), double priorpdf(vector<double> param), vector<double>param) {
			return(
				sumloglik(pdf, param) +
				priorlog(priorpdf, param)
				);
		}

		MetropolisMCMC(vector<double> parameters, mat &data, 
			int MaxIterations, vector<double> startingvalue) :
				Likelihood(parameters,data),
				maxiter(MaxIterations),StartingParam(startingvalue),
				PriorParamNum(startingvalue.size()){
		};

		void normal_proposal(vector<double> param) {
			for (size_t i = 0; i < PriorParamNum; i++) {
				Proposal[i] = rnd.rproposal()*.1 + abs(param[i]);
			}
		}

		mat MH_MCMC(double pdf(vector<double> param, mat x), double priorpdf(vector<double> param)) {
			iter = 0;
			MonteCarloChain.reshape(maxiter+1, PriorParamNum); // 1+Maxiter * Num of Params
			Proposal.resize(PriorParamNum);
			
			for (int i = 0; i < PriorParamNum; i++) { 
				MonteCarloChain(0, i) = StartingParam[i]; 
			}
			vector<double> tmpVec;
			tmpVec.resize(PriorParamNum);
			double ProbDiff;
			for (int i = 1; i <= maxiter; i++) {
				if (i % 100 == 0) cout << "Pct: " << 100.*i / maxiter << "%\t Acc rate: " << 100.*iter/i <<'%'<<endl;
				for (int jj = 0; jj < PriorParamNum; jj++) tmpVec[jj] = MonteCarloChain(i - 1, jj);
				
				normal_proposal(tmpVec); // proposal

				ProbDiff = posterior(pdf, priorpdf, Proposal) - posterior(pdf, priorpdf, tmpVec);
		
				if (log(rnd.runif()) < ProbDiff) {
					for (int ii = 0; ii < PriorParamNum; ii++)  MonteCarloChain(i, ii) = Proposal[ii];
					iter++;
				}
				else {
					for (int jj = 0; jj < PriorParamNum; jj++) MonteCarloChain(i, jj) = MonteCarloChain(i - 1, jj);
				}
				/* check step by step
				cout<<'\n' << i << "-th run:";
				for (double ss : tmpVec)  cout<< ss << '\t';
				system("pause");
				*/
			}
				return(MonteCarloChain);
			}
		};

	}
