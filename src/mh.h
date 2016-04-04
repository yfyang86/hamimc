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
		vector<double> PriorParam;
		vector<double> Proposal;
		mat MonteCarloChain;

		int SampleSizeNum;  // n: # sample size
		int ParamNum;		// p: # of params, redesign need to be dropped.
		int PriorParamNum;		// p: # of params
		mat Xval;
		mh_rnd::rnorm rnd;
		double sumloglik(double pdf(vector<double> param, mat x), vector<double>param) {
			SumOfLogLik = 0.;
			for (int i = 0; i < SampleSizeNum; i++) {
				SumOfLogLik += pdf(param, Xval.row(i));
			}
			return(SumOfLogLik);
		}
        void ProgressBar(int progress_int) {
			double progressPct = double(progress_int) / maxiter;
			int count = 0;
			cout << "[";
			int barWidth =50;
			int pos = barWidth * progressPct;
			for (int i = 0; i < barWidth; ++i) {
				if (i < pos) std::cout << "=";
				else if (i == pos) std::cout << ">";
				else std::cout << " ";
			}
			cout << "] " << int(progressPct * 100.0)<< "%\r";
			cout.flush();
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
				StartingParam(startingvalue)
				{
			SampleSizeNum = data.n_rows;
			ParamNum = data.n_cols;
			Xval = data;
			iter = 0;
			PriorParamNum = startingvalue.size();
			maxiter=MaxIterations;
		};

		void normal_proposal(vector<double> param) { // ~ N(param, .2), i.e.: param+.2*rnorm(1)
			for (int i = 0; i < PriorParamNum; i++) {
				Proposal[i] = rnd.r()*.2 + param[i];
			}
		}

		mat MH_MCMC(double pdf(vector<double> param, mat x), double priorpdf(vector<double> param)) {
			
			MonteCarloChain.reshape(maxiter+1, PriorParamNum); // 1+Maxiter * Num of Params
			Proposal.resize(PriorParamNum);
			std::random_device rda;
			std::mt19937 gen(rda());
			std::uniform_real_distribution<> runif(0, 1);

			for (int i = 0; i < PriorParamNum; i++) { 
				MonteCarloChain(0, i) = StartingParam[i]; 
			}
			vector<double> tmpVec;
			tmpVec.resize(PriorParamNum);
			double ProbDiff;
            int progressSEP=maxiter/50;
			for (int i = 1; i <= maxiter; i++) {
                if (i % progressSEP == 0)ProgressBar(i);
#ifdef DEBUG
				if (i % 1000 == 0) cout << i << "-run\tPct: " << 100.*i / maxiter << "%\t Acc rate: " << 100.*iter/i <<'%'<<endl;
#endif

				for (int jj = 0; jj < PriorParamNum; jj++) tmpVec[jj] = MonteCarloChain(i - 1, jj);
				
				normal_proposal(tmpVec); // proposal

				ProbDiff = posterior(pdf, priorpdf, Proposal) - posterior(pdf, priorpdf, tmpVec);
		
				if (log(runif(gen)) < ProbDiff) {
					for (int ii = 0; ii < PriorParamNum; ii++)  MonteCarloChain(i, ii) = Proposal[ii];
					iter++;
				}
				else {
					for (int jj = 0; jj < PriorParamNum; jj++) MonteCarloChain(i, jj) = MonteCarloChain(i - 1, jj);
				}
			}
				return(MonteCarloChain);
			}
		};



	}
