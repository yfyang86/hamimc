#pragma once
#include <math.h>
#include <vector>
#include <armadillo>
#include <random>
#include <limits>

#include <iostream>
#include <iomanip>
using namespace std;
using namespace arma;

#include <Windows.h>

double DOUBLE_EPSILON_MH = 10 * std::numeric_limits<double>::epsilon();



typedef normal_distribution<> DefaultProposal;

namespace mh_rnd
{
	std::mt19937 make_seeded_engine() {
		std::random_device r;
		std::seed_seq seed{ r(), r(), r(), r(), r(), r(), r(), r() };
		return std::mt19937(seed);
	}

	struct default_rnd {
		std::mt19937 eng = make_seeded_engine();		//initialize: mersenne twister
		std::normal_distribution<> DefaultProposalDist{ 0, 1 }; //Std normal
		std::uniform_real_distribution<> unif{ 0, 1 };
		default_rnd() = default;
		template<typename SeedSeq> default_rnd(SeedSeq &&seed) : eng(seed) {
		}
		double rproposal() {
			return DefaultProposalDist(eng);
		}
		double runif() {
			return unif(eng);
		}
	};

	struct  rnorm{
		std::mt19937 eng = make_seeded_engine();		//initialize: mersenne twister
		std::normal_distribution<> normal_01{ 0, 1 };
		rnorm() = default;
		template<typename SeedSeq> rnorm(SeedSeq &&seed) : eng(seed) {
		}
		double r() {
			return normal_01(eng);
		}
	};

	inline double dnorm(double x, double mean, double sd, bool uselog) {
		double re;
		re = x - mean;
		re = uselog ? -0.918938533204673 - log(sd) - re*re / 2. / sd / sd : exp(-re*re / 2. / sd / sd) / 2.506628274631 / sd;
		return re;
	}

}
