#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd rmse(4);
	rmse.fill(0.0);

	if(estimations.size() == ground_truth.size())
	{
	
		for(int i=0; i<estimations.size(); i++)
		{
			VectorXd error(4);
			error = estimations[i] - ground_truth[i];
			error = error.array() * error.array();
			rmse+= error;
		}
		rmse/=estimations.size();
		rmse = rmse.array().sqrt();
	}
	else
	{
		cout<<"Size of estimations and ground truth matrices do not match."<<endl;
	}
	return rmse;
}

void Tools::restrict_angle(double &value)
{
	int sign = value/abs(value);
	while(value>M_PI || value<-M_PI)
	{
		value-= sign * 2 * M_PI;
	}
}