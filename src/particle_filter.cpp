/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 1000;  // TODO: Set the number of particles

  //Random engine and distributions
  std::default_random_engine gen;
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);


  //Initialize all the particles
  for(int i=0; i<num_particles; i++)
  {
	  Particle particle;
	  particle.weight = 1;
	  particle.id = i;
	  particle.x = dist_x(gen);
	  particle.y = dist_y(gen);
	  particle.theta = dist_theta(gen);
	  particles.push_back(particle);
	  weights.push_back(particle.weight);
  }

  //Raise initialization flag
  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
	  //Random engine and distributions
	  std::default_random_engine gen;
//	  std::normal_distribution<double> dist_velocity(0.0, std_pos[0]);
//	  std::normal_distribution<double> dist_yaw_rate(0.0, std_pos[2]);
	  std::normal_distribution<double> dist_x(0.0, std_pos[0]);
	  std::normal_distribution<double> dist_y(0.0, std_pos[1]);
	  std::normal_distribution<double> dist_theta(0.0, std_pos[1]);

	for(int i=0; i<num_particles; i++)
	{
		Particle & particle = particles[i];

		//Calculate spatial position
		double xFinal, yFinal, thetaFinal;
		if(yaw_rate == 0)
		{
			xFinal = particle.x + velocity * delta_t * cos(particle.theta);
			yFinal = particle.y + velocity * delta_t * sin(particle.theta);
			thetaFinal = particle.theta;
		}
		else
		{
			xFinal = particle.x + velocity/ yaw_rate* (sin(particle.theta + yaw_rate* delta_t) - sin(particle.theta)) + dist_x(gen);
			yFinal = particle.y + velocity/ yaw_rate* (cos(particle.theta) - cos(particle.theta + yaw_rate* delta_t)) + dist_y(gen);
			thetaFinal = particle.theta + yaw_rate* delta_t + dist_theta(gen);
		}

		particle.x = xFinal;
		particle.y = yFinal;
		particle.theta = thetaFinal;
	}
	std::cout << "particle: " << particles[0].x << std::endl;
	int a = 0;
}

//void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted,
//                                     vector<LandmarkObs>& observations) {
//  /**
//   * TODO: Find the predicted measurement that is closest to each
//   *   observed measurement and assign the observed measurement to this
//   *   particular landmark.
//   * NOTE: this method will NOT be called by the grading code. But you will
//   *   probably find it useful to implement this method and use it as a helper
//   *   during the updateWeights phase.
//   */
//
//}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

	//Clear weights from previous iteration..
	weights.clear();

	for(int p=0; p<particles.size(); p++)
	{
		//Get the reference to current particle
		Particle & particle = particles[p];

		//First of all, clear the particle's previous associations
		particle.partial_weights.clear();
		particle.associations.clear();
		particle.sense_x.clear();
		particle.sense_y.clear();

		for(int l=0; l<observations.size(); l++)
		{
			LandmarkObs const & currentObs = observations[l];

			//First translate measurements to map coordinates if they were measured
			//in the reference of frame of particle
			double transformed_obs_x, transformed_obs_y;
			ParticleFilter::transformer((Particle const &)particle,
					currentObs.x, currentObs.y,
					transformed_obs_x, transformed_obs_y);

			//Push transformed coordinates to sense_x and sense_y vectors
			particle.sense_x.push_back(transformed_obs_x);
			particle.sense_y.push_back(transformed_obs_y);

			//Second associate current measured landmark to one of (min distance) landmarks in map
			double min_distance = std::numeric_limits<double>::max();
			int min_matched_landmark_id = -1;
			double min_map_landmark_x;
			double min_map_landmark_y;
			for(int m=0; m<map_landmarks.landmark_list.size(); m++)
			{
				double map_landmark_x = map_landmarks.landmark_list[m].x_f;
				double map_landmark_y = map_landmarks.landmark_list[m].y_f;
				int map_landmark_id = map_landmarks.landmark_list[m].id_i;
				double distance = ParticleFilter::euclidean_distance(transformed_obs_x, transformed_obs_y,
						map_landmark_x, map_landmark_y);
				if (distance < min_distance)
				{
					min_distance = distance;
					min_matched_landmark_id = map_landmark_id;
					min_map_landmark_x = map_landmark_x;
					min_map_landmark_y = map_landmark_y;
				}
			}

			//At this point a map landmark is found that is the nearest to the current measurement
			//Calculate the probability for this observation, this will be used to calculate particle's
			//overall weight in future..
			//First of all push this association to vector
			particle.associations.push_back(min_matched_landmark_id);

			//Calculate probability
			double weight_partial = ParticleFilter::multiv_prob(std_landmark[0], std_landmark[1],
					transformed_obs_x, transformed_obs_y,
					min_map_landmark_x, min_map_landmark_y
					);
			particle.partial_weights.push_back(weight_partial);



		}

		//Third calculate the new weight by using multivariate gaussian distributions..
		double final_weight = 1.0;
		for(int i=0; i<particle.partial_weights.size(); i++)
		{
			final_weight *= particle.partial_weights[i];
		}

		//Update particle's weight
		particle.weight = final_weight;
		weights.push_back(final_weight);
	}
	int a  = 0;

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

	std::vector<Particle> particles_updated;
	std::discrete_distribution<int> dist(weights.begin(), weights.end());
	std::default_random_engine gen;
	for(int i=0; i<num_particles; i++)
	{
		int index = dist(gen);
		particles_updated.push_back(particles[index]);
	}

	//Swap modified particles list and previous particles list
	particles = particles_updated;
}

void ParticleFilter::transformer(Particle const & particle,
		double obs_x, double obs_y,
		double & x_map, double & y_map)
{
	  // define coordinates and theta
	  double x_part, y_part, theta;
	  x_part = particle.x;
	  y_part = particle.y;
	  theta = particle.theta;

	  // transform to map x coordinate
	  x_map = x_part + (cos(theta) * obs_x) - (sin(theta) * obs_y);

	  // transform to map y coordinate
	  y_map = y_part + (sin(theta) * obs_x) + (cos(theta) * obs_y);
}

double ParticleFilter::euclidean_distance(double x1, double y1, double x2, double y2)
{
	return sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2));
}

double ParticleFilter::multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs,
        double mu_x, double mu_y)
{
	  // calculate normalization term
	  double gauss_norm;
	  gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

	  // calculate exponent
	  double exponent;
	  exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
	               + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));

	  double expRes = exp(-exponent);
	  // calculate weight using normalization terms and exponent
	  double weight;
	  weight = gauss_norm * expRes;

	  return weight;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
