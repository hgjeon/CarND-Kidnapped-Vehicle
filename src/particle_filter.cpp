/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	// TODO: Set standard deviations for x, y, and theta

	double std_x, std_y, std_theta;

	num_particles = 100;

	std_x = std[0];
	std_y = std[1];
	std_theta = std[2];

	default_random_engine gen;
	// This line creates a normal (Gaussian) distribution for x
	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);

	Particle part;

	particles.clear();
	for (int i = 0; i < num_particles; ++i) {
		//double sample_x, sample_y, sample_theta;

		// TODO: Sample  and from these normal distrubtions like this:
		//	 sample_x = dist_x(gen);
		//	 where "gen" is the random engine initialized earlier.
		 part.x = dist_x(gen);
		 part.y = dist_y(gen);
		 part.theta = dist_theta(gen);
		 part.id = i;
		 particles.push_back(part);

		 // Print your samples to the terminal.
		 cout << "Init Particle Sample " << i + 1 << " " << part.x << " " << part.y << " " << part.theta << endl;
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;
	// This line creates a normal (Gaussian) distribution for x
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);

	for (int i = 0; i < num_particles; ++i) {
		//double sample_x, sample_y, sample_theta;
		double delta_x, delta_y, theta_o;
		// TODO: Sample  and from these normal distrubtions like this:
		//	 sample_x = dist_x(gen);
		//	 where "gen" is the random engine initialized earlier.

		theta_o = particles[i].theta;
		if (yaw_rate < 0.0001) {
			delta_x = velocity * delta_t * cos(theta_o);
			delta_y = velocity * delta_t * sin(theta_o);
		} else {
			delta_x = velocity / yaw_rate * (sin(theta_o + yaw_rate * delta_t) - sin(theta_o));
			delta_y = velocity / yaw_rate * (cos(theta_o) - cos(theta_o + yaw_rate * delta_t));
		}
		particles[i].x = particles[i].x + delta_x + dist_x(gen);
		particles[i].y = particles[i].y + delta_y + dist_y(gen);
		particles[i].theta = particles[i].theta + (yaw_rate * delta_t) + dist_theta(gen);

		//particles[i].x = particles[i].x + delta_x + 0;
		//particles[i].y = particles[i].y + delta_y + 0;
		//particles[i].theta = particles[i].theta + yaw_rate * delta_t + 0;

		// Print your samples to the terminal.
		//cout << "Init Particle Sample " << i + 1 << " " << sample_x << " " << sample_y << " " << sample_theta << endl;
	}

	// Print out
	//for (int i = 0; i < num_particles; ++i) {
	//	cout << "Predicted Particle:" << i << " " << particles[i].x << " " << particles[i].y << " " << particles[i].theta << endl;
	//}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

	//SetAssociations(id, sense_x, sense_y)
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	std::vector<LandmarkObs> predicted;

	double sig_x, sig_y;
	sig_x = std_landmark[0];
	sig_y = std_landmark[1];

	double weightSum = 0;

	for (int i = 0; i < num_particles; ++i) {
		double gauss_norm, exponent;
		double weightMulti;

		particles[i].associations.clear();
		particles[i].sense_x.clear();
		particles[i].sense_y.clear();

		// Update the closest Landmark ID in Association field, and get mapped coordination in sense_xy
		for (int obs_id = 0; obs_id < observations.size(); ++obs_id) {
			double x_obs, y_obs, theta, mapped_obs_x, mapped_obs_y, mapped_obs_id;
			x_obs = observations[obs_id].x;
			y_obs = observations[obs_id].y;
			theta = particles[i].theta;

			// transform to map x coordinate
			mapped_obs_x = particles[i].x + (cos(theta) * x_obs - sin(theta) * y_obs);
			// transform to map y coordinate
			mapped_obs_y = particles[i].y + (sin(theta) * x_obs + cos(theta) * y_obs);
			mapped_obs_id = observations[obs_id].id;

			//if (obs_id == 0)
			//	cout << "[mapped_obs_xy]:" << i << " " << mapped_obs_x << " " << mapped_obs_y << endl;

			// Find and set association corresponding to the current observation
			double min_dist = sensor_range * 1000;
			double min_land_id = -1;

			for (int land_id = 0; land_id < map_landmarks.landmark_list.size(); ++land_id)
			{
				double l_x, l_y, d_dist;
				l_x = map_landmarks.landmark_list[land_id].x_f;
				l_y = map_landmarks.landmark_list[land_id].y_f;
				d_dist = dist(mapped_obs_x, mapped_obs_y, l_x, l_y);
				if (d_dist < min_dist) {
					min_dist = d_dist;
					min_land_id = land_id;
					//min_land_id = map_landmarks.landmark_list[land_id].id_i;
				}
			}

			if (min_land_id != -1) {
				particles[i].associations.push_back(min_land_id);

				// transform to map x coordinate
				particles[i].sense_x.push_back(mapped_obs_x);
				// transform to map y coordinate
				particles[i].sense_y.push_back(mapped_obs_y);
			}
			else {
				assert(0);
			}
		}

		//dataAssociation(particles[i], observations);

		weightMulti = 1.0;
		for (int obs_id = 0; obs_id < particles[i].associations.size(); ++obs_id) {
			double x_obs, y_obs, mu_x, mu_y;

			x_obs = particles[i].sense_x[obs_id];
			y_obs = particles[i].sense_y[obs_id];
			//mu_x = map_landmarks.landmark_list[particles[i].associations[obs_id]].x_f;
			//mu_y = map_landmarks.landmark_list[particles[i].associations[obs_id]].y_f;
			mu_x = map_landmarks.landmark_list[particles[i].associations[obs_id]].x_f;
			mu_y = map_landmarks.landmark_list[particles[i].associations[obs_id]].y_f;

			// calculate normalization term
			gauss_norm = (1 / (2 * M_PI * sig_x * sig_y));

			// calculate exponent
			exponent = ((x_obs - mu_x) * (x_obs - mu_x)) / (2 * sig_x * sig_x) + ((y_obs - mu_y) * (y_obs - mu_y)) / (2 * sig_y * sig_y);

			// calculate weight using normalization terms and exponent
			weightMulti *= (gauss_norm * exp(-exponent));

			//if (obs_id == 0)
			//	cout << "[weightMulti xy_obs]:" << i << " obs_id " << obs_id << " " << (x_obs )  << " " << (y_obs )
			//			<< " mu_xy " << mu_x << " " << mu_y << endl;
			//cout << "[weightMulti dist]:" << i << obs_id << " " << (x_obs - mu_x)  << " " << (y_obs - mu_y) << endl;
			//cout << "[weightMulti obs]:" << i << obs_id << " " << (gauss_norm * exp(-exponent)) << " " << weightMulti << endl;

			assert(gauss_norm * exp(-exponent));

			// Update association with Landmark ID here
			particles[i].associations[obs_id] = (map_landmarks.landmark_list[particles[i].associations[obs_id]].id_i);

		}

		// Save weight for each particle
		particles[i].weight = weightMulti;
		weightSum += weightMulti;
	}

	// Normalize Weight
	weights.clear();
	for (int i = 0; i < num_particles; ++i) {
		//cout << "[weightSum]:" << i << " " << particles[i].weight << " " << weightSum << endl;
		weights.push_back(particles[i].weight / weightSum);
		particles[i].weight = particles[i].weight / weightSum;
	}

#if 0
	// Print out
	for (int i = 0; i < num_particles; ++i) {
		cout << "NormWeight Particle:" << i << " " << particles[i].x << " " << particles[i].y << " " << particles[i].theta
				<< " Weight" << particles[i].weight
		 		<< " Assoc" << particles[i].associations[0] << " Sense: " << particles[i].sense_x[0] << " " << particles[i].sense_y[0] << endl;
	}
#endif
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	default_random_engine gen;

	discrete_distribution<int> distrubution(weights.begin(), weights.end());
	vector<Particle> resampledParticle(num_particles);
	//vector<Particle> resampledParticle;

	int idx;
	for (int i = 0; i < num_particles; ++i) {
		idx = distrubution(gen);
#if 1
		resampledParticle.at(i) = particles.at(idx);
#else
		resampledParticle[i].x = particles[idx].x;
		resampledParticle[i].y = particles[idx].y;
		resampledParticle[i].theta = particles[idx].theta;
		resampledParticle[i].weight = particles[idx].weight;
		//resampledParticle[i].associations = particles[idx].associations;
		for (int obs_id = 0; obs_id < particles[i].associations.size(); ++obs_id) {
			resampledParticle[i].associations[obs_id] = particles[idx].associations[obs_id];
		}
#endif
		//cout << "[Resampled]:" << idx << " " << particles[idx].x << " " << particles[idx].y << " " << weights[i] << " " << endl;
		//cout << "[Resampled]:" << idx << " " << resampledParticle[i].x << " " << resampledParticle[i].y << " " << resampledParticle[i].theta << endl;

	}

#if 1
	particles = resampledParticle;
#else
	for (int i = 0; i < num_particles; ++i) {
		particles.at(idx) = resampledParticle.at(i);
		particles[i].x = resampledParticle[i].x;
		particles[i].y = resampledParticle[i].y;
		particles[i].theta = resampledParticle[i].theta;
		particles[i].weight = resampledParticle[i].weight;
		//particles[i].associations = resampledParticle[i].associations;
		for (int obs_id = 0; obs_id < particles[i].associations.size(); ++obs_id) {
			particles[idx].associations[obs_id] = resampledParticle[i].associations[obs_id];
		}
	}
#endif


	// Print out
	//for (int i = 0; i < num_particles; ++i) {
	//	cout << "weights:" << i + 1 << " " << weights[i] << endl;
	//}

	// Print out
	for (int i = 0; i < num_particles; ++i) {
		cout << "Resampled Particle:" << i << " " << particles[i].x << " " << particles[i].y << " " << particles[i].theta
				<< " Weight" << particles[i].weight
		 		<< " Assoc" << particles[i].associations[0] << " Sense: " << particles[i].sense_x[0] << " " << particles[i].sense_y[0] << endl;
	}
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
