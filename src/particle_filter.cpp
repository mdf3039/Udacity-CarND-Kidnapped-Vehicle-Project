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

void ParticleFilter::init(double x, double y, double theta, double std_pos[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 1000;

	// This line creates a normal (Gaussian) distribution for x, y, and theta
	default_random_engine gen;
	normal_distribution<double> dist_x(x, std_pos[0]);
	normal_distribution<double> dist_y(y, std_pos[1]);
	normal_distribution<double> dist_theta(theta, std_pos[2]);

	std::vector<double> weights(num_particles);
	std::vector<Particle> particles(num_particles);
	// for each particle, create the Particle and weight
	for (int i = 0; i < num_particles; ++i) {
        weights[i] = 1.0;
        particles[i].x = dist_x(gen);
        particles[i].y = dist_y(gen);
        particles[i].theta = dist_theta(gen);
        particles[i].weight = weights[i];
	}

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// This line creates a normal (Gaussian) distribution for x, y, and theta
	default_random_engine gen;
	normal_distribution<double> dist_x(0, std_pos[0]/2);
	normal_distribution<double> dist_y(0, std_pos[1]/2);
	normal_distribution<double> dist_theta(0, std_pos[2]/2);

	if (abs(yaw_rate)<=0.00001){
        for (int i = 0; i < num_particles; ++i) {
            particles[i].x += velocity/yaw_rate*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta))+dist_x(gen);
            particles[i].y += velocity/yaw_rate*(-1*cos(particles[i].theta+yaw_rate*delta_t)+cos(particles[i].theta))+dist_y(gen);
            particles[i].theta += yaw_rate*delta_t+dist_theta(gen);
        }
	}
    else{
        for (int i = 0; i < num_particles; ++i) {
            particles[i].x += velocity*cos(particles[i].theta)*delta_t+dist_x(gen);
            particles[i].y += velocity*sin(particles[i].theta)*delta_t+dist_y(gen);
            particles[i].theta += yaw_rate*delta_t+dist_theta(gen);
        }
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.


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

	//Update each particle separately
	for (int i = 0; i < num_particles; ++i) {
        //transform each noisy observation to for the specific particle
        //save observation in LandmarkObs vector structure
        std::vector<LandmarkObs> mapped_observations;
        for (int j = 0; j < observations.size(); ++j){
            LandmarkObs mapped_obs;
            mapped_obs.x = particles[i].x + (cos(particles[i].theta)*observations[j].x)-(sin(particles[i].theta)-observations[j].y);
            mapped_obs.y = particles[i].y + (sin(particles[i].theta)*observations[j].x)+(cos(particles[i].theta)-observations[j].y);
            mapped_observations.push_back(mapped_obs);
        }
        //obtain the distance from the first map landmark to all of the transformed
        //observations. save the distances in a vector and change the id in the
        //mapped observations to the first map landmark id
        std::vector<double> landmark_distances;
        for (int j = 0; j < mapped_observations.size(); ++j){
            landmark_distances.push_back(sqrt(pow(map_landmarks.landmark_list[0].x_f-mapped_observations[j].x,2)+pow(map_landmarks.landmark_list[0].y_f-mapped_observations[j].y,2)));
            mapped_observations[j].id = 0;
        }
        //for the other landmarks, obtain the distance from the landmark to each observation
        //if distance is less than the previous distance, change the id to the landmark
        //also change the distance in the landmark_distances vector to match
        double distance_to_landmark;
        for (int k = 1; k < map_landmarks.landmark_list.size(); ++k){
            for (int j = 0; j < mapped_observations.size(); ++j){
                distance_to_landmark = sqrt(pow(map_landmarks.landmark_list[k].x_f-mapped_observations[j].x,2)+pow(map_landmarks.landmark_list[k].y_f-mapped_observations[j].y,2));
                if (distance_to_landmark<landmark_distances[j]){
                    mapped_observations[j].id = k;
                    landmark_distances[j] = distance_to_landmark;
                }
            }
        }

	}


}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
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
