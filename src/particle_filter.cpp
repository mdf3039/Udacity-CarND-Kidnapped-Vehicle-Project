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

#include "Hungarian.h"
#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std_pos[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 300;

	// This line creates a normal (Gaussian) distribution for x, y, and theta
	default_random_engine gen;
	normal_distribution<double> dist_x(x, 3.0);
	normal_distribution<double> dist_y(y, 3.0);
	normal_distribution<double> dist_theta(theta, 0.1);

	//std::vector<double> weights(num_particles);
	//std::vector<Particle> particles(num_particles);
	// for each particle, create the Particle and weight
	for (int i = 0; i < num_particles; ++i) {
        weights.push_back(1.0);
        Particle gen_particle;
        gen_particle.x = dist_x(gen);
        gen_particle.y = dist_y(gen);
        gen_particle.theta = dist_theta(gen);
        //normalize theta to be between -pi and pi
        while (gen_particle.theta>M_PI){
            gen_particle.theta-=M_PI;
        }
        while (gen_particle.theta<-1.0*M_PI){
            gen_particle.theta+=M_PI;
        }
        gen_particle.weight = 1.0;
        particles.push_back(gen_particle);
	}
	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// This line creates a normal (Gaussian) distribution for x, y, and theta
	default_random_engine gen;
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);

	if (abs(yaw_rate)>0.001){
        for (int i = 0; i < num_particles; ++i) {
            particles[i].x += velocity/yaw_rate*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta))+dist_x(gen);
            particles[i].y += velocity/yaw_rate*(-1*cos(particles[i].theta+yaw_rate*delta_t)+cos(particles[i].theta))+dist_y(gen);
            particles[i].theta += yaw_rate*delta_t+dist_theta(gen);
        }
	}
    else{
        for (int i = 0; i < num_particles; ++i) {
            cout<<"X: "<<particles[i].x<<endl;
            particles[i].x += velocity*cos(particles[i].theta)*delta_t+dist_x(gen);
            cout<<"X_new: "<<particles[i].x<<endl;
            cout<<"Y: "<<particles[i].y<<endl;
            particles[i].y += velocity*sin(particles[i].theta)*delta_t+dist_y(gen);
            cout<<"Y_new: "<<particles[i].y<<endl;
            cout<<"Theta: "<<particles[i].theta<<endl;
            particles[i].theta += yaw_rate*delta_t+dist_theta(gen);
            cout<<"Theta_new: "<<particles[i].theta<<endl;
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
	double particle_weight_sum = 0;
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
        //obtain a subset of all landmarks in the vicinity. The sensor_range provides the lidar area
        Map subset_landmarks;
        for (int j=0; j<map_landmarks.landmark_list.size(); ++j){
            //if the distance from the car to the landmark is less than the sensor_range,
            //append to the subset_landmarks list
            if (sqrt(pow(particles[i].x-map_landmarks.landmark_list[j].x_f,2)+pow(particles[i].y-map_landmarks.landmark_list[j].y_f,2))<sensor_range+20){
                subset_landmarks.landmark_list.push_back(map_landmarks.landmark_list[j]);
            }
        }
        //obtain the distance from each landmark to an observation. Put in a vector
        //of vectors with length(vector)=#landmarks for #observations vectors.
        std::vector<std::vector<double>> DistMatrix;
        for (int j=0; j<mapped_observations.size();++j){
            std::vector<double> dist_vector;
            for (int k=0; k<subset_landmarks.landmark_list.size();++k){
                dist_vector.push_back(sqrt(pow(mapped_observations[j].x-subset_landmarks.landmark_list[k].x_f,2)+pow(mapped_observations[j].y-subset_landmarks.landmark_list[k].y_f,2)));
            }
            DistMatrix.push_back(dist_vector);
        }
        //use the Hungarian Algorithm to find the shortest distances connecting observations
        //to landmarks. Create dummy Assignment variable
        //the Solve function will produce a cost. Use this cost as the weight for resampling
        std::vector<int> hung_assignments;
        HungarianAlgorithm hung_alg;
        //double cost;
        hung_assignments = hung_alg.Solve(DistMatrix,hung_assignments);
        //hung_assignments = hung_alg.Assignment;
        /*for (int j=0; j<hung_assignments.size(); ++j){
            cout<<"Particle "<<i<<" assignment "<<j<<": "<<hung_assignments[j]<<endl;
        }*/
        //now that the assignments to each observation are found, find the distance and probabilities
        double particle_prob = 1.0;
        double x2_dist;
        double y2_dist;
        for (int j = 0; j < hung_assignments.size(); ++j){
            if (hung_assignments[j]<0){
                x2_dist = 400;
                y2_dist = 400;
            }
            else{
                x2_dist = pow(subset_landmarks.landmark_list[hung_assignments[j]].x_f-mapped_observations[j].x,2);
                y2_dist = pow(subset_landmarks.landmark_list[hung_assignments[j]].y_f-mapped_observations[j].y,2);
            }
            particle_prob *= exp(-1.0*(x2_dist+y2_dist)/(2*std_landmark[0]*std_landmark[0]));
        }
        //set the inverse cost value as the new weight for this particle
        particles[i].weight = particle_prob;
        //add the weight to the sum of the weights
        particle_weight_sum += particle_prob;
        /*
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
        //now that the closest landmark to each mapped_observation is known, multiply
        //probabilities together to obtain the weight for the particle
        double particle_prob = 1.0;
        double x2_dist;
        double y2_dist;
        for (int j = 0; j < mapped_observations.size(); ++j){
            x2_dist = pow(map_landmarks.landmark_list[mapped_observations[j].id].x_f-mapped_observations[j].x,2);
            y2_dist = pow(map_landmarks.landmark_list[mapped_observations[j].id].y_f-mapped_observations[j].y,2);
            particle_prob *= exp(-1.0*(x2_dist+y2_dist)/(2*std_landmark[0]*std_landmark[0]));
        }
        //set this value as the new weight for this particle
        particles[i].weight = particle_prob;
        //add the weight to the sum of the weights
        particle_weight_sum += particle_prob;
        */
	}
    //Normalize all weight values based on the sum
    for (int i = 0; i < num_particles; ++i) {
        particles[i].weight /= particle_weight_sum;
    }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	//obtain the random number generator and a vector of all of the weights
    std::default_random_engine generator;
    std::vector<double> weight_vector;
    //cout<<"Particle weights in vector: ";
    for (int i = 0; i < num_particles; ++i){
        weight_vector.push_back(particles[i].weight);
        //cout<<particles[i].weight<<" ";
    }
    //cout<<endl;
    //obtain the discrete distribution according to the weight_vector
    std::discrete_distribution<> distribution(weight_vector.begin(), weight_vector.end());
    //Sample with replacement from this distribution, appending the new particles
    //to a list
    //cout<<"Numbers generated: ";
    std::vector<Particle> particles_updated;
    int number_generated;
    Particle generated_particle;
    for (int i = 0; i < num_particles; ++i){
        number_generated = distribution(generator);
        //cout<<number_generated<<" ";
        generated_particle.theta = particles[number_generated].theta;
        //cout<<generated_particle.theta<<" ";
        generated_particle.x = particles[number_generated].x;
        //cout<<generated_particle.x<<" ";
        generated_particle.y = particles[number_generated].y;
        //cout<<generated_particle.y<<" ";
        generated_particle.weight = particles[number_generated].weight;
        particles_updated.push_back(generated_particle);
    }
    /*cout<<endl;
    //replace the old list with the new list and delete the old list
    cout<<"Old particle list: ";
    for (int i = 0; i < num_particles; ++i){
        cout<<"("<<particles[i].x<<","<<particles[i].y<<","
            <<particles[i].weight<<","<<particles[i].theta<<") ";
    }
    cout<<endl;
    cout<<"New particle list: ";
    for (int i = 0; i < num_particles; ++i){
        cout<<"("<<particles_updated[i].x<<","<<particles_updated[i].y<<","
            <<particles_updated[i].weight<<","<<particles_updated[i].theta<<") ";
    }
    cout<<endl;
    */
    particles.swap(particles_updated);
    /*cout<<"New old particle list: ";
    for (int i = 0; i < num_particles; ++i){
        cout<<"("<<particles[i].x<<","<<particles[i].y<<","
            <<particles[i].weight<<","<<particles[i].theta<<") ";
    }
    cout<<endl;
    */

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
