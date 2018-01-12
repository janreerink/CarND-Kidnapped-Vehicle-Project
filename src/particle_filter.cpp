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
#include <map>

#include "particle_filter.h"

using namespace std;

//set random engine
static default_random_engine gen;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    
    //set number of particles and initialize vector
    num_particles = 100; 
    particles.reserve(num_particles);
    weights.reserve(num_particles);
    
    
    // get standard deviations from array
    double std_x, std_y, std_theta;
    std_x = std[0];
    std_y =  std[1];
    std_theta =  std[2];
    
    // create normal distributions for x,y,theta
	normal_distribution<double> dist_x(0, std_x);
	normal_distribution<double> dist_y(0, std_y);
	normal_distribution<double> dist_theta(0, std_theta);

    // initialize particles by sampling from distributions, add them to list of particles
	for (int i = 0; i < num_particles; i++) {

        Particle current_p ;
        current_p.id = i;
        current_p.x = x + dist_x(gen);
        current_p.y = y + dist_y(gen);
        current_p.theta = theta + dist_theta(gen);
        current_p.weight = 1.;

        weights.push_back(1.) ;
        particles.push_back(current_p);
	}
    is_initialized = true;
}


void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    
    // create normal distributions for x,y,theta
	normal_distribution<double> di_x(0, std_pos[0]);
	normal_distribution<double> di_y(0, std_pos[1]);
	normal_distribution<double> di_theta(0, std_pos[2]);
    
    //loop through particles, predict, add noise
    for (int i = 0; i < num_particles; i++) {
        if (fabs(yaw_rate) < 0.0000001) {
        particles[i].x = particles[i].x + (cos(particles[i].theta) * velocity * delta_t) ;
        particles[i].y = particles[i].y + (sin(particles[i].theta) * velocity * delta_t) ;
        //theta remains unchanged
        }
        else { //yaw rate zero
        particles[i].x = particles[i].x + velocity/yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta) ) ;
        particles[i].y = particles[i].y + velocity/yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t) ) ;
        particles[i].theta = (particles[i].theta + ( yaw_rate * delta_t)) ;
        }

        // noise
        particles[i].x = particles[i].x + di_x(gen);
        particles[i].y = particles[i].y + di_y(gen);
        particles[i].theta = particles[i].theta + di_theta(gen);
	}

}


void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
    
    // predicted landmarks: all landmarks within sensor range for a given particle
    // obsrvations: lidar measurements from vehicle
    // assign to lidar sensor measurements the closest landmark id

    
    double current_dist;
    current_dist = 0;
    
    
    for (int n = 0; n < observations.size(); n++) {
        // var for storage of closest landmark dist, initialize to some large value
        double min_dist; 
        min_dist = 1000000.0;
        int mindist_idx;
        mindist_idx = -1 ; 
        
        for (int i = 0; i < predicted.size(); i++) {
            current_dist = dist(predicted[i].x,predicted[i].y, observations[n].x,observations[n].y); 
            if (current_dist < min_dist) {
                min_dist = current_dist;
                mindist_idx = predicted[i].id;
            }
        observations[n].id = mindist_idx;
        }
    }
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
    
    // transform observed landmarks to map coordinates
    
    // for each particle predict measurements to landmarks within range
    
    double sig_x = std_landmark[0];
    double sig_y = std_landmark[1];
    for (int p = 0; p < particles.size(); p++) { // for each particle p
        // init vector of predicted landmarks for this particle
        std::vector<LandmarkObs> p_predicted;
        double theta_p = particles[p].theta;
        double x_p = particles[p].x;
        double y_p = particles[p].y;
        double current_dist;
        for (int o = 0; o < map_landmarks.landmark_list.size(); o++) { // for each lm
            float o_x = map_landmarks.landmark_list[o].x_f;
            float o_y = map_landmarks.landmark_list[o].y_f;
            int o_id = map_landmarks.landmark_list[o].id_i;

        	current_dist = dist(particles[p].x,particles[p].y, o_x,o_y);
            if (current_dist <= sensor_range) {
            //if (fabs(o_x - x_p) <= sensor_range && fabs(o_y - y_p) <= sensor_range) {
                p_predicted.push_back(LandmarkObs{ o_id, o_x, o_y  });

                }
        }
        // transform lidar measurements from vehicle coords to map coord system
        std:vector<LandmarkObs> observations_m;
        observations_m.reserve(observations.size());
        for (int i = 0; i < observations.size(); i++) {
            double x_c = observations[i].x;
            double y_c = observations[i].y;
            int id_c = observations[i].id;
            
            double x_m = x_c * cos(theta_p) - y_c * sin(theta_p) + x_p;
            double y_m = x_c * sin(theta_p) + y_c * cos(theta_p) + y_p;
            //double x_m = cos(theta_p)*observations[i].x - sin(theta_p)*observations[i].y + x_p;
            //double y_m = sin(theta_p)*observations[i].x + cos(theta_p)*observations[i].y + y_p;
            observations_m.push_back(LandmarkObs{id_c, x_m, y_m}); 
        }
        // assign to each sensor measurement ID of closest landmark
        dataAssociation(p_predicted, observations_m);
        // update particle weight
        
        double product_probs = 1.0; 
        double weight = 1.0;
        particles[p].weight = 1.0;
        for (int i = 0; i < observations_m.size(); i++) {
            double x_o = observations_m[i].x;
            double y_o = observations_m[i].y;
            int closest_id = observations_m[i].id;
            // nearest landmark coords
            double mu_x, mu_y; 
            for (int a = 0; a < p_predicted.size(); a++) {
                if (p_predicted[a].id == closest_id) {
                    mu_x = p_predicted[a].x;
                    mu_y = p_predicted[a].y;
                }
            }
            // multi variate gaussian normalization term
            double gauss_norm = (1/(2 * M_PI * sig_x * sig_y));
            // multi variate gaussian exponent
            double exponent = ( pow(mu_x- x_o ,2))/(2 * sig_x*sig_x) + (pow(mu_y-y_o,2))/(2 * sig_y*sig_y);
            weight = gauss_norm * exp(-exponent);
            particles[p].weight *= weight;
        

        }
        // update vector of unnormlaized weights and particle weight
        weights[p] = particles[p].weight;
        
    } // end particle loop
}


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    
    //use discrete distribution to sample weight idx i with probability weights[i]/sum(weights)
    //std::random_device rd;
    //std::mt19937 gen2(rd());
    std::discrete_distribution<int> d(weights.begin(), weights.end());
    
    // create new particles vector and populate
    std::vector<Particle> particles_re(num_particles);
    for(int n = 0; n < num_particles; n++) {
    	int idx = d(gen);
    	particles_re[n] = particles[idx];
    }

    particles = particles_re;
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
